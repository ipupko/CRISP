#!/usr/bin/env python3
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 1: File Validation, MD5 and Input Summary
# Version: 0.2.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
##########################################################################

import argparse
import hashlib
import os
import sys
import datetime
import json
from pathlib import Path


# ─────────────────────────────────────────────
# DEFAULTS (overridden by instruction file)
# ─────────────────────────────────────────────
DEFAULTS = {
    "INPUT_FORMAT"   : None,       # REQUIRED
    "INPUT_PATH"     : None,       # REQUIRED
    "OUTPUT_DIR"     : "./results",
    "PROJECT_NAME"   : "project",
    "RUN_QC"         : "YES",
    "RUN_ANEUPLOIDY" : "YES",
    "RUN_SEX_CHECK"  : "YES",
    "RUN_PCA"        : "YES",
    "RUN_AMEND"      : "YES",
    "MIND"           : "0.05",
    "GENO"           : "0.05",
    "MAF"            : "0.01",
    "HWE"            : "1e-6",
    "HET_SD"         : "3",
    "HOM_Z_HIGH"     : "3",
    "HOM_Z_LOW"      : "-2",
    "XXX_F_THRESHOLD": "-0.15",
    "REMOVE_MONO"    : "YES",
    "PLINK1_PATH"    : "plink",
    "PLINK2_PATH"    : "plink2",
    "RSCRIPT_PATH"   : "Rscript",
}

# Extensions expected per format
FORMAT_EXTENSIONS = {
    "PLINK" : [".bed", ".bim", ".fam"],
    "PED"   : [".ped", ".map"],
    "VCF"   : [".vcf", ".vcf.gz"],
    "BCF"   : [".bcf"],
    "BGEN"  : [".bgen", ".sample"],
}


# ─────────────────────────────────────────────
# INSTRUCTION FILE PARSER
# ─────────────────────────────────────────────
def parse_instruction_file(config_path: str) -> dict:
    """
    Parse key = value instruction file.
    Lines starting with # are comments.
    Returns merged dict of defaults + user overrides.
    """
    config = DEFAULTS.copy()

    if not os.path.isfile(config_path):
        abort(f"Instruction file not found: {config_path}")

    with open(config_path, "r") as fh:
        for line_num, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                warn(f"Line {line_num} skipped (no '=' found): {raw.rstrip()}")
                continue
            key, _, value = line.partition("=")
            key   = key.strip().upper()
            # Strip inline comments (everything after first #)
            value = value.split("#")[0].strip()
            if key in config or key in DEFAULTS:
                config[key] = value
            else:
                warn(f"Line {line_num}: unrecognised key '{key}' — ignored")

    # Validate required fields
    if not config.get("INPUT_FORMAT"):
        abort("INPUT_FORMAT is required in the instruction file.")
    if not config.get("INPUT_PATH"):
        abort("INPUT_PATH is required in the instruction file.")

    config["INPUT_FORMAT"] = config["INPUT_FORMAT"].upper()
    if config["INPUT_FORMAT"] not in FORMAT_EXTENSIONS:
        abort(f"INPUT_FORMAT '{config['INPUT_FORMAT']}' not recognised. "
              f"Valid options: {', '.join(FORMAT_EXTENSIONS.keys())}")

    return config


# ─────────────────────────────────────────────
# MD5 COMPUTATION
# ─────────────────────────────────────────────
def md5_file(filepath: str, chunk_size: int = 1024 * 1024) -> str:
    """Compute MD5 checksum of a file in chunks (memory-safe for large files)."""
    h = hashlib.md5()
    with open(filepath, "rb") as fh:
        while chunk := fh.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


# ─────────────────────────────────────────────
# FILE RESOLUTION
# ─────────────────────────────────────────────
def resolve_input_files(config: dict) -> list:
    """
    Given INPUT_FORMAT and INPUT_PATH, return list of expected file paths.
    For PLINK/PED: INPUT_PATH is the prefix (no extension).
    For VCF/BCF/BGEN: INPUT_PATH is the full file path.
    """
    fmt        = config["INPUT_FORMAT"]
    input_path = config["INPUT_PATH"]
    extensions = FORMAT_EXTENSIONS[fmt]
    files      = []

    if fmt in ("PLINK", "PED"):
        # Prefix-based: append each expected extension
        for ext in extensions:
            files.append(input_path + ext)
    else:
        # Single file (VCF, BCF, BGEN)
        # Try exact path first, then try each known extension
        if os.path.isfile(input_path):
            files.append(input_path)
            # For BGEN, also look for companion .sample file
            if fmt == "BGEN":
                sample_file = os.path.splitext(input_path)[0] + ".sample"
                if os.path.isfile(sample_file):
                    files.append(sample_file)
        else:
            # Try appending extensions (user may have omitted extension)
            for ext in extensions:
                candidate = input_path + ext
                if os.path.isfile(candidate):
                    files.append(candidate)

    return files


# ─────────────────────────────────────────────
# FORMAT-SPECIFIC METADATA EXTRACTION
# ─────────────────────────────────────────────
def count_plink_samples_variants(prefix: str) -> dict:
    """Count samples from .fam and variants from .bim."""
    counts = {"samples": "N/A", "variants": "N/A"}
    fam = prefix + ".fam"
    bim = prefix + ".bim"
    if os.path.isfile(fam):
        with open(fam) as f:
            counts["samples"] = sum(1 for line in f if line.strip())
    if os.path.isfile(bim):
        with open(bim) as f:
            counts["variants"] = sum(1 for line in f if line.strip())
    return counts


def count_ped_samples_variants(prefix: str) -> dict:
    """Count samples from .ped (one per line) and variants from .map."""
    counts = {"samples": "N/A", "variants": "N/A"}
    ped = prefix + ".ped"
    map_ = prefix + ".map"
    if os.path.isfile(ped):
        with open(ped) as f:
            counts["samples"] = sum(1 for line in f if line.strip())
    if os.path.isfile(map_):
        with open(map_) as f:
            counts["variants"] = sum(1 for line in f if line.strip())
    return counts


def count_vcf_samples_variants(vcf_path: str) -> dict:
    """Count samples from VCF header (#CHROM line) and variants (non-header lines)."""
    counts  = {"samples": "N/A", "variants": "N/A"}
    n_vars  = 0
    opener  = open

    # Handle gzipped VCF
    if vcf_path.endswith(".gz"):
        import gzip
        opener = gzip.open

    try:
        with opener(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    # Samples start at column index 9
                    counts["samples"] = max(0, len(cols) - 9)
                elif not line.startswith("#"):
                    n_vars += 1
        counts["variants"] = n_vars
    except Exception as e:
        warn(f"Could not parse VCF metadata: {e}")
    return counts


def count_bgen_samples(bgen_path: str) -> dict:
    """Count samples from companion .sample file."""
    counts      = {"samples": "N/A", "variants": "N/A (BGEN variant count requires bgenix)"}
    sample_file = os.path.splitext(bgen_path)[0] + ".sample"
    if os.path.isfile(sample_file):
        with open(sample_file) as f:
            lines = [l for l in f if l.strip()]
        # .sample format: 2 header lines then one row per sample
        counts["samples"] = max(0, len(lines) - 2)
    return counts


def get_format_metadata(config: dict) -> dict:
    """Dispatch to correct metadata extractor based on format."""
    fmt  = config["INPUT_FORMAT"]
    path = config["INPUT_PATH"]
    if fmt == "PLINK":
        return count_plink_samples_variants(path)
    elif fmt == "PED":
        return count_ped_samples_variants(path)
    elif fmt in ("VCF", "BCF"):
        # Resolve actual file path
        for ext in FORMAT_EXTENSIONS[fmt]:
            candidate = path if os.path.isfile(path) else path + ext
            if os.path.isfile(candidate):
                return count_vcf_samples_variants(candidate)
        return {"samples": "N/A", "variants": "N/A"}
    elif fmt == "BGEN":
        candidate = path if os.path.isfile(path) else path + ".bgen"
        return count_bgen_samples(candidate)
    return {"samples": "N/A", "variants": "N/A"}


# ─────────────────────────────────────────────
# REPORT WRITER
# ─────────────────────────────────────────────
def write_summary_report(config: dict, file_records: list, metadata: dict,
                         missing_files: list, output_dir: str):
    """Write human-readable input summary report and machine-readable JSON."""

    project     = config.get("PROJECT_NAME", "project")
    run_time    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_path = os.path.join(output_dir, f"{project}_step1_input_summary.txt")
    json_path   = os.path.join(output_dir, f"{project}_step1_input_summary.json")

    separator = "=" * 65

    with open(report_path, "w") as out:
        out.write(f"{separator}\n")
        out.write(f"  GENETIC DATA PIPELINE — INPUT SUMMARY REPORT\n")
        out.write(f"{separator}\n")
        out.write(f"  Project      : {project}\n")
        out.write(f"  Run datetime : {run_time}\n")
        out.write(f"  Input format : {config['INPUT_FORMAT']}\n")
        out.write(f"  Input path   : {config['INPUT_PATH']}\n")
        out.write(f"{separator}\n\n")

        out.write(f"FILE INVENTORY\n")
        out.write(f"{'-' * 65}\n")
        if file_records:
            for rec in file_records:
                out.write(f"  File    : {rec['filename']}\n")
                out.write(f"  Path    : {rec['path']}\n")
                out.write(f"  Size    : {rec['size_mb']} MB\n")
                out.write(f"  MD5     : {rec['md5']}\n")
                out.write(f"\n")
        else:
            out.write("  No input files found.\n\n")

        if missing_files:
            out.write(f"MISSING FILES\n")
            out.write(f"{'-' * 65}\n")
            for mf in missing_files:
                out.write(f"  MISSING : {mf}\n")
            out.write("\n")

        out.write(f"DATASET DIMENSIONS\n")
        out.write(f"{'-' * 65}\n")
        out.write(f"  Samples  : {metadata.get('samples', 'N/A')}\n")
        out.write(f"  Variants : {metadata.get('variants', 'N/A')}\n\n")

        out.write(f"PIPELINE PARAMETERS (from instruction file)\n")
        out.write(f"{'-' * 65}\n")
        param_keys = ["RUN_QC", "RUN_ANEUPLOIDY", "RUN_SEX_CHECK", "RUN_PCA",
                      "RUN_AMEND", "MIND", "GENO", "MAF", "HWE", "HET_SD",
                      "REMOVE_MONO", "PLINK1_PATH", "PLINK2_PATH"]
        for k in param_keys:
            out.write(f"  {k:<20} : {config.get(k, 'default')}\n")
        out.write(f"\n{separator}\n")
        out.write(f"  END OF REPORT\n")
        out.write(f"{separator}\n")

    # Machine-readable JSON for downstream scripts to consume
    json_data = {
        "project"       : project,
        "run_datetime"  : run_time,
        "input_format"  : config["INPUT_FORMAT"],
        "input_path"    : config["INPUT_PATH"],
        "output_dir"    : output_dir,
        "files"         : file_records,
        "missing_files" : missing_files,
        "samples"       : metadata.get("samples"),
        "variants"      : metadata.get("variants"),
        "config"        : config,
    }
    with open(json_path, "w") as jf:
        json.dump(json_data, jf, indent=2)

    return report_path, json_path


# ─────────────────────────────────────────────
# LOGGING HELPERS
# ─────────────────────────────────────────────
def log(msg: str):
    print(f"[STEP1] {msg}", flush=True)

def warn(msg: str):
    print(f"[STEP1 WARNING] {msg}", flush=True)

def abort(msg: str):
    print(f"[STEP1 ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Validate input files, compute MD5s, write summary report."
    )
    parser.add_argument("--config", required=True,
                        help="Path to pipeline instruction file")
    args = parser.parse_args()

    # 1. Parse instruction file
    log("Parsing instruction file...")
    config = parse_instruction_file(args.config)
    log(f"  Project      : {config['PROJECT_NAME']}")
    log(f"  Input format : {config['INPUT_FORMAT']}")
    log(f"  Input path   : {config['INPUT_PATH']}")

    # 2. Set up output directory
    output_dir = config["OUTPUT_DIR"]
    os.makedirs(output_dir, exist_ok=True)
    log(f"  Output dir   : {output_dir}")

    # 3. Resolve expected input files
    log("Resolving input files...")
    expected_files = resolve_input_files(config)

    file_records  = []
    missing_files = []

    for filepath in expected_files:
        if not os.path.isfile(filepath):
            warn(f"Expected file not found: {filepath}")
            missing_files.append(filepath)
            continue

        log(f"  Found: {filepath}")
        size_bytes = os.path.getsize(filepath)
        size_mb    = round(size_bytes / (1024 * 1024), 3)

        log(f"    Computing MD5... (size: {size_mb} MB)")
        checksum = md5_file(filepath)
        log(f"    MD5: {checksum}")

        file_records.append({
            "filename" : os.path.basename(filepath),
            "path"     : os.path.abspath(filepath),
            "size_mb"  : size_mb,
            "md5"      : checksum,
        })

    # 4. Abort if critical files are missing
    if missing_files:
        # For PLINK, all three files are required
        if config["INPUT_FORMAT"] in ("PLINK", "PED"):
            abort(f"{len(missing_files)} required input file(s) missing. "
                  f"Cannot proceed. Check INPUT_PATH in instruction file.")
        else:
            warn(f"{len(missing_files)} file(s) missing — proceeding with found files.")

    # 5. Extract format-specific metadata
    log("Extracting dataset dimensions...")
    metadata = get_format_metadata(config)
    log(f"  Samples  : {metadata.get('samples', 'N/A')}")
    log(f"  Variants : {metadata.get('variants', 'N/A')}")

    # 6. Write summary report
    log("Writing summary report...")
    report_path, json_path = write_summary_report(
        config, file_records, metadata, missing_files, output_dir
    )
    log(f"  Text report : {report_path}")
    log(f"  JSON export : {json_path}")

    # 7. Save parsed config as JSON for downstream steps to consume
    config_out = os.path.join(output_dir, f"{config['PROJECT_NAME']}_parsed_config.json")
    with open(config_out, "w") as cf:
        json.dump(config, cf, indent=2)
    log(f"  Parsed config saved : {config_out}")

    log("Step 1 complete.")

    # Exit non-zero if files were missing (SLURM will catch this)
    if missing_files and config["INPUT_FORMAT"] in ("PLINK", "PED"):
        sys.exit(1)


if __name__ == "__main__":
    main()
