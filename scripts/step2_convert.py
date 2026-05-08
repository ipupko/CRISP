#!/usr/bin/env python3
"""
CRISP - Comprehensive Robust Integrated SNP Processing
Step 2: Format Conversion
Version: 0.2.0

Developed by Igor Pupko
https://github.com/ipupko/CRISP
Date Updated : 03/05/2026

Converts input genotype files to PLINK BED/BIM/FAM format.
Supported input formats: PLINK (pass-through), PED/MAP, VCF, BCF, BGEN.

Chromosome cleaning retains 1-22, X, Y, XY, MT only.
Multiallelic variants excluded via --max-alleles 2 (VCF and BGEN).
Intermediate files optionally retained via KEEP_INTERMEDIATE.

Usage:
    python3 step2_convert.py --config <parsed_config.json>
"""

import argparse
import json
import os
import subprocess
import sys
import datetime
import hashlib


# ─────────────────────────────────────────────
# VALID CHROMOSOMES
# ─────────────────────────────────────────────
VALID_CHR = [str(i) for i in range(1, 23)] + ["X", "Y", "XY", "MT"]
VALID_CHR_PREFIXED = ["chr" + c for c in VALID_CHR]
ALL_VALID = set(VALID_CHR + VALID_CHR_PREFIXED)


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
def log(msg):
    print(f"[STEP2] {msg}", flush=True)

def warn(msg):
    print(f"[STEP2 WARNING] {msg}", flush=True)

def abort(msg):
    print(f"[STEP2 ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


# ─────────────────────────────────────────────
# MD5
# ─────────────────────────────────────────────
def md5_file(filepath, chunk_size=1024 * 1024):
    h = hashlib.md5()
    with open(filepath, "rb") as fh:
        while chunk := fh.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


# ─────────────────────────────────────────────
# RUN SHELL COMMAND
# ─────────────────────────────────────────────
def run_cmd(cmd, log_file=None):
    log(f"Running: {' '.join(cmd)}")
    with open(log_file, "w") if log_file else open(os.devnull, "w") as lf:
        result = subprocess.run(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            text=True
        )
    return result.returncode


# ─────────────────────────────────────────────
# COUNT BIM VARIANTS AND CHROMOSOMES
# ─────────────────────────────────────────────
def count_bim(bim_path):
    """Read .bim and return total count, per-chr counts, and non-standard chrs."""
    chr_counts   = {}
    non_standard = {}
    total        = 0

    with open(bim_path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            chrom = parts[0]
            total += 1
            chr_counts[chrom] = chr_counts.get(chrom, 0) + 1
            if chrom not in ALL_VALID:
                non_standard[chrom] = non_standard.get(chrom, 0) + 1

    return total, chr_counts, non_standard


# ─────────────────────────────────────────────
# WRITE CHROMOSOME EXTRACT LIST
# ─────────────────────────────────────────────
def write_chr_list(bim_path, out_path):
    """Write variant IDs on standard chromosomes only."""
    retained = []
    excluded = 0

    with open(bim_path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            chrom  = parts[0]
            var_id = parts[1]
            if chrom in ALL_VALID:
                retained.append(var_id)
            else:
                excluded += 1

    with open(out_path, "w") as f:
        for vid in retained:
            f.write(vid + "\n")

    return len(retained), excluded


# ─────────────────────────────────────────────
# CONVERSION FUNCTIONS
# ─────────────────────────────────────────────

def convert_ped(config, converted_dir, log_dir):
    """PED/MAP to BED/BIM/FAM via PLINK 1.9."""
    plink1     = config.get("PLINK1_PATH", "plink")
    prefix     = config["INPUT_PATH"]
    out_prefix = os.path.join(converted_dir, config["PROJECT_NAME"] + "_converted")
    log_file   = os.path.join(log_dir, "step2_ped_conversion.log")

    cmd = [
        plink1,
        "--ped", prefix + ".ped",
        "--map", prefix + ".map",
        "--make-bed",
        "--allow-extra-chr",
        "--out", out_prefix
    ]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        abort(f"PED/MAP conversion failed. Check log: {log_file}")

    return out_prefix


def convert_vcf(config, converted_dir, log_dir):
    """VCF or BCF to BED/BIM/FAM via PLINK 2."""
    plink2     = config.get("PLINK2_PATH", "plink2")
    fmt        = config["INPUT_FORMAT"].upper()
    path       = config["INPUT_PATH"]
    out_prefix = os.path.join(converted_dir, config["PROJECT_NAME"] + "_converted")
    log_file   = os.path.join(log_dir, "step2_vcf_conversion.log")

    # try plain path first, then common extensions
    if fmt == "VCF":
        candidates = [path, path + ".vcf.gz", path + ".vcf"]
    else:
        candidates = [path, path + ".bcf"]

    vcf_file = next((c for c in candidates if os.path.isfile(c)), None)

    if vcf_file is None:
        abort(
            f"No {'VCF' if fmt == 'VCF' else 'BCF'} file found.\n"
            f"  Tried:\n" +
            "\n".join(f"    {c}" for c in candidates)
        )

    flag = "--vcf" if fmt == "VCF" else "--bcf"

    cmd = [
        plink2,
        flag, vcf_file,
        "--make-bed",
        "--allow-extra-chr",
        "--max-alleles", "2",
        "--out", out_prefix
    ]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        abort(f"VCF/BCF conversion failed. Check log: {log_file}")

    return out_prefix


def convert_bgen(config, converted_dir, log_dir):
    """BGEN to BED/BIM/FAM via PLINK 2."""
    plink2     = config.get("PLINK2_PATH", "plink2")
    path       = config["INPUT_PATH"]
    out_prefix = os.path.join(converted_dir, config["PROJECT_NAME"] + "_converted")
    log_file   = os.path.join(log_dir, "step2_bgen_conversion.log")

    bgen_file = path if os.path.isfile(path) else path + ".bgen"
    if not os.path.isfile(bgen_file):
        abort(f"No BGEN file found at: {path}")

    # need companion .sample file
    sample_file = os.path.splitext(bgen_file)[0] + ".sample"
    if not os.path.isfile(sample_file):
        abort(f"BGEN .sample file not found: {sample_file}")

    cmd = [
        plink2,
        "--bgen", bgen_file, "ref-first",
        "--sample", sample_file,
        "--make-bed",
        "--allow-extra-chr",
        "--max-alleles", "2",
        "--out", out_prefix
    ]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        abort(f"BGEN conversion failed. Check log: {log_file}")

    return out_prefix


def passthrough_plink(config, converted_dir, log_dir):
    """
    BED input, create a working copy so all downstream steps
    always reference the same output directory regardless of format.
    """
    plink1     = config.get("PLINK1_PATH", "plink")
    prefix     = config["INPUT_PATH"]
    out_prefix = os.path.join(converted_dir, config["PROJECT_NAME"] + "_converted")
    log_file   = os.path.join(log_dir, "step2_plink_passthrough.log")

    cmd = [
        plink1,
        "--bfile", prefix,
        "--make-bed",
        "--allow-extra-chr",
        "--out", out_prefix
    ]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        abort(f"PLINK pass-through failed. Check log: {log_file}")

    return out_prefix


# ─────────────────────────────────────────────
# CHROMOSOME CLEANING
# ─────────────────────────────────────────────
def clean_chromosomes(config, converted_prefix, converted_dir, log_dir):
    """
    Drop any non-standard chromosome codes.
    Returns (clean_prefix, variants_after, excluded_count, non_standard_dict).
    If nothing to clean, returns converted_prefix unchanged.
    """
    plink1 = config.get("PLINK1_PATH", "plink")
    bim    = converted_prefix + ".bim"

    log("Checking chromosome composition...")
    total, chr_counts, non_standard = count_bim(bim)

    log(f"  Total variants : {total:,}")

    if non_standard:
        for chrom, count in sorted(non_standard.items()):
            warn(f"  Non-standard chr '{chrom}': {count:,} variants will be excluded")
    else:
        log("  All chromosomes are standard. Skipping chr cleaning.")
        return converted_prefix, total, 0, {}

    # build extract list and re-run PLINK
    extract_list = os.path.join(converted_dir, "chr_extract_list.txt")
    retained, excluded = write_chr_list(bim, extract_list)

    log(f"  Variants retained : {retained:,}")
    log(f"  Variants excluded : {excluded:,}")

    clean_prefix = converted_prefix + "_chrclean"
    log_file     = os.path.join(log_dir, "step2_chr_cleaning.log")

    cmd = [
        plink1,
        "--bfile", converted_prefix,
        "--extract", extract_list,
        "--make-bed",
        "--out", clean_prefix
    ]

    rc = run_cmd(cmd, log_file)
    if rc != 0:
        abort(f"Chromosome cleaning failed. Check log: {log_file}")

    return clean_prefix, retained, excluded, non_standard


# ─────────────────────────────────────────────
# CLEANUP INTERMEDIATE FILES
# ─────────────────────────────────────────────
def cleanup_intermediate(converted_prefix, keep):
    if keep.upper() == "YES":
        log("KEEP_INTERMEDIATE = YES, intermediate files retained.")
        return
    log("KEEP_INTERMEDIATE = NO, removing intermediate files...")
    for ext in [".bed", ".bim", ".fam", ".log", ".nosex"]:
        f = converted_prefix + ext
        if os.path.isfile(f):
            os.remove(f)
            log(f"  Removed: {f}")


# ─────────────────────────────────────────────
# WRITE STEP 2 REPORT
# ─────────────────────────────────────────────
def write_report(config, fmt, converted_prefix, clean_prefix,
                 variants_before, variants_after, excluded_chr,
                 non_standard, output_dir):

    project     = config["PROJECT_NAME"]
    run_time    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_path = os.path.join(output_dir, f"{project}_step2_conversion_report.txt")
    json_path   = os.path.join(output_dir, f"{project}_step2_conversion_report.json")
    sep         = "=" * 65
    dash        = "-" * 65

    # output file MD5s
    file_records = []
    for ext in [".bed", ".bim", ".fam"]:
        fpath = clean_prefix + ext
        if os.path.isfile(fpath):
            size_mb  = round(os.path.getsize(fpath) / (1024 * 1024), 3)
            checksum = md5_file(fpath)
            file_records.append({
                "filename" : os.path.basename(fpath),
                "path"     : os.path.abspath(fpath),
                "size_mb"  : size_mb,
                "md5"      : checksum
            })

    with open(report_path, "w") as out:
        out.write(f"{sep}\n")
        out.write(f"  CRISP: STEP 2 CONVERSION REPORT\n")
        out.write(f"  Comprehensive Robust Integrated SNP Processing\n")
        out.write(f"{sep}\n")
        out.write(f"  Project      : {project}\n")
        out.write(f"  Date         : {run_time}\n")
        out.write(f"  Input format : {fmt}\n")
        out.write(f"  Input path   : {config['INPUT_PATH']}\n")
        out.write(f"{dash}\n")
        out.write(f"  CONVERSION\n")
        tool = "PLINK 1.9" if fmt in ("PED", "PLINK") else "PLINK 2"
        out.write(f"  Tool          : {tool}\n")
        out.write(f"  Output prefix : {os.path.basename(clean_prefix)}\n")
        out.write(f"{dash}\n")
        out.write(f"  VARIANT COUNTS\n")
        out.write(f"  Variants before chr cleaning : {variants_before:,}\n")
        out.write(f"  Variants after chr cleaning  : {variants_after:,}\n")
        out.write(f"  Variants excluded            : {excluded_chr:,}\n")
        if non_standard:
            out.write(f"{dash}\n")
            out.write(f"  NON-STANDARD CHROMOSOMES EXCLUDED\n")
            for chrom, count in sorted(non_standard.items()):
                out.write(f"  chr '{chrom}' : {count:,} variants excluded\n")
        out.write(f"{dash}\n")
        out.write(f"  OUTPUT FILES\n")
        for rec in file_records:
            out.write(f"  {rec['filename']}\n")
            out.write(f"  Size : {rec['size_mb']} MB  |  MD5 : {rec['md5']}\n\n")
        out.write(f"{dash}\n")
        out.write(f"  KEEP_INTERMEDIATE : {config.get('KEEP_INTERMEDIATE', 'YES')}\n")
        out.write(f"{sep}\n")
        out.write(f"  END OF REPORT\n")
        out.write(f"{sep}\n")

    # JSON export for downstream steps
    json_data = {
        "project"           : project,
        "run_datetime"      : run_time,
        "input_format"      : fmt,
        "input_path"        : config["INPUT_PATH"],
        "converted_prefix"  : os.path.abspath(clean_prefix),
        "variants_before"   : variants_before,
        "variants_after"    : variants_after,
        "variants_excluded" : excluded_chr,
        "non_standard_chr"  : non_standard,
        "output_files"      : file_records,
        "keep_intermediate" : config.get("KEEP_INTERMEDIATE", "YES")
    }

    with open(json_path, "w") as jf:
        json.dump(json_data, jf, indent=2)

    return report_path, json_path


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRISP Step 2: Format conversion to BED/BIM/FAM"
    )
    parser.add_argument("--config", required=True,
                        help="Path to parsed config JSON from Step 1")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        abort(f"Config JSON not found: {args.config}")

    with open(args.config) as f:
        config = json.load(f)

    fmt        = config["INPUT_FORMAT"].upper()
    output_dir = config["OUTPUT_DIR"]
    project    = config["PROJECT_NAME"]

    log(f"Project      : {project}")
    log(f"Input format : {fmt}")
    log(f"Input path   : {config['INPUT_PATH']}")

    converted_dir = os.path.join(output_dir, "step2_converted")
    log_dir       = os.path.join(output_dir, "logs")
    os.makedirs(converted_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # conversion
    log("Starting format conversion...")

    if fmt == "PLINK":
        converted_prefix = passthrough_plink(config, converted_dir, log_dir)
    elif fmt == "PED":
        converted_prefix = convert_ped(config, converted_dir, log_dir)
    elif fmt in ("VCF", "BCF"):
        converted_prefix = convert_vcf(config, converted_dir, log_dir)
    elif fmt == "BGEN":
        converted_prefix = convert_bgen(config, converted_dir, log_dir)
    else:
        abort(f"Unrecognised format: {fmt}")

    log("Conversion complete.")

    # chromosome cleaning
    log("Starting chromosome cleaning...")
    clean_prefix, variants_after, excluded_chr, non_standard = clean_chromosomes(
        config, converted_prefix, converted_dir, log_dir
    )

    # count before-cleaning variants from the pre-clean BIM
    bim_before      = converted_prefix + ".bim"
    variants_before = sum(1 for line in open(bim_before) if line.strip())

    log("Chromosome cleaning complete.")

    # remove intermediate files if requested
    if clean_prefix != converted_prefix:
        cleanup_intermediate(
            converted_prefix,
            config.get("KEEP_INTERMEDIATE", "YES")
        )

    # write report and JSON
    log("Writing Step 2 report...")
    report_path, json_path = write_report(
        config, fmt, converted_prefix, clean_prefix,
        variants_before, variants_after, excluded_chr,
        non_standard, output_dir
    )

    log(f"  Report : {report_path}")
    log(f"  JSON   : {json_path}")

    # update config JSON so downstream steps know the converted prefix
    config["CONVERTED_PREFIX"] = os.path.abspath(clean_prefix)
    config_out = os.path.join(output_dir, f"{project}_parsed_config.json")
    with open(config_out, "w") as cf:
        json.dump(config, cf, indent=2)

    log(f"  Updated config : {config_out}")
    log("Step 2 complete.")


if __name__ == "__main__":
    main()
