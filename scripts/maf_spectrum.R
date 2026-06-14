#!/usr/bin/env Rscript
# ##############################################################
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
# ##############################################################
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### MAF Spectrum Binning (Step 4, v0.3.0)
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Computes a MAF spectrum from one or two PLINK .frq files and writes:
#   <PROJECT_NAME>_step4_maf_spectrum.txt  — text table(s)
#   <PROJECT_NAME>_step4_maf_spectrum.pdf  — grouped bar chart
#
# Modes:
#   COMMON  : 6-bin MAF spectrum using default bin edges
#             0, 0.0001, 0.01, 0.05, 0.10, 1
#             Bin 1 : MAF == 0             (actual monomorphic)
#             Bin 2 : 0 < MAF <= 0.0001   (considered monomorphic)
#             Bin 3 : 0.0001 < MAF <= 0.01
#             Bin 4 : 0.01 < MAF <= 0.05
#             Bin 5 : 0.05 < MAF <= 0.10
#             Bin 6 : MAF > 0.10
#
#   CUSTOM  : same 6-bin structure, user-supplied edges via
#             MAF_SPECTRUM_BINS in crisp_instructions.txt
#
# Input frequency file format (PLINK --freq output, .frq):
#   columns: CHR SNP A1 A2 MAF NCHROBS
#
# Colour palette is read from CRISP_PAL_* environment variables exported
# by _crisp_flavour.sh's _crisp_palette() function.
# CRISP_PAL_MODE: STANDARD | COLOURBLIND | NIGHT
#
# Note: this script is self-contained and does not source plot_snprate.R.
# The theme below mirrors crisp_theme in plot_snprate.R, with right-side
# legend per the CRISP plotting standard.
#
# Usage:
#   Rscript maf_spectrum.R \
#       <MODE> \
#       <MAF_SPECTRUM_BINS> \
#       <SAMPLES> \
#       <OUT_DIR> \
#       <PROJECT_NAME> \
#       diagnostic:<freq_file> \
#       [final:<freq_file>]
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
})

##########################################################################
### ARGUMENT PARSING
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    cat("[CRISP-MAFSPEC] ERROR: Insufficient arguments.\n")
    cat("[CRISP-MAFSPEC] Usage: Rscript maf_spectrum.R <MODE> <BINS> <SAMPLES> <OUT_DIR> <PROJECT_NAME> diagnostic:<file> [final:<file>]\n")
    quit(status = 1)
}

MODE         <- toupper(args[1])
BINS_STR     <- args[2]
SAMPLES      <- as.integer(args[3])
OUT_DIR      <- args[4]
PROJECT_NAME <- args[5]

checkpoint_args <- args[6:length(args)]

if (!(MODE %in% c("COMMON", "CUSTOM"))) {
    cat(sprintf("[CRISP-MAFSPEC] ERROR: Unknown MAF_SPECTRUM_MODE '%s'. Valid: COMMON, CUSTOM\n", MODE))
    quit(status = 1)
}

if (is.na(SAMPLES) || SAMPLES <= 0) {
    cat("[CRISP-MAFSPEC] ERROR: SAMPLES must be a positive integer.\n")
    quit(status = 1)
}

checkpoints <- list()
for (ca in checkpoint_args) {
    parts <- strsplit(ca, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
        cat(sprintf("[CRISP-MAFSPEC] ERROR: Malformed checkpoint argument '%s' (expected label:path)\n", ca))
        quit(status = 1)
    }
    label <- parts[1]
    path  <- paste(parts[-1], collapse = ":")
    checkpoints[[label]] <- path
}

if (!"diagnostic" %in% names(checkpoints)) {
    cat("[CRISP-MAFSPEC] ERROR: 'diagnostic' checkpoint is required.\n")
    quit(status = 1)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

##########################################################################
### COLOUR PALETTE
### STANDARD    : default CRISP palette
### COLOURBLIND : Okabe-Ito (Nature Methods recommended)
### NIGHT       : dark-background, Okabe-Ito-derived
###              (provisional pending final NIGHT palette in _crisp_flavour.sh)
##########################################################################

colour_mode <- toupper(Sys.getenv("CRISP_PAL_MODE", unset = "STANDARD"))

if (colour_mode == "COLOURBLIND") {
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",  "#56B4E9")
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS", "#0072B2")
    BG_COLOUR     <- "#FFFFFF"
    FG_COLOUR     <- "#000000"
    GRID_COLOUR   <- "#E0E0E0"
    SUBTITLE_COL  <- "#555555"
    CAPTION_COL   <- "#888888"
} else if (colour_mode == "NIGHT") {
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",  "#56B4E9")
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS", "#009E73")
    BG_COLOUR     <- "#1E1E1E"
    FG_COLOUR     <- "#E0E0E0"
    GRID_COLOUR   <- "#3A3A3A"
    SUBTITLE_COL  <- "#AAAAAA"
    CAPTION_COL   <- "#777777"
} else {
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",  "#7ec8e3")
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS", "#1D9E75")
    BG_COLOUR     <- "#FFFFFF"
    FG_COLOUR     <- "#000000"
    GRID_COLOUR   <- "#E0E0E0"
    SUBTITLE_COL  <- "#555555"
    CAPTION_COL   <- "#888888"
}

COL_DIAGNOSTIC <- PAL_SKY
COL_FINAL      <- PAL_PASS

##########################################################################
### THEME
### Mirrors crisp_theme in plot_snprate.R.
### Right-side legend per CRISP plotting standard.
##########################################################################

crisp_spectrum_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title        = element_text(face = "bold", size = 14, colour = FG_COLOUR),
        plot.subtitle     = element_text(colour = SUBTITLE_COL, size = 10),
        axis.title        = element_text(face = "bold", colour = FG_COLOUR),
        axis.text         = element_text(colour = FG_COLOUR),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = GRID_COLOUR),
        plot.caption      = element_text(colour = CAPTION_COL, size = 8),
        legend.position   = "right",
        legend.title      = element_text(face = "bold", colour = FG_COLOUR),
        legend.text       = element_text(colour = FG_COLOUR),
        plot.background   = element_rect(fill = BG_COLOUR, colour = NA),
        panel.background  = element_rect(fill = BG_COLOUR, colour = NA),
        legend.background = element_rect(fill = BG_COLOUR, colour = NA),
        legend.key        = element_rect(fill = BG_COLOUR, colour = NA)
    )

##########################################################################
### HELPERS
##########################################################################

read_frq <- function(path) {
    if (!file.exists(path)) {
        cat(sprintf("[CRISP-MAFSPEC] ERROR: Frequency file not found: %s\n", path))
        quit(status = 1)
    }
    df <- read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE)
    if (!"MAF" %in% colnames(df)) {
        cat(sprintf("[CRISP-MAFSPEC] ERROR: '%s' missing MAF column. Is it a PLINK --freq .frq file?\n", path))
        quit(status = 1)
    }
    df
}

parse_edges <- function(s) {
    vals <- as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
    if (any(is.na(vals))) {
        cat(sprintf("[CRISP-MAFSPEC] ERROR: Could not parse MAF_SPECTRUM_BINS '%s'\n", s))
        quit(status = 1)
    }
    if (length(vals) != 6) {
        cat(sprintf("[CRISP-MAFSPEC] ERROR: MAF_SPECTRUM_BINS must have exactly 6 values (got %d)\n", length(vals)))
        quit(status = 1)
    }
    if (vals[1] != 0) {
        cat("[CRISP-MAFSPEC] ERROR: First MAF_SPECTRUM_BINS edge must be 0.\n")
        quit(status = 1)
    }
    if (any(diff(vals) <= 0)) {
        cat("[CRISP-MAFSPEC] ERROR: MAF_SPECTRUM_BINS edges must be strictly increasing.\n")
        quit(status = 1)
    }
    vals
}

##########################################################################
### BIN STRUCTURE
### edges = [e1, e2, e3, e4, e5, e6] where e1=0, e6=1 (typical)
###
###   Bin 1 — Actual monomorphic      : MAF == 0
###   Bin 2 — Considered monomorphic  : 0       < MAF <= e2
###   Bin 3 — Low frequency             : e2      < MAF <= e3
###   Bin 4 — Low frequency           : e3      < MAF <= e4
###   Bin 5 — Common                  : e4      < MAF <= e5
###   Bin 6 — Very common             : MAF     > e5
###
### Bins are mutually exclusive and exhaustive over [0, 1].
##########################################################################

build_bins <- function(edges) {
    list(
        edges  = edges,
        labels = c(
            "Actual monomorphic",
            sprintf("Considered mono (0 < MAF <= %.4g%%)", edges[2] * 100),
            sprintf("%.4g%% < MAF <= %.4g%%",  edges[2] * 100, edges[3] * 100),
            sprintf("%.4g%% < MAF <= %.4g%%",  edges[3] * 100, edges[4] * 100),
            sprintf("%.4g%% < MAF <= %.4g%%",  edges[4] * 100, edges[5] * 100),
            sprintf("MAF > %.4g%%",            edges[5] * 100)
        )
    )
}

count_bins <- function(maf, edges) {
    c(
        sum(maf == 0),
        sum(maf >  0        & maf <= edges[2]),
        sum(maf >  edges[2] & maf <= edges[3]),
        sum(maf >  edges[3] & maf <= edges[4]),
        sum(maf >  edges[4] & maf <= edges[5]),
        sum(maf >  edges[5])
    )
}

##########################################################################
### COMPUTE SPECTRA
##########################################################################

edges    <- parse_edges(BINS_STR)
bin_info <- build_bins(edges)
results  <- list()

for (label in names(checkpoints)) {
    frq    <- read_frq(checkpoints[[label]])
    counts <- count_bins(frq$MAF, edges)
    total  <- nrow(frq)

    if (sum(counts) != total) {
        cat(sprintf("[CRISP-MAFSPEC] WARNING: Bin sum (%d) != total variants (%d) for checkpoint '%s'.\n",
                    sum(counts), total, label))
    }

    results[[label]] <- list(counts = counts, total = total)
    cat(sprintf("[CRISP-MAFSPEC] Checkpoint '%s': %d variants processed.\n", label, total))
}

##########################################################################
### TEXT TABLE
##########################################################################

txt_path <- file.path(OUT_DIR, sprintf("%s_step4_maf_spectrum.txt", PROJECT_NAME))
con      <- file(txt_path, open = "wt")

writeLines("==================================================================", con)
writeLines("  CRISP: STEP 4 MAF SPECTRUM", con)
writeLines(sprintf("  Project   : %s", PROJECT_NAME), con)
writeLines(sprintf("  Mode      : %s", MODE), con)
writeLines(sprintf("  Samples   : %d", SAMPLES), con)
writeLines(sprintf("  Bin edges : %s", paste(edges, collapse = ", ")), con)
writeLines("==================================================================", con)
writeLines("", con)

for (label in names(results)) {
    r <- results[[label]]
    writeLines(sprintf("--- %s (n = %d variants) ---", toupper(label), r$total), con)
    writeLines(sprintf("  %-38s %10s  %8s", "Bin", "Count", "%"), con)
    writeLines(sprintf("  %s", strrep("-", 60)), con)
    for (i in seq_along(bin_info$labels)) {
        pct <- if (r$total > 0) 100 * r$counts[i] / r$total else 0
        writeLines(sprintf("  %-38s %10d  %7.2f%%", bin_info$labels[i], r$counts[i], pct), con)
    }
    writeLines("", con)
}

writeLines("==================================================================", con)
writeLines("  END OF MAF SPECTRUM", con)
writeLines("==================================================================", con)
close(con)

cat(sprintf("[CRISP-MAFSPEC] Text table: %s\n", txt_path))

##########################################################################
### GROUPED BAR CHART
##########################################################################

pdf_path <- file.path(OUT_DIR, sprintf("%s_step4_maf_spectrum.pdf", PROJECT_NAME))

plot_df <- do.call(rbind, lapply(names(results), function(label) {
    data.frame(
        bin        = factor(bin_info$labels, levels = bin_info$labels),
        checkpoint = toupper(label),
        count      = results[[label]]$counts,
        stringsAsFactors = FALSE
    )
}))

chk_levels  <- toupper(names(results))
plot_df$checkpoint <- factor(plot_df$checkpoint, levels = chk_levels)

chk_colours <- setNames(
    c(COL_DIAGNOSTIC, COL_FINAL)[seq_along(chk_levels)],
    chk_levels
)

subtitle_text <- if (length(results) > 1) {
    "Diagnostic (post call-rate filter) vs Final (post mono/MAF removal)"
} else {
    "Diagnostic checkpoint (post call-rate filter)"
}

p <- ggplot(plot_df, aes(x = bin, y = count, fill = checkpoint)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = NA) +
    scale_fill_manual(values = chk_colours, name = "Checkpoint") +
    scale_y_continuous(labels = comma) +
    labs(
        title    = sprintf("CRISP — Step 4: MAF spectrum (%s mode)", MODE),
        subtitle = subtitle_text,
        x        = "Minor allele frequency bin",
        y        = "Variant count",
        caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
    ) +
    crisp_spectrum_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(pdf_path, width = 10, height = 6, bg = BG_COLOUR)
print(p)
dev.off()

cat(sprintf("[CRISP-MAFSPEC] Plot:       %s\n", pdf_path))
