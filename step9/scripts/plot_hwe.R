#!/usr/bin/env Rscript
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### Step 8: HWE Filtering, Plotting Script
### Version: 0.1.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Generates HWE p-value distribution plots from PLINK --hardy output.
#
# Plots:
#   Plot 1: Histogram of -log10(HWE p-values) with threshold line
#   Plot 2: QQ plot of observed vs expected -log10(p-values)
#
# Usage:
#   Rscript plot_hwe.R <output_dir> <hwe_file>
#                      <threshold> <mode>
##########################################################################


##########################################################################
### COLOUR PALETTE
### Reads PLOT_COLOUR_MODE from environment (set by crisp_hwe.sh)
### STANDARD: default CRISP palette
### COLOURBLIND: Okabe-Ito palette (Nature Methods recommended)
##########################################################################

colour_mode <- Sys.getenv("CRISP_PAL_MODE", unset = "STANDARD")

if (toupper(colour_mode) == "COLOURBLIND") {
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#0072B2")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#D55E00")
    PAL_WARN      <- Sys.getenv("CRISP_PAL_WARN",      "#E69F00")
    PAL_HIGHLIGHT <- Sys.getenv("CRISP_PAL_HIGHLIGHT",  "#009E73")
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",        "#56B4E9")
    PAL_PINK      <- Sys.getenv("CRISP_PAL_PINK",       "#CC79A7")
    PAL_ORANGE    <- Sys.getenv("CRISP_PAL_YELLOW",     "#F0E442")
} else {
    PAL_PASS      <- "#1D9E75"
    PAL_FAIL      <- "#ff5f57"
    PAL_WARN      <- "#febc2e"
    PAL_HIGHLIGHT <- "#5dcaa5"
    PAL_SKY       <- "#7ec8e3"
    PAL_PINK      <- "#afa9ec"
    PAL_ORANGE    <- "#ff8c00"
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    cat("Usage: Rscript plot_hwe.R <output_dir> <hwe_file> <threshold> <mode> [<hwe_meta>] [<meta_file>]\n")
    quit(status = 1)
}

output_dir  <- args[1]
hwe_file    <- args[2]
threshold   <- as.numeric(args[3])
mode        <- toupper(args[4])
hwe_meta    <- if (length(args) >= 5) toupper(args[5]) else "NO"
meta_file   <- if (length(args) >= 6) args[6] else ""

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[PLOT] CRISP Step 8: HWE plotting\n"))
cat(sprintf("[PLOT] File      : %s\n", hwe_file))
cat(sprintf("[PLOT] Threshold : %g\n", threshold))
cat(sprintf("[PLOT] Mode      : %s\n", mode))

if (!file.exists(hwe_file)) {
    cat(sprintf("[PLOT] ERROR: HWE file not found: %s\n", hwe_file))
    quit(status = 1)
}

df <- fread(hwe_file, header = TRUE)

if (nrow(df) == 0) {
    cat("[PLOT] ERROR: HWE file is empty.\n")
    quit(status = 1)
}

cat(sprintf("[PLOT] Variants loaded: %d\n", nrow(df)))

# use ALL genotype test rows
df <- df[df$TEST == "ALL",]
df <- df[!is.na(df$P) & df$P > 0,]

cat(sprintf("[PLOT] Variants with valid P-values: %d\n", nrow(df)))

df[, log10P := -log10(P)]
threshold_log <- -log10(threshold)
n_fail <- sum(df$P < threshold)

cat(sprintf("[PLOT] Variants failing HWE (p < %g): %d\n", threshold, n_fail))

CAPTION <- "CRISP | Comprehensive Robust Integrated SNP Processing | Compass Genomics"

crisp_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "#555555", size = 10),
        axis.title       = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#f8f9fa", colour = NA),
        plot.caption     = element_text(colour = "#888888", size = 8),
        legend.position  = "right"
    )

##########################################################################
### PLOT 1: -LOG10(P) HISTOGRAM
##########################################################################

p1 <- ggplot(df, aes(x = log10P)) +
    geom_histogram(
        aes(fill = log10P >= threshold_log),
        bins       = 80,
        colour     = "white",
        linewidth  = 0.15
    ) +
    scale_fill_manual(
        values = c("FALSE" = "#1D9E75", "TRUE" = "#ff5f57"),
        labels = c("FALSE" = "Passing HWE",
                   "TRUE"  = sprintf("Failing HWE (p < %g)", threshold)),
        name   = ""
    ) +
    geom_vline(
        xintercept = threshold_log,
        colour     = "#ff5f57",
        linetype   = "dashed",
        linewidth  = 1.0
    ) +
    annotate(
        "text",
        x     = threshold_log + 0.15,
        y     = Inf,
        label = sprintf("-log10(p) = %.1f", threshold_log),
        hjust = 0, vjust = 1.5,
        colour    = "#ff5f57",
        size      = 3.2,
        fontface  = "bold"
    ) +
    scale_y_continuous(labels = comma) +
    labs(
        title    = "CRISP: Step 8 HWE, P-value Distribution",
        subtitle = sprintf(
            "Mode: %s  |  Threshold: p < %g  |  %s failing HWE",
            mode,
            threshold,
            format(n_fail, big.mark = ",")
        ),
        x       = expression(-log[10](italic(p))),
        y       = "Number of variants",
        caption = CAPTION
    ) +
    crisp_theme

out1 <- file.path(output_dir, "HWE_distribution.pdf")
pdf(out1, width = 10, height = 6)
print(p1)
dev.off()
cat(sprintf("[PLOT] Histogram saved: %s\n", out1))

##########################################################################
### PLOT 2: QQ PLOT
##########################################################################

n     <- nrow(df)
obs   <- sort(df$log10P, decreasing = TRUE)
exp   <- -log10((seq_len(n)) / n)
df_qq <- data.frame(expected = exp, observed = obs)

# sample for large datasets to keep plot manageable
if (n > 50000) {
    set.seed(42)
    idx    <- sort(sample(n, 50000))
    df_qq  <- df_qq[idx,]
}

lambda <- median(df$log10P) / (-log10(0.5))

p2 <- ggplot(df_qq, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0,
                colour = "#888", linetype = "dashed", linewidth = 0.8) +
    geom_point(
        colour = "#1D9E75", size = 1.0, alpha = 0.5
    ) +
    geom_hline(
        yintercept = threshold_log,
        colour     = "#ff5f57",
        linetype   = "dashed",
        linewidth  = 0.8
    ) +
    annotate(
        "text",
        x     = 0, y = threshold_log + 0.15,
        label = sprintf("p = %g threshold", threshold),
        hjust = 0, vjust = 0,
        colour   = "#ff5f57",
        size     = 3.0,
        fontface = "bold"
    ) +
    labs(
        title    = "CRISP: Step 8 HWE, QQ Plot",
        subtitle = sprintf(
            "Mode: %s  |  lambda = %.4f  |  %d variants",
            mode, lambda, nrow(df)
        ),
        x       = expression("Expected " ~ -log[10](italic(p))),
        y       = expression("Observed " ~ -log[10](italic(p))),
        caption = CAPTION
    ) +
    crisp_theme

out2 <- file.path(output_dir, "HWE_qq.pdf")
pdf(out2, width = 8, height = 8)
print(p2)
dev.off()
cat(sprintf("[PLOT] QQ plot saved: %s\n", out2))

cat("\n[PLOT] HWE plotting complete.\n")

##########################################################################
### PLOT 3: META-ANALYSIS CLASSIFICATION (OPTIONAL)
##########################################################################

if (hwe_meta == "YES" && nchar(meta_file) > 0 && file.exists(meta_file)) {

    cat(sprintf("[PLOT] Loading meta-analysis results: %s\n", meta_file))
    df_meta <- fread(meta_file, header = TRUE)

    if (nrow(df_meta) > 0 && "CLASS" %in% colnames(df_meta) &&
        "P_META" %in% colnames(df_meta)) {

        df_meta <- df_meta[!is.na(df_meta$P_META) & df_meta$P_META > 0,]
        df_meta[, log10P_meta := -log10(P_META)]

        CLASS_COLORS <- c(
            "PASS"              = "#1D9E75",
            "FAIL_ALL"          = "#ff5f57",
            "FAIL_CASES_ONLY"   = "#febc2e",
            "FAIL_CONTROLS_ONLY"= "#ff8c00"
        )

        p3 <- ggplot(df_meta, aes(x = log10P_meta, fill = CLASS)) +
            geom_histogram(bins = 80, colour = "white", linewidth = 0.15) +
            scale_fill_manual(
                values = CLASS_COLORS,
                labels = c(
                    "PASS"               = "Pass",
                    "FAIL_ALL"           = "Fail all strata (likely error)",
                    "FAIL_CASES_ONLY"    = "Fail cases only (potential signal)",
                    "FAIL_CONTROLS_ONLY" = "Fail controls only"
                ),
                name = "Classification"
            ) +
            geom_vline(xintercept = -log10(threshold),
                       colour = "#ff5f57", linetype = "dashed",
                       linewidth = 1.0) +
            scale_y_continuous(labels = comma) +
            labs(
                title    = "CRISP: Step 8 HWE, Meta-analysis Classification",
                subtitle = sprintf(
                    "Fisher's combined p-value across strata  |  threshold p < %g",
                    threshold),
                x       = expression(-log[10](italic(p[meta]))),
                y       = "Number of variants",
                caption = CAPTION
            ) +
            crisp_theme

        out3 <- file.path(output_dir, "HWE_meta_classification.pdf")
        pdf(out3, width = 11, height = 6)
        print(p3)
        dev.off()
        cat(sprintf("[PLOT] Meta-analysis classification plot saved: %s\n", out3))
    }
}
