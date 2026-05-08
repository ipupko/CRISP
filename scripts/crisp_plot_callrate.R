#!/usr/bin/env Rscript
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 3: Sample Call Rate, Plotting Script
# Version: 0.2.0
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
##########################################################################
# DESCRIPTION
# Generates missingness distribution plots from PLINK .imiss files.
#
# SIMPLE mode:
#   Single histogram with MIND threshold line, failing samples in red.
#
# CASCADE mode:
#   Four individual PDFs plus one faceted 2x2 comparison plot.
#
# CUSTOM mode:
#   One individual PDF per tier plus one faceted plot.
#   Rows scale automatically to tier count.
#   Thresholds parsed from .imiss filenames (mind{value} pattern).
#
# Usage:
#   Rscript plot_callrate.R <mode> <mind> <report_dir> <imiss_file(s)>
#   mode       : SIMPLE, CASCADE or CUSTOM
#   mind       : final MIND threshold (e.g. 0.05)
#   report_dir : output directory for plots
#   imiss_files: one .imiss file (SIMPLE) or one per tier
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    cat("Usage: Rscript plot_callrate.R <mode> <mind> <report_dir> <imiss_file(s)>\n")
    quit(status = 1)
}

mode        <- toupper(args[1])
mind        <- as.numeric(args[2])
report_dir  <- args[3]
imiss_files <- args[4:length(args)]
n_passes    <- length(imiss_files)

dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# shared ggplot2 theme
crisp_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(colour = "#555555", size = 10),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(colour = "#888888", size = 8),
        legend.position  = "bottom",
        strip.text       = element_text(face = "bold", size = 10)
    )

# load a .imiss file, check it has F_MISS
load_imiss <- function(filepath) {
    df <- fread(filepath, header = TRUE)
    if (!"F_MISS" %in% colnames(df)) {
        stop(paste("F_MISS column not found in:", filepath))
    }
    return(df)
}

# extract threshold from filename, e.g. mind0.20.imiss -> 0.20
parse_threshold <- function(filepath, fallback) {
    m <- regmatches(basename(filepath),
                    regexpr("mind([0-9]+\\.?[0-9]*)", basename(filepath)))
    if (length(m) > 0) as.numeric(sub("mind", "", m)) else fallback
}

##########################################################################
# SIMPLE MODE
##########################################################################

plot_simple <- function(imiss_file, mind, report_dir) {

    cat("[PLOT] Loading:", imiss_file, "\n")
    df <- load_imiss(imiss_file)

    n_total   <- nrow(df)
    n_failing <- sum(df$F_MISS > mind)

    cat(sprintf("[PLOT] Samples total   : %d\n", n_total))
    cat(sprintf("[PLOT] Samples failing : %d (F_MISS > %.3f)\n", n_failing, mind))

    df$status <- ifelse(df$F_MISS > mind, "Failing", "Passing")

    p <- ggplot(df, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 50, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = mind, colour = "#ff5f57",
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = mind, y = Inf,
                 label = paste0("MIND = ", mind),
                 hjust = -0.1, vjust = 1.5,
                 colour = "#ff5f57", size = 3.5, fontface = "bold") +
        scale_fill_manual(
            values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
            name   = "Sample status"
        ) +
        scale_x_continuous(
            labels = percent_format(accuracy = 0.1),
            limits = c(0, max(df$F_MISS) * 1.05)
        ) +
        scale_y_continuous(labels = comma) +
        labs(
            title   = "CRISP: Step 3 Sample Call Rate",
            subtitle = sprintf(
                "SIMPLE mode  |  MIND: %.3f  |  %d/%d samples failing",
                mind, n_failing, n_total),
            x       = "Per-sample missingness rate (F_MISS)",
            y       = "Number of samples",
            caption = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    out_file <- file.path(report_dir, "step3_callrate_simple.pdf")
    pdf(out_file, width = 10, height = 6)
    print(p)
    dev.off()

    cat(sprintf("[PLOT] Simple histogram saved: %s\n", out_file))
    return(out_file)
}

##########################################################################
# MULTI-PASS MODE (CASCADE and CUSTOM)
# Handles any number of passes, rows in faceted plot scale to count
##########################################################################

plot_multipass <- function(imiss_files, thresholds, mode, report_dir) {

    n           <- length(imiss_files)
    pass_labels <- sprintf("Pass %d (mind=%.2f)", seq_along(thresholds), thresholds)
    all_data    <- list()
    pdf_files   <- c()

    # individual plot per pass
    for (i in seq_along(imiss_files)) {

        cat(sprintf("\n[PLOT] Pass %d, loading: %s\n", i, imiss_files[i]))
        df        <- load_imiss(imiss_files[i])
        threshold <- thresholds[i]

        n_total   <- nrow(df)
        n_failing <- sum(df$F_MISS > threshold)

        cat(sprintf("[PLOT]   Threshold : %.3f\n", threshold))
        cat(sprintf("[PLOT]   Samples   : %d total, %d failing\n",
                    n_total, n_failing))

        df$status    <- ifelse(df$F_MISS > threshold, "Failing", "Passing")
        df$pass      <- pass_labels[i]
        df$threshold <- threshold

        all_data[[i]] <- df

        p <- ggplot(df, aes(x = F_MISS, fill = status)) +
            geom_histogram(bins = 50, colour = "white", linewidth = 0.2) +
            geom_vline(xintercept = threshold, colour = "#ff5f57",
                       linetype = "dashed", linewidth = 0.8) +
            annotate("text", x = threshold, y = Inf,
                     label = paste0("MIND = ", threshold),
                     hjust = -0.1, vjust = 1.5,
                     colour = "#ff5f57", size = 3.5, fontface = "bold") +
            scale_fill_manual(
                values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
                name   = "Sample status"
            ) +
            scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
            scale_y_continuous(labels = comma) +
            labs(
                title   = sprintf("CRISP: Step 3 Sample Call Rate, %s", pass_labels[i]),
                subtitle = sprintf(
                    "%s mode  |  %d/%d samples failing at this threshold",
                    mode, n_failing, n_total),
                x       = "Per-sample missingness rate (F_MISS)",
                y       = "Number of samples",
                caption = "CRISP | Comprehensive Robust Integrated SNP Processing"
            ) +
            crisp_theme

        out_file <- file.path(report_dir,
                              sprintf("step3_callrate_%s_pass%d.pdf",
                                      tolower(mode), i))
        pdf(out_file, width = 10, height = 6)
        print(p)
        dev.off()

        pdf_files <- c(pdf_files, out_file)
        cat(sprintf("[PLOT]   Individual plot saved: %s\n", out_file))
    }

    # faceted plot, 2 columns, rows scale to pass count
    cat(sprintf("\n[PLOT] Generating faceted %s comparison plot (%d passes)...\n",
                mode, n))

    combined      <- do.call(rbind, all_data)
    combined$pass <- factor(combined$pass, levels = pass_labels)

    thresholds_df <- data.frame(
        pass      = factor(pass_labels, levels = pass_labels),
        threshold = thresholds
    )

    n_cols <- 2
    n_rows <- ceiling(n / n_cols)
    fig_h  <- max(6, n_rows * 4.5)

    p_facet <- ggplot(combined, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 40, colour = "white", linewidth = 0.15) +
        geom_vline(data = thresholds_df,
                   aes(xintercept = threshold),
                   colour = "#ff5f57", linetype = "dashed",
                   linewidth = 0.7) +
        facet_wrap(~ pass, scales = "free", ncol = n_cols) +
        scale_fill_manual(
            values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
            name   = "Sample status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
        scale_y_continuous(labels = comma) +
        labs(
            title   = sprintf("CRISP: Step 3 Sample Call Rate, %s Mode", mode),
            subtitle = sprintf(
                "Missingness distribution across %d threshold passes", n),
            x       = "Per-sample missingness rate (F_MISS)",
            y       = "Number of samples",
            caption = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    facet_file <- file.path(report_dir,
                            sprintf("step3_callrate_%s_faceted.pdf",
                                    tolower(mode)))
    pdf(facet_file, width = 14, height = fig_h)
    print(p_facet)
    dev.off()

    cat(sprintf("[PLOT] Faceted plot saved: %s\n", facet_file))
    return(c(pdf_files, facet_file))
}

##########################################################################
# DISPATCH
##########################################################################

cat("\n[PLOT] CRISP Step 3, Call Rate Plotting\n")
cat(sprintf("[PLOT] Mode       : %s\n", mode))
cat(sprintf("[PLOT] MIND       : %.3f\n", mind))
cat(sprintf("[PLOT] Report dir : %s\n", report_dir))
cat(sprintf("[PLOT] Files      : %d\n\n", n_passes))

if (mode == "SIMPLE") {

    plot_simple(imiss_files[1], mind, report_dir)

} else if (mode == "CASCADE") {

    # fixed cascade thresholds
    thresholds <- c(0.25, 0.20, 0.10, 0.05)
    if (n_passes != 4) {
        cat(sprintf("[PLOT] ERROR: CASCADE expects 4 .imiss files, got %d\n", n_passes))
        quit(status = 1)
    }
    plot_multipass(imiss_files, thresholds, "CASCADE", report_dir)

} else if (mode == "CUSTOM") {

    # infer thresholds from filenames
    thresholds <- sapply(seq_along(imiss_files), function(i) {
        parse_threshold(imiss_files[i],
                        fallback = mind + (n_passes - i) * 0.05)
    })
    plot_multipass(imiss_files, thresholds, "CUSTOM", report_dir)

} else {

    cat(sprintf("[PLOT] ERROR: Unknown mode '%s'. Use SIMPLE, CASCADE or CUSTOM.\n", mode))
    quit(status = 1)
}

cat("\n[PLOT] Plotting complete.\n")### HELPER: LOAD IMISS FILE
##########################################################################

load_imiss <- function(filepath) {
    df <- fread(filepath, header = TRUE)
    # PLINK .imiss columns: FID IID MISS_PHENO N_MISS N_GENO F_MISS
    if (!"F_MISS" %in% colnames(df)) {
        stop(paste("F_MISS column not found in:", filepath))
    }
    return(df)
}

##########################################################################
### SIMPLE MODE -- SINGLE HISTOGRAM
##########################################################################

plot_simple <- function(imiss_file, mind, report_dir) {

    cat("[PLOT] Loading:", imiss_file, "\n")
    df <- load_imiss(imiss_file)

    n_total   <- nrow(df)
    n_failing <- sum(df$F_MISS > mind)
    n_passing <- n_total - n_failing

    cat(sprintf("[PLOT] Samples total   : %d\n", n_total))
    cat(sprintf("[PLOT] Samples failing : %d (F_MISS > %.3f)\n", n_failing, mind))
    cat(sprintf("[PLOT] Samples passing : %d\n", n_passing))

    df$status <- ifelse(df$F_MISS > mind, "Failing", "Passing")

    p <- ggplot(df, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 50, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = mind, colour = "#ff5f57",
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = mind, y = Inf,
                 label = paste0("MIND = ", mind),
                 hjust = -0.1, vjust = 1.5,
                 colour = "#ff5f57", size = 3.5, fontface = "bold") +
        scale_fill_manual(values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
                          name = "Sample status") +
        scale_x_continuous(labels = percent_format(accuracy = 0.1),
                           limits = c(0, max(df$F_MISS) * 1.05)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP -- Step 3: Sample Call Rate",
            subtitle = sprintf("SIMPLE mode  |  MIND threshold: %.3f  |  %d/%d samples failing",
                               mind, n_failing, n_total),
            x        = "Per-sample missingness rate (F_MISS)",
            y        = "Number of samples",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    out_file <- file.path(report_dir, "step3_callrate_simple.pdf")
    pdf(out_file, width = 10, height = 6)
    print(p)
    dev.off()

    cat(sprintf("[PLOT] Simple histogram saved: %s\n", out_file))
    return(out_file)
}

##########################################################################
### CASCADE MODE -- INDIVIDUAL PLOTS + FACETED PLOT
##########################################################################

plot_cascade <- function(imiss_files, mind, report_dir) {

    # Cascade thresholds match the fixed sequence in crisp_callrate.sh
    thresholds <- c(0.25, 0.20, 0.10, 0.05)
    pass_labels <- c("Pass 1 (mind=0.25)", "Pass 2 (mind=0.20)",
                     "Pass 3 (mind=0.10)", "Pass 4 (mind=0.05)")

    if (length(imiss_files) != 4) {
        stop(sprintf("CASCADE mode expects 4 .imiss files, got %d", length(imiss_files)))
    }

    all_data   <- list()
    pdf_files  <- c()

    for (i in seq_along(imiss_files)) {

        cat(sprintf("\n[PLOT] Pass %d -- Loading: %s\n", i, imiss_files[i]))
        df        <- load_imiss(imiss_files[i])
        threshold <- thresholds[i]

        n_total   <- nrow(df)
        n_failing <- sum(df$F_MISS > threshold)

        cat(sprintf("[PLOT]   Threshold : %.2f\n", threshold))
        cat(sprintf("[PLOT]   Samples   : %d total, %d failing\n", n_total, n_failing))

        df$status    <- ifelse(df$F_MISS > threshold, "Failing", "Passing")
        df$pass      <- pass_labels[i]
        df$threshold <- threshold

        all_data[[i]] <- df

        # ── Individual plot per pass ──────────────────────────
        p <- ggplot(df, aes(x = F_MISS, fill = status)) +
            geom_histogram(bins = 50, colour = "white", linewidth = 0.2) +
            geom_vline(xintercept = threshold, colour = "#ff5f57",
                       linetype = "dashed", linewidth = 0.8) +
            annotate("text", x = threshold, y = Inf,
                     label = paste0("MIND = ", threshold),
                     hjust = -0.1, vjust = 1.5,
                     colour = "#ff5f57", size = 3.5, fontface = "bold") +
            scale_fill_manual(
                values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
                name   = "Sample status"
            ) +
            scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
            scale_y_continuous(labels = comma) +
            labs(
                title    = sprintf("CRISP -- Step 3: Sample Call Rate -- %s", pass_labels[i]),
                subtitle = sprintf("CASCADE mode  |  %d/%d samples failing at this threshold",
                                   n_failing, n_total),
                x        = "Per-sample missingness rate (F_MISS)",
                y        = "Number of samples",
                caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
            ) +
            crisp_theme

        out_file <- file.path(report_dir,
                              sprintf("step3_callrate_cascade_pass%d.pdf", i))
        pdf(out_file, width = 10, height = 6)
        print(p)
        dev.off()

        pdf_files <- c(pdf_files, out_file)
        cat(sprintf("[PLOT]   Individual plot saved: %s\n", out_file))
    }

    # ── Faceted plot -- all passes side by side ───────────────
    cat("\n[PLOT] Generating faceted comparison plot...\n")

    combined <- do.call(rbind, all_data)
    combined$pass <- factor(combined$pass, levels = pass_labels)

    # Build per-facet threshold lines
    thresholds_df <- data.frame(
        pass      = factor(pass_labels, levels = pass_labels),
        threshold = thresholds
    )

    p_facet <- ggplot(combined, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 40, colour = "white", linewidth = 0.15) +
        geom_vline(data = thresholds_df,
                   aes(xintercept = threshold),
                   colour = "#ff5f57", linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ pass, scales = "free", ncol = 2) +
        scale_fill_manual(
            values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
            name   = "Sample status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP -- Step 3: Sample Call Rate -- CASCADE Mode",
            subtitle = "Missingness distribution at each cascade threshold pass",
            x        = "Per-sample missingness rate (F_MISS)",
            y        = "Number of samples",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme +
        theme(strip.text = element_text(face = "bold", size = 10))

    facet_file <- file.path(report_dir, "step3_callrate_cascade_faceted.pdf")
    pdf(facet_file, width = 14, height = 10)
    print(p_facet)
    dev.off()

    cat(sprintf("[PLOT] Faceted plot saved: %s\n", facet_file))

    return(c(pdf_files, facet_file))
}

##########################################################################
### DISPATCH
##########################################################################

cat("\n[PLOT] CRISP Step 3 -- Call Rate Plotting\n")
cat(sprintf("[PLOT] Mode       : %s\n", mode))
cat(sprintf("[PLOT] MIND       : %.3f\n", mind))
cat(sprintf("[PLOT] Report dir : %s\n", report_dir))
cat(sprintf("[PLOT] Files      : %d\n\n", length(imiss_files)))

if (mode == "SIMPLE") {
    plot_simple(imiss_files[1], mind, report_dir)
} else if (mode == "CASCADE") {
    plot_cascade(imiss_files, mind, report_dir)
} else {
    cat(sprintf("[PLOT] ERROR: Unknown mode '%s'. Use SIMPLE or CASCADE.\n", mode))
    quit(status = 1)
}

cat("\n[PLOT] Plotting complete.\n")
