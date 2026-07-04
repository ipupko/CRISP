#!/usr/bin/env Rscript
##########################################################################
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 3: Sample Call Rate, Plotting Script
# Version: 0.3.5
# Developed by Igor Pupko
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# DESCRIPTION
# Generates missingness distribution plots from PLINK .imiss files.
#
# SIMPLE mode:
#   Single histogram with MIND threshold line, failing samples highlighted.
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
#
# Colour/background are read from environment variables (set by
# crisp_callrate.sh via _crisp_flavour.sh / _crisp_palette):
#   CRISP_PAL_PASS, CRISP_PAL_FAIL, CRISP_PAL_BG, CRISP_PAL_PANEL,
#   CRISP_PAL_TEXT, CRISP_PAL_SUBTEXT, CRISP_PAL_GRID,
#   CRISP_PAL_MODE, CRISP_PAL_BACKGROUND
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

##########################################################################
# PALETTE (from environment, with safe STANDARD/light fallbacks)
##########################################################################

env_or <- function(name, default) {
    val <- Sys.getenv(name)
    if (nzchar(val)) val else default
}

PAL_PASS    <- env_or("CRISP_PAL_PASS",    "#1D9E75")
PAL_FAIL    <- env_or("CRISP_PAL_FAIL",    "#ff5f57")
PAL_BG      <- env_or("CRISP_PAL_BG",      "#FFFFFF")
PAL_PANEL   <- env_or("CRISP_PAL_PANEL",   "#f8f9fa")
PAL_TEXT    <- env_or("CRISP_PAL_TEXT",    "#1a1a1a")
PAL_SUBTEXT <- env_or("CRISP_PAL_SUBTEXT", "#555555")
PAL_GRID    <- env_or("CRISP_PAL_GRID",    "#e0e0e0")
PAL_MODE    <- env_or("CRISP_PAL_MODE",    "STANDARD")
PAL_BACKGROUND <- env_or("CRISP_PAL_BACKGROUND", "LIGHT")

cat(sprintf("[PLOT] Colour mode : %s\n", PAL_MODE))
cat(sprintf("[PLOT] Background  : %s\n", PAL_BACKGROUND))

# shared ggplot2 theme -- right-side legend (CRISP standard) and full
# light/dark background support
crisp_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(face = "bold", size = 14, colour = PAL_TEXT),
        plot.subtitle    = element_text(colour = PAL_SUBTEXT, size = 10),
        axis.title       = element_text(face = "bold", colour = PAL_TEXT),
        axis.text        = element_text(colour = PAL_SUBTEXT),
        panel.grid.major = element_line(colour = PAL_GRID),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(colour = PAL_SUBTEXT, size = 8),
        legend.position  = "right",
        legend.title     = element_text(colour = PAL_TEXT),
        legend.text      = element_text(colour = PAL_TEXT),
        strip.text       = element_text(face = "bold", size = 10, colour = PAL_TEXT),
        strip.background = element_rect(fill = PAL_PANEL, colour = NA),
        plot.background  = element_rect(fill = PAL_BG, colour = NA),
        panel.background = element_rect(fill = PAL_PANEL, colour = NA),
        legend.background = element_rect(fill = PAL_BG, colour = NA),
        legend.key       = element_rect(fill = PAL_BG, colour = NA)
    )

# load a .imiss file, check it has F_MISS
load_imiss <- function(filepath) {
    if (!file.exists(filepath)) {
        cat(sprintf("[PLOT] ERROR: .imiss file not found: %s\n", filepath))
        quit(status = 1)
    }
    df <- tryCatch(
        fread(filepath, header = TRUE),
        error = function(e) {
            cat(sprintf("[PLOT] ERROR: Failed to read %s: %s\n", filepath, conditionMessage(e)))
            quit(status = 1)
        }
    )
    if (!"F_MISS" %in% colnames(df)) {
        cat(sprintf("[PLOT] ERROR: F_MISS column not found in: %s\n", filepath))
        quit(status = 1)
    }
    return(df)
}

# extract threshold from filename, e.g. mind0.20.imiss -> 0.20
parse_threshold <- function(filepath, fallback) {
    m <- regmatches(basename(filepath),
                    regexpr("mind([0-9]+\\.?[0-9]*)", basename(filepath)))
    if (length(m) > 0) as.numeric(sub("mind", "", m)) else fallback
}

# sensible x-axis upper limit: at least 2x the threshold, but never
# less than the max observed F_MISS * 1.05 (so points aren't clipped)
x_axis_limit <- function(df, threshold) {
    max(max(df$F_MISS) * 1.05, threshold * 2)
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
        geom_histogram(bins = 50, colour = PAL_BG, linewidth = 0.2) +
        geom_vline(xintercept = mind, colour = PAL_FAIL,
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = mind, y = Inf,
                 label = paste0("MIND = ", mind),
                 hjust = -0.1, vjust = 1.5,
                 colour = PAL_FAIL, size = 3.5, fontface = "bold") +
        scale_fill_manual(
            values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
            name   = "Sample status"
        ) +
        scale_x_continuous(
            labels = percent_format(accuracy = 0.1),
            limits = c(0, x_axis_limit(df, mind))
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
    pdf(out_file, width = 10, height = 6, bg = PAL_BG)
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
            geom_histogram(bins = 50, colour = PAL_BG, linewidth = 0.2) +
            geom_vline(xintercept = threshold, colour = PAL_FAIL,
                       linetype = "dashed", linewidth = 0.8) +
            annotate("text", x = threshold, y = Inf,
                     label = paste0("MIND = ", threshold),
                     hjust = -0.1, vjust = 1.5,
                     colour = PAL_FAIL, size = 3.5, fontface = "bold") +
            scale_fill_manual(
                values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
                name   = "Sample status"
            ) +
            scale_x_continuous(
                labels = percent_format(accuracy = 0.1),
                limits = c(0, x_axis_limit(df, threshold))
            ) +
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
        pdf(out_file, width = 10, height = 6, bg = PAL_BG)
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
        geom_histogram(bins = 40, colour = PAL_BG, linewidth = 0.15) +
        geom_vline(data = thresholds_df,
                   aes(xintercept = threshold),
                   colour = PAL_FAIL, linetype = "dashed",
                   linewidth = 0.7) +
        facet_wrap(~ pass, scales = "free", ncol = n_cols) +
        scale_fill_manual(
            values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
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
    pdf(facet_file, width = 14, height = fig_h, bg = PAL_BG)
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

    if (n_passes != 4) {
        cat(sprintf("[PLOT] ERROR: CASCADE expects 4 .imiss files, got %d\n", n_passes))
        quit(status = 1)
    }
    # thresholds parsed from filenames (mind{value} pattern), with the
    # documented CASCADE defaults as fallback if a filename is unparseable
    cascade_defaults <- c(0.25, 0.20, 0.10, 0.05)
    thresholds <- sapply(seq_along(imiss_files), function(i) {
        parse_threshold(imiss_files[i], fallback = cascade_defaults[i])
    })
    plot_multipass(imiss_files, thresholds, "CASCADE", report_dir)

} else if (mode == "CUSTOM") {

    if (n_passes < 2) {
        cat(sprintf("[PLOT] ERROR: CUSTOM mode requires at least 2 .imiss files, got %d\n", n_passes))
        quit(status = 1)
    }
    # infer thresholds from filenames; fallback descends from `mind` in
    # 0.05 steps for any pass whose filename doesn't match mind{value}
    thresholds <- sapply(seq_along(imiss_files), function(i) {
        parse_threshold(imiss_files[i],
                        fallback = mind + (n_passes - i) * 0.05)
    })
    plot_multipass(imiss_files, thresholds, "CUSTOM", report_dir)

} else {

    cat(sprintf("[PLOT] ERROR: Unknown mode '%s'. Use SIMPLE, CASCADE or CUSTOM.\n", mode))
    quit(status = 1)
}

cat("\n[PLOT] Plotting complete.\n")
