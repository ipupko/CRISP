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
### Step 4: Variant Call Rate — Plotting Script
### Version: 0.3.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
### DESCRIPTION
### Generates per-SNP missingness distribution plots from PLINK
### .lmiss files.
###
### SIMPLE mode:
###   Single histogram of per-variant missingness rates
###   Vertical line marking the GENO threshold
###   Variants failing threshold highlighted in red
###
### CASCADE mode:
###   Four individual histograms, one per pass
###   Single faceted 2x2 plot showing all passes
###
### CUSTOM mode:
###   One individual histogram per tier
###   Single faceted plot — 2 columns, rows scale to tier count
###
### Usage:
###   Rscript plot_snprate.R <mode> <geno> <report_dir> <lmiss_files...>
###   mode       : SIMPLE, CASCADE, or CUSTOM
###   geno       : final GENO threshold (e.g. 0.05)
###   report_dir : output directory for plots
###   lmiss_files: one .lmiss file (SIMPLE) or one per tier
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    cat("Usage: Rscript plot_snprate.R <mode> <geno> <report_dir> <lmiss_file(s)>\n")
    quit(status = 1)
}

mode        <- toupper(args[1])
geno        <- as.numeric(args[2])
report_dir  <- args[3]
lmiss_files <- args[4:length(args)]
n_passes    <- length(lmiss_files)

dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

##########################################################################
### COLOUR PALETTE
### Reads CRISP_PAL_* environment variables exported by _crisp_flavour.sh's
### _crisp_palette() function in crisp_snprate.sh.
###
### WHY ENV VARS: R cannot source a bash script. Environment variable
### passing is the only clean cross-language mechanism. The _crisp_palette()
### call in crisp_snprate.sh runs before any Rscript call, so all vars
### are available in the R process environment.
###
### WHY FALLBACK VALUES: Hardcoded fallbacks match the previous hardcoded
### values exactly, so STANDARD mode produces byte-identical output to the
### pre-v0.3.0 script. This makes the palette change a zero-risk upgrade
### for users not using COLOURBLIND or NIGHT mode.
###
### NOTE: CRISP_PAL_SUBTITLE and CRISP_PAL_CAPTION are not yet exported
### by _crisp_flavour.sh. Fallbacks are used until they are added.
##########################################################################

colour_mode <- toupper(Sys.getenv("CRISP_PAL_MODE", unset = "STANDARD"))

if (colour_mode == "COLOURBLIND") {
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#0072B2")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#D55E00")
    PAL_SUBTITLE  <- Sys.getenv("CRISP_PAL_SUBTITLE",  "#555555")
    PAL_CAPTION   <- Sys.getenv("CRISP_PAL_CAPTION",   "#888888")
    BG_COLOUR     <- "#FFFFFF"
    FG_COLOUR     <- "#000000"
    GRID_COLOUR   <- "#E0E0E0"
} else if (colour_mode == "NIGHT") {
    # Dark-background variant, Okabe-Ito-derived.
    # Provisional pending final NIGHT palette in _crisp_flavour.sh.
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#56B4E9")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#E69F00")
    PAL_SUBTITLE  <- Sys.getenv("CRISP_PAL_SUBTITLE",  "#AAAAAA")
    PAL_CAPTION   <- Sys.getenv("CRISP_PAL_CAPTION",   "#777777")
    BG_COLOUR     <- "#1E1E1E"
    FG_COLOUR     <- "#E0E0E0"
    GRID_COLOUR   <- "#3A3A3A"
} else {
    # STANDARD — fallbacks match pre-v0.3.0 hardcoded values exactly
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#1D9E75")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#ff5f57")
    PAL_SUBTITLE  <- Sys.getenv("CRISP_PAL_SUBTITLE",  "#555555")
    PAL_CAPTION   <- Sys.getenv("CRISP_PAL_CAPTION",   "#888888")
    BG_COLOUR     <- "#FFFFFF"
    FG_COLOUR     <- "#000000"
    GRID_COLOUR   <- "#E0E0E0"
}

##########################################################################
### THEME
### WHY legend.position = "right": CRISP plotting standard. Right-side
### legends do not obscure the plot area and scale gracefully with long
### legend labels. The previous "bottom" position was inconsistent with
### the right-side legend standard established in maf_spectrum.R and the
### broader CRISP plotting conventions. Fixed in v0.3.0.
##########################################################################

crisp_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title        = element_text(face = "bold", size = 14,
                                         colour = FG_COLOUR),
        plot.subtitle     = element_text(colour = PAL_SUBTITLE, size = 10),
        axis.title        = element_text(face = "bold", colour = FG_COLOUR),
        axis.text         = element_text(colour = FG_COLOUR),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = GRID_COLOUR),
        plot.caption      = element_text(colour = PAL_CAPTION, size = 8),
        legend.position   = "right",
        legend.title      = element_text(face = "bold", colour = FG_COLOUR),
        legend.text       = element_text(colour = FG_COLOUR),
        strip.text        = element_text(face = "bold", size = 10,
                                         colour = FG_COLOUR),
        plot.background   = element_rect(fill = BG_COLOUR, colour = NA),
        panel.background  = element_rect(fill = BG_COLOUR, colour = NA),
        legend.background = element_rect(fill = BG_COLOUR, colour = NA),
        legend.key        = element_rect(fill = BG_COLOUR, colour = NA)
    )

##########################################################################
### HELPER: LOAD LMISS FILE
##########################################################################

load_lmiss <- function(filepath) {
    df <- fread(filepath, header = TRUE)
    # PLINK .lmiss columns: CHR SNP N_MISS N_GENO F_MISS
    if (!"F_MISS" %in% colnames(df)) {
        stop(paste("F_MISS column not found in:", filepath))
    }
    return(df)
}

##########################################################################
### SIMPLE MODE
##########################################################################

plot_simple <- function(lmiss_file, geno, report_dir) {

    cat("[PLOT] Loading:", lmiss_file, "\n")
    df <- load_lmiss(lmiss_file)

    n_total   <- nrow(df)
    n_failing <- sum(df$F_MISS > geno)
    n_passing <- n_total - n_failing

    cat(sprintf("[PLOT] Variants total   : %d\n", n_total))
    cat(sprintf("[PLOT] Variants failing : %d (F_MISS > %.3f)\n", n_failing, geno))
    cat(sprintf("[PLOT] Variants passing : %d\n", n_passing))

    df$status <- ifelse(df$F_MISS > geno, "Failing", "Passing")

    p <- ggplot(df, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = geno, colour = PAL_FAIL,
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = geno, y = Inf,
                 label = paste0("GENO = ", geno),
                 hjust = -0.1, vjust = 1.5,
                 colour = PAL_FAIL, size = 3.5, fontface = "bold") +
        scale_fill_manual(
            values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
            name   = "Variant status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1),
                           limits = c(0, max(df$F_MISS) * 1.05)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP — Step 4: Variant Call Rate",
            subtitle = sprintf("SIMPLE mode  |  GENO threshold: %.3f  |  %s/%s variants failing",
                               geno,
                               format(n_failing, big.mark = ","),
                               format(n_total, big.mark = ",")),
            x        = "Per-variant missingness rate (F_MISS)",
            y        = "Number of variants",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    out_file <- file.path(report_dir, "step4_snprate_simple.pdf")
    pdf(out_file, width = 10, height = 6)
    print(p)
    dev.off()

    cat(sprintf("[PLOT] Simple histogram saved: %s\n", out_file))
    return(out_file)
}

##########################################################################
### MULTI-PASS MODE (CASCADE and CUSTOM)
### Handles any number of passes — 2 column faceted grid
##########################################################################

plot_multipass <- function(lmiss_files, thresholds, mode, report_dir) {

    n <- length(lmiss_files)

    if (length(thresholds) != n) {
        stop("Number of thresholds must match number of lmiss files")
    }

    pass_labels <- sprintf("Pass %d (geno=%.2f)", seq_along(thresholds), thresholds)
    all_data    <- list()
    pdf_files   <- c()

    for (i in seq_along(lmiss_files)) {

        cat(sprintf("\n[PLOT] Pass %d — Loading: %s\n", i, lmiss_files[i]))
        df        <- load_lmiss(lmiss_files[i])
        threshold <- thresholds[i]

        n_total   <- nrow(df)
        n_failing <- sum(df$F_MISS > threshold)

        cat(sprintf("[PLOT]   Threshold : %.3f\n", threshold))
        cat(sprintf("[PLOT]   Variants  : %s total, %s failing\n",
                    format(n_total, big.mark = ","),
                    format(n_failing, big.mark = ",")))

        df$status    <- ifelse(df$F_MISS > threshold, "Failing", "Passing")
        df$pass      <- pass_labels[i]
        df$threshold <- threshold

        all_data[[i]] <- df

        # Individual plot per pass
        p <- ggplot(df, aes(x = F_MISS, fill = status)) +
            geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
            geom_vline(xintercept = threshold, colour = PAL_FAIL,
                       linetype = "dashed", linewidth = 0.8) +
            annotate("text", x = threshold, y = Inf,
                     label = paste0("GENO = ", threshold),
                     hjust = -0.1, vjust = 1.5,
                     colour = PAL_FAIL, size = 3.5, fontface = "bold") +
            scale_fill_manual(
                values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
                name   = "Variant status"
            ) +
            scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
            scale_y_continuous(labels = comma) +
            labs(
                title    = sprintf("CRISP — Step 4: Variant Call Rate — %s", pass_labels[i]),
                subtitle = sprintf("%s mode  |  %s/%s variants failing at this threshold",
                                   mode,
                                   format(n_failing, big.mark = ","),
                                   format(n_total, big.mark = ",")),
                x        = "Per-variant missingness rate (F_MISS)",
                y        = "Number of variants",
                caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
            ) +
            crisp_theme

        out_file <- file.path(report_dir,
                              sprintf("step4_snprate_%s_pass%d.pdf",
                                      tolower(mode), i))
        pdf(out_file, width = 10, height = 6)
        print(p)
        dev.off()

        pdf_files <- c(pdf_files, out_file)
        cat(sprintf("[PLOT]   Individual plot saved: %s\n", out_file))
    }

    # Faceted plot — 2 columns, rows scale to number of passes
    cat(sprintf("\n[PLOT] Generating faceted %s comparison plot (%d passes)...\n",
                mode, n))

    combined <- do.call(rbind, all_data)
    combined$pass <- factor(combined$pass, levels = pass_labels)

    thresholds_df <- data.frame(
        pass      = factor(pass_labels, levels = pass_labels),
        threshold = thresholds
    )

    n_cols <- 2
    n_rows <- ceiling(n / n_cols)
    fig_h  <- max(6, n_rows * 4.5)
    fig_w  <- 14

    p_facet <- ggplot(combined, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 50, colour = "white", linewidth = 0.15) +
        geom_vline(data = thresholds_df,
                   aes(xintercept = threshold),
                   colour = PAL_FAIL, linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ pass, scales = "free", ncol = n_cols) +
        scale_fill_manual(
            values = c("Failing" = PAL_FAIL, "Passing" = PAL_PASS),
            name   = "Variant status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = sprintf("CRISP — Step 4: Variant Call Rate — %s Mode", mode),
            subtitle = sprintf("Per-variant missingness distribution across %d threshold passes",
                               n),
            x        = "Per-variant missingness rate (F_MISS)",
            y        = "Number of variants",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    facet_file <- file.path(report_dir,
                            sprintf("step4_snprate_%s_faceted.pdf", tolower(mode)))
    pdf(facet_file, width = fig_w, height = fig_h)
    print(p_facet)
    dev.off()

    cat(sprintf("[PLOT] Faceted plot saved: %s\n", facet_file))

    return(c(pdf_files, facet_file))
}

##########################################################################
### DISPATCH
##########################################################################

cat("\n[PLOT] CRISP — Step 4: Variant Call Rate Plotting\n")
cat(sprintf("[PLOT] Mode       : %s\n", mode))
cat(sprintf("[PLOT] GENO       : %.3f\n", geno))
cat(sprintf("[PLOT] Report dir : %s\n", report_dir))
cat(sprintf("[PLOT] Files      : %d\n\n", n_passes))

if (mode == "SIMPLE") {
    plot_simple(lmiss_files[1], geno, report_dir)

} else if (mode == "CASCADE") {
    # Fixed cascade thresholds
    thresholds <- c(0.25, 0.20, 0.10, 0.05)
    if (n_passes != 4) {
        cat(sprintf("[PLOT] ERROR: CASCADE expects 4 .lmiss files, got %d\n", n_passes))
        quit(status = 1)
    }
    plot_multipass(lmiss_files, thresholds, "CASCADE", report_dir)

} else if (mode == "CUSTOM") {
    # Thresholds passed as fifth argument onwards — extract from filenames
    # or passed explicitly. We infer from the number of files.
    # Thresholds are embedded in filenames by the shell script.
    # Extract from file names: step4_custom_pass1_geno0.30.lmiss
    thresholds <- numeric(n_passes)
    for (i in seq_along(lmiss_files)) {
        fname <- basename(lmiss_files[i])
        m     <- regmatches(fname, regexpr("geno([0-9.]+)", fname))
        if (length(m) > 0) {
            thresholds[i] <- as.numeric(sub("geno", "", m))
        } else {
            # Fall back to evenly spaced if cannot parse
            thresholds[i] <- geno + (n_passes - i) * 0.05
        }
    }
    plot_multipass(lmiss_files, thresholds, "CUSTOM", report_dir)

} else {
    cat(sprintf("[PLOT] ERROR: Unknown mode '%s'. Use SIMPLE, CASCADE, or CUSTOM.\n", mode))
    quit(status = 1)
}

cat("\n[PLOT] Plotting complete.\n")        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(colour = "#888888", size = 8),
        legend.position  = "bottom",
        strip.text       = element_text(face = "bold", size = 10)
    )

##########################################################################
# HELPER: LOAD LMISS FILE
##########################################################################

load_lmiss <- function(filepath) {
    df <- fread(filepath, header = TRUE)
    # PLINK .lmiss columns: CHR SNP N_MISS N_GENO F_MISS
    if (!"F_MISS" %in% colnames(df)) {
        stop(paste("F_MISS column not found in:", filepath))
    }
    return(df)
}

##########################################################################
# SIMPLE MODE
##########################################################################

plot_simple <- function(lmiss_file, geno, report_dir) {

    cat("[PLOT] Loading:", lmiss_file, "\n")
    df <- load_lmiss(lmiss_file)

    n_total   <- nrow(df)
    n_failing <- sum(df$F_MISS > geno)
    n_passing <- n_total - n_failing

    cat(sprintf("[PLOT] Variants total   : %d\n", n_total))
    cat(sprintf("[PLOT] Variants failing : %d (F_MISS > %.3f)\n", n_failing, geno))
    cat(sprintf("[PLOT] Variants passing : %d\n", n_passing))

    df$status <- ifelse(df$F_MISS > geno, "Failing", "Passing")

    p <- ggplot(df, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = geno, colour = "#ff5f57",
                   linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = geno, y = Inf,
                 label = paste0("GENO = ", geno),
                 hjust = -0.1, vjust = 1.5,
                 colour = "#ff5f57", size = 3.5, fontface = "bold") +
        scale_fill_manual(
            values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
            name   = "Variant status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1),
                           limits = c(0, max(df$F_MISS) * 1.05)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP -- Step 4: Variant Call Rate",
            subtitle = sprintf("SIMPLE mode  |  GENO threshold: %.3f  |  %s/%s variants failing",
                               geno,
                               format(n_failing, big.mark = ","),
                               format(n_total, big.mark = ",")),
            x        = "Per-variant missingness rate (F_MISS)",
            y        = "Number of variants",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    out_file <- file.path(report_dir, "step4_snprate_simple.pdf")
    pdf(out_file, width = 10, height = 6)
    print(p)
    dev.off()

    cat(sprintf("[PLOT] Simple histogram saved: %s\n", out_file))
    return(out_file)
}

##########################################################################
# MULTI-PASS MODE (CASCADE and CUSTOM)
# Handles any number of passes -- 2 column faceted grid
##########################################################################

plot_multipass <- function(lmiss_files, thresholds, mode, report_dir) {

    n <- length(lmiss_files)

    if (length(thresholds) != n) {
        stop("Number of thresholds must match number of lmiss files")
    }

    pass_labels <- sprintf("Pass %d (geno=%.2f)", seq_along(thresholds), thresholds)
    all_data    <- list()
    pdf_files   <- c()

    for (i in seq_along(lmiss_files)) {

        cat(sprintf("\n[PLOT] Pass %d -- Loading: %s\n", i, lmiss_files[i]))
        df        <- load_lmiss(lmiss_files[i])
        threshold <- thresholds[i]

        n_total   <- nrow(df)
        n_failing <- sum(df$F_MISS > threshold)

        cat(sprintf("[PLOT]   Threshold : %.3f\n", threshold))
        cat(sprintf("[PLOT]   Variants  : %s total, %s failing\n",
                    format(n_total, big.mark = ","),
                    format(n_failing, big.mark = ",")))

        df$status    <- ifelse(df$F_MISS > threshold, "Failing", "Passing")
        df$pass      <- pass_labels[i]
        df$threshold <- threshold

        all_data[[i]] <- df

        # Individual plot per pass
        p <- ggplot(df, aes(x = F_MISS, fill = status)) +
            geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
            geom_vline(xintercept = threshold, colour = "#ff5f57",
                       linetype = "dashed", linewidth = 0.8) +
            annotate("text", x = threshold, y = Inf,
                     label = paste0("GENO = ", threshold),
                     hjust = -0.1, vjust = 1.5,
                     colour = "#ff5f57", size = 3.5, fontface = "bold") +
            scale_fill_manual(
                values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
                name   = "Variant status"
            ) +
            scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
            scale_y_continuous(labels = comma) +
            labs(
                title    = sprintf("CRISP -- Step 4: Variant Call Rate -- %s", pass_labels[i]),
                subtitle = sprintf("%s mode  |  %s/%s variants failing at this threshold",
                                   mode,
                                   format(n_failing, big.mark = ","),
                                   format(n_total, big.mark = ",")),
                x        = "Per-variant missingness rate (F_MISS)",
                y        = "Number of variants",
                caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
            ) +
            crisp_theme

        out_file <- file.path(report_dir,
                              sprintf("step4_snprate_%s_pass%d.pdf",
                                      tolower(mode), i))
        pdf(out_file, width = 10, height = 6)
        print(p)
        dev.off()

        pdf_files <- c(pdf_files, out_file)
        cat(sprintf("[PLOT]   Individual plot saved: %s\n", out_file))
    }

    # Faceted plot -- 2 columns, rows scale to number of passes
    cat(sprintf("\n[PLOT] Generating faceted %s comparison plot (%d passes)...\n",
                mode, n))

    combined <- do.call(rbind, all_data)
    combined$pass <- factor(combined$pass, levels = pass_labels)

    thresholds_df <- data.frame(
        pass      = factor(pass_labels, levels = pass_labels),
        threshold = thresholds
    )

    n_cols <- 2
    n_rows <- ceiling(n / n_cols)
    fig_h  <- max(6, n_rows * 4.5)
    fig_w  <- 14

    p_facet <- ggplot(combined, aes(x = F_MISS, fill = status)) +
        geom_histogram(bins = 50, colour = "white", linewidth = 0.15) +
        geom_vline(data = thresholds_df,
                   aes(xintercept = threshold),
                   colour = "#ff5f57", linetype = "dashed", linewidth = 0.7) +
        facet_wrap(~ pass, scales = "free", ncol = n_cols) +
        scale_fill_manual(
            values = c("Failing" = "#ff5f57", "Passing" = "#1D9E75"),
            name   = "Variant status"
        ) +
        scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = sprintf("CRISP -- Step 4: Variant Call Rate -- %s Mode", mode),
            subtitle = sprintf("Per-variant missingness distribution across %d threshold passes",
                               n),
            x        = "Per-variant missingness rate (F_MISS)",
            y        = "Number of variants",
            caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
        ) +
        crisp_theme

    facet_file <- file.path(report_dir,
                            sprintf("step4_snprate_%s_faceted.pdf", tolower(mode)))
    pdf(facet_file, width = fig_w, height = fig_h)
    print(p_facet)
    dev.off()

    cat(sprintf("[PLOT] Faceted plot saved: %s\n", facet_file))

    return(c(pdf_files, facet_file))
}

##########################################################################
# DISPATCH
##########################################################################

cat("\n[PLOT] CRISP Step 4 -- Variant Call Rate Plotting\n")
cat(sprintf("[PLOT] Mode       : %s\n", mode))
cat(sprintf("[PLOT] GENO       : %.3f\n", geno))
cat(sprintf("[PLOT] Report dir : %s\n", report_dir))
cat(sprintf("[PLOT] Files      : %d\n\n", n_passes))

if (mode == "SIMPLE") {
    plot_simple(lmiss_files[1], geno, report_dir)

} else if (mode == "CASCADE") {
    # Fixed cascade thresholds
    thresholds <- c(0.25, 0.20, 0.10, 0.05)
    if (n_passes != 4) {
        cat(sprintf("[PLOT] ERROR: CASCADE expects 4 .lmiss files, got %d\n", n_passes))
        quit(status = 1)
    }
    plot_multipass(lmiss_files, thresholds, "CASCADE", report_dir)

} else if (mode == "CUSTOM") {
    # Thresholds passed as fifth argument onwards -- extract from filenames
    # or passed explicitly. We infer from the number of files.
    # Thresholds are embedded in filenames by the shell script.
    # Extract from file names: step4_custom_pass1_geno0.30.lmiss
    thresholds <- numeric(n_passes)
    for (i in seq_along(lmiss_files)) {
        fname <- basename(lmiss_files[i])
        m     <- regmatches(fname, regexpr("geno([0-9.]+)", fname))
        if (length(m) > 0) {
            thresholds[i] <- as.numeric(sub("geno", "", m))
        } else {
            # Fall back to evenly spaced if cannot parse
            thresholds[i] <- geno + (n_passes - i) * 0.05
        }
    }
    plot_multipass(lmiss_files, thresholds, "CUSTOM", report_dir)

} else {
    cat(sprintf("[PLOT] ERROR: Unknown mode '%s'. Use SIMPLE, CASCADE, or CUSTOM.\n", mode))
    quit(status = 1)
}

cat("\n[PLOT] Plotting complete.\n")
