#!/usr/bin/env Rscript
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### Step 7: Relatedness, Plotting Script
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Generates relatedness plots from PLINK IBD or KING output.
#
# IBD mode:
#   Plot 1: IBD0 vs IBD1 scatter coloured by relationship class
#           Theoretical relationship lines overlaid
#           (parent-offspring, full siblings, half siblings, etc.)
#   Plot 2: PI_HAT distribution histogram with threshold line
#
# KING mode:
#   Plot 1: Kinship coefficient distribution histogram
#   Plot 2: Kinship vs IBS0 scatter coloured by relationship class
#
# Relationship classification thresholds:
#   IBD  : PI_HAT >= 0.9   duplicate / MZ twin
#          PI_HAT >= 0.4   first degree
#          PI_HAT >= 0.185 second degree
#   KING : kinship >= 0.354  duplicate / MZ twin
#          kinship >= 0.177  first degree
#          kinship >= 0.0884 second degree
#
# Usage:
#   Rscript plot_relatedness.R <output_dir> <relatedness_file>
#                              <method> <pihat_cutoff> <king_cutoff>
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    cat("Usage: Rscript plot_relatedness.R <output_dir> <relatedness_file>",
        "<method> <pihat_cutoff> <king_cutoff>\n")
    quit(status = 1)
}

output_dir    <- args[1]
rel_file      <- args[2]
method        <- toupper(args[3])
pihat_cutoff  <- as.numeric(args[4])
king_cutoff   <- as.numeric(args[5])

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[PLOT] CRISP Step 7: Relatedness plotting\n"))
cat(sprintf("[PLOT] Method  : %s\n", method))
cat(sprintf("[PLOT] File    : %s\n", rel_file))

if (!file.exists(rel_file)) {
    cat(sprintf("[PLOT] ERROR: Relatedness file not found: %s\n", rel_file))
    quit(status = 1)
}

# shared theme
crisp_theme <- theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "#555555", size = 10),
        axis.title       = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#f8f9fa", colour = NA),
        plot.caption     = element_text(colour = "#888888", size = 8),
        legend.position  = "right",
        legend.title     = element_text(face = "bold", size = 9),
        strip.text       = element_text(face = "bold")
    )

CAPTION <- "CRISP | Comprehensive Robust Integrated SNP Processing | Compass Genomics"

classify_ibd <- function(pihat) {
    dplyr::case_when(
        pihat >= 0.9   ~ "Duplicate / MZ twin",
        pihat >= 0.4   ~ "First degree",
        pihat >= 0.185 ~ "Second degree",
        TRUE           ~ "Distant"
    )
}

classify_king <- function(kinship) {
    dplyr::case_when(
        kinship >= 0.354  ~ "Duplicate / MZ twin",
        kinship >= 0.177  ~ "First degree",
        kinship >= 0.0884 ~ "Second degree",
        TRUE              ~ "Distant"
    )
}

REL_COLORS <- c(
    "Duplicate / MZ twin" = "#ff5f57",
    "First degree"        = "#ff8c00",
    "Second degree"       = "#febc2e",
    "Distant"             = "#1D9E75"
)

# label pairs of interest with boundary-aware collision-resolving arrows
# only labels Duplicate, First degree and Second degree
label_pairs <- function(df_flag, rel_label, col,
                        x_col, y_col,
                        iid1_col = "IID1", iid2_col = "IID2",
                        show_labels = TRUE) {
    if (nrow(df_flag) == 0 || !show_labels) return(invisible(NULL))

    usr   <- par("usr")
    x_rng <- usr[2] - usr[1]
    y_rng <- usr[4] - usr[3]

    candidates <- list(
        c( x_rng*0.10,  y_rng*0.10),
        c(-x_rng*0.10,  y_rng*0.10),
        c( x_rng*0.10, -y_rng*0.10),
        c(-x_rng*0.10, -y_rng*0.10),
        c( x_rng*0.16,  y_rng*0.05),
        c(-x_rng*0.16,  y_rng*0.05),
        c( x_rng*0.05,  y_rng*0.16),
        c( x_rng*0.05, -y_rng*0.16)
    )

    margin  <- x_rng * 0.01
    placed  <- list()

    points(df_flag[[x_col]], df_flag[[y_col]],
           col = col, pch = 8, cex = 1.6, lwd = 2.0)

    for (i in seq_len(nrow(df_flag))) {
        lbl <- paste0(df_flag[[iid1_col]][i], ":",
                      df_flag[[iid2_col]][i], " (", rel_label, ")")
        x0  <- df_flag[[x_col]][i]
        y0  <- df_flag[[y_col]][i]
        lbl_w <- strwidth(lbl,  cex = 0.55)
        lbl_h <- strheight(lbl, cex = 0.55)

        best_x <- NA; best_y <- NA

        for (cand in candidates) {
            cx <- x0 + cand[1]; cy <- y0 + cand[2]
            if (cx - lbl_w/2 < usr[1] + margin) next
            if (cx + lbl_w/2 > usr[2] - margin) next
            if (cy - lbl_h   < usr[3] + margin) next
            if (cy           > usr[4] - margin) next
            overlap <- FALSE
            for (p in placed) {
                if (abs(cx - p[1]) < lbl_w * 0.9 &&
                    abs(cy - p[2]) < lbl_h * 1.5) {
                    overlap <- TRUE; break
                }
            }
            if (!overlap) { best_x <- cx; best_y <- cy; break }
        }

        if (is.na(best_x)) {
            best_x <- x0 + candidates[[1]][1]
            best_y <- y0 + candidates[[1]][2]
            best_x <- max(usr[1]+margin+lbl_w/2,
                          min(usr[2]-margin-lbl_w/2, best_x))
            best_y <- max(usr[3]+margin+lbl_h,
                          min(usr[4]-margin, best_y))
        }

        placed <- c(placed, list(c(best_x, best_y)))

        arrows(x0, y0, best_x, best_y,
               col = col, length = 0.07, lwd = 1.0, angle = 20)
        text(best_x, best_y, labels = lbl,
             col = col, cex = 0.55, font = 2, adj = c(0.5, 0.5))
    }
    return(invisible(NULL))
}

##########################################################################
### IBD PLOTS
##########################################################################

if (method == "IBD") {

    cat("[PLOT] Loading IBD .genome file...\n")
    df <- fread(rel_file, header = TRUE)

    if (nrow(df) == 0) {
        cat("[PLOT] No pairs above PI_HAT threshold. Skipping plots.\n")
        quit(status = 0)
    }

    df[, RELATIONSHIP := classify_ibd(PI_HAT)]
    df[, RELATIONSHIP := factor(RELATIONSHIP,
        levels = c("Duplicate / MZ twin","First degree",
                   "Second degree","Distant"))]

    cat(sprintf("[PLOT] Pairs loaded : %d\n", nrow(df)))
    print(table(df$RELATIONSHIP))

    # Plot 1: IBD0 vs IBD1 scatter
    p1 <- ggplot(df, aes(x = Z0, y = Z1, colour = RELATIONSHIP)) +
        geom_point(alpha = 0.65, size = 1.8) +
        # theoretical relationship lines
        annotate("point", x = 0,    y = 1,    size = 4,
                 shape = 3, colour = "#333", stroke = 1.2) +
        annotate("text",  x = 0,    y = 1.03, label = "Parent-offspring",
                 size = 2.8, colour = "#333", hjust = 0) +
        annotate("point", x = 0.25, y = 0.5,  size = 4,
                 shape = 3, colour = "#333", stroke = 1.2) +
        annotate("text",  x = 0.26, y = 0.53, label = "Full siblings",
                 size = 2.8, colour = "#333", hjust = 0) +
        annotate("point", x = 0.5,  y = 0.5,  size = 4,
                 shape = 3, colour = "#555", stroke = 1.0) +
        annotate("text",  x = 0.51, y = 0.53, label = "Half siblings",
                 size = 2.8, colour = "#555", hjust = 0) +
        annotate("point", x = 1,    y = 0,    size = 4,
                 shape = 3, colour = "#888", stroke = 1.0) +
        annotate("text",  x = 0.96, y = 0.03, label = "Unrelated",
                 size = 2.8, colour = "#888", hjust = 1) +
        scale_colour_manual(values = REL_COLORS, name = "Relationship") +
        scale_x_continuous(limits = c(-0.05, 1.05)) +
        scale_y_continuous(limits = c(-0.05, 1.05)) +
        labs(
            title    = "CRISP: Step 7 Relatedness, IBD0 vs IBD1",
            subtitle = sprintf(
                "PI_HAT cutoff: %.3f  |  %d related pairs",
                pihat_cutoff, nrow(df)),
            x       = "IBD0 (proportion of alleles sharing 0 copies)",
            y       = "IBD1 (proportion of alleles sharing 1 copy)",
            caption = CAPTION
        ) +
        crisp_theme

    out1 <- file.path(output_dir, "Relatedness_IBD_scatter.pdf")
    pdf(out1, width = 10, height = 7)
    print(p1)
    # label duplicates and MZ twins only
    label_pairs(df[df$RELATIONSHIP == "Duplicate / MZ twin"],
                "Dup/MZ", "#ff5f57", "Z0", "Z1")
    dev.off()
    cat(sprintf("[PLOT] IBD scatter saved: %s\n", out1))

    # Plot 2: PI_HAT distribution
    p2 <- ggplot(df, aes(x = PI_HAT, fill = RELATIONSHIP)) +
        geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = pihat_cutoff,
                   colour = "#ff5f57", linetype = "dashed", linewidth = 0.9) +
        annotate("text",
                 x     = pihat_cutoff + 0.01,
                 y     = Inf,
                 label = paste0("PI_HAT = ", pihat_cutoff),
                 hjust = 0, vjust = 1.5,
                 colour = "#ff5f57", size = 3.2, fontface = "bold") +
        geom_vline(xintercept = 0.4,   colour = "#ff8c00",
                   linetype = "dotted", linewidth = 0.8) +
        geom_vline(xintercept = 0.9,   colour = "#ff5f57",
                   linetype = "dotted", linewidth = 0.8) +
        scale_fill_manual(values = REL_COLORS, name = "Relationship") +
        scale_x_continuous(limits = c(pihat_cutoff - 0.01, 1.01)) +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP: Step 7 Relatedness, PI_HAT Distribution",
            subtitle = sprintf(
                "%d pairs above PI_HAT >= %.3f",
                nrow(df), pihat_cutoff),
            x       = "PI_HAT (proportion of IBD)",
            y       = "Number of pairs",
            caption = CAPTION
        ) +
        crisp_theme

    out2 <- file.path(output_dir, "Relatedness_IBD_pihat.pdf")
    pdf(out2, width = 10, height = 6)
    print(p2)
    dev.off()
    cat(sprintf("[PLOT] PI_HAT histogram saved: %s\n", out2))

    # text report
    report_file <- file.path(output_dir, "relatedness_report.txt")
    sink(report_file)
    cat("==================================================================\n")
    cat("  CRISP: STEP 7 RELATEDNESS PAIRS REPORT\n")
    cat("  Comprehensive Robust Integrated SNP Processing\n")
    cat("==================================================================\n")
    cat(sprintf("  Method       : IBD\n"))
    cat(sprintf("  PI_HAT cutoff: %.3f\n", pihat_cutoff))
    cat(sprintf("  Total pairs  : %d\n", nrow(df)))
    cat("------------------------------------------------------------------\n")
    for (rel in c("Duplicate / MZ twin","First degree","Second degree","Distant")) {
        grp <- df[df$RELATIONSHIP == rel]
        cat(sprintf("\n  %s (%d pairs)\n", rel, nrow(grp)))
        cat(rep("-", 60), "\n", sep="")
        if (nrow(grp) > 0) {
            for (j in seq_len(nrow(grp))) {
                cat(sprintf("  FID1=%-12s IID1=%-12s  FID2=%-12s IID2=%-12s  PI_HAT=%.4f  Z0=%.3f  Z1=%.3f  Z2=%.3f\n",
                    grp$FID1[j], grp$IID1[j],
                    grp$FID2[j], grp$IID2[j],
                    grp$PI_HAT[j], grp$Z0[j], grp$Z1[j], grp$Z2[j]))
            }
        }
    }
    cat("\n==================================================================\n")
    cat("  END OF REPORT\n")
    cat("==================================================================\n")
    sink()
    cat(sprintf("[PLOT] Relatedness report saved: %s\n", report_file))

##########################################################################
### KING PLOTS
##########################################################################

} else if (method == "KING") {

    cat("[PLOT] Loading KING .kin0 file...\n")
    df <- fread(rel_file, header = TRUE)

    if (nrow(df) == 0) {
        cat("[PLOT] No pairs above KING threshold. Skipping plots.\n")
        quit(status = 0)
    }

    setnames(df, old = "#FID1", new = "FID1", skip_absent = TRUE)

    df[, RELATIONSHIP := classify_king(KINSHIP)]
    df[, RELATIONSHIP := factor(RELATIONSHIP,
        levels = c("Duplicate / MZ twin","First degree",
                   "Second degree","Distant"))]

    cat(sprintf("[PLOT] Pairs loaded : %d\n", nrow(df)))
    print(table(df$RELATIONSHIP))

    # Plot 1: Kinship distribution
    p1 <- ggplot(df, aes(x = KINSHIP, fill = RELATIONSHIP)) +
        geom_histogram(bins = 60, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = king_cutoff,
                   colour = "#ff5f57", linetype = "dashed", linewidth = 0.9) +
        annotate("text",
                 x     = king_cutoff + 0.005,
                 y     = Inf,
                 label = paste0("KING = ", king_cutoff),
                 hjust = 0, vjust = 1.5,
                 colour = "#ff5f57", size = 3.2, fontface = "bold") +
        geom_vline(xintercept = 0.177, colour = "#ff8c00",
                   linetype = "dotted", linewidth = 0.8) +
        geom_vline(xintercept = 0.354, colour = "#ff5f57",
                   linetype = "dotted", linewidth = 0.8) +
        scale_fill_manual(values = REL_COLORS, name = "Relationship") +
        scale_y_continuous(labels = comma) +
        labs(
            title    = "CRISP: Step 7 Relatedness, KING Kinship Distribution",
            subtitle = sprintf(
                "%d pairs above KING >= %.4f",
                nrow(df), king_cutoff),
            x       = "Kinship coefficient",
            y       = "Number of pairs",
            caption = CAPTION
        ) +
        crisp_theme

    out1 <- file.path(output_dir, "Relatedness_KING_kinship.pdf")
    pdf(out1, width = 10, height = 6)
    print(p1)
    dev.off()
    cat(sprintf("[PLOT] KING kinship histogram saved: %s\n", out1))

    # Plot 2: Kinship vs IBS0 if available
    if ("IBS0" %in% colnames(df)) {

        p2 <- ggplot(df, aes(x = IBS0, y = KINSHIP, colour = RELATIONSHIP)) +
            geom_point(alpha = 0.65, size = 1.8) +
            geom_hline(yintercept = c(0.0884, 0.177, 0.354),
                       linetype = "dotted", colour = "#888", linewidth = 0.6) +
            scale_colour_manual(values = REL_COLORS, name = "Relationship") +
            scale_y_continuous(labels = number_format(accuracy = 0.001)) +
            labs(
                title    = "CRISP: Step 7 Relatedness, KING Kinship vs IBS0",
                subtitle = sprintf(
                    "%d related pairs  |  KING cutoff: %.4f",
                    nrow(df), king_cutoff),
                x       = "IBS0 (proportion of loci sharing 0 alleles)",
                y       = "Kinship coefficient",
                caption = CAPTION
            ) +
            crisp_theme

        out2 <- file.path(output_dir, "Relatedness_KING_scatter.pdf")
        pdf(out2, width = 10, height = 7)
        print(p2)
        # label duplicates and MZ twins only
        label_pairs(df[df$RELATIONSHIP == "Duplicate / MZ twin"],
                    "Dup/MZ", "#ff5f57", "IBS0", "KINSHIP")
        dev.off()
        cat(sprintf("[PLOT] KING scatter saved: %s\n", out2))
    }

} else {
    cat(sprintf("[PLOT] ERROR: Unknown method '%s'. Use IBD or KING.\n", method))
    quit(status = 1)
}

cat("\n[PLOT] Relatedness plotting complete.\n")
