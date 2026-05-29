#!/usr/bin/env Rscript
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### Step 5: Sex Check and Aneuploidy Detection, Plotting Script
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Performs sex verification and chromosomal aneuploidy detection
# using PLINK F-statistics and Y chromosome counts.
#
# Detects:
#   Sex mismatches (both directions)
#   Turner syndrome (X0)       : high F, low Y in reported females
#   Klinefelter syndrome (XXY) : low F, high Y in reported males
#   Triple-X syndrome (XXX)    : extreme X heterozygosity in females
#
# Usage:
#   Rscript plot_sexcheck.R <input_dir> <output_dir> <f_stat_y_count_file>
#                           <f_female_max> <f_male_min> <f_turner>
#                           <f_klinefelter> <f_xxx> <y_use_mean> <y_manual>
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
    cat("Usage: Rscript plot_sexcheck.R <input_dir> <output_dir>",
        "<f_stat_y_count_file> <f_female_max> <f_male_min>",
        "<f_turner> <f_klinefelter> <f_xxx> <y_use_mean> <y_manual>\n")
    quit(status = 1)
}

input_dir       <- args[1]
output_dir      <- args[2]
f_stat_y_file   <- args[3]
f_female_max    <- as.numeric(args[4])
f_male_min      <- as.numeric(args[5])
f_turner        <- as.numeric(args[6])
f_klinefelter   <- as.numeric(args[7])
f_xxx           <- as.numeric(args[8])
y_use_mean      <- toupper(args[9])
y_manual        <- as.numeric(args[10])
label_field     <- if (length(args) >= 11) args[11] else "IID"
show_labels     <- if (length(args) >= 12) toupper(args[12]) == "YES" else TRUE

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[SEXCHECK] Parameters:\n"))
cat(sprintf("  F female max    : %.3f\n", f_female_max))
cat(sprintf("  F male min      : %.3f\n", f_male_min))
cat(sprintf("  F Turner        : %.3f\n", f_turner))
cat(sprintf("  F Klinefelter   : %.3f\n", f_klinefelter))
cat(sprintf("  F XXX threshold : %.3f\n", f_xxx))
cat(sprintf("  Y use mean      : %s\n",   y_use_mean))
if (y_use_mean != "YES") {
    cat(sprintf("  Y manual cutoff : %.1f\n", y_manual))
}
cat("\n")

##########################################################################
# LOAD DATA
##########################################################################

cat(sprintf("[SEXCHECK] Loading: %s\n", f_stat_y_file))
df <- fread(f_stat_y_file, header = TRUE)

cat(sprintf("[SEXCHECK] Samples loaded: %d\n\n", nrow(df)))

# Determine Y count threshold
if (y_use_mean == "YES") {
    y_threshold <- mean(df$YCOUNT, na.rm = TRUE)
    cat(sprintf("[SEXCHECK] Y count threshold: mean = %.2f\n\n", y_threshold))
} else {
    y_threshold <- y_manual
    cat(sprintf("[SEXCHECK] Y count threshold: manual = %.2f\n\n", y_threshold))
}

##########################################################################
# SUBSETS
##########################################################################

df_mal   <- df[df$PEDSEX == 1, ]
df_femal <- df[df$PEDSEX == 2, ]

cat(sprintf("[SEXCHECK] Reported males   : %d\n", nrow(df_mal)))
cat(sprintf("[SEXCHECK] Reported females : %d\n\n", nrow(df_femal)))

##########################################################################
# PLOTS
##########################################################################

cat("[SEXCHECK] Generating sex check plots...\n")

pdf_file <- file.path(output_dir, "Sex_check.pdf")
pdf(pdf_file)

# Overall
plot(df$F, df$YCOUNT,
     col  = df$PEDSEX,
     xlab = "F Statistic",
     ylab = "Y Count",
     main = "Overall Sex Check",
     pch  = 19, cex = 0.6)
abline(v = f_female_max, col = "blue",   lty = 2)
abline(v = f_male_min,   col = "red",    lty = 2)
abline(h = y_threshold,  col = "grey50", lty = 2)
legend("topright",
       legend = c("Male", "Female",
                  paste0("F=", f_female_max),
                  paste0("F=", f_male_min),
                  paste0("Y mean=", round(y_threshold, 1))),
       col    = c("red", "black", "blue", "red", "grey50"),
       lty    = c(NA, NA, 2, 2, 2),
       pch    = c(19, 19, NA, NA, NA),
       cex    = 0.8)

# label anomalies with boundary-aware, collision-resolving arrows
# called after each plot is drawn, before dev.off()
# placed_labels is an environment used to track label positions across calls
label_samples <- function(df_flag, type_label, col,
                          label_field  = "IID",
                          show_labels  = TRUE,
                          placed_list  = NULL) {
    if (nrow(df_flag) == 0) return(invisible(NULL))

    # draw star overlay
    points(df_flag$F, df_flag$YCOUNT,
           col = col, pch = 8, cex = 1.6, lwd = 2.0)

    if (!show_labels) return(invisible(NULL))

    usr   <- par("usr")         # c(xmin, xmax, ymin, ymax)
    x_rng <- usr[2] - usr[1]
    y_rng <- usr[4] - usr[3]

    for (i in seq_len(nrow(df_flag))) {
        lbl <- switch(label_field,
            "FID_IID" = paste0(df_flag$FID[i], "_", df_flag$IID[i],
                               " (", type_label, ")"),
            "FID"     = paste0(df_flag$FID[i], " (", type_label, ")"),
                        paste0(df_flag$IID[i], " (", type_label, ")")
        )

        x0 <- df_flag$F[i]
        y0 <- df_flag$YCOUNT[i]

        # candidate offsets in priority order
        candidates <- list(
            c( x_rng*0.10,  y_rng*0.10),
            c(-x_rng*0.10,  y_rng*0.10),
            c( x_rng*0.10, -y_rng*0.10),
            c(-x_rng*0.10, -y_rng*0.10),
            c( x_rng*0.15,  y_rng*0.05),
            c(-x_rng*0.15,  y_rng*0.05),
            c( x_rng*0.05,  y_rng*0.15),
            c( x_rng*0.05, -y_rng*0.15)
        )

        lbl_w <- strwidth(lbl,  cex = 0.62)
        lbl_h <- strheight(lbl, cex = 0.62)
        margin <- x_rng * 0.01

        best_x <- NA
        best_y <- NA

        for (cand in candidates) {
            cx <- x0 + cand[1]
            cy <- y0 + cand[2]

            # boundary check -- keep label inside plot area
            if (cx - lbl_w/2 < usr[1] + margin) next
            if (cx + lbl_w/2 > usr[2] - margin) next
            if (cy - lbl_h/2 < usr[3] + margin) next
            if (cy + lbl_h/2 > usr[4] - margin) next

            # collision check against already placed labels
            overlap <- FALSE
            if (!is.null(placed_list) && length(placed_list$x) > 0) {
                for (j in seq_along(placed_list$x)) {
                    if (abs(cx - placed_list$x[j]) < lbl_w * 0.9 &&
                        abs(cy - placed_list$y[j]) < lbl_h * 1.5) {
                        overlap <- TRUE
                        break
                    }
                }
            }
            if (!overlap) { best_x <- cx; best_y <- cy; break }
        }

        # fallback if no clean position found
        if (is.na(best_x)) {
            best_x <- x0 + candidates[[1]][1]
            best_y <- y0 + candidates[[1]][2]
            best_x <- max(usr[1] + margin + lbl_w/2,
                          min(usr[2] - margin - lbl_w/2, best_x))
            best_y <- max(usr[3] + margin + lbl_h,
                          min(usr[4] - margin, best_y))
        }

        # record placement
        if (!is.null(placed_list)) {
            placed_list$x <- c(placed_list$x, best_x)
            placed_list$y <- c(placed_list$y, best_y)
        }

        # draw curved arrow
        arrows(x0, y0, best_x, best_y,
               col    = col,
               length = 0.07,
               lwd    = 1.0,
               angle  = 20)

        # draw label
        text(best_x, best_y,
             labels = lbl,
             col    = col,
             cex    = 0.62,
             font   = 2,
             adj    = c(0.5, 0.5))
    }
    return(invisible(placed_list))
}

label_samples(sex_mismatch_females, "mismatch",       "#ff8c00", label_field, show_labels)
label_samples(sex_mismatch_males,   "mismatch",       "#9400d3", label_field, show_labels)
label_samples(X0_Turner,            "Turner X0",      "#0000cc", label_field, show_labels)
label_samples(XXY_Klinefelter,      "Klinefelter XXY","#cc00cc", label_field, show_labels)
label_samples(XXX_TripleX,          "Triple-X XXX",   "#00aeae", label_field, show_labels)

# Males
plot(df_mal$F, df_mal$YCOUNT,
     col  = "red",
     xlab = "F Statistic",
     ylab = "Y Count",
     main = "Male Sex Check",
     pch  = 19, cex = 0.6)
abline(v = f_klinefelter, col = "darkred", lty = 2)
abline(h = y_threshold,   col = "grey50",  lty = 2)
legend("topright",
       legend = c("Male",
                  paste0("F=", f_klinefelter),
                  paste0("Y mean=", round(y_threshold, 1))),
       col    = c("red", "darkred", "grey50"),
       lty    = c(NA, 2, 2),
       pch    = c(19, NA, NA),
       cex    = 0.8)

label_samples(sex_mismatch_males, "mismatch",       "#9400d3", label_field, show_labels)
label_samples(XXY_Klinefelter,    "Klinefelter XXY","#cc00cc", label_field, show_labels)

# Females
plot(df_femal$F, df_femal$YCOUNT,
     col  = "black",
     xlab = "F Statistic",
     ylab = "Y Count",
     main = "Female Sex Check",
     pch  = 19, cex = 0.6)
abline(v = f_turner, col = "purple", lty = 2)
abline(v = f_xxx,    col = "blue",   lty = 2)
abline(h = y_threshold, col = "grey50", lty = 2)
legend("topright",
       legend = c("Female",
                  paste0("F=", f_turner, " (Turner)"),
                  paste0("F=", f_xxx,    " (XXX)"),
                  paste0("Y mean=", round(y_threshold, 1))),
       col    = c("black", "purple", "blue", "grey50"),
       lty    = c(NA, 2, 2, 2),
       pch    = c(19, NA, NA, NA),
       cex    = 0.8)

label_samples(sex_mismatch_females, "mismatch",   "#ff8c00", label_field, show_labels)
label_samples(X0_Turner,            "Turner X0",  "#0000cc", label_field, show_labels)
label_samples(XXX_TripleX,          "Triple-X XXX","#00aeae", label_field, show_labels)

dev.off()
cat(sprintf("[SEXCHECK] Plots saved: %s\n\n", pdf_file))

##########################################################################
# DETECTION
##########################################################################

cat("[SEXCHECK] Running detection...\n\n")

# FEMALE SEX MISMATCH
# sex mismatch females: reported females with high F and high Y, likely males
sex_mismatch_females <- subset(df_femal,
    F > f_male_min & YCOUNT > y_threshold)

# MALE SEX MISMATCH
# sex mismatch males: reported males with low F and low Y, likely females
sex_mismatch_males <- subset(df_mal,
    F < f_female_max & YCOUNT < y_threshold)

# TURNER SYNDROME (X0)
# Reported females with high F but low Y
# High F in a female suggests only one X chromosome
X0_Turner <- subset(df_femal,
    F > f_turner & YCOUNT < y_threshold)

# KLINEFELTER SYNDROME (XXY)
# Reported males with low F but high Y
# Low F in a male suggests extra X chromosome
XXY_Klinefelter <- subset(df_mal,
    F < f_klinefelter & YCOUNT > y_threshold)

# TRIPLE-X SYNDROME (XXX)
# Reported females with extreme X heterozygosity and low Y
# Very negative F suggests extra X chromosome
XXX_TripleX <- subset(df_femal,
    F < f_xxx & YCOUNT < y_threshold)

##########################################################################
# REPORT
##########################################################################

report_file <- file.path(output_dir, "report_sex.mismatch_aneuploidies.txt")
sink(report_file)

cat("================================================================\n")
cat("  CRISP: STEP 5 SEX CHECK AND ANEUPLOIDY REPORT\n")
cat("  Comprehensive Robust Integrated SNP Processing\n")
cat("================================================================\n\n")

cat(sprintf("  Date           : %s\n", Sys.time()))
cat(sprintf("  F female max   : %.3f\n", f_female_max))
cat(sprintf("  F male min     : %.3f\n", f_male_min))
cat(sprintf("  Y threshold    : %.2f (%s)\n\n",
            y_threshold,
            ifelse(y_use_mean == "YES", "mean", "manual")))

cat("----------------------------------------------------------------\n")
cat("SAMPLE COUNTS\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  Reported males   : %d\n", nrow(df_mal)))
cat(sprintf("  Reported females : %d\n\n", nrow(df_femal)))

cat("----------------------------------------------------------------\n")
cat("SEX MISMATCHES\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  Sex mismatches in females : %d\n", nrow(sex_mismatch_females)))
if (nrow(sex_mismatch_females) > 0) {
    cat("\n  Details:\n")
    print(as.data.frame(sex_mismatch_females))
    cat("\n")
}

cat(sprintf("  Sex mismatches in males   : %d\n", nrow(sex_mismatch_males)))
if (nrow(sex_mismatch_males) > 0) {
    cat("\n  Details:\n")
    print(as.data.frame(sex_mismatch_males))
    cat("\n")
}

cat("----------------------------------------------------------------\n")
cat("CHROMOSOMAL ANEUPLOIDIES\n")
cat("----------------------------------------------------------------\n")

cat(sprintf("  Turner syndrome (X0)      : %d\n", nrow(X0_Turner)))
if (nrow(X0_Turner) > 0) {
    cat("\n  Details:\n")
    print(as.data.frame(X0_Turner))
    cat("\n")
}

cat(sprintf("  Klinefelter syndrome (XXY): %d\n", nrow(XXY_Klinefelter)))
if (nrow(XXY_Klinefelter) > 0) {
    cat("\n  Details:\n")
    print(as.data.frame(XXY_Klinefelter))
    cat("\n")
}

cat(sprintf("  Triple-X syndrome (XXX)   : %d\n", nrow(XXX_TripleX)))
if (nrow(XXX_TripleX) > 0) {
    cat("\n  Details:\n")
    print(as.data.frame(XXX_TripleX))
    cat("\n")
}

# Combined aneuploidy summary
df_aneuploidies <- unique(rbind(X0_Turner, XXY_Klinefelter, XXX_TripleX))
cat("----------------------------------------------------------------\n")
cat(sprintf("  Total aneuploidies flagged: %d\n", nrow(df_aneuploidies)))
if (nrow(df_aneuploidies) > 0) {
    cat("\n  All aneuploidy samples:\n")
    print(as.data.frame(df_aneuploidies))
}

cat("\n================================================================\n")
cat("  END OF REPORT\n")
cat("================================================================\n")

sink()
cat(sprintf("[SEXCHECK] Report written: %s\n\n", report_file))

##########################################################################
# EXCLUSION LISTS
##########################################################################

# Sex mismatch exclusion list (both directions combined)
df_sex_mismatch <- unique(rbind(sex_mismatch_females, sex_mismatch_males))
mismatch_file   <- file.path(output_dir, "id_list_sex_mismatch.txt")
write.table(df_sex_mismatch[, c(1, 2)],
            mismatch_file,
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            col.names = FALSE)
cat(sprintf("[SEXCHECK] Sex mismatch exclusion list: %s (%d samples)\n",
            mismatch_file, nrow(df_sex_mismatch)))

# Aneuploidy exclusion list
aneuploidy_file <- file.path(output_dir, "id_list_aneuploidies.txt")
write.table(df_aneuploidies[, c(1, 2)],
            aneuploidy_file,
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            col.names = FALSE)
cat(sprintf("[SEXCHECK] Aneuploidy exclusion list  : %s (%d samples)\n",
            aneuploidy_file, nrow(df_aneuploidies)))

cat("\n[SEXCHECK] Sex check complete.\n")
