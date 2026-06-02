#!/usr/bin/env Rscript
##########################################################################
### CRISP - Comprehensive Robust Integrated SNP Processing
### Step 6: Heterozygosity and Homozygosity, Plotting Script
### Version: 0.2.0
### Developed by Igor Pupko
### https://github.com/ipupko/CRISP
### Part of the Compass Genomics suite
##########################################################################
# Performs heterozygosity and homozygosity QC using PLINK output.
# Detects outliers using Z-score thresholds on autosomal F-statistics
# and visualises total runs of homozygosity (ROH) against Z-scores.
#
# Inputs:
#   .sexcheck file   : X chromosome F-statistics (sex check step)
#   .het file        : autosomal heterozygosity
#   .hom.indiv file  : runs of homozygosity per individual
#
# Outputs:
#   Homozygosity_Runs_vs_Zscore.pdf    ROH vs Z-score scatter plot
#   report_homozygosity.txt            Detailed text report
#   outliers_high_het.txt              High heterozygosity exclusions
#   outliers_excess_homo.txt           Excess homozygosity exclusions
#   outliers_combined.txt              Combined exclusion list
#
# Usage:
#   Rscript plot_homozygosity.R <output_dir> <sexcheck_file>
#                               <het_file> <hom_indiv_file>
#                               <z_high> <z_low>
##########################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    cat("Usage: Rscript plot_homozygosity.R",
        "<output_dir> <sexcheck_file> <het_file>",
        "<hom_indiv_file> <z_high> <z_low>\n")
    quit(status = 1)
}

output_dir      <- args[1]
sexcheck_file   <- args[2]
het_file        <- args[3]
hom_indiv_file  <- args[4]
z_high          <- as.numeric(args[5])
z_low           <- as.numeric(args[6])

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[HOMO] Parameters:\n"))
cat(sprintf("  Z-score upper threshold : %.1f (excess homozygosity)\n", z_high))
cat(sprintf("  Z-score lower threshold : %.1f (high heterozygosity)\n", z_low))
cat("\n")

##########################################################################
### LOAD DATA
##########################################################################

cat(sprintf("[HOMO] Loading sexcheck file : %s\n", sexcheck_file))
x_homozygosity <- fread(sexcheck_file, header = TRUE)

cat(sprintf("[HOMO] Loading het file      : %s\n", het_file))
autosome_het <- fread(het_file, header = TRUE)

cat(sprintf("[HOMO] Loading hom.indiv    : %s\n", hom_indiv_file))
hom_indiv <- fread(hom_indiv_file, header = TRUE)

cat(sprintf("[HOMO] Samples loaded: %d\n\n", nrow(autosome_het)))

##########################################################################
### PREPARE DATA
##########################################################################

# Create sex index from sexcheck for merging
sex_index <- x_homozygosity[, c("IID", "PEDSEX"), with = FALSE]

# Merge sex into autosomal het
autosome_het <- merge(autosome_het, sex_index,
                      by = "IID", all = TRUE)

# Add suffixes to avoid column name clashes
colnames(autosome_het) <- paste0(colnames(autosome_het), "_auto")
colnames(x_homozygosity) <- paste0(colnames(x_homozygosity), "_XChr")

##########################################################################
### CALCULATE Z-SCORES
##########################################################################

# X chromosome Z-score
x_homozygosity$Z_Score_XChr <- (
    x_homozygosity$F_XChr - mean(x_homozygosity$F_XChr, na.rm = TRUE)
) / sd(x_homozygosity$F_XChr, na.rm = TRUE)

# Autosomal Z-score
autosome_het$Z_Score_auto <- (
    autosome_het$F_auto - mean(autosome_het$F_auto, na.rm = TRUE)
) / sd(autosome_het$F_auto, na.rm = TRUE)

cat(sprintf("[HOMO] Autosomal F-stat mean : %.4f\n",
            mean(autosome_het$F_auto, na.rm = TRUE)))
cat(sprintf("[HOMO] Autosomal F-stat SD   : %.4f\n",
            sd(autosome_het$F_auto, na.rm = TRUE)))
cat("\n")

##########################################################################
### MERGE WITH ROH DATA
##########################################################################

autosome_het_roh <- merge(autosome_het, hom_indiv,
                           by.x = "IID_auto",
                           by.y = "IID",
                           all  = TRUE)

##########################################################################
### CLASSIFY OUTLIERS
##########################################################################

autosome_het_roh$Category <- ifelse(
    autosome_het_roh$Z_Score_auto > z_high, "Excess Homozygosity",
    ifelse(
        autosome_het_roh$Z_Score_auto < z_low, "High Heterozygosity",
        "Normal"
    )
)

n_excess_homo <- sum(autosome_het_roh$Category == "Excess Homozygosity",
                     na.rm = TRUE)
n_high_het    <- sum(autosome_het_roh$Category == "High Heterozygosity",
                     na.rm = TRUE)
n_normal      <- sum(autosome_het_roh$Category == "Normal", na.rm = TRUE)

cat(sprintf("[HOMO] Excess homozygosity (Z > %.1f)  : %d\n",
            z_high, n_excess_homo))
cat(sprintf("[HOMO] High heterozygosity (Z < %.1f)  : %d\n",
            z_low, n_high_het))
cat(sprintf("[HOMO] Normal                           : %d\n\n",
            n_normal))

##########################################################################
### PLOT: ROH KB VS HOMOZYGOSITY Z-SCORE
##########################################################################

cat("[HOMO] Generating plots...\n")

pdf_file <- file.path(output_dir, "Homozygosity_Runs_vs_Zscore.pdf")
pdf(pdf_file, width = 10, height = 7)

category_colors <- c(
    "Excess Homozygosity" = "#ff5f57",
    "High Heterozygosity" = "#febc2e",
    "Normal"              = "#1D9E75"
)

p <- ggplot(autosome_het_roh,
            aes(x = KB, y = Z_Score_auto, colour = Category)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = 0,      linetype = "dashed",
               colour = "black",    linewidth = 0.5) +
    geom_hline(yintercept = z_high, linetype = "dashed",
               colour = "#ff5f57",  linewidth = 0.8) +
    geom_hline(yintercept = z_low,  linetype = "dashed",
               colour = "#febc2e",  linewidth = 0.8) +
    annotate("text", x = max(autosome_het_roh$KB, na.rm = TRUE) * 0.05,
             y = z_high + 0.15,
             label = paste0("Z = ", z_high, " (excess homozygosity)"),
             colour = "#ff5f57", size = 3, hjust = 0) +
    annotate("text", x = max(autosome_het_roh$KB, na.rm = TRUE) * 0.05,
             y = z_low - 0.15,
             label = paste0("Z = ", z_low, " (high heterozygosity)"),
             colour = "#febc2e", size = 3, hjust = 0) +
    scale_colour_manual(values = category_colors,
                        name   = "Category") +
    scale_x_continuous(labels = label_comma()) +
    labs(
        title    = "CRISP: Step 6: Total Homozygosity Runs vs Homozygosity Z-score",
        subtitle = sprintf(
            "Z > %.1f: %d excess homozygosity  |  Z < %.1f: %d high heterozygosity  |  %d normal",
            z_high, n_excess_homo, z_low, n_high_het, n_normal
        ),
        x        = "Total runs of homozygosity (KB)",
        y        = "Homozygosity Z-score",
        caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "#555555", size = 9),
        axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        plot.caption     = element_text(colour = "#888888", size = 8)
    )

print(p)
dev.off()

cat(sprintf("[HOMO] Plot saved: %s\n\n", pdf_file))

##########################################################################
### WRITE REPORT
##########################################################################

report_file <- file.path(output_dir, "report_homozygosity.txt")
sink(report_file)

cat("================================================================\n")
cat("  CRISP: STEP 6 HETEROZYGOSITY AND HOMOZYGOSITY REPORT\n")
cat("  Comprehensive Robust Integrated SNP Processing\n")
cat("================================================================\n\n")

cat(sprintf("  Date              : %s\n", Sys.time()))
cat(sprintf("  Z-score upper     : %.1f (excess homozygosity)\n", z_high))
cat(sprintf("  Z-score lower     : %.1f (high heterozygosity)\n\n", z_low))

cat("----------------------------------------------------------------\n")
cat("AUTOSOMAL HETEROZYGOSITY\n")
cat("----------------------------------------------------------------\n")
cat(sprintf("  Samples total          : %d\n", nrow(autosome_het_roh)))
cat(sprintf("  F-stat mean            : %.4f\n",
            mean(autosome_het_roh$F_auto, na.rm = TRUE)))
cat(sprintf("  F-stat SD              : %.4f\n",
            sd(autosome_het_roh$F_auto, na.rm = TRUE)))
cat(sprintf("  Excess homozygosity    : %d (Z > %.1f)\n",
            n_excess_homo, z_high))
cat(sprintf("  High heterozygosity    : %d (Z < %.1f)\n\n",
            n_high_het, z_low))

cat("----------------------------------------------------------------\n")
cat("EXCESS HOMOZYGOSITY SAMPLES\n")
cat("----------------------------------------------------------------\n")
outliers_excess <- autosome_het_roh[
    autosome_het_roh$Category == "Excess Homozygosity", ]
if (nrow(outliers_excess) > 0) {
    print(as.data.frame(outliers_excess[,
        c("IID_auto", "F_auto", "Z_Score_auto", "KB")]))
} else {
    cat("  None detected.\n")
}

cat("\n----------------------------------------------------------------\n")
cat("HIGH HETEROZYGOSITY SAMPLES\n")
cat("----------------------------------------------------------------\n")
outliers_het <- autosome_het_roh[
    autosome_het_roh$Category == "High Heterozygosity", ]
if (nrow(outliers_het) > 0) {
    print(as.data.frame(outliers_het[,
        c("IID_auto", "F_auto", "Z_Score_auto", "KB")]))
} else {
    cat("  None detected.\n")
}

cat("\n================================================================\n")
cat("  END OF REPORT\n")
cat("================================================================\n")
sink()

cat(sprintf("[HOMO] Report written: %s\n\n", report_file))

##########################################################################
### WRITE EXCLUSION LISTS
##########################################################################

# High heterozygosity (your original naming convention preserved)
het_excl_file <- file.path(output_dir, "outliers_high_het.txt")
if (nrow(outliers_het) > 0) {
    write.table(outliers_het[, c("FID_auto", "IID_auto")],
                het_excl_file,
                quote     = FALSE,
                row.names = FALSE,
                sep       = "\t",
                col.names = FALSE)
} else {
    file.create(het_excl_file)
}
cat(sprintf("[HOMO] High het exclusion list     : %s (%d samples)\n",
            het_excl_file, nrow(outliers_het)))

# Excess homozygosity
homo_excl_file <- file.path(output_dir, "outliers_excess_homo.txt")
if (nrow(outliers_excess) > 0) {
    write.table(outliers_excess[, c("FID_auto", "IID_auto")],
                homo_excl_file,
                quote     = FALSE,
                row.names = FALSE,
                sep       = "\t",
                col.names = FALSE)
} else {
    file.create(homo_excl_file)
}
cat(sprintf("[HOMO] Excess homo exclusion list  : %s (%d samples)\n",
            homo_excl_file, nrow(outliers_excess)))

# Combined exclusion list
combined_file <- file.path(output_dir, "outliers_combined.txt")
all_outliers  <- unique(rbind(
    outliers_het[,    c("FID_auto", "IID_auto")],
    outliers_excess[, c("FID_auto", "IID_auto")]
))
write.table(all_outliers,
            combined_file,
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            col.names = FALSE)
cat(sprintf("[HOMO] Combined exclusion list     : %s (%d samples)\n\n",
            combined_file, nrow(all_outliers)))

cat("[HOMO] Homozygosity check complete.\n")
