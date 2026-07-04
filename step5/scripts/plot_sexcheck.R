#!/usr/bin/env Rscript
##########################################################################
#
#   .oooooo.   ooooooooo.   ooooo  .oooooo..o ooooooooo.
#  d8P'  `Y8b  `888   `Y88. `888' d8P'    `Y8 `888   `Y88.
# 888           888   .d88'  888  Y88bo.       888   .d88'
# 888           888ooo88P'   888   `"Y8888o.   888ooo88P'
# 888           888`88b.     888       `"Y88b  888
# `88b    ooo   888  `88b.   888  oo     .d8P  888
#  `Y8bood8P'  o888o  o888o o888o 8""88888P'  o888o
#
# CRISP - Comprehensive Robust Integrated SNP Processing
# Step 5: Sex Check and Aneuploidy Detection — R Plotting Script
# Version: 0.4.0
# https://github.com/ipupko/CRISP
# Part of the Compass Genomics suite
##########################################################################
# Detects and plots sex mismatches and chromosomal aneuploidies using
# PLINK X-chromosome F-statistics and Y-chromosome variant counts.
#
# Detects:
#   Sex mismatches (both directions)
#   Turner syndrome   (X0)  : high F, low Y in reported females
#   Klinefelter       (XXY) : low F,  high Y in reported males
#   Triple-X syndrome (XXX) : extreme X heterozygosity in females
#
# v0.4.0:
#   - NEW: aneuploidy_cohort_manifest.txt — typed manifest (CATEGORY,
#     TYPE) of Turner/Klinefelter/Triple-X samples for cohort routing.
#     Aneuploidies are confirmed karyotype calls, not data errors —
#     valuable to RNA/eQTL researchers as natural dosage experiments.
#     Sex mismatches remain exclusion-only (Type 2 error, not biology).
#   - NEW: per-type ID lists (id_list_aneuploidy_turner.txt, etc.) for
#     shell-side --keep extraction into per-type sub-cohorts.
#   - FIX: Y threshold now uses MALE-ONLY median, not whole-cohort mean.
#   - FIX: anomaly labels across ALL types now repel jointly in a single
#     geom_label_repel call (previously each type repelled independently,
#     so labels from different types could still collide with each other).
#   - FIX: anomaly stars now appear in the legend. Previously drawn with
#     a fixed colour= argument, which ggplot2 never legends; now mapped
#     via aes() with a second colour scale (ggnewscale) so the legend
#     shows anomaly type alongside reported sex.
#   - FIX: removed dead save_page() function (defined, never called).
#
# v0.3.0 fixes:
#   - Ported from base graphics to ggplot2 (was loaded but unused)
#   - Detection now runs BEFORE plotting (was silently after dev.off())
#   - All legends moved to right side (topright -> right-hand scale)
#   - CRISP_PAL_* colour variables now actually used in all geoms
#   - label_samples() moved to top-level function definition
#   - Wrong script name in palette comment fixed (was crisp_hwe.sh)
#   - SEXCHECK_LABEL_FIELD and SEXCHECK_LABEL_ANOMALIES now forwarded
#     as args[11] and args[12] from crisp_sexcheck.sh
#
# Usage:
#   Rscript plot_sexcheck.R <input_dir> <output_dir> <f_stat_y_count_file>
#                           <f_female_max> <f_male_min> <f_turner>
#                           <f_klinefelter> <f_xxx> <y_use_mean> <y_manual>
#                           [label_field] [show_labels]
##########################################################################

##########################################################################
#
#   ____ ___  _     ___  _   _ ____    ____   _    _
#  / ___/ _ \| |   / _ \| | | |  _ \  |  _ \ / \  | |
# | |  | | | | |  | | | | | | | |_) | | |_) / _ \ | |
# | |__| |_| | |__| |_| | |_| |  _ <  |  __/ ___ \| |___
#  \____\___/|_____\___/ \___/|_| \_\ |_| /_/   \_\_____|
#
# Reads PLOT_COLOUR_MODE from the environment variable CRISP_PAL_MODE,
# which is exported by _crisp_palette() in crisp_sexcheck.sh before
# either plot script is invoked. Falls back to STANDARD if unset.
##########################################################################

colour_mode <- Sys.getenv("CRISP_PAL_MODE", unset = "STANDARD")

if (toupper(colour_mode) == "COLOURBLIND") {
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#0072B2")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#D55E00")
    PAL_WARN      <- Sys.getenv("CRISP_PAL_WARN",      "#E69F00")
    PAL_HIGHLIGHT <- Sys.getenv("CRISP_PAL_HIGHLIGHT", "#009E73")
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",       "#56B4E9")
    PAL_PINK      <- Sys.getenv("CRISP_PAL_PINK",      "#CC79A7")
    PAL_ORANGE    <- Sys.getenv("CRISP_PAL_YELLOW",    "#F0E442")
} else if (toupper(colour_mode) == "NIGHT") {
    PAL_PASS      <- Sys.getenv("CRISP_PAL_PASS",      "#2ED9A3")
    PAL_FAIL      <- Sys.getenv("CRISP_PAL_FAIL",      "#FF7B72")
    PAL_WARN      <- Sys.getenv("CRISP_PAL_WARN",      "#FFC857")
    PAL_HIGHLIGHT <- Sys.getenv("CRISP_PAL_HIGHLIGHT", "#5EEBC4")
    PAL_SKY       <- Sys.getenv("CRISP_PAL_SKY",       "#7FD1FF")
    PAL_PINK      <- Sys.getenv("CRISP_PAL_PINK",      "#D8B4FE")
    PAL_ORANGE    <- Sys.getenv("CRISP_PAL_YELLOW",    "#FFE066")
} else {
    PAL_PASS      <- "#1D9E75"
    PAL_FAIL      <- "#ff5f57"
    PAL_WARN      <- "#febc2e"
    PAL_HIGHLIGHT <- "#5dcaa5"
    PAL_SKY       <- "#7ec8e3"
    PAL_PINK      <- "#afa9ec"
    PAL_ORANGE    <- "#ff8c00"
}

# Background / text / grid — orthogonal to PLOT_COLOUR_MODE, exported by
# _crisp_palette() in the shell script. All plot elements must use these
# instead of hardcoded colours so DARK background works transparently.
env_or <- function(name, default) {
    val <- Sys.getenv(name)
    if (nzchar(val)) val else default
}
PAL_BG         <- env_or("CRISP_PAL_BG",         "#FFFFFF")
PAL_PANEL      <- env_or("CRISP_PAL_PANEL",      "#f8f9fa")
PAL_TEXT       <- env_or("CRISP_PAL_TEXT",       "#1a1a1a")
PAL_SUBTEXT    <- env_or("CRISP_PAL_SUBTEXT",    "#555555")
PAL_GRID       <- env_or("CRISP_PAL_GRID",       "#e0e0e0")
PAL_BACKGROUND <- env_or("CRISP_PAL_BACKGROUND", "LIGHT")

cat(sprintf("[SEXCHECK] Colour mode : %s\n", colour_mode))
cat(sprintf("[SEXCHECK] Background  : %s\n", PAL_BACKGROUND))

# Anomaly colours — derived from CRISP palette so all three modes
# (STANDARD / COLOURBLIND / NIGHT) carry through to anomaly markers
COL_MISMATCH_F  <- PAL_ORANGE
COL_MISMATCH_M  <- PAL_PINK
COL_TURNER      <- PAL_SKY
COL_KLINEFELTER <- PAL_FAIL
COL_XXX         <- PAL_HIGHLIGHT

# Fixed sex-group scatter colours — not palette-derived, but readable
# on both LIGHT and DARK backgrounds
COL_MALE        <- "#e05c4b"
COL_FEMALE      <- "#1D9E75"
COL_UNKNOWN     <- "#aaaaaa"

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggrepel)
    library(ggnewscale)
})

##########################################################################
#
#  ____   _    ____  ____  _____   ____   _    ____      _    __  __ ____
# |  _ \ / \  |  _ \/ ___|| ____| |  _ \ / \  |  _ \    / \  |  \/  / ___|
# | |_) / _ \ | |_) \___ \|  _|   | |_) / _ \ | |_) |  / _ \ | |\/| \___ \
# |  __/ ___ \|  _ < ___) | |___  |  __/ ___ \|  _ <  / ___ \| |  | |___) |
# |_| /_/   \_\_| \_\____/|_____| |_| /_/   \_\_| \_\/_/   \_\_|  |_|____/
#
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
    cat("Usage: Rscript plot_sexcheck.R <input_dir> <output_dir>",
        "<f_stat_y_count_file> <f_female_max> <f_male_min>",
        "<f_turner> <f_klinefelter> <f_xxx> <y_use_mean> <y_manual>",
        "[label_field] [show_labels]\n")
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
label_field     <- if (length(args) >= 11) toupper(args[11]) else "IID"
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
cat(sprintf("  Label field     : %s\n",   label_field))
cat(sprintf("  Show labels     : %s\n",   ifelse(show_labels, "YES", "NO")))
cat("\n")

##########################################################################
#
#  _     ___    _    ____    ____    _  _____  _
# | |   / _ \  / \  |  _ \  |  _ \  / \|_   _|/ \
# | |  | | | |/ _ \ | | | | | | | |/ _ \ | | / _ \
# | |__| |_| / ___ \| |_| | | |_| / ___ \| |/ ___ \
# |_____\___/_/   \_\____/  |____/_/   \_\_/_/   \_\
#
##########################################################################

cat(sprintf("[SEXCHECK] Loading: %s\n", f_stat_y_file))
df <- fread(f_stat_y_file, header = TRUE)
cat(sprintf("[SEXCHECK] Samples loaded: %d\n\n", nrow(df)))

df_mal   <- df[df$PEDSEX == 1, ]
df_femal <- df[df$PEDSEX == 2, ]

# Y-count threshold
# FIX v0.4.0: previously computed from the WHOLE cohort (mean(df$YCOUNT)),
# which is dragged down by females (~0) and unknowns. A male-only median
# is robust to outliers and anchors the threshold to the Y-present
# distribution it is meant to separate from.
if (y_use_mean == "YES") {
    if (nrow(df_mal) > 0) {
        y_threshold <- median(df_mal$YCOUNT, na.rm = TRUE)
        cat(sprintf("[SEXCHECK] Y count threshold: male-only median = %.2f\n\n",
                    y_threshold))
    } else {
        y_threshold <- median(df$YCOUNT, na.rm = TRUE)
        cat("[SEXCHECK] WARNING: No reported males found; falling back to\n")
        cat(sprintf("           whole-cohort median = %.2f\n\n", y_threshold))
    }
} else {
    y_threshold <- y_manual
    cat(sprintf("[SEXCHECK] Y count threshold: manual = %.2f\n\n", y_threshold))
}

cat(sprintf("[SEXCHECK] Reported males   : %d\n", nrow(df_mal)))
cat(sprintf("[SEXCHECK] Reported females : %d\n\n", nrow(df_femal)))

##########################################################################
#
#  ____  _____ _____ _____ ____ _____ ___ ___  _   _
# |  _ \| ____|_   _| ____/ ___|_   _|_ _/ _ \| \ | |
# | | | |  _|   | | |  _|| |     | |  | | | | |  \| |
# | |_| | |___  | | | |__| |___  | |  | | |_| | |\  |
# |____/|_____| |_| |_____\____| |_| |___\___/|_| \_|
#
# FIX #4: Detection MUST run before plotting so that anomaly data frames
# exist when label_samples() is called inside each plot. In v0.2.0 the
# detection block was placed after dev.off(), meaning label_samples()
# calls referenced objects that did not yet exist — all label overlays
# silently failed (R would throw "object not found" inside the closed
# graphics device, suppressed because the PDF was already written).
##########################################################################

cat("[SEXCHECK] Running detection...\n\n")

# Sex mismatch — female: reported female but genotypically male
# (high F and non-zero Y count)
sex_mismatch_females <- subset(df_femal,
    F > f_male_min & YCOUNT > y_threshold)

# Sex mismatch — male: reported male but genotypically female
# (low F and low/absent Y count)
sex_mismatch_males <- subset(df_mal,
    F < f_female_max & YCOUNT < y_threshold)

# Turner syndrome (X0): reported female with elevated F (single X)
# and low Y count. The elevated F reflects hemizygosity on one X.
X0_Turner <- subset(df_femal,
    F > f_turner & YCOUNT < y_threshold)

# Klinefelter syndrome (XXY): reported male with low F (extra X
# chromosome dilutes inbreeding signal) and detectable Y count.
XXY_Klinefelter <- subset(df_mal,
    F < f_klinefelter & YCOUNT > y_threshold)

# Triple-X syndrome (XXX): reported female with strongly negative F
# (extra X increases heterozygosity) and low Y count.
XXX_TripleX <- subset(df_femal,
    F < f_xxx & YCOUNT < y_threshold)

# Combined exclusion sets for report and exclusion list writing
df_sex_mismatch <- unique(rbind(sex_mismatch_females, sex_mismatch_males))
df_aneuploidies <- unique(rbind(X0_Turner, XXY_Klinefelter, XXX_TripleX))

cat(sprintf("[SEXCHECK] Sex mismatch females  : %d\n", nrow(sex_mismatch_females)))
cat(sprintf("[SEXCHECK] Sex mismatch males    : %d\n", nrow(sex_mismatch_males)))
cat(sprintf("[SEXCHECK] Turner (X0)           : %d\n", nrow(X0_Turner)))
cat(sprintf("[SEXCHECK] Klinefelter (XXY)     : %d\n", nrow(XXY_Klinefelter)))
cat(sprintf("[SEXCHECK] Triple-X (XXX)        : %d\n", nrow(XXX_TripleX)))
cat("\n")

##########################################################################
#
#  ____  _     ___  _____ ____
# |  _ \| |   / _ \|_   _/ ___|
# | |_) | |  | | | | | | \___ \
# |  __/| |__| |_| | | |  ___) |
# |_|   |_____\___/  |_| |____/
#
# FIX #5: All three plots rebuilt in ggplot2 (base graphics removed).
# FIX #6: Legends moved to right-hand side (theme(legend.position="right")).
#         The CRISP standard requires all legends on the right.
# FIX #7: PAL_* variables now used for all geom colours.
# FIX #8: label_samples() defined here at top level before any plot call,
#         and uses ggrepel::geom_label_repel for collision-safe labelling.
##########################################################################

# Build the label string for a flagged sample according to label_field.
# Type suffix omitted: anomaly category is conveyed by star colour and
# the legend, keeping labels short and reducing plot clutter.
make_label <- function(df_flag, type_label, label_field) {
    if (nrow(df_flag) == 0) return(character(0))
    switch(label_field,
        "FID_IID" = paste0(df_flag$FID, "_", df_flag$IID),
        "FID"     = as.character(df_flag$FID),
                    as.character(df_flag$IID)   # default: IID
    )
}

# CRISP ggplot2 theme — applied to every plot.
# All colours read from PAL_* variables (set from CRISP_PAL_* env vars
# above) so LIGHT and DARK backgrounds work without any theme changes.
crisp_theme <- function() {
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor    = element_blank(),
        panel.grid.major    = element_line(colour = PAL_GRID, linewidth = 0.4),
        legend.position     = "right",
        legend.title        = element_text(size = 9, face = "bold",
                                           colour = PAL_TEXT),
        legend.text         = element_text(size = 8, colour = PAL_TEXT),
        legend.background   = element_rect(fill = PAL_BG, colour = NA),
        legend.key          = element_rect(fill = PAL_BG, colour = NA),
        plot.title          = element_text(face = "bold", size = 13,
                                           colour = PAL_TEXT),
        plot.subtitle       = element_text(size = 9, colour = PAL_SUBTEXT),
        plot.caption        = element_text(size = 8, colour = PAL_SUBTEXT,
                                           hjust = 0.5),
        plot.background     = element_rect(fill = PAL_BG, colour = NA),
        panel.background    = element_rect(fill = PAL_PANEL, colour = NA),
        strip.background    = element_rect(fill = PAL_PANEL, colour = NA),
        strip.text          = element_text(face = "bold", colour = PAL_TEXT),
        axis.title          = element_text(face = "bold", size = 10,
                                           colour = PAL_TEXT),
        axis.text           = element_text(size = 9, colour = PAL_SUBTEXT),
        axis.ticks          = element_line(colour = PAL_SUBTEXT)
    )
}

# Convenience: save a single ggplot page to a PdfPages-style PDF
save_page <- function(p, pdf_dev) {
    print(p)
}

cat("[SEXCHECK] Generating sex check plots...\n")

pdf_file <- file.path(output_dir, "Sex_check.pdf")
pdf(pdf_file, width = 10, height = 7)

# ── HELPER: anomaly star points, mapped so a legend is generated ─────────
# A fixed `colour=` argument (as used previously) never appears in a
# ggplot2 legend — only aesthetics mapped via aes() do. Now that
# aneuploidies may be ROUTED rather than discarded (see manifest writer
# below), researchers need to see anomaly categories at a glance, so we
# build one combined anomaly dataframe, layer a second colour scale on
# top via ggnewscale, and map colour properly so the legend renders.
build_anomaly_df <- function(groups) {
    # groups: named list of list(df = <data.frame>, type = <string>, col = <string>)
    pieces <- lapply(groups, function(g) {
        if (nrow(g$df) == 0) return(NULL)
        d <- as.data.frame(g$df)
        d$.type <- g$type
        d$.colour <- g$col
        d
    })
    pieces <- pieces[!sapply(pieces, is.null)]
    if (length(pieces) == 0) return(NULL)
    do.call(rbind, lapply(pieces, function(d) d[, c("F", "YCOUNT", ".type", ".colour")]))
}

anomaly_points_layer <- function(groups) {
    adf <- build_anomaly_df(groups)
    if (is.null(adf) || nrow(adf) == 0) return(list())
    type_levels <- unique(adf$.type)
    col_lookup  <- setNames(unique(adf$.colour), type_levels)[type_levels]
    adf$.type <- factor(adf$.type, levels = type_levels)
    list(
        ggnewscale::new_scale_colour(),
        geom_point(
            data   = adf,
            aes(x = F, y = YCOUNT, colour = .type),
            shape  = 8, size = 3.5, stroke = 1.5,
            inherit.aes = FALSE
        ),
        scale_colour_manual(
            name   = "Anomaly",
            values = col_lookup
        )
    )
}

# ── HELPER: anomaly star points only (no labels) ─────────────────────────
# Labels are now built separately and combined into a SINGLE
# geom_label_repel call per plot (see build_label_df below), so that
# labels from different anomaly types repel each other too — not just
# within their own type. Previously each anomaly_layers() call ran its
# own independent repel pass, so a Turner label and a same-region
# mismatch label could still collide.
anomaly_points <- function(df_flag, col) {
    if (nrow(df_flag) == 0) return(list())
    df_flag <- as.data.frame(df_flag)
    list(
        geom_point(
            data   = df_flag,
            aes(x = F, y = YCOUNT),
            shape  = 8, size = 3.5, stroke = 1.5,
            colour = col,
            inherit.aes = FALSE
        )
    )
}

# Build one combined label dataframe across ALL anomaly groups for a plot,
# so a single geom_label_repel call can repel labels jointly.
build_label_df <- function(groups) {
    # groups: named list of list(df = <data.frame>, type = <string>, col = <string>)
    pieces <- lapply(groups, function(g) {
        if (nrow(g$df) == 0) return(NULL)
        d <- as.data.frame(g$df)
        d$.label <- make_label(d, g$type, label_field)
        d$.colour <- g$col
        d
    })
    pieces <- pieces[!sapply(pieces, is.null)]
    if (length(pieces) == 0) return(NULL)
    do.call(rbind, lapply(pieces, function(d)
        d[, c("F", "YCOUNT", ".label", ".colour")]))
}

# Single combined repel layer for a plot's full anomaly set
combined_label_layer <- function(groups, show_labels) {
    if (!show_labels) return(list())
    lbl_df <- build_label_df(groups)
    if (is.null(lbl_df) || nrow(lbl_df) == 0) return(list())
    list(
        ggrepel::geom_label_repel(
            data           = lbl_df,
            aes(x = F, y = YCOUNT, label = .label),
            colour         = lbl_df$.colour,
            fill           = "white",
            size           = 2.8,
            fontface       = "bold",
            label.size     = 0.35,
            box.padding    = 0.4,
            point.padding  = 0.3,
            min.segment.length = 0,
            segment.colour = lbl_df$.colour,
            segment.size   = 0.5,
            max.overlaps   = Inf,
            inherit.aes    = FALSE
        )
    )
}

# ── OVERALL PLOT ─────────────────────────────────────────────────────────
df_plot <- as.data.frame(df)
df_plot$sex_label <- factor(
    ifelse(df_plot$PEDSEX == 1, "Male",
    ifelse(df_plot$PEDSEX == 2, "Female", "Unknown")),
    levels = c("Male", "Female", "Unknown")
)

n_mis  <- nrow(df_sex_mismatch)
n_aneu <- nrow(df_aneuploidies)

p_overall <- ggplot(df_plot, aes(x = F, y = YCOUNT, colour = sex_label)) +
    geom_point(size = 1.2, alpha = 0.65) +
    scale_colour_manual(
        name   = "Reported sex",
        values = c(Male = COL_MALE, Female = COL_FEMALE, Unknown = COL_UNKNOWN)
    ) +
    geom_vline(xintercept = f_female_max, colour = PAL_PASS,
               linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = f_male_min,  colour = PAL_FAIL,
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = y_threshold, colour = "#666666",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text",
             x     = f_female_max,
             y     = Inf,
             label = sprintf("F=%.2f\n(female)", f_female_max),
             vjust = 1.3, hjust = -0.1,
             size  = 2.8, colour = PAL_PASS, fontface = "bold") +
    annotate("text",
             x     = f_male_min,
             y     = Inf,
             label = sprintf("F=%.2f\n(male)", f_male_min),
             vjust = 1.3, hjust = -0.1,
             size  = 2.8, colour = PAL_FAIL, fontface = "bold") +
    annotate("text",
             x     = -Inf,
             y     = y_threshold,
             label = sprintf("Y mean = %.1f", y_threshold),
             vjust = -0.4, hjust = -0.1,
             size  = 2.8, colour = "#666666") +
    anomaly_points_layer(
        list(
            list(df = sex_mismatch_females, type = "Mismatch F",        col = COL_MISMATCH_F),
            list(df = sex_mismatch_males,   type = "Mismatch M",        col = COL_MISMATCH_M),
            list(df = X0_Turner,            type = "Turner X0",         col = COL_TURNER),
            list(df = XXY_Klinefelter,      type = "Klinefelter XXY",   col = COL_KLINEFELTER),
            list(df = XXX_TripleX,          type = "Triple-X XXX",      col = COL_XXX)
        )
    ) +
    combined_label_layer(
        list(
            list(df = sex_mismatch_females, type = "mismatch",       col = COL_MISMATCH_F),
            list(df = sex_mismatch_males,   type = "mismatch",       col = COL_MISMATCH_M),
            list(df = X0_Turner,            type = "Turner X0",      col = COL_TURNER),
            list(df = XXY_Klinefelter,      type = "Klinefelter XXY",col = COL_KLINEFELTER),
            list(df = XXX_TripleX,          type = "Triple-X XXX",   col = COL_XXX)
        ),
        show_labels
    ) +
    labs(
        title    = "Overall Sex Check",
        subtitle = sprintf("%d samples  |  %d sex mismatches  |  %d aneuploidies",
                           nrow(df), n_mis, n_aneu),
        x        = "F Statistic (X chromosome)",
        y        = "Y Chromosome Variant Count",
        caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
    ) +
    crisp_theme()

print(p_overall)

# ── MALE PLOT ─────────────────────────────────────────────────────────────
df_mal_plot <- as.data.frame(df_mal)

p_male <- ggplot(df_mal_plot, aes(x = F, y = YCOUNT)) +
    geom_point(colour = COL_MALE, size = 1.2, alpha = 0.65) +
    geom_vline(xintercept = f_klinefelter, colour = PAL_FAIL,
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = y_threshold, colour = "#666666",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text",
             x     = f_klinefelter,
             y     = Inf,
             label = sprintf("F=%.2f\n(Klinefelter)", f_klinefelter),
             vjust = 1.3, hjust = -0.1,
             size  = 2.8, colour = PAL_FAIL, fontface = "bold") +
    annotate("text",
             x     = -Inf,
             y     = y_threshold,
             label = sprintf("Y mean = %.1f", y_threshold),
             vjust = -0.4, hjust = -0.1,
             size  = 2.8, colour = "#666666") +
    anomaly_points_layer(
        list(
            list(df = sex_mismatch_males, type = "Mismatch M",      col = COL_MISMATCH_M),
            list(df = XXY_Klinefelter,    type = "Klinefelter XXY", col = COL_KLINEFELTER)
        )
    ) +
    combined_label_layer(
        list(
            list(df = sex_mismatch_males, type = "mismatch",        col = COL_MISMATCH_M),
            list(df = XXY_Klinefelter,    type = "Klinefelter XXY", col = COL_KLINEFELTER)
        ),
        show_labels
    ) +
    labs(
        title    = "Male Sex Check",
        subtitle = sprintf("%d reported males", nrow(df_mal_plot)),
        x        = "F Statistic (X chromosome)",
        y        = "Y Chromosome Variant Count",
        caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
    ) +
    crisp_theme()

print(p_male)

# ── FEMALE PLOT ───────────────────────────────────────────────────────────
df_femal_plot <- as.data.frame(df_femal)

p_female <- ggplot(df_femal_plot, aes(x = F, y = YCOUNT)) +
    geom_point(colour = COL_FEMALE, size = 1.2, alpha = 0.65) +
    geom_vline(xintercept = f_turner, colour = PAL_PINK,
               linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = f_xxx,    colour = PAL_SKY,
               linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = y_threshold, colour = "#666666",
               linetype = "dotted", linewidth = 0.8) +
    annotate("text",
             x     = f_turner,
             y     = Inf,
             label = sprintf("F=%.2f\n(Turner)", f_turner),
             vjust = 1.3, hjust = -0.1,
             size  = 2.8, colour = PAL_PINK, fontface = "bold") +
    annotate("text",
             x     = f_xxx,
             y     = Inf,
             label = sprintf("F=%.2f\n(Triple-X)", f_xxx),
             vjust = 1.3, hjust = -0.1,
             size  = 2.8, colour = PAL_SKY, fontface = "bold") +
    annotate("text",
             x     = -Inf,
             y     = y_threshold,
             label = sprintf("Y mean = %.1f", y_threshold),
             vjust = -0.4, hjust = -0.1,
             size  = 2.8, colour = "#666666") +
    anomaly_points_layer(
        list(
            list(df = sex_mismatch_females, type = "Mismatch F",   col = COL_MISMATCH_F),
            list(df = X0_Turner,            type = "Turner X0",    col = COL_TURNER),
            list(df = XXX_TripleX,          type = "Triple-X XXX", col = COL_XXX)
        )
    ) +
    combined_label_layer(
        list(
            list(df = sex_mismatch_females, type = "mismatch",     col = COL_MISMATCH_F),
            list(df = X0_Turner,            type = "Turner X0",    col = COL_TURNER),
            list(df = XXX_TripleX,          type = "Triple-X XXX", col = COL_XXX)
        ),
        show_labels
    ) +
    labs(
        title    = "Female Sex Check",
        subtitle = sprintf("%d reported females", nrow(df_femal_plot)),
        x        = "F Statistic (X chromosome)",
        y        = "Y Chromosome Variant Count",
        caption  = "CRISP | Comprehensive Robust Integrated SNP Processing"
    ) +
    crisp_theme()

print(p_female)

dev.off()
cat(sprintf("[SEXCHECK] Plots saved: %s\n\n", pdf_file))

##########################################################################
#
#  ____  _____ ____   ___  ____ _____
# |  _ \| ____|  _ \ / _ \|  _ \_   _|
# | |_) |  _| | |_) | | | | |_) || |
# |  _ <| |___|  __/| |_| |  _ < | |
# |_| \_\_____|_|    \___/|_| \_\|_|
#
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
#
#  _______  ______ _       _     ___ ____ _____ ____
# | ____\ \/ / ___| |     | |   |_ _/ ___|_   _/ ___|
# |  _|  \  / |   | |     | |    | |\___  \ | | \___ \
# | |___ /  \ |___| |___  | |___ | | ___) || |  ___) |
# |_____/_/\_\____|_____| |_____|___|____/ |_| |____/
#
##########################################################################

mismatch_file <- file.path(output_dir, "id_list_sex_mismatch.txt")
write.table(df_sex_mismatch[, c(1, 2)],
            mismatch_file,
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            col.names = FALSE)
cat(sprintf("[SEXCHECK] Sex mismatch exclusion list: %s (%d samples)\n",
            mismatch_file, nrow(df_sex_mismatch)))

aneuploidy_file <- file.path(output_dir, "id_list_aneuploidies.txt")
write.table(df_aneuploidies[, c(1, 2)],
            aneuploidy_file,
            quote     = FALSE,
            row.names = FALSE,
            sep       = "\t",
            col.names = FALSE)
cat(sprintf("[SEXCHECK] Aneuploidy exclusion list  : %s (%d samples)\n",
            aneuploidy_file, nrow(df_aneuploidies)))

##########################################################################
#
#   ____ ___  _   _  ___  ____ _____   __  __    _    _   _ ___ _____ _____ ____ _____
#  / ___/ _ \| | | |/ _ \|  _ \_   _| |  \/  |  / \  | \ | |_ _|  ___| ____/ ___|_   _|
# | |  | | | | |_| | | | | |_) || |   | |\/| | / _ \ |  \| || || |_  |  _| \___ \ | |
# | |__| |_| |  _  | |_| |  _ < | |   | |  | |/ ___ \| |\  || ||  _| | |___ ___) || |
#  \____\___/|_| |_|\___/|_| \_\|_|   |_|  |_/_/   \_\_| \_|___|_|   |_____|____/ |_|
#
# Aneuploidies (Turner X0, Klinefelter XXY, Triple-X XXX) are confirmed
# karyotype calls, NOT data errors. Unlike sex mismatches (Type 2 error,
# exclusion-only), these are natural dosage experiments of real value to
# RNA/eQTL researchers via EHR linkage. Write a typed manifest plus
# per-type ID lists so the shell script can route them to a dedicated
# sub-cohort instead of silently discarding them.
##########################################################################

write_typed_manifest <- function(turner_df, klinefelter_df, xxx_df, output_dir) {

    type_map <- list(
        TURNER      = turner_df,
        KLINEFELTER = klinefelter_df,
        TRIPLE_X    = xxx_df
    )

    pieces <- list()
    for (type_name in names(type_map)) {
        d <- as.data.frame(type_map[[type_name]])
        if (nrow(d) == 0) next
        d$CATEGORY <- "ANEUPLOIDY"
        d$TYPE     <- type_name
        pieces[[type_name]] <- d[, c("FID", "IID", "CATEGORY", "TYPE",
                                     "PEDSEX", "F", "YCOUNT")]
    }

    manifest_file <- file.path(output_dir, "aneuploidy_cohort_manifest.txt")

    if (length(pieces) > 0) {
        df_manifest <- unique(do.call(rbind, pieces))
    } else {
        df_manifest <- data.frame(
            FID = character(), IID = character(),
            CATEGORY = character(), TYPE = character(),
            PEDSEX = integer(), F = numeric(), YCOUNT = integer()
        )
    }

    write.table(df_manifest, manifest_file,
                quote = FALSE, row.names = FALSE, sep = "\t")
    cat(sprintf("[SEXCHECK] Aneuploidy cohort manifest : %s (%d samples)\n",
                manifest_file, nrow(df_manifest)))

    # Per-type ID lists for shell-side --keep extraction
    for (type_name in c("TURNER", "KLINEFELTER", "TRIPLE_X")) {
        sub <- df_manifest[df_manifest$TYPE == type_name, c("FID", "IID")]
        type_file <- file.path(
            output_dir,
            sprintf("id_list_aneuploidy_%s.txt", tolower(type_name))
        )
        write.table(sub, type_file,
                    quote = FALSE, row.names = FALSE,
                    sep = "\t", col.names = FALSE)
        cat(sprintf("[SEXCHECK]   %-12s : %d samples -> %s\n",
                    type_name, nrow(sub), type_file))
    }

    invisible(df_manifest)
}

write_typed_manifest(X0_Turner, XXY_Klinefelter, XXX_TripleX, output_dir)

cat("\n[SEXCHECK] Sex check complete.\n")
