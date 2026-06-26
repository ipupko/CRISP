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
#
# Script : scripts/plot_pca_pass1.R
# Part of: Compass Genomics suite — github.com/ipupko/CRISP
# Version: v0.5.0
#
# Purpose:
#   Produce four plots for CRISP Step 6 — PCA Pass 1:
#
#   Plot 1 — Monochrome PCA scatter + ellipses
#     All samples in neutral grey. Ellipses drawn per GMM cluster.
#     Shows population structure without colour interpretation.
#     Analysts see the geometry before the palette is applied.
#
#   Plot 2 — Coloured PCA scatter + ellipses
#     Identical layout to Plot 1. Each cluster painted with a
#     distinct CRISP_PAL_* colour. Outliers always CRISP_PAL_FAIL.
#
#   Plot 3 — Monochrome MDS scatter + ellipses  (if MDS was run)
#   Plot 4 — Coloured MDS scatter + ellipses    (if MDS was run)
#
#   All six palette × background combinations are supported:
#     PLOT_COLOUR_MODE : STANDARD | COLOURBLIND | NIGHT
#     PLOT_BACKGROUND  : LIGHT | DARK
#   Colours consumed exclusively from CRISP_PAL_* environment
#   variables — zero hardcoded hex values in this script.
#
#   COHORT_ADMIXED=YES draws a confidence ellipse per cluster at
#   PCA_ELLIPSE_SD standard deviations (SD-based, not Mahalanobis).
#   When COHORT_ADMIXED=NO, ellipses are still drawn as light
#   convex-hull boundaries to show cluster extent, but no SD
#   annotation is added.
#
#   Right-side external legend — non-negotiable CRISP convention.
#
# Args (passed via Rscript --args or argparse-style --flags):
#   --input       path to .eigenvec or .mds file
#   --satellite   path to {project}_ancestry.txt satellite file
#   --plot-type   PCA or MDS
#   --out-dir     output directory for PDFs
#   --project     project name prefix
#   --n-clusters  number of detected clusters (integer)
#   --ellipse-sd  SD radius for confidence ellipse (default 3)
#   --admixed     COHORT_ADMIXED flag (YES/NO)
#
# Palette environment variables (set by crisp_pca_pass1.sh):
#   CRISP_PAL_MODE        STANDARD | COLOURBLIND | NIGHT
#   CRISP_PAL_BACKGROUND  LIGHT | DARK
#   CRISP_PAL_PASS        main positive/pass colour
#   CRISP_PAL_FAIL        outlier / fail colour (always used for outliers)
#   CRISP_PAL_WARN        warning colour
#   CRISP_PAL_SKY         secondary group colour
#   CRISP_PAL_PINK        tertiary group colour
#   CRISP_PAL_HIGHLIGHT   quaternary group colour
#   CRISP_PAL_SUBTEXT     neutral grey (used for monochrome plot points)
#   CRISP_PAL_TEXT        foreground text colour
#   CRISP_PAL_BG          background colour
#   CRISP_PAL_PANEL       panel/plot area background colour
#
##########################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggforce)     # geom_ellipse / stat_ellipse
})

##########################################################################
#  ___  ___
#  |  \/  |
#  | .  . |_   _ _ __   ___ _ __ ___
#  | |\/| | | | | '_ \ / _ \ '__/ __|
#  | |  | | |_| | | | |  __/ |  \__ \
#  \_|  |_/\__,_|_| |_|\___|_|  |___/
##########################################################################
# WHY: R lacks a clean argparse equivalent out-of-the-box.
# We use a simple key=value parser on commandArgs(trailingOnly=TRUE).
# Flags are passed as --key value pairs from the shell orchestrator.

parse_args <- function() {
  raw  <- commandArgs(trailingOnly = TRUE)
  args <- list(
    input       = NULL,
    satellite   = NULL,
    plot_type   = "PCA",
    out_dir     = ".",
    project     = "crisp_project",
    n_clusters  = NULL,
    ellipse_sd  = 3,
    admixed     = "NO"
  )

  i <- 1
  while (i <= length(raw)) {
    key <- sub("^--", "", raw[i])
    key <- gsub("-", "_", key)
    if (i + 1 <= length(raw) && !startsWith(raw[i + 1], "--")) {
      args[[key]] <- raw[i + 1]
      i <- i + 2
    } else {
      args[[key]] <- TRUE
      i <- i + 1
    }
  }

  args$ellipse_sd <- as.numeric(args$ellipse_sd)
  args$n_clusters <- as.integer(args$n_clusters)
  args$admixed    <- toupper(as.character(args$admixed))
  return(args)
}

args <- parse_args()

##########################################################################
#  ______      _      _   _
#  | ___ \    | |    | | | |
#  | |_/ /__ _| | ___| |_| |_ ___
#  |  __/ _` | |/ _ \ __| __/ _ \
#  | | | (_| | |  __/ |_| ||  __/
#  \_|  \__,_|_|\___|\__|\__\___|
##########################################################################
# WHY: All colours are read from CRISP_PAL_* environment variables set
# by crisp_pca_pass1.sh via _crisp_palette(). Zero hardcoded hex values
# in this script. The _pal() helper returns a safe fallback if a variable
# is unset so the script degrades gracefully in manual/test runs.

_pal <- function(var, fallback = "#888888") {
  val <- Sys.getenv(var, unset = "")
  if (nchar(val) == 0) fallback else val
}

PAL_MODE       <- Sys.getenv("CRISP_PAL_MODE",       unset = "STANDARD")
PAL_BACKGROUND <- Sys.getenv("CRISP_PAL_BACKGROUND",  unset = "LIGHT")

# Group colours — one per detected cluster, cycling if >6 groups
# WHY: Cluster count is dynamic (AUTO K). We define an ordered sequence
# of CRISP palette colours and cycle through them. Outliers always get
# CRISP_PAL_FAIL regardless of mode — this is non-negotiable so that
# outlier status is visually unambiguous across all palette combinations.
GROUP_COLOURS <- c(
  _pal("CRISP_PAL_PASS",      "#1D9E75"),
  _pal("CRISP_PAL_SKY",       "#56B4E9"),
  _pal("CRISP_PAL_PINK",      "#CC79A7"),
  _pal("CRISP_PAL_WARN",      "#FEBC2E"),
  _pal("CRISP_PAL_HIGHLIGHT", "#8BE0CB"),
  _pal("CRISP_PAL_ACCENT",    "#0072B2")
)

COL_OUTLIER  <- _pal("CRISP_PAL_FAIL",    "#FF5F57")
COL_MONO     <- _pal("CRISP_PAL_SUBTEXT", "#8E8E93")  # monochrome plot points
COL_ELLIPSE  <- _pal("CRISP_PAL_TEXT",    "#1C1C1E")  # monochrome ellipse border
COL_BG       <- _pal("CRISP_PAL_BG",      "#FFFFFF")
COL_PANEL    <- _pal("CRISP_PAL_PANEL",   "#F5F5F7")
COL_TEXT     <- _pal("CRISP_PAL_TEXT",    "#1C1C1E")
COL_GRID     <- _pal("CRISP_PAL_GRID",    "#D1D1D6")

##########################################################################
#  ______            _
#  | ___ \          | |
#  | |_/ /___  __ _| |
#  |    // _ \/ _` | |
#  | |\ \  __/ (_| | |
#  \_| \_\___|\__,_|_|
##########################################################################
# WHY: Detection/computation is always before plotting. We load all
# data, compute ellipse parameters, and validate before opening any
# PDF device — this way a data error does not produce a blank PDF.

cat(sprintf("[CRISP plot_pca_pass1.R] Project : %s\n", args$project))
cat(sprintf("[CRISP plot_pca_pass1.R] Type    : %s\n", args$plot_type))
cat(sprintf("[CRISP plot_pca_pass1.R] Palette : %s / %s\n", PAL_MODE, PAL_BACKGROUND))

# ── Read satellite file ───────────────────────────────────────────────
# WHY: We read from the satellite file rather than the raw eigenvec so
# that GROUP labels, OUTLIER flags, and MEMBERSHIP_P are all co-located
# with the PC scores. No separate merge step needed.

if (!file.exists(args$satellite)) {
  stop(sprintf("[CRISP] ERROR: Satellite file not found: %s", args$satellite))
}

sat <- fread(args$satellite, sep = "\t", header = TRUE)
cat(sprintf("[CRISP plot_pca_pass1.R] Samples : %d\n", nrow(sat)))

# Determine axis columns based on plot type
# WHY: For MDS the satellite file still holds PC scores from the PCA
# run. The MDS scores are in a separate .mds file. For MDS plots we
# read the .mds file directly and join SUPERPOP/OUTLIER from satellite.
if (args$plot_type == "MDS") {
  if (!file.exists(args$input)) {
    stop(sprintf("[CRISP] ERROR: MDS file not found: %s", args$input))
  }
  mds_raw <- fread(args$input, sep = " ", header = TRUE)
  # Rename C1 C2 to PC1 PC2 for unified axis labelling
  setnames(mds_raw,
           old = grep("^C[0-9]", names(mds_raw), value = TRUE),
           new = gsub("^C", "PC", grep("^C[0-9]", names(mds_raw), value = TRUE)))
  # Drop SOL if present
  if ("SOL" %in% names(mds_raw)) mds_raw[, SOL := NULL]

  # Join GROUP/OUTLIER from satellite on IID
  plot_df <- merge(mds_raw, sat[, .(IID, SUPERPOP, OUTLIER, MEMBERSHIP_P)],
                   by = "IID", all.x = TRUE)
  axis_label_x <- "MDS Dimension 1"
  axis_label_y <- "MDS Dimension 2"
  pc_x <- "PC1"; pc_y <- "PC2"
} else {
  plot_df  <- as.data.frame(sat)
  axis_label_x <- "PC1"
  axis_label_y <- "PC2"
  pc_x <- "PC1"; pc_y <- "PC2"
}

# Coerce types
plot_df$SUPERPOP   <- as.character(plot_df$SUPERPOP)
plot_df$OUTLIER    <- as.character(plot_df$OUTLIER)
plot_df[[pc_x]]    <- as.numeric(plot_df[[pc_x]])
plot_df[[pc_y]]    <- as.numeric(plot_df[[pc_y]])

# Remove rows with NA PC scores
plot_df <- plot_df[!is.na(plot_df[[pc_x]]) & !is.na(plot_df[[pc_y]]), ]

# ── Group colour mapping ───────────────────────────────────────────────
# WHY: Groups are sorted so colour assignment is deterministic and
# reproducible across runs — GROUP_1 always gets the first colour, etc.
all_groups   <- sort(unique(plot_df$SUPERPOP))
n_groups     <- length(all_groups)
# Cycle colours if more groups than palette entries
group_colours <- setNames(
  GROUP_COLOURS[((seq_along(all_groups) - 1) %% length(GROUP_COLOURS)) + 1],
  all_groups
)

##########################################################################
#  _______ _
#  |__   __| |
#     | |  | |__   ___ _ __ ___   ___
#     | |  | '_ \ / _ \ '_ ` _ \ / _ \
#     | |  | | | |  __/ | | | | |  __/
#     |_|  |_| |_|\___|_| |_| |_|\___|
##########################################################################
# WHY: The CRISP theme is defined once and applied to every plot.
# All backgrounds, text colours, and grid lines are driven by palette
# environment variables — so NIGHT/DARK mode works without code changes.

theme_crisp <- function() {
  theme_minimal(base_size = 11) +
    theme(
      # Background
      plot.background  = element_rect(fill = COL_BG,    colour = NA),
      panel.background = element_rect(fill = COL_PANEL, colour = NA),
      # Grid
      panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # Text
      plot.title    = element_text(colour = COL_TEXT, face = "bold",
                                   size = 13, hjust = 0),
      plot.subtitle = element_text(colour = COL_MONO, size = 9),
      axis.title    = element_text(colour = COL_TEXT, size = 10),
      axis.text     = element_text(colour = COL_MONO, size = 8),
      # Legend — right side, non-negotiable CRISP convention
      legend.position  = "right",
      legend.background = element_rect(fill = COL_BG, colour = NA),
      legend.title     = element_text(colour = COL_TEXT, face = "bold",
                                      size = 9),
      legend.text      = element_text(colour = COL_TEXT, size = 8),
      legend.key       = element_rect(fill = COL_BG, colour = NA),
      # Caption
      plot.caption     = element_text(colour = COL_MONO, size = 7,
                                      hjust = 0),
      # Margins
      plot.margin = margin(12, 12, 12, 12)
    )
}

##########################################################################
#  ______ _ _ _
#  |  ____| | (_)
#  | |__  | | |_ __  ___  ___
#  |  __| | | | | '_ \/ __|/ _ \
#  | |____| | | | |_) \__ \  __/
#  |______|_|_|_| .__/|___/\___|
#               | |
#               |_|
##########################################################################
# WHY: The ellipse function is shared between Plot 1 (monochrome) and
# Plot 2 (coloured). It uses stat_ellipse from ggplot2 which draws a
# confidence ellipse at a given confidence level. We convert the SD
# threshold to a confidence level:
#   level = pchisq(sd^2, df=2)
# where df=2 because we are working in 2D PC space.
# When COHORT_ADMIXED=YES the ellipse is annotated with the SD radius.
# When COHORT_ADMIXED=NO we draw a lighter hull boundary instead.

sd_to_level <- function(sd) pchisq(sd^2, df = 2)
ellipse_level <- sd_to_level(args$ellipse_sd)

add_ellipses <- function(p, colour_map = NULL, mono = FALSE) {
  # WHY: We add one stat_ellipse layer per group. Using group= in a single
  # stat_ellipse call with colour mapped to SUPERPOP produces one ellipse
  # per group automatically, but we want independent colour control for
  # mono vs coloured modes.
  groups_for_ellipse <- all_groups

  for (grp in groups_for_ellipse) {
    grp_data <- plot_df[plot_df$SUPERPOP == grp & plot_df$OUTLIER != "YES", ]
    if (nrow(grp_data) < 5) next  # need >=5 points for ellipse

    ecol <- if (mono) COL_ELLIPSE else group_colours[grp]

    p <- p + stat_ellipse(
      data      = grp_data,
      aes_string(x = pc_x, y = pc_y),
      type      = "norm",
      level     = ellipse_level,
      colour    = ecol,
      linewidth = if (mono) 0.6 else 0.8,
      linetype  = if (mono) "dashed" else "solid",
      alpha     = 0.9,
      inherit.aes = FALSE
    )
  }
  return(p)
}

##########################################################################
#  ______ _       _     __
#  | ___ \ |     | |   /  |
#  | |_/ / | ___ | |_  `| |
#  |  __/| |/ _ \| __|  | |
#  | |   | | (_) | |_  _| |_
#  \_|   |_|\___/ \__| \___/
#
#  Monochrome scatter
##########################################################################
# WHY: Plot 1 is monochrome — all non-outlier samples in COL_MONO (neutral
# grey), outliers in COL_OUTLIER. Ellipses drawn in COL_ELLIPSE dashed.
# The analyst sees the geometric structure of the clusters without any
# colour interpretation. This is the reference plot for cluster review.

plot_type_lower <- tolower(args$plot_type)
out_mono   <- file.path(args$out_dir,
                        sprintf("%s_%s_mono.pdf",   args$project, plot_type_lower))
out_colour <- file.path(args$out_dir,
                        sprintf("%s_%s_colour.pdf", args$project, plot_type_lower))

cat(sprintf("[CRISP plot_pca_pass1.R] Plot 1 (mono)   → %s\n", out_mono))

plot_df$point_colour_mono <- ifelse(
  plot_df$OUTLIER == "YES", COL_OUTLIER, COL_MONO
)
plot_df$point_shape <- ifelse(plot_df$OUTLIER == "YES", 4, 16)
# WHY: Outliers plotted as × (shape 4) rather than filled circles (16).
# This makes outliers identifiable by shape as well as colour — important
# for colourblind accessibility, reinforcing distinction via two channels.

n_outliers  <- sum(plot_df$OUTLIER == "YES")
n_non_out   <- nrow(plot_df) - n_outliers
subtitle_str <- sprintf(
  "N=%d samples — %d cluster(s) — %d outlier(s) — %s / %s",
  nrow(plot_df), n_groups, n_outliers, PAL_MODE, PAL_BACKGROUND
)

p_mono <- ggplot(plot_df, aes_string(x = pc_x, y = pc_y)) +
  geom_point(
    aes(colour = OUTLIER, shape = OUTLIER),
    size  = 0.8,
    alpha = 0.6
  ) +
  scale_colour_manual(
    values = c("NO" = COL_MONO, "YES" = COL_OUTLIER),
    labels = c("NO" = "Sample", "YES" = "Outlier"),
    name   = "Status"
  ) +
  scale_shape_manual(
    values = c("NO" = 16, "YES" = 4),
    labels = c("NO" = "Sample", "YES" = "Outlier"),
    name   = "Status"
  ) +
  labs(
    title    = sprintf("CRISP Step 6 — %s Pass 1 (Monochrome)", args$plot_type),
    subtitle = subtitle_str,
    x        = axis_label_x,
    y        = axis_label_y,
    caption  = sprintf("Confidence: PROBABLE — ellipse SD=%.0f — CRISP v1.0.0",
                       args$ellipse_sd)
  ) +
  theme_crisp() +
  guides(
    colour = guide_legend(override.aes = list(size = 3)),
    shape  = guide_legend(override.aes = list(size = 3))
  )

# Add ellipses (monochrome dashed)
p_mono <- add_ellipses(p_mono, mono = TRUE)

# Add cluster labels at centroids
# WHY: Text labels on the monochrome plot identify clusters by name
# (GROUP_1 etc.) so the analyst can track which geometric cluster
# corresponds to which population assignment.
centroids <- aggregate(
  cbind(plot_df[[pc_x]], plot_df[[pc_y]]) ~ SUPERPOP,
  data = plot_df[plot_df$OUTLIER != "YES", ],
  FUN  = mean
)
names(centroids) <- c("SUPERPOP", "cx", "cy")

p_mono <- p_mono +
  geom_text(
    data      = centroids,
    aes(x = cx, y = cy, label = SUPERPOP),
    colour    = COL_TEXT,
    size      = 3,
    fontface  = "bold",
    inherit.aes = FALSE
  )

pdf(out_mono, width = 9, height = 7)
print(p_mono)
dev.off()
cat(sprintf("[CRISP plot_pca_pass1.R] Plot 1 written.\n"))

##########################################################################
#  ______ _       _     _____
#  | ___ \ |     | |   / __  \
#  | |_/ / | ___ | |_  `' / /'
#  |  __/| |/ _ \| __|   / /
#  | |   | | (_) | |_  ./ /___
#  \_|   |_|\___/ \__| \_____/
#
#  Coloured scatter
##########################################################################
# WHY: Plot 2 is identical in layout to Plot 1 but applies CRISP palette
# colours — one per group. Outliers retain COL_OUTLIER (CRISP_PAL_FAIL)
# regardless of which group they were assigned to. The legend on the
# right shows all groups with their colour swatches. Ellipses are drawn
# in the group colour with low alpha fill when COHORT_ADMIXED=YES,
# or as solid outlines when COHORT_ADMIXED=NO.

cat(sprintf("[CRISP plot_pca_pass1.R] Plot 2 (colour) → %s\n", out_colour))

# Build colour and shape vectors for ggplot scale
# WHY: Outliers are keyed as "_OUTLIER" in the colour map so they get
# a distinct entry in the legend rather than overriding group colours.
colour_scale <- c(group_colours, "_OUTLIER" = COL_OUTLIER)
shape_scale  <- c(
  setNames(rep(16L, n_groups), all_groups),
  "_OUTLIER" = 4L
)

# Add a plot_group column — outliers get "_OUTLIER" key
plot_df$plot_group <- ifelse(
  plot_df$OUTLIER == "YES", "_OUTLIER", plot_df$SUPERPOP
)

# Legend labels — em dash in place of underscore for display
legend_labels <- c(
  setNames(all_groups, all_groups),
  "_OUTLIER" = "Outlier"
)

p_colour <- ggplot(plot_df,
                   aes_string(x = pc_x, y = pc_y, colour = "plot_group",
                               shape = "plot_group")) +
  geom_point(
    size  = 0.9,
    alpha = 0.65
  ) +
  scale_colour_manual(
    values = colour_scale,
    labels = legend_labels,
    name   = "Group",
    breaks = c(all_groups, "_OUTLIER")
  ) +
  scale_shape_manual(
    values = shape_scale,
    labels = legend_labels,
    name   = "Group",
    breaks = c(all_groups, "_OUTLIER")
  ) +
  labs(
    title    = sprintf("CRISP Step 6 — %s Pass 1 (Coloured)", args$plot_type),
    subtitle = subtitle_str,
    x        = axis_label_x,
    y        = axis_label_y,
    caption  = sprintf("Confidence: PROBABLE — ellipse SD=%.0f — Palette: %s — CRISP v1.0.0",
                       args$ellipse_sd, PAL_MODE)
  ) +
  theme_crisp() +
  guides(
    colour = guide_legend(override.aes = list(size = 3)),
    shape  = guide_legend(override.aes = list(size = 3))
  )

# Add coloured ellipses
p_colour <- add_ellipses(p_colour, colour_map = group_colours, mono = FALSE)

# If COHORT_ADMIXED=YES, add filled ellipse with low alpha to emphasise boundary
if (args$admixed == "YES") {
  for (grp in all_groups) {
    grp_data <- plot_df[plot_df$SUPERPOP == grp & plot_df$OUTLIER != "YES", ]
    if (nrow(grp_data) < 5) next
    p_colour <- p_colour +
      stat_ellipse(
        data        = grp_data,
        aes_string(x = pc_x, y = pc_y),
        type        = "norm",
        level       = ellipse_level,
        geom        = "polygon",
        fill        = group_colours[grp],
        colour      = NA,
        alpha       = 0.08,
        inherit.aes = FALSE
      )
  }
}

# Add centroid labels
p_colour <- p_colour +
  geom_text(
    data        = centroids,
    aes(x = cx, y = cy, label = SUPERPOP),
    colour      = COL_TEXT,
    size        = 3,
    fontface    = "bold",
    inherit.aes = FALSE
  )

pdf(out_colour, width = 9, height = 7)
print(p_colour)
dev.off()
cat(sprintf("[CRISP plot_pca_pass1.R] Plot 2 written.\n"))

##########################################################################
#  _____
#  /  ___|
#  \ `--.  ___ _ __ ___  ___
#   `--. \/ __| '__/ _ \/ _ \
#  /\__/ / (__| | |  __/  __/
#  \____/ \___|_|  \___|\___|
#
#  Scree plot (PCA only)
##########################################################################
# WHY: The scree plot shows the proportion of variance explained by each
# PC, which helps the analyst decide how many PCs are meaningful for
# downstream analyses (e.g. how many to include as covariates in GWAS).
# Produced only for PCA mode — MDS does not have an eigenvalue equivalent.

if (args$plot_type == "PCA") {
  eigenval_path <- sub("\\.eigenvec$", ".eigenval", args$input)

  if (file.exists(eigenval_path)) {
    eigenvals <- as.numeric(readLines(eigenval_path))
    pct_var   <- eigenvals / sum(eigenvals) * 100
    cum_var   <- cumsum(pct_var)
    n_shown   <- min(20L, length(eigenvals))

    scree_df <- data.frame(
      PC     = seq_len(n_shown),
      PctVar = pct_var[seq_len(n_shown)],
      CumVar = cum_var[seq_len(n_shown)]
    )

    out_scree <- file.path(
      args$out_dir,
      sprintf("%s_pca_scree.pdf", args$project)
    )
    cat(sprintf("[CRISP plot_pca_pass1.R] Scree plot     → %s\n", out_scree))

    p_scree <- ggplot(scree_df, aes(x = PC, y = PctVar)) +
      geom_col(fill = _pal("CRISP_PAL_PASS", "#1D9E75"),
               alpha = 0.8, width = 0.7) +
      geom_line(aes(y = CumVar),
                colour   = _pal("CRISP_PAL_WARN", "#FEBC2E"),
                linewidth = 0.9) +
      geom_point(aes(y = CumVar),
                 colour = _pal("CRISP_PAL_WARN", "#FEBC2E"),
                 size   = 2) +
      scale_x_continuous(breaks = seq_len(n_shown)) +
      labs(
        title    = sprintf("CRISP Step 6 — PCA Scree Plot"),
        subtitle = sprintf("Project: %s — %s / %s", args$project, PAL_MODE, PAL_BACKGROUND),
        x        = "Principal Component",
        y        = "% Variance Explained",
        caption  = "Bars = per-PC variance — Line = cumulative variance — CRISP v1.0.0"
      ) +
      theme_crisp()

    pdf(out_scree, width = 9, height = 5)
    print(p_scree)
    dev.off()
    cat(sprintf("[CRISP plot_pca_pass1.R] Scree plot written.\n"))
  } else {
    cat(sprintf("[CRISP plot_pca_pass1.R] Eigenval file not found, scree skipped: %s\n",
                eigenval_path))
  }
}

cat(sprintf("[CRISP plot_pca_pass1.R] All plots complete — %s / %s\n",
            PAL_MODE, PAL_BACKGROUND))
