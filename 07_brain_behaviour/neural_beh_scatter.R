## neural_beh_scatter.R
##
## Scatter panels of Figure 8B and Supplementary Figure S1B: each
## participant's mean MPH-minus-PBO Cue contrast estimate in the principal
## parietal cluster (extract_neural_beh_link.m), plotted against their
## MPH-minus-PBO change in hit rate (Figure 8B) or premature responses
## (Figure S1B), with the least-squares fit and its 95% confidence interval.
##
## The clusters were defined by their association with these behavioural
## measures in the whole-brain covariate models, so the panels display the
## form of that association and are not an independent test of it. No
## statistics are computed here.
##
## Inputs:
##   - neural_roi_values.csv   (extract_neural_beh_link.m)
##   - MID_beh_covariates.csv  (02_behaviour/MID_analysis.Rmd: per-subject
##     drug-effect differences, after the behavioural session exclusions)
##
## Output: MID_fig_cue_scatter_hitrate.{png,pdf} and
##         MID_fig_cue_scatter_premature.{png,pdf}

rm(list = ls())

## --- Paths -----------------------------------------------------------------
## Optional local copy of the data (same relative layout as the project
## folder), set via the MID_LOCAL_DATA environment variable; otherwise the
## project storage is used.
localCandidates <- Sys.getenv("MID_LOCAL_DATA")
localHit <- localCandidates[nzchar(localCandidates) & dir.exists(localCandidates)]

if (length(localHit) > 0) {
  projectDir <- localHit[1]
} else if (Sys.info()["sysname"] == "Windows") {
  projectDir <- "P:/3017048.01"
} else {
  projectDir <- "/project/3017048.01"
}
behDir <- file.path(projectDir, "bids", "derivatives", "beh", "MID", "stats")
outDir <- file.path(behDir, "MID_output")
dir.create(outDir, showWarnings = FALSE)

neural_file  <- file.path(behDir, "neural_roi_values.csv")
beh_cov_file <- file.path(outDir, "MID_beh_covariates.csv")

## --- Packages ----------------------------------------------------------
packages <- c("dplyr", "tidyr", "readr", "ggplot2", "tibble")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))

## --- Figure settings ---------------------------------------------------
fig_width  <- 3.6    # inches; sized to sit beside the slice-display image
fig_height <- 2.8

panel_map <- tibble::tribble(
  ~cluster,                        ~beh_var,        ~x_lab, ~y_lab, ~pub_base,
  "parietal_MPHvsPBO_Cue_x_hit",   "MPH_PBO_hit",
    "Methylphenidate effect in hit rate",
    "Parietal cluster \n (MPH \u2212 PBO contrast estimate)",
    "MID_fig_cue_scatter_hitrate",
  "parietal_MPHvsPBO_Cue_x_early", "MPH_PBO_early",
    "Methylphenidate effect in premature responses",
    "Parietal cluster \n (MPH \u2212 PBO contrast estimate)",
    "MID_fig_cue_scatter_premature"
)

## --- Load ----------------------------------------------------------------
neural  <- readr::read_csv(neural_file, show_col_types = FALSE)
beh_cov <- readr::read_csv(beh_cov_file, show_col_types = FALSE)

## --- Panels --------------------------------------------------------------
for (i in seq_len(nrow(panel_map))) {
  cl <- panel_map$cluster[i]
  bv <- panel_map$beh_var[i]

  df <- neural %>%
    dplyr::filter(cluster == !!cl) %>%
    dplyr::select(subject, roiMean) %>%
    dplyr::left_join(dplyr::select(beh_cov, subject, beh_val = dplyr::all_of(bv)),
                     by = "subject") %>%
    tidyr::drop_na()

  if (nrow(df) == 0) {
    message("No data for ", cl, " -- skipped (run extract_neural_beh_link.m first).")
    next
  }
  message(sprintf("%s: n = %d", cl, nrow(df)))

  pub <- ggplot2::ggplot(df, ggplot2::aes(x = beh_val, y = roiMean)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = .3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                         fill = "grey75", alpha = .35, linewidth = .7) +
    ggplot2::geom_point(shape = 21, size = 2.1, stroke = .4,
                        colour = "grey25", fill = "#a8bef7", alpha = .85) +
    ggplot2::labs(x = panel_map$x_lab[i], y = panel_map$y_lab[i]) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 10),
      axis.text  = ggplot2::element_text(size = 9, colour = "black"),
      axis.line  = ggplot2::element_line(linewidth = .4),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(file.path(outDir, paste0(panel_map$pub_base[i], ".", ext)),
                     pub, width = fig_width, height = fig_height, dpi = 300)
  }
}

message("Saved scatter panels to: ", outDir)
