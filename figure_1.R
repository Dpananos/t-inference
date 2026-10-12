#!/usr/bin/env Rapp
#| name: figure_1
#| description: Build the normal-mixture figure and write it to assets/.

# ---------------------------------------------------------------------------
# Panels
#   A  log density of the Student's t target and its mixture approximation
#   B  posterior-mean z-score against the observed z-score
#   C  implied shrinkage factor
#   D  simulation-based calibration of the approximation
#
# Run from the project root; the helper sources below are relative to it.
# ---------------------------------------------------------------------------

#| description: File to write. The extension picks the device (.png, .pdf, ...).
#| short: "o"
output <- "assets/figure_1.png"

#| description: Scale of the target Student's t prior.
scale <- 0.05

#| description: Degrees of freedom of the target Student's t prior.
t_df <- 3

#| description: Number of zero-mean normal components in the dictionary.
n_tau0 <- 80L

#| description: Gauss-Legendre nodes used to build the Pearson chi-square matrix.
n_quad <- 2048L

#| description: Simulated datasets for simulation-based calibration.
n_sims <- 10000L

#| description: Posterior draws per SBC simulation (the rank denominator).
n_draws <- 999L

#| description: Standard error at which the shrinkage panels are drawn.
s <- 0.1

#| description: Largest observed z-score shown in the shrinkage panels.
z_max <- 10

#| description: Figure width in inches.
width <- 12

#| description: Figure height in inches.
height <- 3.4

#| description: Resolution for raster devices.
dpi <- 300

#| description: Seed for the SBC simulations.
seed <- 80902991L

#| description: Suppress the fit and SBC diagnostic tables.
quiet <- FALSE

# The body lives in local() so intermediate objects do not become CLI options.
local({
  if (!file.exists("R/priors.R")) {
    stop(
      "run figure_1.R from the project root; R/priors.R was not found in ",
      getwd(),
      call. = FALSE
    )
  }

  library(ggplot2)
  library(patchwork)

  source("R/priors.R")
  source("R/figure_helpers.R")
  source("R/sbc_helpers.R")

  # Approximate the target by minimizing Pearson chi-square over a fixed normal
  # dictionary. make_priors() constructs the C matrix once.
  priors <- make_priors(
    sc = scale,
    t_df = t_df,
    n_tau0 = n_tau0,
    n_quad = n_quad
  )

  sbc_result <- run_sbc(
    priors,
    n_sims = n_sims,
    L = n_draws,
    seed = seed
  )

  figure <- paper_figure(priors, sbc_result, s = s, z_max = z_max)

  directory <- dirname(output)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  ggsave(output, figure, width = width, height = height, dpi = dpi)

  if (!quiet) {
    cat("\nMixture fit\n")
    print(fit_diagnostic_table(priors))
    cat("\nSimulation-based calibration\n")
    print(sbc_result$summary)
    cat("\n")
  }

  cat("wrote ", output, " (", width, " x ", height, " in)\n", sep = "")
})
