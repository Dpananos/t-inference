# ---------------------------------------------------------------------------
# Simulation-based calibration for approximation fidelity.
#
# Effects are generated from the exact target prior; inference uses the fitted
# normal mixture. Ranks are drawn directly from Binomial(L, F_post(theta_true)),
# which is distributionally identical to generating L independent posterior
# draws and counting those below the true effect.
# ---------------------------------------------------------------------------

suppressMessages({
  library(SBC)
  library(ggplot2)
  library(patchwork)
})

sbc_ranks <- function(prior, n_sims = 10000L, L = 999L,
                      s_min = 0.001, s_max = 0.1) {
  theta_true <- prior$r(n_sims)
  s <- runif(n_sims, s_min, s_max)
  theta_hat <- rnorm(n_sims, theta_true, s)

  mixture <- active_mixture(prior$fit)
  s2 <- s^2
  log_weight <- vapply(
    seq_along(mixture$weights),
    function(k) {
      log(mixture$weights[k]) +
        dnorm(
          theta_hat,
          mixture$mu[k],
          sqrt(s2 + mixture$tau[k]^2),
          log = TRUE
        )
    },
    numeric(n_sims)
  )
  log_weight <- matrix(
    log_weight,
    nrow = n_sims,
    ncol = length(mixture$weights)
  )
  row_max <- apply(log_weight, 1, max)
  weight <- exp(log_weight - row_max)

  component_cdf <- vapply(
    seq_along(mixture$weights),
    function(k) {
      tau_post2 <- 1 / (1 / s2 + 1 / mixture$tau[k]^2)
      mu_post <- tau_post2 * (
        theta_hat / s2 +
          mixture$mu[k] / mixture$tau[k]^2
      )
      pnorm(theta_true, mu_post, sqrt(tau_post2))
    },
    numeric(n_sims)
  )
  component_cdf <- matrix(
    component_cdf,
    nrow = n_sims,
    ncol = length(mixture$weights)
  )

  p_below <- rowSums(weight * component_cdf) / rowSums(weight)
  if (any(!is.finite(p_below))) {
    stop("non-finite posterior CDF encountered during SBC")
  }
  p_below <- pmin(pmax(p_below, 0), 1)

  data.frame(
    theta_true = theta_true,
    s = s,
    theta_hat = theta_hat,
    posterior_cdf_at_truth = p_below,
    rank = rbinom(n_sims, L, p_below)
  )
}

run_sbc <- function(priors, n_sims = 10000L, L = 999L,
                    s_min = 0.001, s_max = 0.1,
                    seed = 80902991L) {
  set.seed(seed)
  simulations <- lapply(names(priors), function(name) {
    result <- sbc_ranks(
      priors[[name]],
      n_sims = n_sims,
      L = L,
      s_min = s_min,
      s_max = s_max
    )
    result$target <- priors[[name]]$label
    result
  })
  names(simulations) <- names(priors)

  ranks <- do.call(rbind, simulations)
  row.names(ranks) <- NULL
  summary <- do.call(rbind, lapply(simulations, function(result) {
    n <- nrow(result)
    observed_mean <- mean(result$rank)
    expected_mean <- L / 2
    expected_variance <- L * (L + 2) / 12
    mean_standard_error <- sqrt(expected_variance / n)
    rank_counts <- tabulate(result$rank + 1L, nbins = L + 1L)
    empirical_cdf <- cumsum(rank_counts) / n
    uniform_cdf <- seq_len(L + 1L) / (L + 1L)

    data.frame(
      target = result$target[1],
      simulations = n,
      posterior_draws = L,
      mean_rank = observed_mean,
      expected_mean_rank = expected_mean,
      mean_rank_z = (observed_mean - expected_mean) / mean_standard_error,
      rank_variance = var(result$rank),
      expected_rank_variance = expected_variance,
      maximum_ecdf_difference = max(abs(empirical_cdf - uniform_cdf)),
      mean_posterior_cdf_at_truth = mean(result$posterior_cdf_at_truth)
    )
  }))
  row.names(summary) <- NULL

  panels <- lapply(names(simulations), function(name) {
    sbc_ecdf_panel(
      simulations[[name]]$rank,
      L = L,
      title = priors[[name]]$label
    )
  })
  names(panels) <- names(simulations)

  plot <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = sprintf(
        "SBC of the normal-mixture approximation (%d simulations per target)",
        n_sims
      )
    )

  list(
    ranks = ranks,
    summary = summary,
    panels = panels,
    plot = plot,
    configuration = list(
      n_sims = n_sims,
      L = L,
      s_min = s_min,
      s_max = s_max,
      seed = seed
    )
  )
}

sbc_ecdf_panel <- function(ranks, L, title = NULL, base_size = 11,
                           show_parameter_strip = TRUE) {
  ranks_matrix <- matrix(
    ranks,
    ncol = 1,
    dimnames = list(NULL, "theta")
  )

  plot <- SBC::plot_ecdf_diff(ranks_matrix, max_rank = L) +
    labs(
      title = title,
      x = "Normalized rank",
      y = "ECDF difference"
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  if (!show_parameter_strip) {
    plot <- plot + theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  }
  plot
}
