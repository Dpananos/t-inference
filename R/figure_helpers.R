# ---------------------------------------------------------------------------
# Panels and assembly for figure_1.R.
# ---------------------------------------------------------------------------

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

mixture_mean_variance <- function(fit) {
  mixture <- active_mixture(fit)
  mean <- sum(mixture$weights * mixture$mu)
  variance <- sum(
    mixture$weights * (mixture$tau^2 + mixture$mu^2)
  ) - mean^2
  c(mean = mean, variance = variance)
}

log_density_plot <- function(fit, name, target_colour = "black",
                             mixture_colour = "red") {
  line_df <- data.frame(theta = fit$xg, density = fit$hg)
  mixture_df <- data.frame(theta = fit$xg, density = fit$gg)

  ggplot(mapping = aes(theta, density)) +
    geom_line(
      data = line_df,
      aes(colour = "Target"),
      linewidth = 1
    ) +
    geom_line(
      data = mixture_df,
      aes(colour = "Mixture"),
      linewidth = 1.5,
      alpha = 0.55
    ) +
    scale_colour_manual(
      values = c(Target = target_colour, Mixture = mixture_colour),
      labels = c(Target = expression(h(theta)), Mixture = "Mixture"),
      name = NULL
    ) +
    scale_y_log10() +
    labs(x = expression(theta), y = "log density", title = name) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.95),
      legend.justification = c(1, 1),
      legend.background = element_blank()
    )
}

student_t_shrinkage_data <- function(fit, s = 0.1, z_max = 10,
                                     n_points = 600L) {
  moments <- mixture_mean_variance(fit)
  normal_shrinkage <- moments[["variance"]] / (s^2 + moments[["variance"]])
  z <- seq(-z_max, z_max, length.out = n_points)
  theta_hat <- z * s

  list(
    zscore = rbind(
      data.frame(
        z = z,
        value = posterior_mean_mix(theta_hat, s, fit) / s,
        method = "Mixture"
      ),
      data.frame(
        z = z,
        value = normal_shrinkage * z,
        method = "Normal-normal"
      )
    ),
    factor = rbind(
      data.frame(
        z = z[z >= 0],
        value = shrinkage_factor_mix(theta_hat, s, fit)[z >= 0],
        method = "Mixture"
      ),
      data.frame(
        z = z[z >= 0],
        value = normal_shrinkage,
        method = "Normal-normal"
      )
    ),
    normal_shrinkage = normal_shrinkage,
    s = s,
    z_max = z_max
  )
}

shrinkage_z_plot <- function(shrinkage, mixture_colour = "red",
                             normal_colour = "black") {
  ggplot(
    shrinkage$zscore,
    aes(z, value, colour = method, linetype = method)
  ) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_line(linewidth = 1) +
    scale_colour_manual(
      values = c(Mixture = mixture_colour, "Normal-normal" = normal_colour),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(Mixture = 1, "Normal-normal" = 2),
      name = NULL
    ) +
    coord_cartesian(
      xlim = c(-shrinkage$z_max, shrinkage$z_max),
      ylim = c(-shrinkage$z_max, shrinkage$z_max)
    ) +
    labs(
      x = expression(z == hat(theta) / s),
      y = expression(tilde(z) == E(theta ~ "|" ~ hat(theta), s) / s),
      title = "Shrinkage of the z-score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "inside",
      legend.position.inside = c(0.28, 0.85),
      legend.background = element_blank()
    )
}

shrinkage_factor_plot <- function(shrinkage, mixture_colour = "red",
                                  normal_colour = "black") {
  ggplot(
    shrinkage$factor,
    aes(z, value, colour = method, linetype = method)
  ) +
    geom_line(linewidth = 1) +
    scale_colour_manual(
      values = c(Mixture = mixture_colour, "Normal-normal" = normal_colour),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(Mixture = 1, "Normal-normal" = 2),
      name = NULL
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = expression(z == hat(theta) / s),
      y = expression(tilde(z) / z),
      title = "Shrinkage factor"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "inside",
      legend.position.inside = c(0.72, 0.22),
      legend.background = element_blank()
    )
}

paper_figure <- function(priors, sbc_result, s = 0.1, z_max = 10) {
  shrinkage <- student_t_shrinkage_data(
    priors$studentt$fit,
    s = s,
    z_max = z_max
  )
  student_t_ranks <- sbc_result$ranks$rank[
    sbc_result$ranks$target == priors$studentt$label
  ]
  if (length(student_t_ranks) == 0L) {
    stop("sbc_result does not contain Student-t ranks")
  }
  sbc_panel <- sbc_ecdf_panel(
    student_t_ranks,
    L = sbc_result$configuration$L,
    title = "Simulation-based calibration",
    base_size = 10,
    show_parameter_strip = FALSE
  )

  (
    log_density_plot(priors$studentt$fit, "Student's t") |
      shrinkage_z_plot(shrinkage) |
      shrinkage_factor_plot(shrinkage) |
      sbc_panel
  ) +
    plot_layout(widths = c(1, 1, 1, 1.05)) +
    plot_annotation(tag_levels = "A")
}

trapezoid_integral <- function(x, y) {
  sum((y[-1] + y[-length(y)]) / 2 * diff(x))
}

fit_diagnostic_table <- function(priors) {
  do.call(rbind, lapply(names(priors), function(name) {
    fit <- priors[[name]]$fit
    absolute_error <- abs(fit$gg - fit$hg)
    data.frame(
      target = priors[[name]]$label,
      candidate_components = length(fit$pi),
      active_components = fit$solution$active_components,
      chisq = fit$solution$chisq,
      tv_upper_bound = 0.5 * sqrt(max(fit$solution$chisq, 0)),
      fitting_interval_tv = 0.5 * trapezoid_integral(
        fit$xg,
        absolute_error
      ),
      fitting_interval_max_abs_error = max(absolute_error),
      kkt_residual = unname(fit$solution$kkt[["maximum"]]),
      minimum_C_diagonal = min(diag(fit$C)),
      C_symmetry_error = max(abs(fit$C - t(fit$C))),
      quadrature_nodes = nrow(fit$quadrature),
      row.names = NULL
    )
  }))
}
