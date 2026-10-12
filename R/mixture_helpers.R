# ---------------------------------------------------------------------------
# Normal-mixture approximation with Pearson chi-square fitting.
#
# For a target density pi and fixed normal atoms phi_j, the fitting criterion is
#
#   chi^2(q_w || pi) = w' C w - 1,
#   C_jk = integral phi_j(theta) phi_k(theta) / pi(theta) d theta.
#
# build_chisq_problem() constructs C exactly once for a dictionary. The solver,
# diagnostics, fitted-density evaluation, and downstream posterior calculations
# all reuse that stored matrix.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(quadprog))

# --- Deterministic quadrature ----------------------------------------------

# Dependency-free Gauss-Legendre rule on [0, 1]. The probability-scale
# transformation theta = F_pi^{-1}(u) turns
#
#   C_jk = integral phi_j phi_k / pi d theta
#
# into E_pi[(phi_j / pi) (phi_k / pi)], which this rule evaluates
# deterministically without truncating the fitting interval at q_alpha.
gauss_legendre_rule <- function(n = 512L, tolerance = 1e-14) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("n must be a single integer of at least 2")
  }

  nodes <- numeric(n)
  weights <- numeric(n)
  n_roots <- (n + 1L) %/% 2L

  for (i in seq_len(n_roots)) {
    z <- cos(pi * (i - 0.25) / (n + 0.5))
    repeat {
      p_nm1 <- 1
      p_n <- z
      if (n > 1L) {
        for (j in 2:n) {
          p_np1 <- ((2 * j - 1) * z * p_n - (j - 1) * p_nm1) / j
          p_nm1 <- p_n
          p_n <- p_np1
        }
      }
      derivative <- n * (z * p_n - p_nm1) / (z^2 - 1)
      z_next <- z - p_n / derivative
      if (abs(z_next - z) <= tolerance) {
        z <- z_next
        break
      }
      z <- z_next
    }

    # Re-evaluate the derivative at the converged root.
    p_nm1 <- 1
    p_n <- z
    if (n > 1L) {
      for (j in 2:n) {
        p_np1 <- ((2 * j - 1) * z * p_n - (j - 1) * p_nm1) / j
        p_nm1 <- p_n
        p_n <- p_np1
      }
    }
    derivative <- n * (z * p_n - p_nm1) / (z^2 - 1)
    weight <- 2 / ((1 - z^2) * derivative^2)

    nodes[i] <- -z
    nodes[n + 1L - i] <- z
    weights[i] <- weight
    weights[n + 1L - i] <- weight
  }

  list(
    probability = (nodes + 1) / 2,
    weight = weights / 2
  )
}

normal_dictionary <- function(theta, atoms) {
  if (!all(c("mu", "tau") %in% names(atoms))) {
    stop("atoms must contain mu and tau columns")
  }
  if (any(!is.finite(atoms$mu)) || any(!is.finite(atoms$tau)) ||
      any(atoms$tau <= 0)) {
    stop("atom means must be finite and atom standard deviations positive")
  }

  values <- vapply(
    seq_len(nrow(atoms)),
    function(j) dnorm(theta, atoms$mu[j], atoms$tau[j]),
    numeric(length(theta))
  )
  matrix(values, nrow = length(theta), ncol = nrow(atoms))
}

# Construct the fixed chi-square quadratic program for one target/dictionary.
# C is formed as a weighted Gram matrix, so the numerical approximation is
# symmetric and positive semidefinite by construction.
build_chisq_problem <- function(td, tq, atoms, n_quad = 2048L,
                                diagonal_tolerance = 5e-3) {
  rule <- gauss_legendre_rule(n_quad)
  theta <- tq(rule$probability)
  target <- td(theta)

  if (any(!is.finite(theta))) {
    stop("target quantile function returned non-finite quadrature nodes")
  }
  if (any(!is.finite(target)) || any(target <= 0)) {
    stop(
      "target density must be finite and positive at every quadrature node; ",
      "check support and tail compatibility"
    )
  }

  atom_density <- normal_dictionary(theta, atoms)
  density_ratio <- atom_density / target
  if (any(!is.finite(density_ratio))) {
    stop(
      "non-finite component-to-target density ratio; ",
      "the dictionary is not numerically tail-compatible with the target"
    )
  }

  weighted_ratio <- density_ratio * sqrt(rule$weight)
  C <- crossprod(weighted_ratio)
  C <- (C + t(C)) / 2

  if (any(!is.finite(C))) {
    stop("C contains non-finite entries; check support and tail compatibility")
  }
  if (max(abs(C - t(C))) > 100 * .Machine$double.eps * max(1, max(abs(C)))) {
    stop("C is not symmetric within numerical tolerance")
  }
  if (min(diag(C)) < 1 - diagonal_tolerance) {
    stop(
      "C has a diagonal entry below 1 beyond quadrature tolerance; ",
      "increase n_quad or revise the dictionary"
    )
  }

  structure(
    list(
      C = C,
      atoms = atoms,
      quadrature = data.frame(
        probability = rule$probability,
        theta = theta,
        weight = rule$weight,
        target = target
      ),
      n_quad = as.integer(n_quad)
    ),
    class = "chisq_mixture_problem"
  )
}

# --- Simplex quadratic program and diagnostics ------------------------------
#
# The ridge is an absolute addition to the diagonal of C. It is deliberately
# not scaled by C so the numerical regularization is fixed across dictionaries.

chisq_kkt_diagnostics <- function(weights, C, ridge = 0,
                                  active_tolerance = 1e-10) {
  gradient <- as.vector(2 * (C %*% weights + ridge * weights))
  active <- weights > active_tolerance
  if (!any(active)) {
    stop("KKT diagnostics require at least one active weight")
  }

  multiplier <- mean(gradient[active])
  inactive <- !active
  active_residual <- max(abs(gradient[active] - multiplier))
  inactive_residual <- if (any(inactive)) {
    max(pmax(multiplier - gradient[inactive], 0))
  } else {
    0
  }

  c(
    equality = abs(sum(weights) - 1),
    nonnegativity = max(pmax(-weights, 0)),
    stationarity = active_residual,
    dual_feasibility = inactive_residual,
    maximum = max(
      abs(sum(weights) - 1),
      max(pmax(-weights, 0)),
      active_residual,
      inactive_residual
    )
  )
}

solve_chisq_weights <- function(C, ridge = 1e-8,
                                feasibility_tolerance = 1e-8,
                                kkt_tolerance = 1e-6,
                                active_tolerance = 1e-10) {
  if (!is.matrix(C) || nrow(C) != ncol(C) || nrow(C) < 1L) {
    stop("C must be a non-empty square matrix")
  }
  if (any(!is.finite(C))) {
    stop("C must contain only finite values")
  }

  C <- (C + t(C)) / 2
  n_atoms <- nrow(C)

  solution <- quadprog::solve.QP(
    Dmat = 2 * (C + ridge * diag(n_atoms)),
    dvec = numeric(n_atoms),
    Amat = cbind(rep(1, n_atoms), diag(n_atoms)),
    bvec = c(1, numeric(n_atoms)),
    meq = 1
  )
  weights <- solution$solution

  if (abs(sum(weights) - 1) > feasibility_tolerance ||
      min(weights) < -feasibility_tolerance) {
    stop("quadprog returned weights outside the simplex tolerance")
  }

  kkt <- chisq_kkt_diagnostics(
    weights,
    C,
    ridge = ridge,
    active_tolerance = active_tolerance
  )
  kkt_scale <- max(1, max(abs(2 * C %*% weights)))
  if (unname(kkt["maximum"]) > kkt_tolerance * kkt_scale) {
    stop(sprintf(
      "quadprog KKT residual %.3e exceeds tolerance %.3e",
      unname(kkt["maximum"]),
      kkt_tolerance * kkt_scale
    ))
  }

  structure(
    list(
      weights = weights,
      chisq = as.numeric(crossprod(weights, C %*% weights)) - 1,
      regularized_objective = as.numeric(
        crossprod(weights, C %*% weights) +
          ridge * crossprod(weights)
      ),
      ridge = ridge,
      active_tolerance = active_tolerance,
      active_components = sum(weights > active_tolerance),
      kkt = kkt,
      iterations = solution$iterations,
      lagrangian = solution$Lagrangian,
      unconstrained_solution = solution$unconstrained.solution
    ),
    class = "chisq_mixture_solution"
  )
}

fit_chisq_problem <- function(problem, td, tq, alpha = 1e-4, n_grid = 2000L,
                              ridge = 1e-8, active_tolerance = 1e-10) {
  if (!inherits(problem, "chisq_mixture_problem")) {
    stop("problem must be constructed by build_chisq_problem()")
  }

  solution <- solve_chisq_weights(
    problem$C,
    ridge = ridge,
    active_tolerance = active_tolerance
  )
  xg <- seq(tq(alpha), tq(1 - alpha), length.out = n_grid)
  hg <- td(xg)
  dictionary <- normal_dictionary(xg, problem$atoms)
  gg <- as.vector(dictionary %*% solution$weights)

  list(
    atoms = problem$atoms,
    pi = solution$weights,
    sd = problem$atoms$tau,
    mu_grid = unique(problem$atoms$mu),
    tau_grid = unique(problem$atoms$tau),
    Z = matrix(
      solution$weights,
      length(unique(problem$atoms$mu)),
      length(unique(problem$atoms$tau))
    ),
    xg = xg,
    hg = hg,
    gg = gg,
    alpha = alpha,
    criterion = "Pearson chi-square",
    C = problem$C,
    quadrature = problem$quadrature,
    solution = solution,
    active_tolerance = active_tolerance
  )
}

# User-facing convenience wrapper. It builds the atoms, constructs C once, and
# passes that same C through the solver and fit object.
fit_dictionary <- function(td, tq, scale, zero_mean = FALSE, alpha = 1e-4,
                           n_mu = 60L, n_tau = 24L, n_tau0 = 80L,
                           n_grid = 2000L, n_quad = 2048L,
                           mu_grid = NULL, tau_grid = NULL, ridge = 1e-8,
                           active_tolerance = 1e-10) {
  lo <- tq(alpha)
  hi <- tq(1 - alpha)
  if (is.null(mu_grid)) {
    mu_grid <- if (zero_mean) 0 else seq(lo, hi, length.out = n_mu)
  }
  if (is.null(tau_grid)) {
    tau_grid <- exp(seq(
      log(scale * 0.04),
      log(scale * 30),
      length.out = if (zero_mean) n_tau0 else n_tau
    ))
  }
  atoms <- expand.grid(mu = mu_grid, tau = tau_grid)

  problem <- build_chisq_problem(td, tq, atoms, n_quad = n_quad)
  fit_chisq_problem(
    problem,
    td = td,
    tq = tq,
    alpha = alpha,
    n_grid = n_grid,
    ridge = ridge,
    active_tolerance = active_tolerance
  )
}

# --- Mixture prior and closed-form posterior --------------------------------

active_mixture <- function(fit) {
  tolerance <- if (is.null(fit$active_tolerance)) 0 else fit$active_tolerance
  keep <- fit$pi > tolerance
  if (!any(keep)) {
    stop("fit has no positive mixture weights")
  }
  weights <- fit$pi[keep]
  weights <- weights / sum(weights)
  list(
    weights = weights,
    mu = fit$atoms$mu[keep],
    tau = fit$atoms$tau[keep],
    keep = keep
  )
}

dprior_mix <- function(theta, fit) {
  mixture <- active_mixture(fit)
  density <- vapply(
    seq_along(mixture$weights),
    function(k) {
      mixture$weights[k] * dnorm(theta, mixture$mu[k], mixture$tau[k])
    },
    numeric(length(theta))
  )
  rowSums(matrix(
    density,
    nrow = length(theta),
    ncol = length(mixture$weights)
  ))
}

mixture_posterior <- function(theta_hat, s, fit) {
  mixture <- active_mixture(fit)
  pi_k <- mixture$weights
  mu_k <- mixture$mu
  tau_k <- mixture$tau

  tau_post2 <- 1 / (1 / s^2 + 1 / tau_k^2)
  mu_post <- tau_post2 * (theta_hat / s^2 + mu_k / tau_k^2)

  logw <- log(pi_k) +
    dnorm(theta_hat, mu_k, sqrt(s^2 + tau_k^2), log = TRUE)
  weights <- exp(logw - max(logw))
  weights <- weights / sum(weights)

  list(w = weights, mu = mu_post, sd = sqrt(tau_post2))
}

rposterior_mix <- function(n, theta_hat, s, fit) {
  post <- mixture_posterior(theta_hat, s, fit)
  component <- sample.int(length(post$w), n, replace = TRUE, prob = post$w)
  rnorm(n, post$mu[component], post$sd[component])
}

dposterior_mix <- function(theta_grid, theta_hat, s, fit) {
  post <- mixture_posterior(theta_hat, s, fit)
  density <- vapply(
    seq_along(post$w),
    function(k) post$w[k] * dnorm(theta_grid, post$mu[k], post$sd[k]),
    numeric(length(theta_grid))
  )
  rowSums(matrix(
    density,
    nrow = length(theta_grid),
    ncol = length(post$w)
  ))
}

# Zero-mean specialization used by the Student-t shrinkage figures.
posterior_weights <- function(theta_hat, s, fit) {
  mixture <- active_mixture(fit)
  if (any(abs(mixture$mu) > sqrt(.Machine$double.eps))) {
    stop("posterior_weights() is only valid for zero-mean mixtures")
  }
  component_likelihood <- vapply(
    seq_along(mixture$weights),
    function(k) {
      mixture$weights[k] *
        dnorm(theta_hat, 0, sqrt(s^2 + mixture$tau[k]^2))
    },
    numeric(length(theta_hat))
  )
  component_likelihood <- matrix(
    component_likelihood,
    nrow = length(theta_hat),
    ncol = length(mixture$weights)
  )
  component_likelihood / rowSums(component_likelihood)
}

posterior_mean_mix <- function(theta_hat, s, fit) {
  mixture <- active_mixture(fit)
  shrinkage <- mixture$tau^2 / (s^2 + mixture$tau^2)
  rowSums(
    posterior_weights(theta_hat, s, fit) *
      outer(theta_hat, shrinkage)
  )
}

shrinkage_factor_mix <- function(theta_hat, s, fit) {
  mixture <- active_mixture(fit)
  shrinkage <- mixture$tau^2 / (s^2 + mixture$tau^2)
  as.vector(posterior_weights(theta_hat, s, fit) %*% shrinkage)
}
