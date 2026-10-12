# ---------------------------------------------------------------------------
# Target prior and its Pearson-chi-square normal-mixture fit.
#
# The target has full support and is paired with a tail-compatible zero-mean
# normal dictionary. A Gamma target is intentionally not supported: a
# full-support normal mixture is not absolutely continuous with respect to a
# Gamma target, so Pearson chi-square divergence would be infinite.
# ---------------------------------------------------------------------------

source("R/mixture_helpers.R")

# Returns list(studentt=). The entry contains the exact target
# density/quantile/sampler and one fitted normal mixture.
make_priors <- function(sc = 0.05, t_df = 3, n_tau0 = 80L,
                        alpha = 1e-4, n_quad = 2048L, ridge = 1e-8) {
  d_t <- function(x) dt(x / sc, df = t_df) / sc
  q_t <- function(p) sc * qt(p, df = t_df)

  fit_t <- fit_dictionary(
    d_t,
    q_t,
    scale = sc,
    zero_mean = TRUE,
    n_tau0 = n_tau0,
    alpha = alpha,
    n_quad = n_quad,
    ridge = ridge
  )

  list(
    studentt = list(
      label = "Student-t",
      d = d_t,
      q = q_t,
      r = function(n) sc * rt(n, df = t_df),
      fit = fit_t
    )
  )
}
