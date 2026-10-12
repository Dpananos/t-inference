# t-inference

Code and manuscript for **Approximating Heavy-Tailed Priors With Normal Mixtures**
by A. Demetri Pananos and Tyler Buffington.

This project approximates heavy-tailed priors with finite normal mixtures for
Bayesian inference on A/B test effects under a normal likelihood. It fits mixture
weights by minimizing Pearson's chi-square divergence, with a small ridge penalty,
through convex quadratic programming. The fitted mixture supports fast posterior
inference through weighted conjugate normal updates.

The example uses a Student's t prior with 3 degrees of freedom and scale 0.05,
with 80 candidate normal components. The figure compares prior densities and
shrinkage, and checks simulation-based calibration on 10,000 datasets from the
exact target prior. The approximation applies over a practical range: a finite
normal mixture still has Gaussian tails asymptotically.

- [abstract.tex](abstract.tex): manuscript, formulas, results, and limitations.
- [figure_1.R](figure_1.R): entry point for the figure and diagnostic tables.
- [R/](R/README.md): mixture fitting, posterior calculations, plotting, and calibration helpers.
- [assets/figure_1.png](assets/figure_1.png): figure included in the manuscript.

To generate the figure, install R, the `quadprog`, `ggplot2`, `patchwork`, and `SBC`
packages, and the [Rapp](https://github.com/r-lib/Rapp) command-line tool.
Run from the project root:

```sh
Rapp figure_1.R --help
Rapp figure_1.R # writes assets/figure_1.png
```

The script also prints fit and calibration diagnostics. Use `--n-sims 1000` to
change the simulation count, `-o` to change the output path, or `--quiet` to
suppress the diagnostic tables.
