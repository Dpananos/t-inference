# R helpers

`figure_1.R` in the project root is the entry point. It is a
[Rapp](https://github.com/r-lib/Rapp) CLI script, so it must be run from the
project root (the `source()` calls below are relative to it):

```sh
Rapp figure_1.R --help          # every knob, with defaults
Rapp figure_1.R                 # writes assets/figure_1.png
Rapp figure_1.R -o assets/figure_1.pdf --n-sims 1000 --quiet
```

It constructs the Pearson chi-square matrix `C` once for the Student's t
target, solves the mixture, runs simulation-based calibration against the exact
target, saves the four-panel figure, and prints the fit and SBC diagnostic
tables unless `--quiet` is passed.

The modules have separate responsibilities:

- `mixture_helpers.R`: quadrature, one-time `C` construction, simplex QP,
  KKT diagnostics, and mixture posterior calculations.
- `priors.R`: target definition and the tail-compatible normal dictionary.
- `figure_helpers.R`: figure panels and the fit diagnostic table.
- `sbc_helpers.R`: exact posterior-rank simulation and SBC plots.

