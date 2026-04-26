# Tweedie Kernel Density Estimation for Semicontinuous Data: Simulation Code

This repository contains reproducible simulation code for "Tweedie Kernel Density Estimation for Semicontinuous Data". The code supports the simulation studies reported in the accompanying manuscript and is designed to evaluate finite-sample performance, bandwidth selection behavior, and integrated error properties of kernel density estimators for non-negative, zero-inflated outcomes.

## Environment

The script uses `parallel::mclapply` and is therefore Unix-only (Linux / macOS / WSL2 on Windows).

## Dependencies

```r
install.packages(c("tweedie", "ks", "kdensity", "parallel",
                   "ggplot2", "dplyr", "tidyr", "purrr", "stringr"))
```

## Usage

Run `tkde_simulation.R`. Adjust `n`, `p0`, and `ncores` at the top of each scenario as needed.
