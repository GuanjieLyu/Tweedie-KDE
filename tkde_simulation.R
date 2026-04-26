# =============================================================================
# Tweedie Kernel Density Estimation (TKDE): Simulation Study
# -----------------------------------------------------------------------------
# Companion code for the TKDE manuscript (TKDE).
#
# This script reproduces the Monte Carlo simulation comparing the proposed
# Tweedie Kernel Density Estimator (TKDE) against three benchmarks for
# semicontinuous (zero-inflated, non-negative) data:
#
#   1. TKDE                      -- profile LSCV over (p, h)
#   2. Two-part Gaussian KDE     -- plug-in bandwidth (`ks::kde`)
#   3. Two-part Gaussian KDE     -- LSCV bandwidth   (`ks::hlscv`)
#   4. Two-part Gamma KDE        -- Chen (2000), via `kdensity`
#
# Four data-generating processes (DGPs) are considered:
#   M.1  Tweedie compound Poisson-Gamma
#   M.2  Zero-Inflated Gamma
#   M.3  Zero-Inflated Gamma Mixture (boundary spike with heavy tail)
#   M.4  Zero-Inflated Gamma Mixture (bimodal with separated modes)
#
# Each scenario is run for B = 500 replicates at sample sizes
# n in {100, 200, 500} and zero-inflation levels p0 in {0.15, 0.30, 0.45}.
# Performance is measured by Integrated Squared Error (ISE) and
# Integrated Absolute Error (IAE).
#
# -----------------------------------------------------------------------------
# ENVIRONMENT
# -----------------------------------------------------------------------------
# This script is designed to run on a Unix-like system (Linux / macOS /
# WSL2 on Windows). It relies on `parallel::mclapply` and `parallel::mcparallel`
# for forked-process parallelism, which are NOT supported on native Windows R.
# Windows users should run this script under WSL2 (tested on Ubuntu 22.04).
#
# Adjust `ncores` below to match the number of physical cores available.
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
#   tweedie, ks, kdensity, parallel
#   ggplot2, dplyr, tidyr, purrr, stringr   (post-processing / plots)
# =============================================================================

rm(list = ls())

library(tweedie)
library(ks)
library(kdensity)
library(parallel)
library(ggplot2)


# =============================================================================
# 1. GLOBAL SIMULATION SETTINGS
# =============================================================================

ncores <- 28        # number of cores for mclapply (set to detectCores() - 1 if unsure)
B      <- 500       # number of Monte Carlo replicates
n      <- 100       # sample size per replicate (run separately for n = 100, 200, 500)


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

#' Least-Squares Cross-Validation (LSCV) score for TKDE
#'
#' Computes
#'     LSCV(h) = integral fhat^2 dx  -  (2/n) * sum_i fhat_{-i}(y_i)
#' which is to be MINIMISED over h (and, in profile mode, over p).
#'
#' Note: leave-one-out log-likelihood CV is unsuitable here because the
#' Tweedie kernel integrates to 1 in y (the data argument), not in mu (the
#' evaluation point), so fhat is not a proper density in mu and its total
#' mass varies with h. LSCV targets ISE directly and avoids this pitfall.
#'
#' @param y       Numeric vector of observations (>= 0).
#' @param h       Tweedie dispersion (bandwidth) parameter.
#' @param power   Tweedie power index p in (1, 2).
#' @param x_grid  Numeric grid for numerical integration of fhat^2.
#' @param dx      Grid spacing of x_grid.
#' @return        Scalar LSCV score.
tkde_lscv <- function(y, h, power, x_grid, dx) {
  n <- length(y)

  # Term 1: integral of fhat^2 over x_grid
  fhat <- sapply(x_grid, function(x) mean(dtweedie(y, mu = x, phi = h, power = power)))
  int_fhat2 <- sum(fhat^2) * dx

  # Term 2: (2/n) * sum_i fhat_{-i}(y_i), summing only over y_i > 0
  y_pos <- y[y > 0]
  loo_sum <- sum(sapply(y_pos, function(xi) {
    full <- mean(dtweedie(y, mu = xi, phi = h, power = power))
    self <- dtweedie(xi, mu = xi, phi = h, power = power)
    (n * full - self) / (n - 1)
  }))

  int_fhat2 - (2 / n) * loo_sum
}


#' Profile LSCV bandwidth and power selection for TKDE
#'
#' For each candidate p in a grid, finds the LSCV-optimal h, then selects
#' the (p, h) pair with the smallest LSCV score.
#'
#' @param y        Numeric vector of observations (>= 0).
#' @param x_grid   Numeric grid for ISE integration.
#' @param dx       Grid spacing.
#' @param p_range  2-vector giving lower/upper bounds for the power index p.
#' @param n_p      Number of p candidates.
#' @param h_range  2-vector giving lower/upper bounds for the bandwidth h.
#' @param n_grid   Number of h candidates.
#' @return         List with `p_opt`, `h_opt`, `cv_opt`.
tkde_select_ph <- function(y, x_grid, dx,
                           p_range = c(1.1, 1.9), n_p   = 18,
                           h_range = c(0.005, 0.5), n_grid = 20) {
  p_vec <- seq(p_range[1], p_range[2], length.out = n_p)
  h_vec <- seq(h_range[1], h_range[2], length.out = n_grid)

  profile <- sapply(p_vec, function(p) {
    scores <- sapply(h_vec, function(h)
      tkde_lscv(y, h = h, power = p, x_grid = x_grid, dx = dx))
    i_best <- which.min(scores)
    c(h = h_vec[i_best], lscv = scores[i_best])
  })

  i_best <- which.min(profile["lscv", ])
  list(p_opt  = p_vec[i_best],
       h_opt  = profile["h",   i_best],
       cv_opt = profile["lscv", i_best])
}


#' TKDE evaluator on a fixed grid
tkde_grid <- function(y, x_grid, h, power) {
  sapply(x_grid, function(x) mean(dtweedie(y, mu = x, phi = h, power = power)))
}


#' Two-part Gaussian KDE with plug-in bandwidth (`ks::kde`)
ks_twopart_grid <- function(y, x_grid) {
  p0_hat <- mean(y == 0)
  y_pos  <- as.numeric(y[y > 0])
  if (length(y_pos) < 2) return(rep(0, length(x_grid)))
  kde_obj <- kde(x = y_pos, eval.points = x_grid)
  (1 - p0_hat) * as.vector(kde_obj$estimate)
}


#' Two-part Gaussian KDE with LSCV bandwidth (`ks::hlscv`)
ks_lscv_twopart_grid <- function(y, x_grid) {
  p0_hat <- mean(y == 0)
  y_pos  <- as.numeric(y[y > 0])
  if (length(y_pos) < 2) return(rep(0, length(x_grid)))
  h_lscv  <- hlscv(y_pos)
  kde_obj <- kde(x = y_pos, h = h_lscv, eval.points = x_grid)
  (1 - p0_hat) * as.vector(kde_obj$estimate)
}


#' Two-part Gamma KDE (Chen 2000) via `kdensity`
gamma_twopart_grid <- function(y, x_grid) {
  p0_hat <- mean(y == 0)
  y_pos  <- as.numeric(y[y > 0])
  if (length(y_pos) < 2) return(rep(0, length(x_grid)))
  fit <- kdensity(y_pos, kernel = "gamma")
  (1 - p0_hat) * fit(x_grid)
}


#' Parallel `lapply` with a lightweight progress indicator
#'
#' Each worker appends a token to a shared temp file; the master process
#' polls the file and prints percent-complete to the console.
#'
#' Unix-only (uses `mcparallel` / `mccollect`).
mclapply_progress <- function(X, FUN, mc.cores, label = "Running") {
  n_total   <- length(X)
  prog_file <- tempfile("progress_")
  writeLines("0", prog_file)

  wrapped_FUN <- function(x) {
    res <- FUN(x)
    con <- file(prog_file, open = "a")
    writeLines("x", con)
    close(con)
    res
  }

  job <- mcparallel({
    mclapply(X, wrapped_FUN, mc.cores = mc.cores)
  })

  while (TRUE) {
    Sys.sleep(2)
    lines <- readLines(prog_file, warn = FALSE)
    done  <- length(lines) - 1
    pct   <- min(100, round(100 * done / n_total))
    cat(sprintf("\r  [%s] %d / %d  (%d%%)", label, done, n_total, pct))
    flush.console()
    if (done >= n_total) break
    if (!is.null(mccollect(job, wait = FALSE))) break
  }
  cat("\n")
  result <- mccollect(job)[[1]]
  unlink(prog_file)
  result
}


# =============================================================================
# 3. TRUE DENSITY CURVES (FIGURE: scenario overview)
# =============================================================================
# Self-contained block: parameters defined inline so the figure can be
# produced without running any of the simulation code below.

# --- M.1: Tweedie ------------------------------------------------------------
.p0_1 <- 0.30; .p1 <- 1.1; .mu1 <- 2
.phi1 <- .mu1^(2 - .p1) / ((2 - .p1) * (-log(.p0_1)))
.xg1  <- seq(1e-4, qtweedie(0.995, power = .p1, mu = .mu1, phi = .phi1),
             length.out = 400)
.ft1  <- dtweedie(.xg1, mu = .mu1, phi = .phi1, power = .p1)

# --- M.2: Zero-Inflated Gamma ------------------------------------------------
.p0_2 <- 0.30; .sh2 <- 1.3; .ra2 <- 6
.xg2  <- seq(1e-4, qgamma(0.995, shape = .sh2, rate = .ra2), length.out = 400)
.ft2  <- (1 - .p0_2) * dgamma(.xg2, shape = .sh2, rate = .ra2)

# --- M.3: Boundary spike with heavy right tail -------------------------------
.p0_3 <- 0.30
.xg3  <- seq(1e-4, qgamma(0.995, shape = 15, rate = 1), length.out = 400)
.ft3  <- (1 - .p0_3) * (0.55 * dgamma(.xg3, shape =  2, rate = 6) +
                        0.45 * dgamma(.xg3, shape = 15, rate = 1))

# --- M.4: Bimodal with separated modes ---------------------------------------
# Component 1: Gamma(4, 6)   mode = 0.5  (small, low peak)
# Component 2: Gamma(20, 3)  mode = 6.3  (dominant tall peak)
.p0_4 <- 0.30
.xg4  <- seq(1e-4, qgamma(0.995, shape = 20, rate = 3), length.out = 400)
.ft4  <- (1 - .p0_4) * (0.35 * dgamma(.xg4, shape =  4, rate = 6) +
                        0.65 * dgamma(.xg4, shape = 20, rate = 3))

true_curves <- rbind(
  data.frame(x = .xg1, density = .ft1, Scenario = "M.1: Tweedie"),
  data.frame(x = .xg2, density = .ft2, Scenario = "M.2: ZI-Gamma"),
  data.frame(x = .xg3, density = .ft3,
             Scenario = "M.3: ZI-Gamma Mixture\n(boundary spike with heavy tail)"),
  data.frame(x = .xg4, density = .ft4,
             Scenario = "M.4: ZI-Gamma Mixture\n(bimodal with separated modes)")
)
true_curves$Scenario <- factor(
  true_curves$Scenario,
  levels = c("M.1: Tweedie",
             "M.2: ZI-Gamma",
             "M.3: ZI-Gamma Mixture\n(boundary spike with heavy tail)",
             "M.4: ZI-Gamma Mixture\n(bimodal with separated modes)")
)

ggplot(true_curves, aes(x = x, y = density, colour = Scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Scenario, nrow = 2, ncol = 2, scales = "free") +
  labs(x = "x", y = "Density") +
  guides(colour = "none") +
  theme_classic(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text       = element_text(size = 11),
        panel.border     = element_rect(color = "black", fill = NA))

# Clean up the temporary objects used only for the figure
rm(.p1, .mu1, .p0_1, .phi1, .xg1, .ft1,
   .p0_2, .sh2, .ra2, .xg2, .ft2,
   .p0_3, .xg3, .ft3,
   .p0_4, .xg4, .ft4)


# =============================================================================
# 4. SCENARIO M.1 -- TWEEDIE DGP
# =============================================================================

# ---- DGP parameters ---------------------------------------------------------
# Run separately for each zero-inflation level by uncommenting one of:
# p0_target <- 0.15
# p0_target <- 0.30
p0_target <- 0.45

p_true   <- 1.1
mu_true  <- 2
phi_true <- mu_true^(2 - p_true) / ((2 - p_true) * (-log(p0_target)))

r_tweed <- function(n) rtweedie(n, mu = mu_true, phi = phi_true, power = p_true)

# ---- Data generation --------------------------------------------------------
set.seed(123)
Y_list <- replicate(B, r_tweed(n), simplify = FALSE)

x_max  <- qtweedie(0.995, power = p_true, mu = mu_true, phi = phi_true)
x_grid <- seq(1e-4, x_max, length.out = 200)
dx     <- x_grid[2] - x_grid[1]

f_true <- function(x) dtweedie(x, mu = mu_true, phi = phi_true, power = p_true)
ft     <- f_true(x_grid)

# ---- Worker ----------------------------------------------------------------
run_one_rep <- function(b) {
  y <- Y_list[[b]]

  cv_res    <- tkde_select_ph(y, x_grid = x_grid, dx = dx)
  fhat_tkde <- tkde_grid(y, x_grid, h = cv_res$h_opt, power = cv_res$p_opt)
  ise_t     <- sum((fhat_tkde - ft)^2) * dx

  fhat_ks      <- ks_twopart_grid(y, x_grid)
  ise_k        <- sum((fhat_ks - ft)^2) * dx

  fhat_ks_lscv <- ks_lscv_twopart_grid(y, x_grid)
  ise_k_lscv   <- sum((fhat_ks_lscv - ft)^2) * dx

  fhat_gam <- gamma_twopart_grid(y, x_grid)
  ise_gam  <- sum((fhat_gam - ft)^2) * dx

  iae_t      <- sum(abs(fhat_tkde    - ft)) * dx
  iae_k      <- sum(abs(fhat_ks      - ft)) * dx
  iae_k_lscv <- sum(abs(fhat_ks_lscv - ft)) * dx
  iae_gam    <- sum(abs(fhat_gam     - ft)) * dx

  c(p_cv = cv_res$p_opt, h_cv = cv_res$h_opt,
    ise_tkde = ise_t, ise_ks = ise_k, ise_ks_lscv = ise_k_lscv, ise_gamma = ise_gam,
    iae_tkde = iae_t, iae_ks = iae_k, iae_ks_lscv = iae_k_lscv, iae_gamma = iae_gam)
}

# ---- Run --------------------------------------------------------------------
cat(sprintf("\n=== Scenario M.1: Tweedie DGP (p0 = %.2f) ===\n", p0_target))
cat(sprintf("Running %d replicates on %d cores (n = %d) ...\n", B, ncores, n))
t0      <- proc.time()
results <- mclapply_progress(seq_len(B), run_one_rep,
                             mc.cores = ncores, label = "Scenario M.1")
res_mat <- do.call(rbind, results)
saveRDS(res_mat, file = sprintf("S1_%d_%s.rds", n, format(p0_target, nsmall = 2)))
elapsed <- (proc.time() - t0)[3]

# ---- Summary ----------------------------------------------------------------
res_tab <- data.frame(
  Method = c("TKDE (LSCV)", "Gauss KDE plug-in (ks)", "Gauss KDE LSCV (ks)",
             "Gamma KDE (kdensity)"),
  MISE   = c(mean(res_mat[, "ise_tkde"]),    mean(res_mat[, "ise_ks"]),
             mean(res_mat[, "ise_ks_lscv"]), mean(res_mat[, "ise_gamma"])),
  SD     = c(sd(res_mat[, "ise_tkde"]),    sd(res_mat[, "ise_ks"]),
             sd(res_mat[, "ise_ks_lscv"]), sd(res_mat[, "ise_gamma"])),
  MIAE   = c(mean(res_mat[, "iae_tkde"]),    mean(res_mat[, "iae_ks"]),
             mean(res_mat[, "iae_ks_lscv"]), mean(res_mat[, "iae_gamma"])),
  SD_2   = c(sd(res_mat[, "iae_tkde"]),    sd(res_mat[, "iae_ks"]),
             sd(res_mat[, "iae_ks_lscv"]), sd(res_mat[, "iae_gamma"]))
)
cat("\n"); print(res_tab, digits = 4)
cat(sprintf("\nTotal time: %.1f seconds\n", elapsed))
cat("\nSelected p distribution:\n"); print(table(res_mat[, "p_cv"]))
cat("\nSelected h distribution:\n"); print(summary(res_mat[, "h_cv.h"]))


# =============================================================================
# 5. SCENARIO M.2 -- ZERO-INFLATED GAMMA
# =============================================================================

# ---- DGP parameters ---------------------------------------------------------
# p0_zi <- 0.15
# p0_zi <- 0.30
p0_zi    <- 0.45
shape_zi <- 1.3
rate_zi  <- 6.0

f_true_zi <- function(x) (1 - p0_zi) * dgamma(x, shape = shape_zi, rate = rate_zi)

r_zigamma <- function(n) {
  u <- rbinom(n, 1, 1 - p0_zi)
  u * rgamma(n, shape = shape_zi, rate = rate_zi)
}

# ---- Data generation --------------------------------------------------------
set.seed(456)
Y_list_zi <- replicate(B, r_zigamma(n), simplify = FALSE)

x_max_zi  <- qgamma(0.995, shape = shape_zi, rate = rate_zi)
x_grid_zi <- seq(1e-4, x_max_zi, length.out = 200)
dx_zi     <- x_grid_zi[2] - x_grid_zi[1]
ft_zi     <- f_true_zi(x_grid_zi)

# ---- Worker -----------------------------------------------------------------
run_one_rep_zi <- function(b) {
  y <- Y_list_zi[[b]]

  ph        <- tkde_select_ph(y, x_grid = x_grid_zi, dx = dx_zi)
  fhat_tkde <- tkde_grid(y, x_grid_zi, h = ph$h_opt, power = ph$p_opt)
  ise_t     <- sum((fhat_tkde - ft_zi)^2) * dx_zi

  fhat_ks      <- ks_twopart_grid(y, x_grid_zi)
  ise_k        <- sum((fhat_ks - ft_zi)^2) * dx_zi

  fhat_ks_lscv <- ks_lscv_twopart_grid(y, x_grid_zi)
  ise_k_lscv   <- sum((fhat_ks_lscv - ft_zi)^2) * dx_zi

  fhat_gam <- gamma_twopart_grid(y, x_grid_zi)
  ise_gam  <- sum((fhat_gam - ft_zi)^2) * dx_zi

  iae_t      <- sum(abs(fhat_tkde    - ft_zi)) * dx_zi
  iae_k      <- sum(abs(fhat_ks      - ft_zi)) * dx_zi
  iae_k_lscv <- sum(abs(fhat_ks_lscv - ft_zi)) * dx_zi
  iae_gam    <- sum(abs(fhat_gam     - ft_zi)) * dx_zi

  c(p_cv = ph$p_opt, h_cv = ph$h_opt,
    ise_tkde = ise_t, ise_ks = ise_k, ise_ks_lscv = ise_k_lscv, ise_gamma = ise_gam,
    iae_tkde = iae_t, iae_ks = iae_k, iae_ks_lscv = iae_k_lscv, iae_gamma = iae_gam)
}

# ---- Run --------------------------------------------------------------------
cat(sprintf("\n\n=== Scenario M.2: Zero-Inflated Gamma (p0=%.2f, shape=%.1f, rate=%.1f) ===\n",
            p0_zi, shape_zi, rate_zi))
cat(sprintf("Running %d replicates on %d cores (n = %d) ...\n", B, ncores, n))
t0         <- proc.time()
results_zi <- mclapply_progress(seq_len(B), run_one_rep_zi,
                                mc.cores = ncores, label = "Scenario M.2")
res_mat_zi <- do.call(rbind, results_zi)
saveRDS(res_mat_zi, file = sprintf("S2_%d_%s.rds", n, format(p0_zi, nsmall = 2)))
elapsed_zi <- (proc.time() - t0)[3]

# ---- Summary ----------------------------------------------------------------
res_tab_zi <- data.frame(
  Method = c("TKDE (profile LSCV)", "Gauss KDE plug-in (ks)", "Gauss KDE LSCV (ks)",
             "Gamma KDE (kdensity)"),
  MISE   = c(mean(res_mat_zi[, "ise_tkde"]),    mean(res_mat_zi[, "ise_ks"]),
             mean(res_mat_zi[, "ise_ks_lscv"]), mean(res_mat_zi[, "ise_gamma"])),
  SD     = c(sd(res_mat_zi[, "ise_tkde"]),    sd(res_mat_zi[, "ise_ks"]),
             sd(res_mat_zi[, "ise_ks_lscv"]), sd(res_mat_zi[, "ise_gamma"])),
  MIAE   = c(mean(res_mat_zi[, "iae_tkde"]),    mean(res_mat_zi[, "iae_ks"]),
             mean(res_mat_zi[, "iae_ks_lscv"]), mean(res_mat_zi[, "iae_gamma"]))
)
cat("\n"); print(res_tab_zi, digits = 4)
cat("\nSelected p distribution:\n"); print(table(res_mat_zi[, "p_cv"]))
cat("\nSelected h distribution:\n"); print(summary(res_mat_zi[, "h_cv.h"]))
cat(sprintf("Total time: %.1f seconds\n", elapsed_zi))


# =============================================================================
# 6. SCENARIO M.3 -- ZI GAMMA MIXTURE (BOUNDARY SPIKE + HEAVY TAIL)
# =============================================================================

# ---- DGP parameters ---------------------------------------------------------
# p0_mix <- 0.15
# p0_mix <- 0.30
p0_mix <- 0.45

dmix <- function(x) 0.55 * dgamma(x, shape =  2, rate = 6) +
                    0.45 * dgamma(x, shape = 15, rate = 1)

f_true_mix <- function(x) (1 - p0_mix) * dmix(x)

r_zimix <- function(n) {
  u    <- rbinom(n, 1, 1 - p0_mix)
  comp <- rbinom(n, 1, 0.55)
  u * ifelse(comp,
             rgamma(n, shape =  2, rate = 6),
             rgamma(n, shape = 15, rate = 1))
}

# ---- Data generation --------------------------------------------------------
set.seed(789)
Y_list_mix <- replicate(B, r_zimix(n), simplify = FALSE)

x_max_mix  <- qgamma(0.995, shape = 15, rate = 1)
x_grid_mix <- seq(1e-4, x_max_mix, length.out = 200)
dx_mix     <- x_grid_mix[2] - x_grid_mix[1]
ft_mix     <- f_true_mix(x_grid_mix)

# ---- Worker -----------------------------------------------------------------
run_one_rep_mix <- function(b) {
  y <- Y_list_mix[[b]]

  ph        <- tkde_select_ph(y, x_grid = x_grid_mix, dx = dx_mix)
  fhat_tkde <- tkde_grid(y, x_grid_mix, h = ph$h_opt, power = ph$p_opt)
  ise_t     <- sum((fhat_tkde - ft_mix)^2) * dx_mix

  fhat_ks      <- ks_twopart_grid(y, x_grid_mix)
  ise_k        <- sum((fhat_ks - ft_mix)^2) * dx_mix

  fhat_ks_lscv <- ks_lscv_twopart_grid(y, x_grid_mix)
  ise_k_lscv   <- sum((fhat_ks_lscv - ft_mix)^2) * dx_mix

  fhat_gam <- gamma_twopart_grid(y, x_grid_mix)
  ise_gam  <- sum((fhat_gam - ft_mix)^2) * dx_mix

  iae_t      <- sum(abs(fhat_tkde    - ft_mix)) * dx_mix
  iae_k      <- sum(abs(fhat_ks      - ft_mix)) * dx_mix
  iae_k_lscv <- sum(abs(fhat_ks_lscv - ft_mix)) * dx_mix
  iae_gam    <- sum(abs(fhat_gam     - ft_mix)) * dx_mix

  c(p_cv = ph$p_opt, h_cv = ph$h_opt,
    ise_tkde = ise_t, ise_ks = ise_k, ise_ks_lscv = ise_k_lscv, ise_gamma = ise_gam,
    iae_tkde = iae_t, iae_ks = iae_k, iae_ks_lscv = iae_k_lscv, iae_gamma = iae_gam)
}

# ---- Run --------------------------------------------------------------------
cat(sprintf("\n\n=== Scenario M.3: ZI Gamma Mixture (p0 = %.2f, boundary + heavy tail) ===\n",
            p0_mix))
cat(sprintf("Running %d replicates on %d cores (n = %d) ...\n", B, ncores, n))
t0          <- proc.time()
results_mix <- mclapply_progress(seq_len(B), run_one_rep_mix,
                                 mc.cores = ncores, label = "Scenario M.3")
res_mat_mix <- do.call(rbind, results_mix)
saveRDS(results_mix, file = sprintf("S3_%d_%s.rds", n, format(p0_mix, nsmall = 2)))
elapsed_mix <- (proc.time() - t0)[3]

# ---- Summary ----------------------------------------------------------------
res_tab_mix <- data.frame(
  Method = c("TKDE (profile LSCV)", "Gauss KDE plug-in (ks)", "Gauss KDE LSCV (ks)",
             "Gamma KDE (kdensity)"),
  MISE   = c(mean(res_mat_mix[, "ise_tkde"]),    mean(res_mat_mix[, "ise_ks"]),
             mean(res_mat_mix[, "ise_ks_lscv"]), mean(res_mat_mix[, "ise_gamma"])),
  SD     = c(sd(res_mat_mix[, "ise_tkde"]),    sd(res_mat_mix[, "ise_ks"]),
             sd(res_mat_mix[, "ise_ks_lscv"]), sd(res_mat_mix[, "ise_gamma"])),
  MIAE   = c(mean(res_mat_mix[, "iae_tkde"]),    mean(res_mat_mix[, "iae_ks"]),
             mean(res_mat_mix[, "iae_ks_lscv"]), mean(res_mat_mix[, "iae_gamma"]))
)
cat("\n"); print(res_tab_mix, digits = 4)
cat("\nSelected p distribution:\n"); print(table(res_mat_mix[, "p_cv"]))
cat("\nSelected h distribution:\n"); print(summary(res_mat_mix[, "h_cv.h"]))
cat(sprintf("Total time: %.1f seconds\n", elapsed_mix))


# =============================================================================
# 7. SCENARIO M.4 -- ZI GAMMA MIXTURE (BIMODAL, MODES AWAY FROM ZERO)
# =============================================================================
#
# This scenario has two well-separated modes, neither of which sits at the
# boundary. Empirically we found the LSCV surface to be more multimodal in
# (p, h) here, so we override `tkde_select_ph` to use a finer grid.
# Larger (n_p, n_grid) values increase compute cost roughly multiplicatively.
# =============================================================================

# ---- Finer profile-LSCV grid for this scenario ------------------------------
tkde_select_ph <- function(y, x_grid, dx,
                           p_range = c(1.1, 1.9), n_p   = 27,   # 9, 18, or 27
                           h_range = c(0.005, 0.5), n_grid = 40) {  # 20, 40, or 60
  p_vec <- seq(p_range[1], p_range[2], length.out = n_p)
  h_vec <- seq(h_range[1], h_range[2], length.out = n_grid)

  profile <- sapply(p_vec, function(p) {
    scores <- sapply(h_vec, function(h)
      tkde_lscv(y, h = h, power = p, x_grid = x_grid, dx = dx))
    i_best <- which.min(scores)
    c(h = h_vec[i_best], lscv = scores[i_best])
  })

  i_best <- which.min(profile["lscv", ])
  list(p_opt  = p_vec[i_best],
       h_opt  = profile["h",   i_best],
       cv_opt = profile["lscv", i_best])
}

# ---- DGP parameters ---------------------------------------------------------
# Component 1: Gamma(4, 6)   -> mode = 0.5  (small peak, wide spread)
# Component 2: Gamma(20, 3)  -> mode = 6.3  (tall dominant peak)
# p0_mix2 <- 0.15
# p0_mix2 <- 0.30
p0_mix2 <- 0.45

dmix2 <- function(x) 0.35 * dgamma(x, shape =  4, rate = 6) +
                     0.65 * dgamma(x, shape = 20, rate = 3)

f_true_mix2 <- function(x) (1 - p0_mix2) * dmix2(x)

r_zimix2 <- function(n) {
  u    <- rbinom(n, 1, 1 - p0_mix2)
  comp <- rbinom(n, 1, 0.35)
  u * ifelse(comp,
             rgamma(n, shape =  4, rate = 6),
             rgamma(n, shape = 20, rate = 3))
}

# ---- Data generation --------------------------------------------------------
set.seed(555)
Y_list_mix2 <- replicate(B, r_zimix2(n), simplify = FALSE)

x_max_mix2  <- qgamma(0.995, shape = 20, rate = 3)
x_grid_mix2 <- seq(1e-4, x_max_mix2, length.out = 200)
dx_mix2     <- x_grid_mix2[2] - x_grid_mix2[1]
ft_mix2     <- f_true_mix2(x_grid_mix2)

# ---- Worker -----------------------------------------------------------------
run_one_rep_mix2 <- function(b) {
  y <- Y_list_mix2[[b]]

  ph        <- tkde_select_ph(y, x_grid = x_grid_mix2, dx = dx_mix2)
  fhat_tkde <- tkde_grid(y, x_grid_mix2, h = ph$h_opt, power = ph$p_opt)
  ise_t     <- sum((fhat_tkde - ft_mix2)^2) * dx_mix2

  fhat_ks      <- ks_twopart_grid(y, x_grid_mix2)
  ise_k        <- sum((fhat_ks - ft_mix2)^2) * dx_mix2

  fhat_ks_lscv <- ks_lscv_twopart_grid(y, x_grid_mix2)
  ise_k_lscv   <- sum((fhat_ks_lscv - ft_mix2)^2) * dx_mix2

  fhat_gam <- gamma_twopart_grid(y, x_grid_mix2)
  ise_gam  <- sum((fhat_gam - ft_mix2)^2) * dx_mix2

  iae_t      <- sum(abs(fhat_tkde    - ft_mix2)) * dx_mix2
  iae_k      <- sum(abs(fhat_ks      - ft_mix2)) * dx_mix2
  iae_k_lscv <- sum(abs(fhat_ks_lscv - ft_mix2)) * dx_mix2
  iae_gam    <- sum(abs(fhat_gam     - ft_mix2)) * dx_mix2

  c(p_cv = ph$p_opt, h_cv = ph$h_opt,
    ise_tkde = ise_t, ise_ks = ise_k, ise_ks_lscv = ise_k_lscv, ise_gamma = ise_gam,
    iae_tkde = iae_t, iae_ks = iae_k, iae_ks_lscv = iae_k_lscv, iae_gamma = iae_gam)
}

# ---- Run --------------------------------------------------------------------
cat(sprintf("\n\n=== Scenario M.4: ZI Gamma Mixture (p0 = %.2f, separated modes) ===\n",
            p0_mix2))
cat(sprintf("Running %d replicates on %d cores (n = %d) ...\n", B, ncores, n))
t0           <- proc.time()
results_mix2 <- mclapply_progress(seq_len(B), run_one_rep_mix2,
                                  mc.cores = ncores, label = "Scenario M.4")
res_mat_mix2 <- do.call(rbind, results_mix2)
saveRDS(results_mix2,
        file = sprintf("S4_%d_%s_p27h40.rds", n, format(p0_mix2, nsmall = 2)))
elapsed_mix2 <- (proc.time() - t0)[3]

# ---- Summary ----------------------------------------------------------------
res_tab_mix2 <- data.frame(
  Method = c("TKDE (profile LSCV)", "Gauss KDE plug-in (ks)", "Gauss KDE LSCV (ks)",
             "Gamma KDE (kdensity)"),
  MISE   = c(mean(res_mat_mix2[, "ise_tkde"]),    mean(res_mat_mix2[, "ise_ks"]),
             mean(res_mat_mix2[, "ise_ks_lscv"]), mean(res_mat_mix2[, "ise_gamma"])),
  SD     = c(sd(res_mat_mix2[, "ise_tkde"]),    sd(res_mat_mix2[, "ise_ks"]),
             sd(res_mat_mix2[, "ise_ks_lscv"]), sd(res_mat_mix2[, "ise_gamma"])),
  MIAE   = c(mean(res_mat_mix2[, "iae_tkde"]),    mean(res_mat_mix2[, "iae_ks"]),
             mean(res_mat_mix2[, "iae_ks_lscv"]), mean(res_mat_mix2[, "iae_gamma"]))
)
cat("\n"); print(res_tab_mix2, digits = 4)
cat("\nSelected p distribution:\n"); print(table(res_mat_mix2[, "p_cv"]))
cat("\nSelected h distribution:\n"); print(summary(res_mat_mix2[, "h_cv.h"]))
cat(sprintf("Total time: %.1f seconds\n", elapsed_mix2))


# =============================================================================
# 8. POST-PROCESSING -- BOXPLOTS OF ISE / IAE
# =============================================================================
# Run AFTER all .rds result files have been produced (across n and p0).
# Recommended export size: ~700 x 500 px.
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# -----------------------------------------------------------------------------
# 8.1  Helper: load all result files for a given scenario tag (S1/S2/S3/S4)
# -----------------------------------------------------------------------------
# Scenarios M.1 / M.2 stored as a matrix; M.3 / M.4 stored as a list -- handle both.
load_scenario_results <- function(scenario_tag, dir_path = ".",
                                  is_list_format = FALSE) {
  pattern <- sprintf("^%s_\\d+_0\\.(15|3|45)\\.rds$", scenario_tag)
  files   <- list.files(dir_path, pattern = pattern, full.names = TRUE)

  purrr::map_dfr(files, function(f) {
    x <- readRDS(f)
    if (is_list_format) x <- as.data.frame(do.call(rbind, x)) else x <- as.data.frame(x)

    bname  <- basename(f)
    m      <- str_match(bname, sprintf("^%s_(\\d+)_([0-9.]+)\\.rds$", scenario_tag))
    x$SampleSize <- as.integer(m[2])
    x$p0         <- as.numeric(m[3])
    x$Scenario   <- sub("S", "M.", scenario_tag)
    x
  })
}

# -----------------------------------------------------------------------------
# 8.2  Helper: build long plot data and shared ggplot theme
# -----------------------------------------------------------------------------
build_plot_df <- function(sim_data, metric = c("ISE", "IAE")) {
  metric <- match.arg(metric)
  prefix <- if (metric == "ISE") "ise_" else "iae_"

  sim_data %>%
    transmute(
      SampleSize = factor(SampleSize, levels = c(100, 200, 500)),
      p0         = factor(p0),
      Tweedie       = .data[[paste0(prefix, "tkde")]],
      `KDE plug-in` = .data[[paste0(prefix, "ks")]],
      `KDE LSCV`    = .data[[paste0(prefix, "ks_lscv")]],
      `Gamma KDE`   = .data[[paste0(prefix, "gamma")]]
    ) %>%
    pivot_longer(cols = c(Tweedie, `KDE plug-in`, `KDE LSCV`, `Gamma KDE`),
                 names_to = "Model", values_to = metric) %>%
    mutate(
      Model = factor(Model,
                     levels = c("Tweedie", "KDE plug-in", "KDE LSCV", "Gamma KDE")),
      p0    = factor(paste0("p[0]==", p0),
                     levels = paste0("p[0]==", sort(unique(as.numeric(as.character(p0))))))
    )
}

boxplot_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title        = element_text(size = 18, hjust = 0.5, face = "bold"),
    legend.text       = element_text(size = 19),
    legend.title      = element_text(size = 19),
    legend.key.size   = unit(0.6, "cm"),
    panel.spacing     = unit(1.1, "lines"),
    strip.text.x      = element_text(size = 21),
    axis.text         = element_text(size = 21),
    axis.title        = element_text(size = 21),
    axis.ticks        = element_line(linewidth = 0.7),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border      = element_rect(color = "black", fill = NA),
    axis.line         = element_line(color = "black"),
    strip.background  = element_rect(fill = "grey90", color = "black")
  )

# -----------------------------------------------------------------------------
# 8.3  ISE boxplot, M.1 (Tweedie)
# -----------------------------------------------------------------------------
sim_data <- load_scenario_results("S1", is_list_format = FALSE)
plot_df  <- build_plot_df(sim_data, metric = "ISE")

ggplot(plot_df, aes(x = SampleSize, y = ISE, fill = Model)) +
  geom_boxplot(outlier.size = 0.6, width = 0.75,
               position = position_dodge(width = 0.8)) +
  facet_grid(. ~ p0, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 0.025)) +
  labs(title = "M.1: Tweedie", x = "Sample size", y = "ISE", fill = "Model") +
  boxplot_theme + theme(legend.position = "none")

# -----------------------------------------------------------------------------
# 8.4  ISE boxplot, M.3 / M.4 (mixture scenarios stored in list format)
# -----------------------------------------------------------------------------
# Switch tag to "S4" for M.4.
sim_data <- load_scenario_results("S3", is_list_format = TRUE)
plot_df  <- build_plot_df(sim_data, metric = "ISE")

ggplot(plot_df, aes(x = SampleSize, y = ISE, fill = Model)) +
  geom_boxplot(outlier.size = 0.6, width = 0.75,
               position = position_dodge(width = 0.8)) +
  facet_grid(. ~ p0, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 0.3)) +
  labs(title = "M.3: ZI-Gamma Mixture\n(bimodal with dominant mode away from zero)",
       x = "Sample size", y = "ISE", fill = "Model") +
  boxplot_theme + theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# 8.5  IAE boxplot, M.2 (ZI-Gamma)
# -----------------------------------------------------------------------------
sim_data <- load_scenario_results("S2", is_list_format = FALSE)
plot_df  <- build_plot_df(sim_data, metric = "IAE")

ggplot(plot_df, aes(x = SampleSize, y = IAE, fill = Model)) +
  geom_boxplot(outlier.size = 0.6, width = 0.75,
               position = position_dodge(width = 0.8)) +
  facet_grid(. ~ p0, labeller = label_parsed) +
  coord_cartesian(ylim = c(0, 0.4)) +
  labs(title = "M.2: ZI-Gamma", x = "Sample size", y = "IAE", fill = "Model") +
  boxplot_theme + theme(legend.position = "none")

# -----------------------------------------------------------------------------
# 8.6  IAE boxplot, M.4 (separated modes)
# -----------------------------------------------------------------------------
sim_data <- load_scenario_results("S4", is_list_format = TRUE)
plot_df  <- build_plot_df(sim_data, metric = "IAE")

ggplot(plot_df, aes(x = SampleSize, y = IAE, fill = Model)) +
  geom_boxplot(outlier.size = 0.6, width = 0.75,
               position = position_dodge(width = 0.8)) +
  facet_grid(. ~ p0, labeller = label_parsed) +
  labs(title = "M.4: ZI-Gamma Mixture\n(bimodal with separated modes)",
       x = "Sample size", y = "IAE", fill = "Model") +
  boxplot_theme + theme(legend.position = "bottom")

# =============================================================================
# End of script
# =============================================================================
