# =============================================================================
# Monte Carlo Study: Parametric Bootstrap Estimation
# Standard Normal Distribution N(mu = 0, sigma = 1)
#
# Nested simulation structure:
#
#   for each sample size n in {20, 100, 500}:
#     for each of R = 1000 Monte Carlo replications:
#       1. Take the paired original sample x of size n from N(true_mu,
#          true_sigma) -- the same x used by mom.R / MLE.R for this
#          replication (common random numbers, see section 0b).
#       2. Compute the plug-in estimates (mu_hat, sigma_hat) from x -- the
#          MLE/MOM estimator, which coincides for the Normal distribution
#          (see mom.R / MLE.R).
#       3. Parametric bootstrap: treating (mu_hat, sigma_hat) as the "true"
#          parameters, draw B = 200 bootstrap samples of size n from
#          N(mu_hat, sigma_hat) and compute the estimator on each of them.
#       4. Bootstrap bias-corrected point estimates:
#            theta_BC = 2 * theta_hat - mean(theta_hat*)
#       5. Build a 95% percentile-method CI for mu and sigma from the B
#          bootstrap estimates: [quantile(2.5%), quantile(97.5%)].
#       6. Record whether the *true* parameter falls inside that CI
#          (coverage) and the CI length.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

rm(list = ls())
set.seed(2026)                     # full reproducibility, as specified

# Output directory for plots / tables (created if it does not exist yet)
out_dir <- "boots_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# True parameter values (Standard Normal)
true_mu    <- 0
true_sigma <- 1

# Simulation design
n_values  <- c(20, 100, 500)       # sample sizes
R         <- 1000                  # Monte Carlo replications per sample size
B         <- 200                   # parametric bootstrap resamples per replication
conf_level <- 0.95
alpha      <- 1 - conf_level

# -----------------------------------------------------------------------------
# 0b. Paired Monte Carlo samples (common random numbers)
# -----------------------------------------------------------------------------

set.seed(2026)
sample_bank <- lapply(n_values, function(n) {
  matrix(rnorm(n * R, mean = true_mu, sd = true_sigma), nrow = R, ncol = n, byrow = TRUE)
})
names(sample_bank) <- paste0("n_", n_values)

# -----------------------------------------------------------------------------
# 1. Point estimator (MLE = MOM for the Normal distribution)
# -----------------------------------------------------------------------------

#' Vectorized point estimates of mu and sigma for one or many samples at once
#'
#' @param mat a matrix with samples in columns (n rows x k columns)
#' @return list(mu_hat, sigma_hat), each of length k
estimate_params_mat <- function(mat) {
  m1 <- colMeans(mat)
  m2 <- colMeans(mat^2)
  mu_hat    <- m1
  # pmax guards against tiny negative values from floating-point round-off
  sigma_hat <- sqrt(pmax(m2 - m1^2, 0))
  list(mu_hat = mu_hat, sigma_hat = sigma_hat)
}

# -----------------------------------------------------------------------------
# 2. Parametric bootstrap Monte Carlo simulation
# -----------------------------------------------------------------------------

#' Run the nested Monte Carlo / parametric bootstrap simulation for one n
#'
#' @param samples matrix of paired samples (R rows = replications, n columns)
#' @param B number of parametric bootstrap resamples per replication
#' @return list of pre-allocated result vectors (one entry per replication)
run_bootstrap_mc <- function(samples, B) {

  n <- ncol(samples)
  R <- nrow(samples)

  # Pre-allocate all result containers up front
  mu_hat_vec        <- numeric(R)
  sigma_hat_vec     <- numeric(R)
  ci_mu_lower       <- numeric(R)
  ci_mu_upper       <- numeric(R)
  ci_sigma_lower    <- numeric(R)
  ci_sigma_upper    <- numeric(R)
  coverage_mu       <- logical(R)
  coverage_sigma    <- logical(R)

  for (m in seq_len(R)) {

    # --- 1. Paired original sample for this Monte Carlo replication --------
    x <- samples[m, ]

    m1 <- mean(x)
    m2 <- mean(x^2)
    mu_hat    <- m1
    sigma_hat <- sqrt(max(m2 - m1^2, 0))

    # --- 2. Parametric bootstrap: B resamples drawn in one shot -------------
    boot_mat <- matrix(rnorm(n * B, mean = mu_hat, sd = sigma_hat),
                        nrow = n, ncol = B)

    boot_est <- estimate_params_mat(boot_mat)   # vectorized over all B columns

    # --- 3. Bootstrap bias-corrected point estimates -------------------------
    mu_bc    <- 2 * mu_hat    - mean(boot_est$mu_hat,    na.rm = TRUE)
    sigma_bc <- 2 * sigma_hat - mean(boot_est$sigma_hat, na.rm = TRUE)
    if (!is.finite(sigma_bc) || sigma_bc <= 0) sigma_bc <- sigma_hat

    mu_hat_vec[m]    <- mu_bc
    sigma_hat_vec[m] <- sigma_bc

    # --- 4. 95% percentile-method confidence intervals -----------------------
    ci_mu    <- quantile(boot_est$mu_hat,    probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    ci_sigma <- quantile(boot_est$sigma_hat, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

    ci_mu_lower[m]    <- ci_mu[1]
    ci_mu_upper[m]    <- ci_mu[2]
    ci_sigma_lower[m] <- ci_sigma[1]
    ci_sigma_upper[m] <- ci_sigma[2]

    # --- 5. Coverage: does the TRUE parameter fall inside this CI? -----------
    coverage_mu[m]    <- (true_mu    >= ci_mu_lower[m])    && (true_mu    <= ci_mu_upper[m])
    coverage_sigma[m] <- (true_sigma >= ci_sigma_lower[m]) && (true_sigma <= ci_sigma_upper[m])
  }

  list(
    mu_hat = mu_hat_vec, sigma_hat = sigma_hat_vec,
    ci_mu_lower = ci_mu_lower, ci_mu_upper = ci_mu_upper,
    ci_sigma_lower = ci_sigma_lower, ci_sigma_upper = ci_sigma_upper,
    coverage_mu = coverage_mu, coverage_sigma = coverage_sigma
  )
}

# Run the simulation for every sample size and store results in a list
sim_results <- setNames(
  lapply(sample_bank, run_bootstrap_mc, B = B),
  paste0("n_", n_values)
)

# Combine point estimates into one long data frame for plotting
sim_long <- do.call(rbind, lapply(n_values, function(n) {
  res <- sim_results[[paste0("n_", n)]]
  data.frame(n = n, mu_hat = res$mu_hat, sigma_hat = res$sigma_hat)
}))
sim_long$n <- factor(sim_long$n, levels = n_values)

# -----------------------------------------------------------------------------
# 3. Performance metrics
# -----------------------------------------------------------------------------

performance_metrics <- function(estimates, truth) {
  bias     <- mean(estimates, na.rm = TRUE) - truth
  rel_bias <- if (truth != 0) bias / truth else NA   # undefined when truth = 0
  variance <- var(estimates, na.rm = TRUE)
  mse      <- mean((estimates - truth)^2, na.rm = TRUE)
  rel_mse  <- if (truth != 0) mse / truth^2 else NA  # undefined when truth = 0

  c(Bias = bias, RelBias = rel_bias, Variance = variance,
    MSE = mse, RelMSE = rel_mse)
}

interval_metrics <- function(coverage, lower, upper) {
  c(CoverageProbability = mean(coverage, na.rm = TRUE),
    AvgCILength         = mean(upper - lower, na.rm = TRUE))
}

results_table <- do.call(rbind, lapply(n_values, function(n) {
  res <- sim_results[[paste0("n_", n)]]

  mu_point <- performance_metrics(res$mu_hat, true_mu)
  mu_int   <- interval_metrics(res$coverage_mu, res$ci_mu_lower, res$ci_mu_upper)

  sigma_point <- performance_metrics(res$sigma_hat, true_sigma)
  sigma_int   <- interval_metrics(res$coverage_sigma, res$ci_sigma_lower, res$ci_sigma_upper)

  rbind(
    data.frame(n = n, Parameter = "mu",    t(mu_point),    t(mu_int)),
    data.frame(n = n, Parameter = "sigma", t(sigma_point), t(sigma_int))
  )
}))
rownames(results_table) <- NULL

cat("\n========== Monte Carlo Results: Parametric Bootstrap Estimation ==========\n")
print(results_table, row.names = FALSE, digits = 5)
cat("============================================================================\n\n")

write.csv(results_table,
          file.path(out_dir, "bootstrap_normal_performance_metrics.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Visualization
# -----------------------------------------------------------------------------

## 4.1 Boxplots of the bootstrap point estimates across sample sizes ---------

png(file.path(out_dir, "boxplots_mu_sigma.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))

boxplot(mu_hat ~ n, data = sim_long,
        main = expression(paste("Bootstrap estimates of ", mu)),
        xlab = "Sample size (n)", ylab = expression(hat(mu)),
        col = "lightblue")
abline(h = true_mu, col = "red", lty = 2, lwd = 2)

boxplot(sigma_hat ~ n, data = sim_long,
        main = expression(paste("Bootstrap estimates of ", sigma)),
        xlab = "Sample size (n)", ylab = expression(hat(sigma)),
        col = "lightgreen")
abline(h = true_sigma, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))
dev.off()

## 4.2 Density plots of the bootstrap point estimates across sample sizes ----

colors <- c("steelblue", "darkorange", "forestgreen")

dens_mu_20     <- density(sim_results$n_20$mu_hat,     na.rm = TRUE)
dens_mu_100    <- density(sim_results$n_100$mu_hat,    na.rm = TRUE)
dens_mu_500    <- density(sim_results$n_500$mu_hat,    na.rm = TRUE)
dens_sigma_20  <- density(sim_results$n_20$sigma_hat,  na.rm = TRUE)
dens_sigma_100 <- density(sim_results$n_100$sigma_hat, na.rm = TRUE)
dens_sigma_500 <- density(sim_results$n_500$sigma_hat, na.rm = TRUE)

ymax_mu    <- max(dens_mu_20$y, dens_mu_100$y, dens_mu_500$y, na.rm = TRUE) * 1.15
ymax_sigma <- max(dens_sigma_20$y, dens_sigma_100$y, dens_sigma_500$y, na.rm = TRUE) * 1.15

png(file.path(out_dir, "densities_mu_sigma.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))

plot(dens_mu_20,
     main = expression(paste("Sampling distribution of ", hat(mu))),
     xlab = expression(hat(mu)), col = colors[1], lwd = 2,
     xlim = range(sim_long$mu_hat, na.rm = TRUE), ylim = c(0, ymax_mu))
lines(dens_mu_100, col = colors[2], lwd = 2)
lines(dens_mu_500, col = colors[3], lwd = 2)
abline(v = true_mu, col = "red", lty = 2)
legend("topright", legend = paste0("n = ", n_values),
       col = colors, lwd = 2, bty = "n")

plot(dens_sigma_20,
     main = expression(paste("Sampling distribution of ", hat(sigma))),
     xlab = expression(hat(sigma)), col = colors[1], lwd = 2,
     xlim = range(sim_long$sigma_hat, na.rm = TRUE), ylim = c(0, ymax_sigma))
lines(dens_sigma_100, col = colors[2], lwd = 2)
lines(dens_sigma_500, col = colors[3], lwd = 2)
abline(v = true_sigma, col = "red", lty = 2)
legend("topright", legend = paste0("n = ", n_values),
       col = colors, lwd = 2, bty = "n")

par(mfrow = c(1, 1))
dev.off()

cat("Plots saved to:", normalizePath(out_dir), "\n")
cat("Results table saved to:",
    normalizePath(file.path(out_dir, "bootstrap_normal_performance_metrics.csv")), "\n")
