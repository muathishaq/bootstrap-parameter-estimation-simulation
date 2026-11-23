# =============================================================================
# Monte Carlo Study: Maximum Likelihood Estimation (MLE)
# Standard Normal Distribution N(mu = 0, sigma = 1)
#
# For each Monte Carlo replication, a sample of size n is drawn from
# N(0, 1). The MLEs are obtained by maximizing the Normal log-likelihood
#
#   l(mu, sigma) = -n/2 * log(2*pi) - n * log(sigma)
#                  - 1/(2*sigma^2) * sum_i (x_i - mu)^2
#
# Setting the partial derivatives w.r.t. mu and sigma^2 to zero gives the
# well-known closed-form MLEs:
#
#   mu_hat      = xbar             = mean(x)
#   sigma_hat^2 = (1/n) sum (x_i - xbar)^2   (population form, divide by n)
#
# Note: for the Normal distribution the MLEs of mu and sigma coincide
# exactly with the Method of Moments (MOM) estimators used in mom.R -- this
# is a known theoretical result, not a coding error, and is worth remarking
# on in the Results chapter when comparing the two methods.
#
# Paired design: the Monte Carlo samples are generated up front from the same
# seed with identical code in the MOM, MLE and Bootstrap scripts, so every
# replication r of size n uses exactly the same sample x in all three methods.
#
# Performance is assessed via Bias, Relative Bias, Variance, MSE and
# Relative MSE, computed across 1000 replications for each of the sample
# sizes n = 20, 100, 500.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

rm(list = ls())
set.seed(2026)                     # full reproducibility, as specified

# Output directory for plots / tables (created if it does not exist yet)
out_dir <- "mle_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# True parameter values (Standard Normal)
true_mu    <- 0
true_sigma <- 1

# Simulation design
n_values <- c(20, 100, 500)        # sample sizes
R        <- 1000                   # Monte Carlo replications per sample size

# -----------------------------------------------------------------------------
# 0b. Paired Monte Carlo samples (common random numbers)
# -----------------------------------------------------------------------------
# Exactly one sample x per replication: row r of sample_bank[["n_<n>"]] is the
# sample of replication r. Identical across the MOM / MLE / Bootstrap scripts
# (same seed, same generation code, generated before anything else).

set.seed(2026)
sample_bank <- lapply(n_values, function(n) {
  matrix(rnorm(n * R, mean = true_mu, sd = true_sigma), nrow = R, ncol = n, byrow = TRUE)
})
names(sample_bank) <- paste0("n_", n_values)

# -----------------------------------------------------------------------------
# 1. Maximum Likelihood estimator
# -----------------------------------------------------------------------------

#' Compute the MLE estimates of mu and sigma for a single sample
#'
#' @param x numeric vector, a random sample
#' @return named numeric vector c(mu_hat, sigma_hat)
mle_normal <- function(x) {
  n <- length(x)

  mu_hat     <- mean(x)                        # argmax of the log-likelihood
  sigma2_hat <- mean((x - mu_hat)^2)            # MLE variance (divide by n)
  sigma_hat  <- sqrt(sigma2_hat)

  c(mu_hat = mu_hat, sigma_hat = sigma_hat)
}

# -----------------------------------------------------------------------------
# 2. Monte Carlo simulation
# -----------------------------------------------------------------------------

#' Run the Monte Carlo simulation for one sample size on the paired samples
#'
#' @param samples matrix of paired samples (R rows = replications, n columns)
#' @return data.frame with R rows and columns mu_hat, sigma_hat
run_simulation <- function(samples) {
  estimates <- t(apply(samples, 1, mle_normal))
  as.data.frame(estimates)
}

# Run the simulation for every sample size and store results in a list
sim_results <- setNames(
  lapply(sample_bank, run_simulation),
  paste0("n_", n_values)
)

# Combine into one long data frame: columns n, mu_hat, sigma_hat
sim_long <- do.call(rbind, lapply(n_values, function(n) {
  df   <- sim_results[[paste0("n_", n)]]
  df$n <- n
  df
}))
sim_long$n <- factor(sim_long$n, levels = n_values)

# -----------------------------------------------------------------------------
# 3. Performance metrics: Bias, Relative Bias, Variance, MSE, Relative MSE
# -----------------------------------------------------------------------------

#' Compute performance metrics for a vector of estimates
#'
#' @param estimates numeric vector of parameter estimates across replications
#' @param truth     true parameter value
#' @return named numeric vector of metrics
performance_metrics <- function(estimates, truth) {
  bias     <- mean(estimates, na.rm = TRUE) - truth
  rel_bias <- if (truth != 0) bias / truth else NA   # undefined when truth = 0
  variance <- var(estimates, na.rm = TRUE)
  mse      <- mean((estimates - truth)^2, na.rm = TRUE)
  rel_mse  <- if (truth != 0) mse / truth^2 else NA  # undefined when truth = 0

  c(Bias = bias, RelBias = rel_bias, Variance = variance,
    MSE = mse, RelMSE = rel_mse)
}

# Build a results table: one row per (sample size, parameter) combination
results_table <- do.call(rbind, lapply(n_values, function(n) {
  df <- sim_results[[paste0("n_", n)]]

  mu_metrics    <- performance_metrics(df$mu_hat,    true_mu)
  sigma_metrics <- performance_metrics(df$sigma_hat, true_sigma)

  rbind(
    data.frame(n = n, Parameter = "mu",    t(mu_metrics)),
    data.frame(n = n, Parameter = "sigma", t(sigma_metrics))
  )
}))
rownames(results_table) <- NULL

# Note on Relative Bias / Relative MSE for mu:
# Since the true value of mu is 0, dividing by it is mathematically
# undefined, so these columns are reported as NA for mu (this is expected
# and should simply be noted in the Results chapter, not treated as an error).

cat("\n================ Monte Carlo Results: MLE Estimation ================\n")
print(results_table, row.names = FALSE, digits = 5)
cat("=======================================================================\n\n")

# Save the results table for direct inclusion in a thesis/report
write.csv(results_table,
          file.path(out_dir, "mle_normal_performance_metrics.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Visualization
# -----------------------------------------------------------------------------

## 4.1 Boxplots of the estimates across sample sizes -------------------------

png(file.path(out_dir, "boxplots_mu_sigma.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))

boxplot(mu_hat ~ n, data = sim_long,
        main = expression(paste("MLE estimates of ", mu)),
        xlab = "Sample size (n)", ylab = expression(hat(mu)),
        col = "lightblue")
abline(h = true_mu, col = "red", lty = 2, lwd = 2)

boxplot(sigma_hat ~ n, data = sim_long,
        main = expression(paste("MLE estimates of ", sigma)),
        xlab = "Sample size (n)", ylab = expression(hat(sigma)),
        col = "lightgreen")
abline(h = true_sigma, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1))
dev.off()

## 4.2 Density plots of the estimates across sample sizes --------------------

colors <- c("steelblue", "darkorange", "forestgreen")

# Pre-compute densities once so they can be reused for both the PNG export
# and the interactive redraw below, without calling density() twice.
dens_mu_20     <- density(sim_results$n_20$mu_hat,     na.rm = TRUE)
dens_mu_100    <- density(sim_results$n_100$mu_hat,    na.rm = TRUE)
dens_mu_500    <- density(sim_results$n_500$mu_hat,    na.rm = TRUE)
dens_sigma_20  <- density(sim_results$n_20$sigma_hat,  na.rm = TRUE)
dens_sigma_100 <- density(sim_results$n_100$sigma_hat, na.rm = TRUE)
dens_sigma_500 <- density(sim_results$n_500$sigma_hat, na.rm = TRUE)

# Dynamic y-axis limits: use the tallest peak across all three sample sizes
# (rather than assuming which one it will be), plus 15% headroom on top.
ymax_mu    <- max(dens_mu_20$y, dens_mu_100$y, dens_mu_500$y, na.rm = TRUE) * 1.15
ymax_sigma <- max(dens_sigma_20$y, dens_sigma_100$y, dens_sigma_500$y, na.rm = TRUE) * 1.15

png(file.path(out_dir, "densities_mu_sigma.png"), width = 1000, height = 500)
par(mfrow = c(1, 2))

# Density of mu_hat for each sample size
plot(dens_mu_20,
     main = expression(paste("Sampling distribution of ", hat(mu))),
     xlab = expression(hat(mu)), col = colors[1], lwd = 2,
     xlim = range(sim_long$mu_hat, na.rm = TRUE), ylim = c(0, ymax_mu))
lines(dens_mu_100, col = colors[2], lwd = 2)
lines(dens_mu_500, col = colors[3], lwd = 2)
abline(v = true_mu, col = "red", lty = 2)
legend("topright", legend = paste0("n = ", n_values),
       col = colors, lwd = 2, bty = "n")

# Density of sigma_hat for each sample size
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

# Re-draw both figures to the active graphics device as well, so they are
# visible immediately in an interactive R / RStudio session.
par(mfrow = c(1, 2))
boxplot(mu_hat ~ n, data = sim_long,
        main = expression(paste("MLE estimates of ", mu)),
        xlab = "Sample size (n)", ylab = expression(hat(mu)), col = "lightblue")
abline(h = true_mu, col = "red", lty = 2, lwd = 2)
boxplot(sigma_hat ~ n, data = sim_long,
        main = expression(paste("MLE estimates of ", sigma)),
        xlab = "Sample size (n)", ylab = expression(hat(sigma)), col = "lightgreen")
abline(h = true_sigma, col = "red", lty = 2, lwd = 2)
par(mfrow = c(1, 1))

# Re-draw the density plots to the active graphics device as well, so they
# are visible immediately in an interactive R / RStudio session.
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

cat("Plots saved to:", normalizePath(out_dir), "\n")
cat("Results table saved to:",
    normalizePath(file.path(out_dir, "mle_normal_performance_metrics.csv")), "\n")
