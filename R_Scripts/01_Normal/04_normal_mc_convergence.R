# =============================================================================
# Monte Carlo Convergence Study
# Standard Normal Distribution N(mu = 0, sigma = 1)
#
# Rather than looking at the sampling distribution of a single estimate (as
# in mom.R / MLE.R), this script tracks how the *running (cumulative)
# average* of the estimates behaves as more Monte Carlo replications are
# accumulated:
#
#   running_avg[k] = (1/k) * sum_{i=1}^{k} estimate_i,   k = 1, ..., R
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

rm(list = ls())
set.seed(2026)                     # full reproducibility, as specified

out_dir <- "convergence_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

true_mu    <- 0
true_sigma <- 1

n_values <- c(20, 100, 500)        # sample sizes
R        <- 1000                   # Monte Carlo replications per sample size

# -----------------------------------------------------------------------------
# 1. Estimator (MLE = MOM for the Normal distribution)
# -----------------------------------------------------------------------------

estimate_params <- function(x) {
  mu_hat     <- mean(x)
  sigma_hat  <- sqrt(mean((x - mu_hat)^2))   # population-form (MLE/MOM) sd
  c(mu_hat = mu_hat, sigma_hat = sigma_hat)
}

# -----------------------------------------------------------------------------
# 2. Monte Carlo simulation + running (cumulative) average
# -----------------------------------------------------------------------------

run_convergence <- function(n, R) {
  estimates <- t(replicate(R, {
    x <- rnorm(n, mean = true_mu, sd = true_sigma)
    estimate_params(x)
  }))

  # Cumulative (running) average across iterations 1, 2, ..., R
  running_mu    <- cumsum(estimates[, "mu_hat"])    / seq_len(R)
  running_sigma <- cumsum(estimates[, "sigma_hat"]) / seq_len(R)

  data.frame(iteration = seq_len(R),
             running_mu = running_mu,
             running_sigma = running_sigma)
}

conv_results <- setNames(
  lapply(n_values, run_convergence, R = R),
  paste0("n_", n_values)
)

# -----------------------------------------------------------------------------
# 3. Visualization: running-average convergence plots
# -----------------------------------------------------------------------------

colors <- c("steelblue", "darkorange", "forestgreen")

padded_ylim <- function(values, truth, pad = 0.15) {
  rng   <- range(c(values, truth), na.rm = TRUE)
  width <- diff(rng)
  if (width == 0) width <- 1
  c(rng[1] - pad * width, rng[2] + pad * width)
}
ylim_mu    <- padded_ylim(unlist(lapply(conv_results, function(df) df$running_mu)),    true_mu)
ylim_sigma <- padded_ylim(unlist(lapply(conv_results, function(df) df$running_sigma)), true_sigma)

## 3.1 Convergence plot for mu ------------------------------------------------

png(file.path(out_dir, "convergence_mu.png"), width = 900, height = 600)

plot(conv_results$n_20$iteration, conv_results$n_20$running_mu,
     type = "l", col = colors[1], lwd = 2,
     xlim = c(0, R),
     ylim = ylim_mu,
     main = expression(paste("Monte Carlo Convergence of ", hat(mu),
                              " (Running Average)")),
     xlab = "Monte Carlo Iteration", ylab = "Running Average of Estimate")
lines(conv_results$n_100$iteration, conv_results$n_100$running_mu,
      col = colors[2], lwd = 2)
lines(conv_results$n_500$iteration, conv_results$n_500$running_mu,
      col = colors[3], lwd = 2)

abline(h = true_mu, col = "black", lty = 2, lwd = 3)

legend("bottomright",
       legend = c(paste0("n = ", n_values), "True mu"),
       col = c(colors, "black"),
       lty = c(1, 1, 1, 2),
       lwd = c(2, 2, 2, 3),
       bty = "n")

dev.off()

## 3.2 Convergence plot for sigma ---------------------------------------------

png(file.path(out_dir, "convergence_sigma.png"), width = 900, height = 600)

plot(conv_results$n_20$iteration, conv_results$n_20$running_sigma,
     type = "l", col = colors[1], lwd = 2,
     xlim = c(0, R),
     ylim = ylim_sigma,
     main = expression(paste("Monte Carlo Convergence of ", hat(sigma),
                              " (Running Average)")),
     xlab = "Monte Carlo Iteration", ylab = "Running Average of Estimate")
lines(conv_results$n_100$iteration, conv_results$n_100$running_sigma,
      col = colors[2], lwd = 2)
lines(conv_results$n_500$iteration, conv_results$n_500$running_sigma,
      col = colors[3], lwd = 2)

abline(h = true_sigma, col = "black", lty = 2, lwd = 3)

legend("bottomright",
       legend = c(paste0("n = ", n_values), "True sigma"),
       col = c(colors, "black"),
       lty = c(1, 1, 1, 2),
       lwd = c(2, 2, 2, 3),
       bty = "n")

dev.off()

cat("Convergence plots saved to:", normalizePath(out_dir), "\n")
