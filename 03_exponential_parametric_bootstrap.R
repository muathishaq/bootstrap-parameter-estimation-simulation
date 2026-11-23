# ============================================================================
# Parametric Bootstrap Estimation for the Exponential Distribution
#
# Procedure (per Monte Carlo replication, for a given n):
#   1. Take the paired original sample x_1,...,x_n ~ Exp(rate = true_lambda)
#   2. Compute the plug-in estimate lambda_hat = 1 / xbar (MLE/MOM).
#   3. Parametric bootstrap: draw B bootstrap samples of size n from
#      Exp(rate = lambda_hat) and compute lambda_hat*_b = 1 / xbar*_b
#   4. Bootstrap bias-corrected point estimate:
#        lambda_hat_BC = 2 * lambda_hat - mean(lambda_hat*)
#   5. Build a 95% CI for lambda via the percentile method:
#        CI = [ quantile(lambda_hat*, 0.025), quantile(lambda_hat*, 0.975) ]
# ============================================================================

## ---- 0. Setup -------------------------------------------------------------

set.seed(2026)                     # reproducibility

true_lambda <- 0.5                 # true rate parameter
n_values    <- c(20, 100, 500)     # sample sizes under study
mc_reps     <- 1000                # Monte Carlo replications
B           <- 200                 # parametric bootstrap resamples per MC replication
conf_level  <- 0.95                # confidence level for the percentile CI

out_dir <- "boots_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

alpha <- 1 - conf_level
probs <- c(alpha / 2, 1 - alpha / 2)

## ---- 1. Paired Monte Carlo samples (common random numbers) ----------------

set.seed(2026)
sample_bank <- lapply(n_values, function(n) {
  matrix(rexp(n * mc_reps, rate = true_lambda), nrow = mc_reps, ncol = n, byrow = TRUE)
})
names(sample_bank) <- paste0("n_", n_values)

## ---- 1b. Monte Carlo + parametric bootstrap simulation --------------------

simulate_bootstrap <- function(samples, B, probs) {
  n         <- ncol(samples)
  reps      <- nrow(samples)
  point_est <- numeric(reps)
  ci_lower  <- numeric(reps)
  ci_upper  <- numeric(reps)

  for (i in seq_len(reps)) {
    x          <- samples[i, ]
    lambda_hat <- 1 / mean(x)

    # Parametric bootstrap: resample from the fitted Exp(lambda_hat)
    boot_mat    <- matrix(rexp(n * B, rate = lambda_hat), nrow = n, ncol = B)
    boot_means  <- colMeans(boot_mat)
    boot_lambda <- 1 / boot_means

    # Bootstrap bias-corrected point estimate
    lambda_bc <- 2 * lambda_hat - mean(boot_lambda, na.rm = TRUE)
    if (!is.finite(lambda_bc) || lambda_bc <= 0) lambda_bc <- lambda_hat
    point_est[i] <- lambda_bc

    ci          <- quantile(boot_lambda, probs = probs, names = FALSE, na.rm = TRUE)
    ci_lower[i] <- ci[1]
    ci_upper[i] <- ci[2]
  }

  list(point_est = point_est, ci_lower = ci_lower, ci_upper = ci_upper)
}

sim_results <- lapply(sample_bank, simulate_bootstrap, B = B, probs = probs)
names(sim_results) <- paste0("n_", n_values)

sim_data <- do.call(rbind, lapply(seq_along(n_values), function(i) {
  data.frame(
    n          = n_values[i],
    lambda_hat = sim_results[[i]]$point_est
  )
}))
sim_data$n <- factor(sim_data$n, levels = n_values)

## ---- 2. Point-estimate metrics ---------------------------------------------

compute_point_metrics <- function(estimates, true_value) {
  bias     <- mean(estimates, na.rm = TRUE) - true_value
  rel_bias <- bias / true_value
  variance <- var(estimates, na.rm = TRUE)
  mse      <- mean((estimates - true_value)^2, na.rm = TRUE)
  rel_mse  <- mse / true_value^2

  data.frame(
    Bias          = bias,
    Relative_Bias = rel_bias,
    Variance      = variance,
    MSE           = mse,
    Relative_MSE  = rel_mse
  )
}

## ---- 3. Confidence-interval metrics ----------------------------------------

compute_ci_metrics <- function(ci_lower, ci_upper, true_value) {
  covered <- (ci_lower <= true_value) & (ci_upper >= true_value)
  data.frame(
    Coverage_Probability = mean(covered, na.rm = TRUE),
    Average_CI_Length    = mean(ci_upper - ci_lower, na.rm = TRUE)
  )
}

## ---- 4. Assemble full metrics table ----------------------------------------

metrics_df <- do.call(rbind, lapply(seq_along(n_values), function(i) {
  res      <- sim_results[[i]]
  pt_m     <- compute_point_metrics(res$point_est, true_lambda)
  ci_m     <- compute_ci_metrics(res$ci_lower, res$ci_upper, true_lambda)

  cbind(
    n             = n_values[i],
    True_Lambda   = true_lambda,
    Mean_Estimate = mean(res$point_est, na.rm = TRUE),
    pt_m,
    ci_m
  )
}))

print(metrics_df)

## ---- 5. Export metrics to CSV ----------------------------------------------

write.csv(
  metrics_df,
  file      = file.path(out_dir, "boots_performance_metrics.csv"),
  row.names = FALSE
)

## ---- 6. Boxplots: lambda_hat across sample sizes ---------------------------

png(file.path(out_dir, "boots_boxplot.png"), width = 900, height = 650, res = 120)

boxplot(
  lambda_hat ~ n, data = sim_data,
  main = expression(paste("Parametric Bootstrap Estimates of ", lambda, " by Sample Size")),
  xlab = "Sample Size (n)",
  ylab = expression(hat(lambda)),
  col  = "lightsalmon"
)
abline(h = true_lambda, col = "red", lwd = 2, lty = 2)
legend(
  "topright",
  legend = expression(paste("True ", lambda, " = 0.5")),
  col = "red", lty = 2, lwd = 2, bty = "n"
)

dev.off()

## ---- 7. Density plots: sampling distribution of lambda_hat -----------------

png(file.path(out_dir, "boots_density.png"), width = 900, height = 650, res = 120)

colors    <- c("blue", "darkgreen", "purple")
dens_list <- lapply(sim_results, function(res) density(res$point_est, na.rm = TRUE))

xlim_range <- range(unlist(lapply(dens_list, function(d) d$x)), na.rm = TRUE)
ymax       <- max(unlist(lapply(dens_list, function(d) d$y)), na.rm = TRUE) * 1.15

plot(
  dens_list[[1]], col = colors[1], lwd = 2,
  xlim = xlim_range, ylim = c(0, ymax),
  main = expression(paste("Sampling Distribution of ", hat(lambda), " (Parametric Bootstrap)")),
  xlab = expression(hat(lambda))
)
for (i in 2:length(dens_list)) {
  lines(dens_list[[i]], col = colors[i], lwd = 2)
}
abline(v = true_lambda, col = "red", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c(paste("n =", n_values), "True lambda"),
  col    = c(colors, "red"),
  lwd    = 2,
  lty    = c(1, 1, 1, 2)
)

dev.off()

## ---- 8. Wrap-up -------------------------------------------------------------

cat("\nSimulation complete.\n")
cat("Performance metrics:\n")
print(metrics_df)
cat("\nOutputs saved to:", normalizePath(out_dir), "\n")
