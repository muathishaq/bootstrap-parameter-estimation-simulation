# ============================================================================
# Method of Moments (MOM) Estimation for the Exponential Distribution
#
# Simulation study assessing the finite-sample performance of the MOM
# estimator for the rate parameter (lambda) of an Exponential distribution.
#
#   MOM estimator:  lambda_hat = 1 / X_bar
#
# Paired design: the Monte Carlo samples are generated up front from the same
# seed with identical code in the MOM, MLE and Bootstrap scripts, so every
# replication r of size n uses exactly the same sample x in all three methods.
# ============================================================================

## ---- 0. Setup -------------------------------------------------------------

set.seed(2026)                     # reproducibility

true_lambda <- 0.5                 # true rate parameter
n_values    <- c(20, 100, 500)     # sample sizes under study
n_reps      <- 1000                # Monte Carlo replications

out_dir <- "mom_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

## ---- 1. Paired Monte Carlo samples (common random numbers) ----------------

set.seed(2026)
sample_bank <- lapply(n_values, function(n) {
  matrix(rexp(n * n_reps, rate = true_lambda), nrow = n_reps, ncol = n, byrow = TRUE)
})
names(sample_bank) <- paste0("n_", n_values)

## ---- 1b. Monte Carlo simulation -------------------------------------------

simulate_mom <- function(samples) {
  estimates <- numeric(nrow(samples))
  for (i in seq_len(nrow(samples))) {
    x        <- samples[i, ]
    xbar     <- mean(x)
    estimates[i] <- 1 / xbar
  }
  estimates
}

estimates_list <- lapply(sample_bank, simulate_mom)
names(estimates_list) <- paste0("n_", n_values)

sim_data <- do.call(rbind, lapply(seq_along(n_values), function(i) {
  data.frame(
    n          = n_values[i],
    lambda_hat = estimates_list[[i]]
  )
}))
sim_data$n <- factor(sim_data$n, levels = n_values)

## ---- 2. Performance metrics ------------------------------------------------

compute_metrics <- function(estimates, true_value) {
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

metrics_df <- do.call(rbind, lapply(seq_along(n_values), function(i) {
  m <- compute_metrics(estimates_list[[i]], true_lambda)
  cbind(
    n             = n_values[i],
    True_Lambda   = true_lambda,
    Mean_Estimate = mean(estimates_list[[i]], na.rm = TRUE),
    m
  )
}))

print(metrics_df)

## ---- 3. Export metrics to CSV ----------------------------------------------

write.csv(
  metrics_df,
  file      = file.path(out_dir, "mom_performance_metrics.csv"),
  row.names = FALSE
)

## ---- 4. Boxplots: lambda_hat across sample sizes ---------------------------

png(file.path(out_dir, "mom_boxplot.png"), width = 900, height = 650, res = 120)

boxplot(
  lambda_hat ~ n, data = sim_data,
  main = expression(paste("MOM Estimates of ", lambda, " by Sample Size")),
  xlab = "Sample Size (n)",
  ylab = expression(hat(lambda)),
  col  = "lightblue"
)
abline(h = true_lambda, col = "red", lwd = 2, lty = 2)
legend(
  "topright",
  legend = expression(paste("True ", lambda, " = 0.5")),
  col = "red", lty = 2, lwd = 2, bty = "n"
)

dev.off()

## ---- 5. Density plots: sampling distribution of lambda_hat -----------------

png(file.path(out_dir, "mom_density.png"), width = 900, height = 650, res = 120)

colors    <- c("blue", "darkgreen", "purple")
dens_list <- lapply(estimates_list, density, na.rm = TRUE)

xlim_range <- range(unlist(lapply(dens_list, function(d) d$x)), na.rm = TRUE)
ymax       <- max(unlist(lapply(dens_list, function(d) d$y)), na.rm = TRUE) * 1.15

plot(
  dens_list[[1]], col = colors[1], lwd = 2,
  xlim = xlim_range, ylim = c(0, ymax),
  main = expression(paste("Sampling Distribution of ", hat(lambda), " (MOM)")),
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

## ---- 6. Wrap-up -------------------------------------------------------------

cat("\nSimulation complete.\n")
cat("Performance metrics:\n")
print(metrics_df)
cat("\nOutputs saved to:", normalizePath(out_dir), "\n")
