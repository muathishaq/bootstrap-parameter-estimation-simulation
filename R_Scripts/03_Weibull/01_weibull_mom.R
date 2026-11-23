## ============================================================
## Weibull Distribution Parameter Estimation - Monte Carlo Study
## Method: Method of Moments (MoM) via root-finding (uniroot)
## Target: Weibull(shape = 0.5, scale = 4)
## Paired design: the same sample x per replication is used by
## the MoM, MLE and Bootstrap scripts (common random numbers).
## ============================================================

## ---- 0. Setup -------------------------------------------------
set.seed(2026)

required_pkgs <- c("ggplot2", "tidyr", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(ggplot2)
library(tidyr)
library(dplyr)

out_dir <- "weibull_mom_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

## ---- 1. Simulation parameters ----------------------------------
alpha_true <- 0.5 # true shape
beta_true  <- 4   # true scale
n_vals     <- c(20, 100, 500)
n_rep      <- 1000

## ---- 1b. Paired Monte Carlo samples (common random numbers) -----
set.seed(2026)
sample_bank <- lapply(n_vals, function(n) {
  matrix(stats::rweibull(n * n_rep, shape = alpha_true, scale = beta_true),
         nrow = n_rep, ncol = n, byrow = TRUE)
})
names(sample_bank) <- as.character(n_vals)

## ---- 2. Method of Moments estimator for Weibull ----------------
## For X ~ Weibull(shape = k, scale = lambda):
##   E[X]   = lambda * Gamma(1 + 1/k)
##   Var[X] = lambda^2 * [Gamma(1 + 2/k) - Gamma(1 + 1/k)^2]
## => CV^2 = Var/E[X]^2 = Gamma(1+2/k)/Gamma(1+1/k)^2 - 1
## Solve the above for k numerically (uniroot), then recover lambda.
weibull_mom <- function(x) {
  m1  <- mean(x)
  v   <- stats::var(x)
  cv2 <- v / m1^2

  cv2_eq <- function(k) {
    gamma(1 + 2 / k) / (gamma(1 + 1 / k))^2 - 1 - cv2
  }

  k_hat <- tryCatch(
    {
      root <- stats::uniroot(cv2_eq, interval = c(0.05, 100),
                              extendInt = "yes", tol = .Machine$double.eps^0.25)
      root$root
    },
    error = function(e) NA_real_
  )

  lambda_hat <- if (!is.na(k_hat)) {
    m1 / gamma(1 + 1 / k_hat)
  } else {
    NA_real_
  }

  c(shape = k_hat, scale = lambda_hat)
}

## ---- 3. Monte Carlo simulation ----------------------------------
sim_results <- list()

for (n in n_vals) {

  shape_hat <- numeric(n_rep)
  scale_hat <- numeric(n_rep)
  samples   <- sample_bank[[as.character(n)]]

  for (r in seq_len(n_rep)) {
    x   <- samples[r, ]
    est <- tryCatch(weibull_mom(x),
                    error = function(e) c(shape = NA_real_, scale = NA_real_))
    shape_hat[r] <- est["shape"]
    scale_hat[r] <- est["scale"]
  }

  sim_results[[as.character(n)]] <- data.frame(
    n          = n,
    rep        = seq_len(n_rep),
    shape_hat  = shape_hat,
    scale_hat  = scale_hat
  )
}

estimates_df <- dplyr::bind_rows(sim_results)

# Drop failed (NA) replications, report how many were dropped
n_na <- sum(!is.finite(estimates_df$shape_hat) | !is.finite(estimates_df$scale_hat))
if (n_na > 0) {
  message(sprintf("Dropping %d replication(s) with failed MoM estimation.", n_na))
}
estimates_df <- estimates_df[is.finite(estimates_df$shape_hat) &
                               is.finite(estimates_df$scale_hat), ]

write.csv(estimates_df,
          file = file.path(out_dir, "weibull_mom_raw_estimates.csv"),
          row.names = FALSE)

## ---- 4. Performance metrics -------------------------------------
compute_metrics <- function(est, true_val) {
  bias     <- mean(est, na.rm = TRUE) - true_val
  rel_bias <- bias / true_val
  variance <- stats::var(est, na.rm = TRUE)
  mse      <- mean((est - true_val)^2, na.rm = TRUE)
  rel_mse  <- mse / true_val^2
  c(Bias = bias, Relative_Bias = rel_bias, Variance = variance,
    MSE = mse, Relative_MSE = rel_mse)
}

metrics_list <- lapply(n_vals, function(n) {
  sub <- estimates_df[estimates_df$n == n, ]

  shape_metrics <- compute_metrics(sub$shape_hat, alpha_true)
  scale_metrics <- compute_metrics(sub$scale_hat, beta_true)

  data.frame(
    n         = n,
    Parameter = c("shape (alpha)", "scale (beta)"),
    True_Value = c(alpha_true, beta_true),
    Mean_Estimate = c(mean(sub$shape_hat, na.rm = TRUE), mean(sub$scale_hat, na.rm = TRUE)),
    Bias          = c(shape_metrics["Bias"], scale_metrics["Bias"]),
    Relative_Bias = c(shape_metrics["Relative_Bias"], scale_metrics["Relative_Bias"]),
    Variance      = c(shape_metrics["Variance"], scale_metrics["Variance"]),
    MSE           = c(shape_metrics["MSE"], scale_metrics["MSE"]),
    Relative_MSE  = c(shape_metrics["Relative_MSE"], scale_metrics["Relative_MSE"])
  )
})

metrics_df <- dplyr::bind_rows(metrics_list)

## ---- 5. Export metrics to CSV -----------------------------------
write.csv(metrics_df,
          file = file.path(out_dir, "weibull_mom_metrics.csv"),
          row.names = FALSE)

print(metrics_df)

## ---- 6. Visualizations -------------------------------------------

plot_df <- estimates_df %>%
  tidyr::pivot_longer(
    cols = c(shape_hat, scale_hat),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  dplyr::mutate(
    parameter = factor(parameter,
                        levels = c("shape_hat", "scale_hat"),
                        labels = c("Shape (alpha)", "Scale (beta)")),
    n = factor(n, levels = n_vals)
  )

true_vals_df <- data.frame(
  parameter = factor(c("Shape (alpha)", "Scale (beta)"),
                      levels = c("Shape (alpha)", "Scale (beta)")),
  true_value = c(alpha_true, beta_true)
)

ymax_df <- plot_df %>%
  dplyr::filter(is.finite(estimate)) %>%
  dplyr::group_by(parameter, n) %>%
  dplyr::summarise(peak = if (dplyr::n() >= 2) max(stats::density(estimate, na.rm = TRUE)$y) else NA_real_,
                   .groups = "drop") %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(ymax = max(peak, na.rm = TRUE) * 1.15, .groups = "drop")
ylim_df <- rbind(data.frame(parameter = ymax_df$parameter, y = 0),
                 data.frame(parameter = ymax_df$parameter, y = ymax_df$ymax))

boxplot_gg <- ggplot(plot_df, aes(x = n, y = estimate, fill = n)) +
  geom_boxplot(outlier.alpha = 0.3, na.rm = TRUE) +
  geom_hline(data = true_vals_df, aes(yintercept = true_value),
             linetype = "dashed", color = "red") +
  facet_wrap(~ parameter, scales = "free_y") +
  labs(
    title = "Weibull MoM Estimates by Sample Size",
    subtitle = "Dashed red line = true parameter value",
    x = "Sample size (n)",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "weibull_mom_boxplots.png"),
       plot = boxplot_gg, width = 9, height = 5, dpi = 300)

density_gg <- ggplot(plot_df, aes(x = estimate, fill = n, color = n)) +
  geom_density(alpha = 0.35, na.rm = TRUE) +
  geom_blank(data = ylim_df, aes(y = y), inherit.aes = FALSE) +
  geom_vline(data = true_vals_df, aes(xintercept = true_value),
             linetype = "dashed", color = "red") +
  facet_wrap(~ parameter, scales = "free") +
  labs(
    title = "Sampling Distribution of Weibull MoM Estimates",
    subtitle = "Dashed red line = true parameter value",
    x = "Estimate",
    y = "Density",
    fill = "n", color = "n"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "weibull_mom_density.png"),
       plot = density_gg, width = 9, height = 5, dpi = 300)

## ---- 7. Session info (for reproducibility record) -----------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "session_info.txt"))

message("Done. All outputs written to: ", normalizePath(out_dir))
