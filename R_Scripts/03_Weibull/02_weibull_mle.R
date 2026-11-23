## ============================================================
## Weibull Distribution Parameter Estimation - Monte Carlo Study
## Method: Maximum Likelihood Estimation (MLE)
## Target: Weibull(shape = 0.5, scale = 4)
## Paired design: the same sample x per replication is used by
## the MoM, MLE and Bootstrap scripts (common random numbers).
## ============================================================

## ---- 0. Setup -------------------------------------------------
set.seed(2026)

required_pkgs <- c("fitdistrplus", "ggplot2", "tidyr", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(fitdistrplus)
library(ggplot2)
library(tidyr)
library(dplyr)

out_dir <- "weibull_mle_output"
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

## ---- 2. Robust starting values (Method of Moments) --------------
weibull_mom_start <- function(x) {
  m1  <- mean(x)
  v   <- stats::var(x)
  cv2 <- v / m1^2

  cv2_eq <- function(k) gamma(1 + 2 / k) / (gamma(1 + 1 / k))^2 - 1 - cv2

  k0 <- tryCatch(
    stats::uniroot(cv2_eq, interval = c(0.05, 100), extendInt = "yes",
                    tol = .Machine$double.eps^0.25)$root,
    error = function(e) 1   # safe fallback if MoM root-finding fails
  )
  if (!is.finite(k0) || k0 <= 0) k0 <- 1

  lambda0 <- m1 / gamma(1 + 1 / k0)
  if (!is.finite(lambda0) || lambda0 <= 0) lambda0 <- m1  # safe default

  c(shape = k0, scale = lambda0)
}

## ---- 3. MLE estimator for Weibull -------------------------------
weibull_nll <- function(par, x) {
  shape <- par[1]
  scale <- par[2]
  if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) {
    return(1e10)
  }
  nll <- -sum(stats::dweibull(x, shape = shape, scale = scale, log = TRUE))
  if (!is.finite(nll)) return(1e10)
  nll
}

weibull_mle_optim <- function(x, start) {
  fit <- tryCatch(
    stats::optim(par = start, fn = weibull_nll, x = x,
                 method = "L-BFGS-B", lower = c(1e-4, 1e-4),
                 upper = c(1e4, 1e4)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0 || !all(is.finite(fit$par))) {
    return(c(shape = NA_real_, scale = NA_real_))
  }
  c(shape = unname(fit$par[1]), scale = unname(fit$par[2]))
}

weibull_mle <- function(x) {
  start <- as.list(weibull_mom_start(x))

  fit <- tryCatch(
    suppressWarnings(
      fitdistrplus::fitdist(x, distr = "weibull", method = "mle",
                             start = start)
    ),
    error   = function(e) NULL,
    warning = function(w) NULL
  )

  if (!is.null(fit) &&
      isTRUE(fit$convergence == 0) &&
      all(is.finite(fit$estimate))) {
    return(c(shape = unname(fit$estimate["shape"]),
             scale = unname(fit$estimate["scale"])))
  }

  # Fallback: base-R optimization, seeded with the same MoM start
  weibull_mle_optim(x, start = unlist(start))
}

## ---- 4. Monte Carlo simulation -----------------------------------
sim_results   <- list()
fail_counts   <- setNames(integer(length(n_vals)), as.character(n_vals))

for (n in n_vals) {

  shape_hat <- rep(NA_real_, n_rep)
  scale_hat <- rep(NA_real_, n_rep)
  samples   <- sample_bank[[as.character(n)]]

  for (r in seq_len(n_rep)) {
    x   <- samples[r, ]
    est <- tryCatch(
      weibull_mle(x),
      error = function(e) c(shape = NA_real_, scale = NA_real_)
    )

    if (all(is.finite(est))) {
      shape_hat[r] <- unname(est["shape"])
      scale_hat[r] <- unname(est["scale"])
    }
  }

  fail_counts[as.character(n)] <- sum(!is.finite(shape_hat) | !is.finite(scale_hat))

  sim_results[[as.character(n)]] <- data.frame(
    n         = n,
    rep       = seq_len(n_rep),
    shape_hat = shape_hat,
    scale_hat = scale_hat
  )
}

estimates_df <- dplyr::bind_rows(sim_results)

message("Non-converging replications (recorded as NA) per n:")
print(fail_counts)

write.csv(estimates_df,
          file = file.path(out_dir, "weibull_mle_raw_estimates.csv"),
          row.names = FALSE)

## ---- 5. Performance metrics ---------------------------------------
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
    n              = n,
    Converged_Reps = sum(is.finite(sub$shape_hat) & is.finite(sub$scale_hat)),
    Failed_Reps    = unname(fail_counts[as.character(n)]),
    Parameter      = c("shape (alpha)", "scale (beta)"),
    True_Value     = c(alpha_true, beta_true),
    Mean_Estimate  = c(mean(sub$shape_hat, na.rm = TRUE), mean(sub$scale_hat, na.rm = TRUE)),
    Bias           = c(shape_metrics["Bias"], scale_metrics["Bias"]),
    Relative_Bias  = c(shape_metrics["Relative_Bias"], scale_metrics["Relative_Bias"]),
    Variance       = c(shape_metrics["Variance"], scale_metrics["Variance"]),
    MSE            = c(shape_metrics["MSE"], scale_metrics["MSE"]),
    Relative_MSE   = c(shape_metrics["Relative_MSE"], scale_metrics["Relative_MSE"])
  )
})

metrics_df <- dplyr::bind_rows(metrics_list)

## ---- 6. Export metrics to CSV -------------------------------------
write.csv(metrics_df,
          file = file.path(out_dir, "weibull_mle_metrics.csv"),
          row.names = FALSE)

print(metrics_df)

## ---- 7. Visualizations ---------------------------------------------

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
    title = "Weibull MLE Estimates by Sample Size",
    subtitle = "Dashed red line = true parameter value",
    x = "Sample size (n)",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "weibull_mle_boxplots.png"),
       plot = boxplot_gg, width = 9, height = 5, dpi = 300)

density_gg <- ggplot(plot_df, aes(x = estimate, fill = n, color = n)) +
  geom_density(alpha = 0.35, na.rm = TRUE) +
  geom_blank(data = ylim_df, aes(y = y), inherit.aes = FALSE) +
  geom_vline(data = true_vals_df, aes(xintercept = true_value),
             linetype = "dashed", color = "red") +
  facet_wrap(~ parameter, scales = "free") +
  labs(
    title = "Sampling Distribution of Weibull MLE Estimates",
    subtitle = "Dashed red line = true parameter value",
    x = "Estimate",
    y = "Density",
    fill = "n", color = "n"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "weibull_mle_density.png"),
       plot = density_gg, width = 9, height = 5, dpi = 300)

## ---- 8. Session info (for reproducibility record) -------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "session_info.txt"))

message("Done. All outputs written to: ", normalizePath(out_dir))
