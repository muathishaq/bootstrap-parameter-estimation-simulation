## ============================================================
## Weibull Distribution Parameter Estimation - Monte Carlo Study
## Method: Parametric Bootstrap (bias-corrected point estimates,
##         percentile-method 95% CIs)
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

out_dir <- "weibull_boots_output"
if (!dir.exists(out_dir)) dir.create(out_dir)

## ---- 1. Simulation parameters ----------------------------------
alpha_true <- 0.5   # true shape
beta_true  <- 4     # true scale
n_vals     <- c(20, 100, 500)
n_rep      <- 1000  # Monte Carlo replications
B          <- 200   # parametric bootstrap resamples per replication
conf_level <- 0.95
alpha_ci   <- 1 - conf_level

max_attempt_factor <- 5  # safety cap on bootstrap resample retry attempts

## ---- 1b. Paired Monte Carlo samples (common random numbers) -----
set.seed(2026)
sample_bank <- lapply(n_vals, function(n) {
  matrix(stats::rweibull(n * n_rep, shape = alpha_true, scale = beta_true),
         nrow = n_rep, ncol = n, byrow = TRUE)
})
names(sample_bank) <- as.character(n_vals)

## ---- 2. Robust MoM starting values ------------------------------
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

## ---- 4. Parametric bootstrap for one Monte Carlo replication ------
parametric_bootstrap_ci <- function(n, shape_hat, scale_hat, B, max_attempts) {
  boot_shape <- numeric(0)
  boot_scale <- numeric(0)
  attempts   <- 0L

  while (length(boot_shape) < B && attempts < max_attempts) {
    attempts <- attempts + 1L

    est <- tryCatch(
      {
        xb <- stats::rweibull(n, shape = shape_hat, scale = scale_hat)
        weibull_mle(xb)
      },
      error = function(e) c(shape = NA_real_, scale = NA_real_)
    )

    if (all(is.finite(est))) {
      boot_shape <- c(boot_shape, unname(est["shape"]))
      boot_scale <- c(boot_scale, unname(est["scale"]))
    }
  }

  if (length(boot_shape) < B) {
    warning(sprintf(
      "Bootstrap: only %d / %d resamples converged after %d attempts.",
      length(boot_shape), B, attempts))
  }

  if (length(boot_shape) < 2) {
    return(c(shape_lo = NA_real_, shape_hi = NA_real_,
             scale_lo = NA_real_, scale_hi = NA_real_,
             shape_mean = NA_real_, scale_mean = NA_real_,
             n_boot = length(boot_shape)))
  }

  shape_ci <- stats::quantile(boot_shape, probs = c(alpha_ci / 2, 1 - alpha_ci / 2),
                               names = FALSE, type = 7, na.rm = TRUE)
  scale_ci <- stats::quantile(boot_scale, probs = c(alpha_ci / 2, 1 - alpha_ci / 2),
                               names = FALSE, type = 7, na.rm = TRUE)

  c(shape_lo = shape_ci[1], shape_hi = shape_ci[2],
    scale_lo = scale_ci[1], scale_hi = scale_ci[2],
    shape_mean = mean(boot_shape, na.rm = TRUE),
    scale_mean = mean(boot_scale, na.rm = TRUE),
    n_boot = length(boot_shape))
}

## ---- 5. Monte Carlo simulation ------------------------------------
sim_results <- list()
fail_counts <- setNames(integer(length(n_vals)), as.character(n_vals))
mle_fallback_counts <- setNames(integer(length(n_vals)), as.character(n_vals))

overall_start <- Sys.time()

for (n in n_vals) {

  n_start <- Sys.time()

  shape_hat_v <- rep(NA_real_, n_rep)
  scale_hat_v <- rep(NA_real_, n_rep)
  shape_lo_v  <- rep(NA_real_, n_rep)
  shape_hi_v  <- rep(NA_real_, n_rep)
  scale_lo_v  <- rep(NA_real_, n_rep)
  scale_hi_v  <- rep(NA_real_, n_rep)
  n_fail      <- 0L
  n_mle_fallback    <- 0L
  max_boot_attempts <- B * max_attempt_factor
  samples     <- sample_bank[[as.character(n)]]

  for (r in seq_len(n_rep)) {

    x <- samples[r, ]

    theta <- tryCatch(weibull_mle(x),
                      error = function(e) c(shape = NA_real_, scale = NA_real_))
    if (!all(is.finite(theta))) {
      theta <- tryCatch(weibull_mom_start(x),
                        error = function(e) c(shape = NA_real_, scale = NA_real_))
      n_mle_fallback <- n_mle_fallback + 1L
    }

    rep_result <- tryCatch(
      {
        if (!all(is.finite(theta))) {
          NULL
        } else {
          ci <- parametric_bootstrap_ci(n, theta["shape"], theta["scale"],
                                         B = B, max_attempts = max_boot_attempts)
          if (any(!is.finite(ci[c("shape_lo", "shape_hi", "scale_lo", "scale_hi",
                                  "shape_mean", "scale_mean")]))) {
            NULL
          } else {
            shape_bc <- 2 * unname(theta["shape"]) - unname(ci["shape_mean"])
            scale_bc <- 2 * unname(theta["scale"]) - unname(ci["scale_mean"])
            if (!is.finite(shape_bc) || shape_bc <= 0) shape_bc <- unname(theta["shape"])
            if (!is.finite(scale_bc) || scale_bc <= 0) scale_bc <- unname(theta["scale"])

            c(shape_hat = shape_bc, scale_hat = scale_bc,
              shape_lo = unname(ci["shape_lo"]), shape_hi = unname(ci["shape_hi"]),
              scale_lo = unname(ci["scale_lo"]), scale_hi = unname(ci["scale_hi"]))
          }
        }
      },
      error = function(e) NULL
    )

    if (is.null(rep_result)) {
      n_fail <- n_fail + 1L
      next
    }

    shape_hat_v[r] <- rep_result["shape_hat"]
    scale_hat_v[r] <- rep_result["scale_hat"]
    shape_lo_v[r]  <- rep_result["shape_lo"]
    shape_hi_v[r]  <- rep_result["shape_hi"]
    scale_lo_v[r]  <- rep_result["scale_lo"]
    scale_hi_v[r]  <- rep_result["scale_hi"]

    if (r %% 100 == 0) {
      message(sprintf("  n = %d: %d / %d replications done (%.1f sec elapsed)",
                       n, r, n_rep,
                       as.numeric(Sys.time() - n_start, units = "secs")))
    }
  }

  if (n_fail > 0) {
    warning(sprintf(
      "n = %d: %d / %d replications failed and are recorded as NA.",
      n, n_fail, n_rep))
  }

  fail_counts[as.character(n)] <- n_fail
  mle_fallback_counts[as.character(n)] <- n_mle_fallback

  sim_results[[as.character(n)]] <- data.frame(
    n         = n,
    rep       = seq_len(n_rep),
    shape_hat = shape_hat_v,
    scale_hat = scale_hat_v,
    shape_ci_lower = shape_lo_v,
    shape_ci_upper = shape_hi_v,
    scale_ci_lower = scale_lo_v,
    scale_ci_upper = scale_hi_v
  )

  message(sprintf("Finished n = %d in %.1f sec", n,
                   as.numeric(Sys.time() - n_start, units = "secs")))
}

message(sprintf("Total simulation time: %.1f sec",
                 as.numeric(Sys.time() - overall_start, units = "secs")))

estimates_df <- dplyr::bind_rows(sim_results)

estimates_df <- estimates_df %>%
  dplyr::mutate(
    shape_covered = (shape_ci_lower <= alpha_true) & (alpha_true <= shape_ci_upper),
    scale_covered = (scale_ci_lower <= beta_true)  & (beta_true  <= scale_ci_upper),
    shape_ci_length = shape_ci_upper - shape_ci_lower,
    scale_ci_length = scale_ci_upper - scale_ci_lower
  )

message("Failed replications (recorded as NA) per n:")
print(fail_counts)
message("Replications where MLE did not converge and the MoM fit was used as the bootstrap plug-in, per n:")
print(mle_fallback_counts)

write.csv(estimates_df,
          file = file.path(out_dir, "weibull_boots_raw_estimates.csv"),
          row.names = FALSE)

## ---- 6. Performance metrics ----------------------------------------
compute_point_metrics <- function(est, true_val) {
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

  shape_pm <- compute_point_metrics(sub$shape_hat, alpha_true)
  scale_pm <- compute_point_metrics(sub$scale_hat, beta_true)

  data.frame(
    n                   = n,
    Converged_Reps      = sum(is.finite(sub$shape_hat) & is.finite(sub$scale_hat)),
    Failed_Reps         = unname(fail_counts[as.character(n)]),
    Parameter           = c("shape (alpha)", "scale (beta)"),
    True_Value          = c(alpha_true, beta_true),
    Mean_Estimate       = c(mean(sub$shape_hat, na.rm = TRUE), mean(sub$scale_hat, na.rm = TRUE)),
    Bias                = c(shape_pm["Bias"], scale_pm["Bias"]),
    Relative_Bias       = c(shape_pm["Relative_Bias"], scale_pm["Relative_Bias"]),
    Variance            = c(shape_pm["Variance"], scale_pm["Variance"]),
    MSE                 = c(shape_pm["MSE"], scale_pm["MSE"]),
    Relative_MSE        = c(shape_pm["Relative_MSE"], scale_pm["Relative_MSE"]),
    Coverage_Probability = c(mean(sub$shape_covered, na.rm = TRUE),
                             mean(sub$scale_covered, na.rm = TRUE)),
    Average_CI_Length    = c(mean(sub$shape_ci_length, na.rm = TRUE),
                             mean(sub$scale_ci_length, na.rm = TRUE))
  )
})

metrics_df <- dplyr::bind_rows(metrics_list)

## ---- 7. Export metrics to CSV ---------------------------------------
write.csv(metrics_df,
          file = file.path(out_dir, "weibull_boots_metrics.csv"),
          row.names = FALSE)

print(metrics_df)

## ---- 8. Session info (for reproducibility record) --------------------
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "session_info.txt"))

message("Done. All outputs written to: ", normalizePath(out_dir))

## ============================================================
## ---- 9. Visualizations (Boxplots & Density Plots) -----------
## ============================================================
library(ggplot2)
library(tidyr)
library(dplyr)

out_dir    <- "weibull_boots_output"
alpha_true <- 0.5
beta_true  <- 4
n_vals     <- c(20, 100, 500)

estimates_df <- read.csv(file.path(out_dir, "weibull_boots_raw_estimates.csv"))

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
    title = "Weibull Parametric Bootstrap Estimates by Sample Size",
    subtitle = "Dashed red line = true parameter value",
    x = "Sample size (n)",
    y = "Estimate"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "weibull_boots_boxplots.png"),
       plot = boxplot_gg, width = 9, height = 5, dpi = 300)

density_gg <- ggplot(plot_df, aes(x = estimate, fill = n, color = n)) +
  geom_density(alpha = 0.35, na.rm = TRUE) +
  geom_blank(data = ylim_df, aes(y = y), inherit.aes = FALSE) +
  geom_vline(data = true_vals_df, aes(xintercept = true_value),
             linetype = "dashed", color = "red") +
  facet_wrap(~ parameter, scales = "free") +
  labs(
    title = "Sampling Distribution of Weibull Parametric Bootstrap Estimates",
    subtitle = "Dashed red line = true parameter value",
    x = "Estimate",
    y = "Density",
    fill = "n", color = "n"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "weibull_boots_density.png"),
       plot = density_gg, width = 9, height = 5, dpi = 300)

message("Plots saved to: ", normalizePath(out_dir))
