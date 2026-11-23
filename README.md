# Parameter Estimation of Probability Distributions Using Bootstrap Resampling: A Simulation Study

**Department of Statistics, Faculty of Science — King Abdulaziz University (2026)**  
**Authors:** Muath Mohammed Ishaq, Ibrahim Ayad Alharbi  
**Supervisor:** Prof. Bander Al-Zahrani  

---

## Overview
Classical parameter estimation methods such as **Maximum Likelihood Estimation (MLE)** and the **Method of Moments (MoM)** achieve strong asymptotic efficiency in large samples, but can become biased, unstable, or lack closed-form solutions when sample sizes are small. This study evaluates the finite-sample performance of **Bias-Corrected Parametric Bootstrap** resampling ($\hat{\theta}_{BC} = 2\hat{\theta} - \bar{\theta}^*$) against MLE and MoM across symmetric, right-skewed, and complex lifetime distributions using a paired **Monte Carlo simulation** framework in **R**.

## Simulation Design
- **Distributions Evaluated:**
  - **Standard Normal:** $N(\mu = 0, \sigma = 1)$
  - **Exponential:** $\text{Exp}(\lambda = 0.5)$
  - **Weibull:** $\text{Weibull}(\alpha = 0.5, \beta = 4)$
- **Sample Sizes:** Restricted ($n = 20$), Moderate ($n = 100$), and Large ($n = 500$)
- **Replications & Reproducibility:** $M = 1,000$ paired Monte Carlo replications (Common Random Numbers) with $B = 200$ parametric bootstrap resamples per replication (`set.seed(2026)`)
- **Evaluation Metrics:** Bias, Relative Bias, Variance, Mean Squared Error (MSE), Relative MSE, and 95% Percentile Confidence Interval Coverage Probability
- **R Packages Used:** Base R (`stats`, `graphics`), `fitdistrplus`, `ggplot2`, `dplyr`, `tidyr`

## Key Findings (Restricted Sample Size, $n = 20$)
| Distribution | Parameter | Classical Bias (MoM / MLE) | Parametric Bootstrap Bias | Key Takeaway |
| :--- | :---: | :---: | :---: | :--- |
| **Normal** | $\sigma$ | $-0.0348$ | **$0.0020$** | Eliminated downward finite-sample scale bias |
| **Exponential** | $\lambda$ | $0.0223$ | **$-0.0048$** | Reduced bias and lowered MSE ($0.0152$ vs. $0.0174$) |
| **Weibull** | $\alpha$ (Shape) | $0.1383$ / $0.0372$ | **$-0.0041$** | Lowest bias & MSE ($0.0086$); MoM was severely unstable |
| **Weibull** | $\beta$ (Scale) | $1.3311$ / $0.3288$ | **$0.0441$** | Reduced scale MSE from $8.9797$ (MoM) and $4.6461$ (MLE) to **$4.2179$** |

- **Asymptotic Convergence:** At $n = 500$, MLE and Parametric Bootstrap converge with near-zero bias across all three distributions, while Weibull MoM remains less efficient.
- **Interval Coverage:** Percentile bootstrap intervals exhibit under-coverage at $n = 20$ (0.871–0.922) and approach the nominal 95% level as $n$ grows.
├── README.md
├── Report/
│   └── Graduation_Project_Report.pdf        # Full thesis report, derivations, tables & plots
└── R_Scripts/
    ├── 01_Normal/
    │   ├── 01_normal_mom.R                  # Method of Moments simulation for N(0, 1)
    │   ├── 02_normal_mle.R                  # Maximum Likelihood simulation for N(0, 1)
    │   ├── 03_normal_parametric_bootstrap.R # Parametric Bootstrap & 95% CI simulation
    │   └── 04_normal_mc_convergence.R       # Monte Carlo running-average convergence diagnostics
    ├── 02_Exponential/
    │   ├── 01_exponential_mom.R             # Method of Moments simulation for Exp(0.5)
    │   ├── 02_exponential_mle.R             # Maximum Likelihood simulation for Exp(0.5)
    │   └── 03_exponential_parametric_bootstrap.R # Parametric Bootstrap & 95% CI simulation
    └── 03_Weibull/
        ├── 01_weibull_mom.R                 # Numerical MoM (uniroot) simulation for Weibull(0.5, 4)
        ├── 02_weibull_mle.R                 # Numerical MLE (fitdistrplus / L-BFGS-B) simulation
        └── 03_weibull_parametric_bootstrap.R # Parametric Bootstrap & 95% CI simulation
```
