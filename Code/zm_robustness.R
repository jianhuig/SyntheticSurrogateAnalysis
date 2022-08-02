# Purpose: Demonstrate robustness to the choice of surrogate and power 
# advantage with increasing missingness.
# Updated: 2022-07-30
setwd("~/Documents/GitHub/SyntheticSurrogateAnalysis/Code")

#' Simulation
#' 
#' Evaluate performance of the synthetic surrogate estimator, compared
#' with oracle and marginal, as a function of the target-surrogate
#' correlation.
#' 
#' @param n0 Subjects with observed target outcomes.
#' @param reps Simulation replicates.
#' @param rho Target-surrogate correlation.
#' @return Data.frame.
Sim <- function(
    n0 = 1e3,
    reps = 500,
    rho = 0.00
) {
  
  
  # Proportion of variance explained by genotype.
  PVE_G <- 0.01
  
  # Missingness levels. miss = 0 included by default.
  MISS <- c(0.25, 0.50, 0.75)
  
  # The outer loop runs over simulation replicates.
  # The inner loop runs over target-surrogate correlations.
  # Because the performance of the oracle and marginal estimators
  #   does not depend on correlation, these are estimated using
  #   the rho = 0 setting. 
  
  outer_loop <- lapply(seq_len(reps), function(i) {
    
    # Generate data with miss = 0.
    data <- DGP(n0 = n0, miss = 0, rho = rho, pve_g = PVE_G)
    
    # Oracle estimator.
    oracle_est <- OracleEst(data)
    
    # Observed data estimator.
    obs_est <- ObsEst(data)
    
    # Synthetic surrogate estimator.
    ss_est <- SynSurrEst(data)
    
    # Collect oracle and observed.
    ests <- rbind(oracle_est, obs_est, ss_est)
    
    out <- data.frame(
      miss = 0,
      rho = rho,
      method = c("oracle", "obs", "ss"),
      est = ests[, 1],
      se = ests[, 2],
      row.names = NULL
    )
    
    # Inner loop over non-zero correlation levels.
    ss_est <- lapply(MISS, function(miss) {
      
      data <- DGP(n0 = n0, miss = miss, rho = rho, pve_g = PVE_G)
      ss_est <- SynSurrEst(data)
      out <- c(
        miss = miss,
        rho = rho,
        est = as.numeric(ss_est[1]),
        se = as.numeric(ss_est[2])
      )
      
      return(out)
    })
    ss_est <- do.call(rbind, ss_est)
    ss_est <- data.frame(ss_est) %>%
      tibble::add_column(method = "ss", .after = "rho")
    
    out <- rbind(out, ss_est)
    out$bias <- out$est - sqrt(PVE_G)
    return(out)
  })
  
  results <- do.call(rbind, outer_loop)
  return(results)
}


# -----------------------------------------------------------------------------


#' Data Generating Process
#' 
#' Generate data. "g" is genotype, "x" is a covariate, "y" is the outcome, and
#' "yobs" is the outcome after introduction of missingness, "s" is the
#' surrogate.
#'
#' The number of subjects with observed target outcomes is always n. If
#' missingness is > 0, then subjects with target outcomes missing are added to
#' achieve the desired level of missingness.
#' 
#' @param n0 Subjects with observed target outcomes.
#' @param miss Missingness rate.
#' @param rho Target-surrogate correlation.
#' @param maf Minor allele frequency.
#' @param pve_g Proportion of variation explained by genotype.
#' @param pve_x Proportion of variation explained by covariate.
#' @return Data.frame.
DGP <- function(
    n0,
    miss,
    rho,
    maf = 0.25,
    pve_g = 0.005,
    pve_x = 0.10
) {
  
  # Subjects with observed target outcomes.
  n1 <- round(n0 * miss / (1 - miss))
  n <- n0 + n1
  
  # Genotype and covariate.
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  g <- scale(g)
  x <- stats::rnorm(n = n)
  x <- scale(x)
  
  # Linear predictor.
  beta <- c(sqrt(pve_g), sqrt(pve_x))
  eta <- as.numeric(cbind(g, x) %*% beta)
  
  # Residuals.
  sigma <- 1 - sum(c(pve_g, pve_x))
  eps_s <- stats::rnorm(n, sd = sqrt(sigma))
  eps_y <- rho * eps_s + sqrt(1 - rho^2) * stats::rnorm(n, sd = sqrt(sigma))
  
  # Oracle outcome.
  s <- eta + eps_s
  y <- eta + eps_y
  
  # Outcome with missingness.
  draw <- sample(seq_len(n), size = round(miss * n), replace = FALSE)
  yobs <- y
  yobs[draw] <- NA
  
  out <- data.frame(
    g = g,
    x = x,
    s = s,
    y = y,
    yobs = yobs
  )
  return(out)
}


#' Estimate Genetic Effect
#' 
#' Estimates the genetic effect and standard error.
#' 
#' @param covar Data.frame of covariates. Should not include "g" or "y".
#' @param geno Genotype vector.
#' @param outcome Outcome vector.
#' @return Numeric vector containing the estimate "est" and standard error "se".
EstBetaG <- function(covar, geno, outcome) {
  df <- data.frame(
    y = outcome,
    g = geno,
    x = covar
  ) %>% na.omit()
  
  fit <- lm(y ~ g + x, data = df)
  fit_summary <- summary(fit)$coefficients
  
  out <- fit_summary[
    rownames(fit_summary) == "g", c("Estimate", "Std. Error")]
  names(out) <- c("est", "se")
  return(out)
}


#' Oracle Estimator
#' 
#' @param data Data.frame.
#' @return Numeric vector containing "oracle_est" and "oracle_se".
OracleEst <- function(data) {
  oracle_est <- EstBetaG(
    covar = data$x,
    geno = data$g,
    outcome = data$y
  )
  names(oracle_est) <- paste0("oracle_", names(oracle_est))
  return(oracle_est)
}


#' Observed Data Estimator
#' 
#' @param data Data.frame.
#' @return Numeric vector containing "obs_est" and "obs_se".
ObsEst <- function(data) {
  obs_est <- EstBetaG(
    covar = data$x,
    geno = data$g,
    outcome = data$yobs
  )
  names(obs_est) <- paste0("obs_", names(obs_est))
  return(obs_est)
}


#' SynSurr Estimator
#' 
#' @param data Data.frame.
#' @param imp_param Imputation parameters.
#' @return Numeric vector containing "ss_est" and "ss_se".
SynSurrEst <- function(data) {
  
  # SynSurr estimator.
  fit <- SurrogateRegression::Fit.BNLS(
    t = data$yobs,
    s = data$s, 
    X = data[, c("g", "x")],
  )
  param <- coef(fit, type = "Target")
  
  out <- c(
    ss_est = param$Point[1],
    ss_se = param$SE[1]
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Tabulate Simulation Results
#' 
#' @param sim Data.frame of simulation results.
#' @return Data.frame.
TabSim <- function(sim) {
  tab <- sim %>%
    dplyr::group_by(method, rho, miss) %>%
    dplyr::summarise(
      bias = mean(bias),
      se = sqrt(mean(se^2)),
      .groups = "drop"
    ) 
  return(tab)
}


#' Plot Simulation Results
#' 
#' @param sim Data.frame of simulation results.
#' @return ggplot.
PlotSim <- function(sim) {
  
  df <- sim %>%
    dplyr::filter(method != "oracle")
  df$miss[df$method == "obs"] <- -0.01
  df$miss <- factor(
    df$miss,
    levels = c(-0.01, 0, 0.25, 0.5, 0.75),
    ordered = TRUE
  )
  df$method <- paste0(df$method, "_", df$miss)
  df$method <- factor(
    x = df$method,
    levels = c("obs_-0.01", "ss_0", "ss_0.25", "ss_0.5", "ss_0.75"),
  )
  
  q <- ggplot(data = df) + 
    theme_bw() + 
    theme(legend.position = "bottom") + 
    geom_hline(
      yintercept = 0.0,
      linetype = "dashed",
      color = "gray",
      size = 1.2
    ) + 
    geom_boxplot(
      aes(x = method, y = bias, fill = miss),
    ) +
    scale_fill_brewer(
      name = "Missingness",
      palette = "Blues",
      labels = c("NA", "0%", "25%", "50%", "75%")
    ) + 
    scale_x_discrete(
      name = "Method",
      labels = c("Marginal", rep("SynSurr", times = 4))
    ) + 
    scale_y_continuous(
      name = "Bias"
    )
  return(q)
}


# -----------------------------------------------------------------------------


# Case 1: Rho = 0.0.
set.seed(10101)
sim <- Sim(rho = 0.0, reps = 1000)
tab <- TabSim(sim)
show(tab)
q0 <- PlotSim(sim) + 
  ggtitle(expression(Synthetic~Surrogate~Correlation~rho==0.00))


# Case 2: Rho = 0.75
set.seed(10101)
sim <- Sim(rho = 0.75, reps = 1000)
tab <- TabSim(sim)
show(tab)
q1 <- PlotSim(sim) + 
  ggtitle(expression(Synthetic~Surrogate~Correlation~rho==0.75)) +
  theme(legend.position = "none")


# Figure.
q <- cowplot::plot_grid(
  plotlist = list(q0, q1),
  nrow = 2,
  rel_heights = c(5, 4),
  labels = c("A", "B")
)

ggsave(
  file = "results/robust_sim.png",
  device = "png",
  width = 8.0,
  height = 8.0,
  units = "in",
  dpi = 480
)
