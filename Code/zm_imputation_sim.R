#' Purpose: Examine bias and variance of multiple imputation procedure
#' under correct and incorrect model specification.
#' Updated: 2022-07-25
setwd("~/Documents/Lab/Projects/Synthetic Surrogates/")

library(dplyr)
library(ggplot2)
library(ggsci)
library(RNOmni)
library(SurrogateRegression)

# -----------------------------------------------------------------------------


#' Simulation
#' 
#' @param n Sample size.
#' @param reps Simulation replicates.
#' @return Data.frame.
Sim <- function(n = 1e3, reps = 1e3) {
  
  # Imputation scenarios, listing which covariates to include.
  scenarios <- list(
    c("g", "x"),
    c("g"),
    c("x")
  )
  
  # Outer loop iterates over realizations of the data.
  # For each outer loop, the following two estimators are calculated:
  # 1. Oracle, based on the true Y (without missing values).
  # 2. Observed data, based on only the observed values of Y.
  # Note that these estimators do not depend on the imputation procedure.
  
  # Next, an inner loop iterates over imputation models:
  # 1. Correclty specified, including G and X.
  # 2. Misspecified: G only.
  # 3. Misspecified: X only.

  # For each inner loop, the following estimates are calculated:
  # 1. Single imputation.
  # 2. Multiple imputation.
  # 3. Synthetic surrogate.
  
  outer_loop <- lapply(seq_len(reps), function(i) {
    
    # Generate data. 
    # The same data sets are used for each imputation scenario
    #  to make the results comparable.
    train_data <- DGP(n)
    eval_data <- DGP(n)
    
    # Oracle estimator.
    oracle_est <- OracleEst(eval_data)
    
    # Observed data estimator.
    obs_est <- ObsEst(eval_data)
    
    # Collect oracle and observed.
    out <- data.frame(
      config = c(NA, NA),
      method = c("oracle", "obs"),
      est = c(oracle_est[1], obs_est[1]),
      se = c(oracle_est[2], obs_est[2]),
      row.names = NULL
    )
    
    # Inner loop over imputation models.
    imp_est <- lapply(scenarios, function(config) {
      
      y_config <- c("y", config)
      
      # Train imputation model.
      imp_data <- train_data %>% 
        dplyr::select(dplyr::all_of(y_config))
      imp_param <- FitImpModel(imp_data)
      
      # Single imputation.
      imp_eval_data <- GenSingleImp(
        eval_data,
        imp_param,
        add_error = FALSE
      )
      si_est <- ImpEst(imp_eval_data)
      
      # Multiple imputation. 
      mi_est <- MI(eval_data, imp_param)
      
      # SynSurr.
      ss_est <- SynSurrEst(eval_data, imp_param)
      tmp <- rbind(si_est, mi_est, ss_est)
      
      # Collect results.
      out <- data.frame(
        config = paste(config, collapse = ""),
        method = c("si", "mi", "ss"),
        est = tmp[, 1],
        se = tmp[, 2],
        row.names = NULL
      )
      return(out)
    })
    imp_est <- do.call(rbind, imp_est)
    
    out <- rbind(out, imp_est)
    return(out)
  })
  results <- do.call(rbind, outer_loop)
  return(results)
}


# -----------------------------------------------------------------------------

#' Data Generating Process
#' 
#' Generate data. "g" is genotype, "x" is a covariates, "y"
#' is the outcome, and "yobs" is the outcome after introduction of
#' missingness.
#' 
#' @param n Sample size.
#' @param bg Genetic effect size.
#' @param bx Effect of X. 
#' @param maf Minor allele frequency.
#' @param miss Missingness rate.
#' @param rho Correlation of G with X
#' @param ve Residual variance.
#' @return Data.frame.
DGP <- function(
    n,
    bg = 1.0,
    bx = -0.50,
    maf = 0.25,
    miss = 0.25,
    rho = sqrt(0.5),
    ve = 1.0
) {
  
  # Genotype and covariates. G and Z have correlation rho_gz.
  g <- stats::rnorm(n)
  x <- sqrt(1 - rho^2) * stats::rnorm(n) + rho * g
  
  design <- cbind(g, x)
  eta <- design %*% c(bg, bx)
  
  # Oracle outcome.
  y <- eta + stats::rnorm(n, sd = sqrt(ve))
  
  # Outcome with missingness.
  draw <- sample(seq_len(n), size = round(miss * n), replace = FALSE)
  yobs <- y
  yobs[draw] <- NA
  
  out <- data.frame(
    g = g,
    x = x,
    y = y,
    yobs = yobs
  )
  return(out)
}


# -----------------------------------------------------------------------------


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
    covar
  ) %>% na.omit()
  
  fit <- lm(y ~ ., data = df)
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


#' Single Imputation Estimator
#' 
#' @param data Data.frame.
#' @return Numeric vector containing "si_est" and "si_se".
ImpEst <- function(data) {
  si_est <- EstBetaG(
    covar = data$x,
    geno = data$g,
    outcome = data$yhat
  )
  names(si_est) <- paste0("si_", names(si_est))
  return(si_est)
}


# -----------------------------------------------------------------------------


#' Fit Imputation Model
#' 
#' @param data Data.frame.
#' @return List containing the imputation coefficients and residual SD.
FitImpModel <- function(data) {
  fit <- lm(y ~ ., data = data)
  fit_summary <- summary(fit)
  out <- list()
  out$coef <- coef(fit)
  out$sigma <- fit_summary$sigma
  return(out)
}


#' Generate Imputation
#' 
#' Generate a single imputation of y.
#' 
#' @param data Data.frame.
#' @param imp_param Imputation model parameters.
#' @param add_error Add residual error? TRUE for MI.
#' @return Data.frame.
GenSingleImp <- function(data, imp_param, add_error = TRUE) {
  
  data_miss <- data %>%
    dplyr::filter(is.na(yobs)) 
  
  # Determine which covariates are in the imputation model.
  covar_in_imp <- setdiff(names(imp_param$coef), "(Intercept)")
  
  design_miss <- data_miss %>%
    dplyr::select(dplyr::all_of(covar_in_imp)) %>%
    as.matrix()
  design_miss <- cbind(1, design_miss)
  
  eta <- design_miss %*% imp_param$coef
  if (add_error) {
    yhat <- eta + stats::rnorm(length(eta), sd = imp_param$sigma)
  } else {
    yhat <- eta
  }
  
  out <- data
  out$yhat <- data$yobs
  out$yhat[is.na(out$yhat)] <- yhat
  
  return(out)
}


#' Multiple Imputation
#' 
#' Impute the data n_imp times. For each imputation, estimate the genetic effect
#' and standard error. Combine results across the n_imp imputations, using the
#' mean to estimate the genetic effect, and Rubin's rules to estimate the
#' standard error.
#' 
#' @param data Data.frame.
#' @param imp_param Imputation parameters.
#' @param n_imp Number of imputations
#' @return Numeric vector containing "mi_est" and "mi_se".
MI <- function(
    data,
    imp_param,
    n_imp = 10
) {
  
  # Run n_imp imputations.
  results <- lapply(seq_len(n_imp), function(i) {
    single_imp <- GenSingleImp(data, imp_param)
    est <- EstBetaG(
      covar = single_imp$x,
      geno = single_imp$g,
      outcome = single_imp$yhat
    )
    return(est)
  })
  results <- do.call(rbind, results)
  
  # Final estimate and SE.
  bg_bar <- mean(results[, "est"])
  v_w <- mean(results[, "se"]^2)  # Variance within.
  v_b <- var(results[, "est"])  # Variance between.
  v <- v_w + (1 + 1 / n_imp) * v_b  # Total variance.
  bg_se <- sqrt(v)
  
  out <- c(mi_est = bg_bar, mi_se = bg_se)
  return(out)
}


#' SynSurr Estimator
#' 
#' @param data Data.frame.
#' @param imp_param Imputation parameters.
#' @return Numeric vector containing "ss_est" and "ss_se".
SynSurrEst <- function(
  data,
  imp_param
) {
  
  # Determine which covariates are in the imputation model.
  covar_in_imp <- setdiff(names(imp_param$coef), "(Intercept)")
  
  # Generate synthetic surrogate for all subjects.
  design <- data %>%
    dplyr::select(dplyr::all_of(covar_in_imp)) %>%
    as.matrix()
  design <- cbind(1, design)
  yhat <- design %*% imp_param$coef
  yhat <- RankNorm(as.numeric(yhat))
  
  # SynSurr estimator.
  fit <- SurrogateRegression::Fit.BNLS(
    t = data$yobs,
    s = yhat, 
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


#' Plot Simulation Results
#' 
#' @param sim Data.frame of simulation results.
#' @return ggplot.
PlotSim <- function(sim) {
  df <- sim %>%
    dplyr::group_by(method, config) %>%
    dplyr::summarise(
      est = mean(est),
      se = sqrt(mean(se^2)),
      lower = est - 2 * se,
      upper = est + 2 * se
    )
  df$x <- c(
    "mi_g",
    "mi_gx",
    "mi_x",
    "obs",
    "oracle",
    "si_g",
    "si_gx",
    "si_x",
    "ss_g",
    "ss_gx",
    "ss_x"
  )
  df$x <- factor(
    df$x,
    levels = c(
      "oracle",
      "obs",
      "mi_gx",
      "si_gx",
      "ss_gx",
      "mi_g",
      "si_g",
      "ss_g",
      "mi_x",
      "si_x",
      "ss_x"
    ),
    labels = c(
      "Oracle",
      "Marginal",
      "MI\n(G,X)",
      "SI\n(G,X)",
      "SS\n(G,X)",
      "MI\n(G)",
      "SI\n(G)",
      "SS\n(G)",
      "MI\n(X)",
      "SI\n(X)",
      "SS\n(X)"
    ),
    ordered = TRUE
  )
  
  q <- ggplot(data = df) + 
    theme_bw() + 
    theme(legend.position = "top") + 
    geom_col(
      aes(x = x, y = est, fill = method),
    ) + 
    scale_fill_jco(
      name = NULL,
      labels = c(
        "Multiple\nImputation",
        "Marginal",
        "Oracle",
        "Single\nImputation",
        "Synthetic\nSurrogate"
      )
    ) +
    geom_errorbar(
      aes(x = x, ymin = lower, ymax = upper),
      width = 0.5
    ) + 
    xlab("Estimator") +
    ylab("Estimate") + 
    geom_hline(
      yintercept = 1.0,
      linetype = "dashed",
      color = "gray"
    )
  return(q)
}


#' Tabulate Simulation Results
#' 
#' Tabulates the mean estimate and root-mean-square standard error.
#' 
#' @param sim Data.frame of simulation results.
#' @return Data.frame.
TabulateSim <- function(sim) {
  tab <- sim %>%
    dplyr::group_by(method, config) %>%
    dplyr::summarise(
      est = mean(est),
      se = sqrt(mean(se^2)),
      lower = est - 2 * se,
      upper = est + 2 * se
    )
  return(tab)
}


# -----------------------------------------------------------------------------

# Simulation with correctly specified imputation model.
seed <- 110
set.seed(seed)
sim <- Sim(n = 1000, reps = 200)

# Bar plots.
q <- PlotSim(sim)
show(q)
ggsave(
  plot = q,
  file = "results/imputation_sim_s110.png",
  device = "png",
  width = 8.0,
  height = 4.0,
  units = "in",
  dpi = 480
)

# Tabulate results.
tab <- TabulateSim(sim)
show(tab)
file <- paste0("results/imputation_sim_s", seed, ".tsv")
data.table::fwrite(
  x = tab,
  file = file,
  sep = "\t"
)
