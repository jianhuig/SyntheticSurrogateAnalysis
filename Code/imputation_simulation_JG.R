#' Purpose: Examine bias and variance of multiple imputation procedure
#' under correct and incorrect model specification.
#' Updated: 2022-06-25
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------

#' Data Generating Process
#'
#' Generate data. "g" is genotype, "x" and "z" are covariates, "y"
#' is the outcome, and "yobs" is the outcome after introduction of
#' missingness.
#'
#' @param n Sample size.
#' @param bg Genetic effect size.
#' @param bx Effect of X.
#' @param bz Effect of Z.
#' @param maf Minor allele frequency.
#' @param miss Missingness rate.
#' @param rho_xz Correlation of (X, Z).
#' @param ve Residual variance.
#' @return Data.frame.
DGP <- function(n,
                bg = 1.0,
                bx = -0.50,
                bz = 0.50,
                maf = 0.25,
                miss = 0.20,
                rho_xz = 0.5,
                ve = 1.0) {
  
  # Genotype and covariates. X and Z have correlation rho_xz.
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  x <- stats::rnorm(n)
  z <- sqrt(1 - rho_xz^2) * stats::rnorm(n) + rho_xz * x
  
  design <- cbind(g, x, z)
  eta <- design %*% c(bg, bx, bz)
  
  # Oracle outcome.
  y <- eta + stats::rnorm(n, sd = sqrt(ve))
  
  # Outcome with missingness.
  draw <- sample(seq_len(n), size = round(miss * n), replace = FALSE)
  yobs <- y
  yobs[draw] <- NA
  
  out <- data.frame(
    g = g,
    x = x,
    z = z,
    y = y,
    yobs = yobs
  )
  return(out)
}


#' Fit Imputation Model
#'
#' Here the imputation model is correctly specified or mispecified.
#'
#' @param data Data.frame.
#' @return List containing the imputation coefficients and residual SD.
FitImpModel <- function(data, cov = c("g", "x", "z")) {
  fit <- lm(as.formula(paste("y", paste(cov, collapse = "+"), sep = "~")), data = data)
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
#' @return Data.frame.
GenSingleImp <- function(data, imp_param = NULL, cov = c("g", "x", "z"), beta = NULL) {
  data_miss <- data %>%
    dplyr::filter(is.na(yobs))
  
  design_miss <- data_miss %>%
    dplyr::select(cov) %>%
    as.matrix()
  
  if (is.null(beta)) {
    design_miss <- cbind(1, design_miss)
    eta <- design_miss %*% imp_param$coef
    yhat <- eta + stats::rnorm(length(eta), sd = imp_param$sigma)
  } else {
    eta <- design_miss %*% beta
    yhat <- eta + stats::rnorm(length(eta), sd = 1)
  }
  
  out <- data
  out$yhat <- data$yobs
  out$yhat[is.na(out$yhat)] <- yhat
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
    covar
  ) %>% na.omit()
  
  fit <- lm(y ~ ., data = df)
  fit_summary <- summary(fit)$coefficients
  
  out <- fit_summary[
    rownames(fit_summary) == "g", c("Estimate", "Std. Error")
  ]
  names(out) <- c("est", "se")
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
#' @return Numeric vector containing the estimate "est" and standard error "se".
MI <- function(data,
               imp_param = NULL,
               n_imp = 10,
               cov = c("g", "x", "z"),
               beta = NULL) {
  
  # Run n_imp imputations.
  results <- lapply(seq_len(n_imp), function(i) {
    single_imp <- GenSingleImp(data, imp_param, cov = cov, beta = beta)
    est <- EstBetaG(
      covar = single_imp[, c("x", "z")],
      geno = single_imp$g,
      outcome = single_imp$yhat
    )
    return(est)
  })
  results <- do.call(rbind, results)
  
  # Final estimate and SE.
  bg_bar <- mean(results[, "est"])
  v_w <- mean(results[, "se"]^2) # Variance within.
  v_b <- var(results[, "est"]) # Variance between.
  v <- v_w + (1 + 1 / n_imp) * v_b # Total variance.
  bg_se <- sqrt(v)
  
  out <- c(est = bg_bar, se = bg_se)
  return(out)
}


#' Simulation
#'
#' Case of correctly specified imputation model.
#'
#' @param n Sample size.
#' @param reps Simulation replicates.
#' @return Data.frame.
Sim <- function(n = 1e3, reps = 1e3, cov = c("g", "x", "z"), beta = NULL, single = FALSE) {
  
  # Run inference procedures:
  # 1. Oracle
  # 2. Observed data
  # 3. Multiple imputation
  results <- lapply(seq_len(reps), function(i) {
    
    # Generate independent inference data.
    data <- DGP(n)
    
    # Oracle estimator.
    oracle_est <- EstBetaG(
      covar = data[, c("x")], # Only use x for the inference model.
      geno = data$g,
      outcome = data$y
    )
    names(oracle_est) <- paste0("oracle_", names(oracle_est))
    
    # Observed data estimator.
    obs_est <- EstBetaG(
      covar = data[, c("x", "z")],
      geno = data$g,
      outcome = data$yobs
    )
    names(obs_est) <- paste0("obs_", names(obs_est))
    
    # Fit imputation model on training data.
    if (is.null(beta)) {
      train_data <- DGP(n)
      imp_param <- FitImpModel(train_data, cov = cov)
      beta <- imp_param$coef
     
    } 
    
      # Single imputation estimator.
      if(isTRUE(single)){
        temp <- GenSingleImp(data = data, imp_param = NULL, cov = cov, beta = beta)
        mi_est <- EstBetaG(
          cov = temp[, c("x","z")],
          geno = temp$g,
          outcome = temp$yhat
        )
      } else {
        # Multiple imputation estimator.
        mi_est <- MI(data, imp_param = NULL, cov = cov, beta = beta)
      }
      names(mi_est) <- paste0("mi_", names(mi_est))
    
    
    out <- c(oracle_est, obs_est, mi_est)
    return(out)
  })
  results <- do.call(rbind, results)
  return(data.frame(results))
}


#' Plot Simulation Results
#'
#' @param sim Data.frame of simulation results.
#' @return ggplot.
PlotSim <- function(sim, single = FALSE) {
  if(isTRUE(single)){
    labels <- c("Oracle", "Single\nImputation", "Observed\nData")
  } else {
    labels <- c("Oracle", "Multiple\nImputation", "Observed\nData")
  }
  
  df <- sim %>%
    dplyr::select(oracle_est, obs_est, mi_est) %>%
    tidyr::pivot_longer(
      cols = tidyr::everything(),
      names_to = "x",
      values_to = "Estimate"
    ) %>%
    dplyr::mutate(
      x = gsub("_est", "", x),
      Estimator = factor(
        x,
        levels = c("oracle", "mi", "obs"),
        labels = labels
      )
    )
  
  q <- ggplot(data = df) +
    theme_bw() +
    geom_boxplot(
      aes(x = Estimator, y = Estimate, fill = Estimator),
      show.legend = FALSE
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
  mu <- sim %>%
    dplyr::select(oracle_est, obs_est, mi_est) %>%
    tidyr::pivot_longer(
      cols = tidyr::everything(),
      names_to = "x",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      estimator = gsub("_est", "", x)
    ) %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(est = mean(y))
  
  se <- sim %>%
    dplyr::select(oracle_se, obs_se, mi_se) %>%
    tidyr::pivot_longer(
      cols = tidyr::everything(),
      names_to = "x",
      values_to = "y"
    ) %>%
    dplyr::mutate(
      estimator = gsub("_se", "", x)
    ) %>%
    dplyr::group_by(estimator) %>%
    dplyr::summarise(se = sqrt(mean(y^2)))
  
  tab <- mu %>%
    dplyr::inner_join(se, by = "estimator")
  return(tab)
}



