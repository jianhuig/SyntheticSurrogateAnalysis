#' Purpose: Examine bias and variance of multiple imputation procedure
#' under correct and incorrect model specification.
#' Updated: 2023-07-25
#' 
setwd("~/Documents/GitHub/SyntheticSurrogateAnalysis")
# Lab/Projects/Synthetic Surrogates/")


# devtools::install_github(repo = 'zrmacc/SurrogateRegression')

library(dplyr)
library(ggplot2)
library(ggsci)
library(RNOmni)
library(SurrogateRegression)

# -----------------------------------------------------------------------------


#' Simulation
#' 
#' @param n0 Size of the model-building set.
#' @param n1 Size of the inference set.
#' @param reps Simulation replicates.
#' @param train_indep Train on independent data? 
#' @return Data.frame.
Sim <- function(
  n0 = 1e3,
  n1 = 1e4,
  reps = 500,
  train_indep = TRUE,
  ycut_u,
  ycut_l
) {
  
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
    if (train_indep) {
      train_data <- DGP(n0, ycut_u = ycut_u, ycut_l = ycut_l)
      eval_data <- DGP(n1, ycut_u = ycut_u, ycut_l = ycut_l)
    } else {
      eval_data <- DGP(n1, ycut_u = ycut_u, ycut_l = ycut_l)
      train_data <- eval_data
    }
    
    
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
#' @param miss Missingness rate.
#' @param rho Correlation of G with X
#' @return Data.frame.
DGP <- function(
  n,
  bg = sqrt(0.01),
  bx = -0.50,
  miss = 0.25,
  rho = sqrt(0.5),
  ycut_u = 0.9,
  ycut_l = 0.1
) {
  
  # Genotype and covariates. G and Z have correlation rho_gz.
  g <- stats::rnorm(n)
  x <- sqrt(1 - rho^2) * stats::rnorm(n) + rho * g
  
  design <- cbind(g, x)
  eta <- design %*% c(bg, bx)
  ve <- as.numeric(var(eta))
  
  # Oracle outcome.
  y <- eta + stats::rnorm(n, sd = sqrt(1 - ve))
  
  # Outcome with missingness based on the value of y.
  # draw <- sample(seq_len(n), size = round(miss * n), replace = FALSE)
  draw_u <-  which(y >= quantile(y, ycut_u))
  draw_l <-  which(y <= quantile(y, ycut_l))
  draw <- c(draw_l, draw_u)
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
  fit <- SurrogateRegression::FitBNR(
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
      ase = sqrt(mean(se^2)),
      ese = sd(est),
      est = mean(est),
      alower = est - 2 * ase,
      aupper = est + 2 * ase,
      elower = est - 2 * ese,
      eupper = est + 2 * ese
    )
  return(tab)
}


#' Plot Simulation Results
#' 
#' @param sim Data.frame of simulation results.
#' @param bg Target effect size.
#' @return ggplot.
PlotSim <- function(sim, ana_se = TRUE, bg = 0.10) {
  
  df <- TabulateSim(sim)
  
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
      "Standard",
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
  
  df2 <- df %>%
    dplyr::select(x, method, config, est, ase, ese) %>%
    tidyr::pivot_longer(
      cols = c("ase", "ese"),
      names_to = "se_method",
      values_to = "se"
    ) %>%
    dplyr::mutate(
      lower = est - 2 * se,
      upper = est + 2 * se
    )
  df2$se_method <- factor(
    df2$se_method,
    levels = c("ase", "ese"),
    labels = c("Analytical", "Empirical")
  )
  
  pal <- ggsci::pal_d3()(5)
  q <- ggplot(data = df2) + 
    theme_bw() + 
    theme(legend.position = "top") + 
    geom_point(
      aes(x = x, y = est, color = method, group = se_method),
      size = 4,
      position = position_dodge(width = 0.8)
    ) + 
    geom_linerange(
      aes(x = x, ymin = lower, ymax = upper, linetype = se_method),
      position = position_dodge(width = 0.8)
    ) +
    scale_color_manual(
      name = NULL,
      labels = c(
        "Multiple\nImputation",
        "Standard",
        "Oracle",
        "Single\nImputation",
        "SynSurr"
      ),
      values = pal
    ) + 
    scale_linetype_manual(
      name = "CI",
      values = c("dotted", "solid")
    ) +
    xlab("Estimator") +
    ylab("Estimate") + 
    geom_hline(
      yintercept = bg,
      linetype = "dashed",
      color = "gray"
    )
  return(q)
}


# -----------------------------------------------------------------------------

# Simulation with correctly specified imputation model.
set.seed(110)
sim <- Sim(n0 = 1e3, n1 = 1e4, train_indep = TRUE, 
           ycut_u = ycut_u, ycut_l = ycut_l)

q <- PlotSim(sim, ana_se = TRUE) +
  scale_y_continuous(
    breaks = seq(from = 0.0, to = 0.15, by = 0.05),
    limits = c(-0.05, 0.15)
  )
show(q)

ggsave(
  plot = q,
  file = paste0("results/imputation_sim_mnar_", setting, ".png"),
  device = "png",
  width = 8.0,
  height = 4.0,
  units = "in",
  dpi = 480
)


# Tabulate results.
tab <- TabulateSim(sim)
show(tab)
file <- paste0("Data/imputation_sim_tab_mnar_", setting, ".tsv")
data.table::fwrite(
  x = tab,
  file = file,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Tabulate result.
# -----------------------------------------------------------------------------

sim <- data.table::fread(file = paste0("Data/imputation_sim_tab_mnar_", 
                                       setting, ".tsv"))
tab <- TabulateSim(sim)

order <- c(5, 4, 7, 2, 10, 6, 1, 9, 8, 3, 11)
ordered_tab <- tab[order, ]

ordered_tab <- ordered_tab %>%
  dplyr::select(method, config, est, ase, ese) %>%
  dplyr::rename(
    Method = method,
    "Imputation Model" = config,
    Estimate = est,
    "Analytical SE" = ase,
    "Empirical SE" = ese
  )
print(xtable::xtable(ordered_tab, digits = 3), include.rownames = FALSE)
