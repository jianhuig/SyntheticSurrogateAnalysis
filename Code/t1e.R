# Update to style of imputation file.

#' Simulation of Type I error and power 
#'
#'
#' @param n Sample size for the number of complete cases.
#' @param reps Simulation replicates.
#' @param missing_rate Missingness rates to consider.
#' @param rho Synthetic surrogate-target correlations.
#' @param pve_g Proportion of variance explained by genetic cov G
#' @return Data.frame.
Sim <- function(n = 1e3, reps = 1e3, missing_rate = 0,
                rho = 0, pve_g) {
  
  outer_loop <- lapply(seq_len(reps), function(i) {
    
    # Generate data under the null
    data <- DGP(n0 = n, miss = missing_rate, rho = rho, pve_g = pve_g)
    
    # Oracle estimator.
    oracle_est <- OracleEst(data)
    
    # Observed data estimator.
    obs_est <- ObsEst(data)
    
    # Synthetic surrogate estimator.
    ss_est <- SynSurrEst(data)
    
    # Collect oracle and observed.
    ests <- rbind(oracle_est, obs_est, ss_est)
    
    out <- data.frame(
      miss = missing_rate,
      rho = rho,
      method = c("oracle", "obs", "ss"),
      est = ests[, 1],
      se = ests[, 2],
      row.names = NULL
    )
    
    return(out)
    })
  
  do.call(rbind, outer_loop)
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

# ---------------------------------------------------------------------------
# Run the type I error simulation.
# pve_g = 0, under the null
# Loop over all missing rate and rho combinations
# Outer loop missing_rate over (0, 0.25, 0.5, 0.75)
# Inner loop rho over (0, 0.25, 0.5, 0.75)
outer_loop <- lapply(c(0, 0.25, 0.5, 0.75),function(RHO){
  inner_loop <- lapply(c(0, 0.25, 0.5, 0.75),function(MISS){
    Sim(n = 1e3, reps = 1e3, missing_rate = MISS, rho = RHO, pve_g = 0)
  })
  do.call(rbind, inner_loop)}
)

t1e_result <- do.call(rbind, outer_loop)

# Tabulate result
t1e_result %>% filter(method == "ss") %>% mutate(chisq2 = (est/se)^2) %>% 
  mutate(t1e = chisq2 > 3.841) %>% group_by(miss, rho) %>% 
  summarise_at(c("t1e","chisq2"), mean)

# Run the power simulation.
# pve_g = 0.12, under the alternative
# Loop over all missing rate and rho combinations
# Outer loop missing_rate over (0, 0.25, 0.5, 0.75)
# Inner loop rho over (0, 0.25, 0.5, 0.75)
outer_loop <- lapply(c(0, 0.25, 0.5, 0.75),function(RHO){
  inner_loop <- lapply(c(0, 0.25, 0.5, 0.75),function(MISS){
    Sim(n = 1e3, reps = 1e3, missing_rate = MISS, rho = RHO, pve_g = 0.005)
  })
  do.call(rbind, inner_loop)}
)

power_result <- do.call(rbind, outer_loop)

# Tabulate result
power_result %>% filter(method == "ss") %>% mutate(chisq2 = (est/se)^2) %>% 
  mutate(power = chisq2 > 3.841) %>% group_by(miss, rho) %>% 
  summarise_at(c("power","chisq2"), mean)









s