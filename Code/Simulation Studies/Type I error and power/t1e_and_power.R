# Update to style of imputation file.
library(dplyr)
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
                rho = 0, pve_g=0) {
  outer_loop <- pbapply::pblapply(seq_len(reps), function(i) {
    
    # Generate data under the null
    data <- DGP(n0 = n, miss = missing_rate, rho = rho, pve_g = 0)
    
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
      pve_g = 0,
      est = ests[, 1],
      se = ests[, 2],
      row.names = NULL
    )
    
    return(out)
    }, cl = cl)
  
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
CMplot::CMplot(try, plot.type = 'q')

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
    rownames(fit_summary) == "g", c("Estimate", "Pr(>|t|)", "Std. Error")]
  names(out) <- c("est", "p", "se")
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
  fit <- SurrogateRegression::Fit.BNR(
    t = data$yobs,
    s = data$s, 
    X = data[, c("g", "x")],
  )
  param <- coef(fit, type = "Target")
  
  out <- c(
    ss_est = param$Point[1],
    ss_p = param$p[1],
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
outer_loop <- pbapply::pblapply(c(0, 0.25, 0.5, 0.75),function(RHO){
  inner_loop <- lapply(c(0, 0.25, 0.5, 0.75, 0.9),function(MISS){
    Sim(n = 1e3, reps = 1e5, missing_rate = MISS, rho = RHO, pve_g = 0)
  })
  do.call(rbind, inner_loop)}
)

t1e_result <- do.call(rbind, outer_loop)

t1e_result <- readRDS("../Data/t1e_result.rds") %>% filter(method == "ss")
library(ggplot2)
color_ramp <- colorRampPalette(colors = c("white", "#0073C2FF"))
biv_palette <- color_ramp(n = 9)[c(3, 5, 7, 9)]
missing_names <- c(`0` = "Missing Rate = 0%",
                   `0.25` = "Missing Rate = 25%",
                   `0.5` = "Missing Rate = 50%",
                   `0.75` = "Missing Rate = 75%",
                   `0.9` = "Missing Rate = 90%"
)
p1 <- t1e_result %>% 
  group_by(miss, rho) %>%
  mutate(expected = -log10(ppoints(n())), observed = -log10(sort(se)),
         clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n(), shape2 = n():1)),
         cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n(), shape2 = n():1))) %>%
  ggplot()+
  geom_point(aes(expected, observed, color = factor(rho)))+
  xlab(expression(paste("Expected -log"[10], plain(P))))+
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  guides(color = guide_legend(title = expression(rho))) +
  scale_color_manual(name = expression(rho ~ "(%)"), values = biv_palette)+ 
  geom_abline(intercept = 0, slope = 1, color = "red")+
  geom_line(aes(expected, cupper), linetype = 2) +
  geom_line(aes(expected, clower), linetype = 2) +
  theme_bw()+
  facet_wrap(~miss, labeller = as_labeller(missing_names))

ggsave(
  plot = p1,
  filename = "qq_sim.png",
  device = "png",
  height = 6,
  width = 7.5,
  units = "in",
  dpi = 360
)

# Run the power simulation for different SNP-heritability level
# Unde the alternative
snp_heritability <- seq(0.001, 0.01, 0.001)
# Loop over all SNP-heritability, missing rate and rho combinations
# Outer Outer loop SNP heritability over seq(0.001, 0.01, 0.001)
# Inner loop missing_rate over (0, 0.25, 0.5, 0.75)
# Inner inner loop rho over (0, 0.25, 0.5, 0.75)
parallel <- FALSE
if (!parallel) {
  outer_loop <- pbapply::pblapply(snp_heritability, function(H) {
    inner_loop <- lapply(c(0, 0.25, 0.5, 0.75), function(RHO) {
      inner_inner_loop <- lapply(c(0, 0.25, 0.5, 0.75, 0.9), function(MISS) {
        Sim(n = 1e3, reps = 1e4, missing_rate = MISS, rho = RHO, pve_g = H)
      })
      do.call(rbind, inner_inner_loop)
    })
    do.call(rbind, inner_loop)
  })
} else {
  cl <- parallel::makeCluster(parallel::detectCores())
  parallel::clusterExport(cl, varlist = as.vector(lsf.str()))
  parallel::clusterEvalQ(cl, library(dplyr))
  outer_loop <- pbapply::pblapply(snp_heritability, function(H) {
    inner_loop <- lapply(c(0, 0.25, 0.5, 0.75), function(RHO) {
      inner_inner_loop <- lapply(c(0, 0.25, 0.5, 0.75, 0.9), function(MISS) {
        Sim(n = 1e3, reps = 1e4, missing_rate = MISS, rho = RHO, pve_g = H)
      })
      do.call(rbind, inner_inner_loop)
    })
    do.call(rbind, inner_loop)
  }, cl = cl)
  closeAllConnections()
}
power_result <- do.call(rbind, outer_loop)

# -----------------------------------------------------------------------------
# Tables and Figures for the Manuscript

library(xtable)
# Table for main text  chi2, t1e + chi2, power for snp heritability 0.005
t1e_table <- t1e_result %>% filter(method == "ss") %>% mutate(chisq2 = (est/se)^2) %>% 
  mutate(t1e = chisq2 > qchisq(0.95, df = 1)) %>% group_by(miss, rho) %>% 
  summarise_at(c("t1e","chisq2"), mean)

power_table <- power_result %>% filter(method == "ss" & pve_g==0.005) %>% mutate(chisq2 = (est/se)^2) %>% 
  mutate(power = chisq2 > qchisq(0.95, df = 1)) %>% group_by(miss, rho) %>% summarise_at(c("power","chisq2"), mean)

print(xtable::xtable(t1e_table %>% 
                 inner_join(power_table, by = c("miss", "rho")) %>% 
                   mutate(miss = miss*100) %>%
                   setNames(c("Missing Rate (\\%)", "$\\rho$",  "Type I Error",
                              "$\\chi^2$",   "Power", "$\\chi^2$")),
                 digits = c(0, 0, rep(2, 5))),
      sanitize.colnames.function = identity,
      include.rownames = FALSE)


# Figure for main text with power curves
color_ramp <- colorRampPalette(colors = c("white", "#0073C2FF"))
biv_palette <- color_ramp(n = 9)[c(3, 5, 7, 9)]

power_plot <- power_result %>%
  filter(method == "ss") %>%
  mutate(chisq2 = (est / se)^2) %>%
  mutate(power = chisq2 > qchisq(0.95, df = 1)) %>%
  group_by(miss, rho, pve_g) %>%
  summarise_at("power", mean) %>%
  ggplot(aes(
    x = pve_g*100, y = power, group = rho,
    color = factor(rho)
  )) +
  geom_point() +
  geom_line() +
  facet_wrap(.~miss, labeller = as_labeller(missing_names))+
  xlab("SNP Heritability (%)") +
  ylab("Power") +
  ylim(0, 1) +
  guides(color = guide_legend(title = expression(rho))) +
  scale_color_manual(name = expression(rho ~ "(%)"), values = biv_palette)+ theme_bw()

ggsave(
  plot = power_plot,
  filename = "t1e_power_power_curves.png",
  device = "png",
  height = 6,
  width = 7.5,
  units = "in",
  dpi = 360
)


# Figure 1 for efficiency analysis
missing_names <- c(`0` = "Missing Rate = 0%",
                   `0.25` = "Missing Rate = 25%",
                   `0.5` = "Missing Rate = 50%",
                   `0.75` = "Missing Rate = 75%",
                   `0.9` = "Missing Rate = 90%"
)

re_plot <- power_result %>%
  filter(!method == "oracle") %>%
  filter(pve_g == 0.005) %>%
  mutate(id = rep(1:(n()/2), each = 2)) %>%
  select(-est) %>%
  tidyr::pivot_wider(names_from = method, values_from = se) %>%
  mutate(re = (obs/ss)^2) %>%
  group_by(miss, rho) %>%
  summarize_at("re", mean) %>%
  ggplot(aes(
    x = 100*miss, y = re, group = rho,
    color = factor(rho)
  )) +
  geom_point() +
  geom_line() +
  xlab("Missing Rate (%)") +
  ylab("Relative Efficiency") +
  guides(color = guide_legend(title = expression(rho))) +
  scale_x_continuous(breaks = c(0,25,50,75,90))+
  scale_color_manual(name = expression(rho ~ "(%)"), values = biv_palette)+theme_bw()

ggsave(
  plot = re_plot,
  filename = "t1e_power_RE.png",
  device = "png",
  height = 6,
  width = 7.5,
  units = "in",
  dpi = 360
)

# Figure 2 for efficiency analysis
re_h2_plot <- power_result %>%
  filter(!method == "oracle") %>%
  mutate(id = rep(1:(n()/2), each = 2)) %>%
  select(-est) %>%
  tidyr::pivot_wider(names_from = method, values_from = se) %>%
  mutate(re = (obs/ss)^2) %>%
  group_by(miss, rho, pve_g) %>%
  summarize_at("re", mean) %>%
  ggplot(aes(x = pve_g*100, y = re, group = factor(rho), color = factor(rho))) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ miss, labeller = as_labeller(missing_names)) +
  scale_color_manual(name = expression(rho), values = biv_palette) +
  xlab("SNP Heritability (%)") +
  scale_x_continuous(breaks = c(0.1,0.3,0.5,0.7,0.9,1.1))+
  ylab("Relative Efficiency")+
  theme_bw()

ggsave(
  plot = re_h2_plot,
  filename = "t1e_power_RE_h2.png",
  device = "png",
  height = 6,
  width = 7.5,
  units = "in",
  dpi = 360
)
