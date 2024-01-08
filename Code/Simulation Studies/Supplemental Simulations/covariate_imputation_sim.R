#' Purpose: Examine the effect of imputing an input to the surrogate model.
#' Updated: 2023-08-22
setwd("~/Documents/Lab/Projects/Synthetic Surrogates/")

library(dplyr)
library(ggplot2)
library(ggsci)
library(RNOmni)
library(SurrogateRegression)


# -----------------------------------------------------------------------------


#' Data Generating Process
#' 
#' X and Z are correlated. Z is utilized to generate the synthetic surrogate. 
#' 
#' @param n0 Number of subjects with observed target outcomes.
#' @param bg Genetic effect size. 
#' @param miss_y Missingness of Y.
#' @param miss_z Missingness of Z. 
#' @param rho_gx Correlation between G and X.
#' @param rho_xz Correlation between X and Z.
#' @return Data.frame.
DGP <- function(
  n0,
  bg = 0,
  miss_y = 0.5,
  miss_z = 0.0,
  rho_gx = 0.5,
  rho_xz = 0.5
) {
  
  # Total sample size.
  n <- n0 / (1 - miss_y)
  
  # Genotype and covariates.
  g <- stats::rnorm(n)
  # x <- stats::rnorm(n)
  x <- rho_gx * g + sqrt(1 - rho_gx^2) * stats::rnorm(n) 
  z <- rho_xz * x + sqrt(1 - rho_xz^2) * stats::rnorm(n)
  
  # Generate polynomial basis for x. 
  z_deg <- 3
  z_basis <- poly(z, degree = z_deg)
  z_basis <- scale(z_basis, center = FALSE)
  
  # Coefficients for x and z. 
  vx <- 0.1
  bx <- sqrt(vx)
  
  vz <- 0.3
  bz <- rep(sqrt(vz / z_deg), z_deg)

  design <- cbind(g, x, z_basis)
  eta <- design %*% c(bg, bx, bz)
  
  # Residual variation.
  ve <- 1 - var(as.numeric(eta))
  
  # Oracle outcome.
  y <- eta + stats::rnorm(n, sd = sqrt(ve))
  
  # Predictor missingness.
  draw <- sample(seq_len(n), size = round(miss_z * n), replace = FALSE)
  zobs <- z
  zobs[draw] <- NA
  
  # Outcome missingness.
  draw <- sample(seq_len(n), size = round(miss_y * n), replace = FALSE)
  yobs <- y
  yobs[draw] <- NA
  
  out <- data.frame(
    g = g,
    x = x,
    y = y,
    yobs = yobs,
    z = z,
    zobs = zobs
  )
  return(out)
}


#' Impute Missing Covariate
#'
#' @param data Data.frame.
#' @return Data.frame.
ImputeZ <- function(data) {
  
  fit <- lm(zobs ~ x, data = data)
  df_impute <- data %>%
    dplyr::filter(is.na(zobs)) %>%
    dplyr::select(x)
  
  zhat <- predict(fit, newdata = df_impute)
  
  data$zhat <- data$zobs
  data$zhat[is.na(data$zhat)] <- zhat
  return(data)
  
}


#' Fit SynSurr Model
#' 
#' @param data Data.frame.
#' @return Prediction function.
FitSynSurr <- function(data) {
  
  df <- data %>% 
    dplyr::select(zhat, y) %>%
    dplyr::filter(complete.cases(.))
  zhat <- df %>% dplyr::pull(zhat)
  y <- df %>% dplyr::pull(y)
  
  # Fit.
  fit <- lm(y ~ splines::bs(zhat, Boundary.knots = c(-5, 5)))
  
  # Define prediction function.
  Out <- function(data) {

    yhat <- predict(fit, newdata = data)
    return(yhat)
    
  }
  return(Out)
}


#' SynSurr Inference
#' 
#' @param data Data.frame.
#' @return SynSurr results.
SynSurrGWAS <- function(data) {
  
  covar <- cbind(1, data$g, data$x)
  colnames(covar) <- c("int", "g", "x")
  fit <- SurrogateRegression::FitBNR(
    t = data$yobs,
    s = data$yhat,
    X = covar,
    is_zero = c(FALSE, TRUE, FALSE)
  )
  
  # Format results.
  coef <- fit@Regression.tab
  out <- data.frame(
    bg = coef$Point[2],
    bg_p = coef$p[2]
  )
  
  return(out)
}


# -----------------------------------------------------------------------------


#' Simulation Instance
#' 
#' @param miss_y Target missingness.
#' @param miss_z Predictor missingness.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set. 
#' @return Data.frame.
SimInstance <- function(
  miss_y,  
  miss_z,
  bg = 0.0,
  n = 1e3
) {
  
  # Training and GWAS data. 
  data_train <- DGP(n = n, bg = bg, miss_y = 0.0, miss_z = 0.0)
  data_gwas <- DGP(n = n, bg = bg, miss_y = miss_y, miss_z = miss_z)
  
  # Impute Z.
  data_train <- ImputeZ(data_train)
  data_gwas <- ImputeZ(data_gwas)
  
  # Fit SynSurr model. 
  synsurr_model <- FitSynSurr(data_train)
  
  # Generate surrogate.
  data_gwas$yhat <- synsurr_model(data_gwas)
  
  # Run SynSurr.
  results <- SynSurrGWAS(data_gwas)
  results$miss_y <- miss_y
  results$miss_z <- miss_z
  return(results)
}


#' Simulation
#' 
#' @param miss_y Target missingness.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set. 
#' @param reps Simulation replicates.
#' @return Data.frame.
Sim <- function(
    miss_y,
    bg = 0.0,
    n = 1e3,
    reps = 500
) {
  
  levels_miss_z <- seq(from = 0.0, to = 0.50, length.out = 5)
  outer <- lapply(levels_miss_z, function(miss_z) {
    
    inner <- lapply(seq_len(reps), function(i) {
      results <- SimInstance(
        miss_y = miss_y, 
        miss_z = miss_z,
        bg = bg, 
        n = n
      )
    })
    inner <- do.call(rbind, inner)
    
  })
  outer <- do.call(rbind, outer)

  return(outer)
}


#' Compile Files
#' 
#' @param file_dir File directory.
#' @return Data.frame.
Compile <- function(file_dir) {
  
  init_dir <- getwd()
  setwd(file_dir)
  
  files <- dir()
  out <- lapply(files, function(f) {data.table::fread(file = f)})
  out <- do.call(rbind, out)
  
  setwd(init_dir)
  return(out)
}


# -----------------------------------------------------------------------------
# # Type I error simulation.
# -----------------------------------------------------------------------------


# Run simulation.
if (FALSE) {
  
  LOOPS <- 5
  REPS <- 1e4
  stem <- "results/covar_impute/t1e"
  levels_miss_y <- c(0.00, 0.25, 0.50, 0.75, 0.90)
  
  sink <- lapply(seq_len(LOOPS), function(i) {
    
    results <- lapply(levels_miss_y, function(miss_y) {
      
      sim <- Sim(miss_y = miss_y, reps = REPS)
      
      date <- Sys.Date()
      time <- format(Sys.time(), "%H%M%S")
      fout <- glue::glue("{stem}/m{miss_y}-{date}-{time}.tsv")
      
      data.table::fwrite(
        x = sim,
        file = fout,
        sep = "\t"
      )
      
      return(NULL)
    })
    
    return(NULL)
  })
  
}


results_file <- "results/covar_impute/covar_impute_sim_t1e.tsv"

if (file.exists(results_file)) {
  
  t1e_sim <- data.table::fread(results_file)
  
} else {
  
  file_dir <- "results/covar_impute/t1e/"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(miss_y, miss_z, bg_p) %>%
    dplyr::mutate(chi2 = qchisq(bg_p, df = 1)) %>%
    dplyr::group_by(miss_y, miss_z) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mu = mean(chi2),
      se = sqrt(var(chi2) / n),
      lower = mu - 2 * se,
      upper = mu + 2 * se
    )
  data.table::fwrite(
    x = agg,
    file = results_file,
    sep = "\t"
  )
  
  t1e_sim <- agg
  rm(unagg, agg)

}


# -----------------------------------------------------------------------------
# Bias simulation.
# -----------------------------------------------------------------------------


if (FALSE) {
  
  LOOPS <- 5
  REPS <- 2e3
  stem <- "results/covar_impute/bias"
  levels_miss_y <- c(0.00, 0.25, 0.50, 0.75, 0.90)
  
  sink <- lapply(seq_len(LOOPS), function(i) {
    
    results <- lapply(levels_miss_y, function(miss_y) {
      
      sim <- Sim(miss_y = miss_y, bg = 0.1, reps = REPS)
      
      date <- Sys.Date()
      time <- format(Sys.time(), "%H%M%S")
      fout <- glue::glue("{stem}/m{miss_y}-{date}-{time}.tsv")
      
      data.table::fwrite(
        x = sim,
        file = fout,
        sep = "\t"
      )
      
      return(NULL)
    })
    
    return(NULL)
  })
  
}


results_file <- "results/covar_impute/covar_impute_sim_bias.tsv"

if (file.exists(results_file)) {
  
  bias_sim <- data.table::fread(results_file)
  
} else {
  
  file_dir <- "results/covar_impute/bias/"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(miss_y, miss_z, bg) %>%
    dplyr::group_by(miss_y, miss_z) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mu = mean(bg),
      se = sqrt(var(bg) / n),
      lower = mu - 2 * se,
      upper = mu + 2 * se
    )
  data.table::fwrite(
    x = agg,
    file = results_file,
    sep = "\t"
  )
  
  bias_sim <- agg
  rm(unagg, agg)
  
}


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

if (!is.factor(t1e_sim$miss_y)) {
  t1e_sim$miss_y <- factor(
    t1e_sim$miss_y,
    levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
    labels = c("0%", "25%", "50%", "75%", "90%")
  )
  
  t1e_sim$miss_z <- factor(
    t1e_sim$miss_z,
    levels = c(0.00, 0.125, 0.250, 0.375, 0.50),
    labels = c("0%", "12.5%", "25%", "37.5%", "50%")
  )
}

if (!is.factor(bias_sim$miss_y)) {
  bias_sim$miss_y <- factor(
    bias_sim$miss_y,
    levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
    labels = c("0%", "25%", "50%", "75%", "90%")
  )
  
  bias_sim$miss_z <- factor(
    bias_sim$miss_z,
    levels = c(0.00, 0.125, 0.250, 0.375, 0.50),
    labels = c("0%", "12.5%", "25%", "37.5%", "50%")
  )
}


# Expected chi2.
q1 <- ggplot() + 
  theme_bw() +
  theme(legend.position = "top") +
  geom_hline(
    yintercept = 1.0,
    linetype = "dashed",
    color = "gray",
    linewidth = 1.2
  ) +
  geom_col(
    data = t1e_sim,
    aes(x = miss_z, y = mu, fill = miss_y),
    position = position_dodge()
  ) +
  geom_errorbar(
    data = t1e_sim,
    aes(x = miss_z, ymin = lower, ymax = upper, group = miss_y),
    position = position_dodge(0.9),
    width = 0.5
  ) + 
  scale_y_continuous(
    name = expression(Expected~X^2),
    limits = c(0.0, 1.4),
    breaks = seq(from = 0.0, to = 1.4, by = 0.2)
  ) +
  scale_x_discrete(
    name = "Predictor Missingness"
  ) + 
  scale_fill_brewer(
    name = "Target\nMissingness",
    palette = "Blues"
  )
  

# Bias.
q2 <- ggplot() + 
  theme_bw() +
  theme(legend.position = "top") +
  geom_hline(
    yintercept = 0.1,
    linetype = "dashed",
    color = "gray",
    linewidth = 1.2
  ) +
  geom_col(
    data = bias_sim,
    aes(x = miss_z, y = mu, fill = miss_y),
    position = position_dodge()
  ) +
  geom_errorbar(
    data = bias_sim,
    aes(x = miss_z, ymin = lower, ymax = upper, group = miss_y),
    position = position_dodge(0.9),
    width = 0.5
  ) + 
  scale_y_continuous(
    name = expression(Estimated~beta[G]),
    limits = c(0.0, 0.14),
    breaks = seq(from = 0.0, to = 0.14, by = 0.02)
  ) +
  scale_x_discrete(
    name = "Predictor Missingness"
  ) + 
  scale_fill_brewer(
    name = "Target\nMissingness",
    palette = "Blues"
  )

legend <- ggpubr::get_legend(q1)
q_legend <- ggpubr::as_ggplot(legend)

q1 <- q1 + theme(legend.position = "none")
q2 <- q2 + theme(legend.position = "none")

q_panel <- cowplot::plot_grid(
  plotlist = list(q_legend, q2, q1),
  labels = c("", "A", "B"),
  axis = "left",
  ncol = 1,
  rel_heights = c(1, 5, 5)
)
show(q_panel)

if (FALSE) {
  ggsave(
    plot = q_panel,
    file = "results/covariate_missingness.pdf",
    width = 9.0,
    height = 6.0
  )
}
