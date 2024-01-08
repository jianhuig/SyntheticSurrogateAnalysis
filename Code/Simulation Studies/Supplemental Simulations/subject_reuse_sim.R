#' Purpose: Examine the effect of reusing the GWAS data set for model training and GWAS. 
#' Updated: 2023-08-21
setwd("~/Documents/Lab/Projects/Synthetic Surrogates/")

library(dplyr)
library(ggplot2)
library(ggsci)
library(RNOmni)
library(SurrogateRegression)


# -----------------------------------------------------------------------------


#' Data Generating Process
#' 
#' @param n0 Number of subjects with observed target outcomes.
#' @param bg Genetic effect size. 
#' @param vx Proportion of variation explained by x. 
#' @param miss Missingness rate.
#' @param rho Genotype-covariate correlation.
#' @return Data.frame.
DGP <- function(
  n0,
  bg = 0,
  vx = 0.3, 
  miss = 0,
  rho = 0.9
) {
  
  # Total sample size.
  n <- n0 / (1 - miss)
  
  # Genotype and covariates. G and X have correlation rho.
  g <- stats::rnorm(n)
  x <- sqrt(1 - rho^2) * stats::rnorm(n) + rho * g
  
  # Generate polynomial basis for x. 
  x_deg <- 3
  x_basis <- poly(x, degree = x_deg)
  x_basis <- scale(x_basis, center = FALSE)
  
  # Coefficient for x. 
  bx <- rep(sqrt(vx / x_deg), x_deg)
  
  design <- cbind(g, x_basis)
  eta <- design %*% c(bg, bx)
  
  # Residual variation.
  ve <- 1 - var(as.numeric(eta))
  
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


#' Fit SynSurr Model
#' 
#' @param data Data.frame.
#' @param degree Model degree.
#' @param spline_model Uses a spline model if TRUE, a scatter plot smooth if FLASE.
#' @return Prediction function.
FitSynSurr <- function(
    data, 
    degree = 2,
    spline_model = TRUE
  ) {
  
  df <- data %>% dplyr::filter(complete.cases(.))
  x <- df %>% dplyr::pull(x)
  y <- df %>% dplyr::pull(y)
  
  # Fit.
  if (spline_model) {
    fit <- lm(y ~ splines::bs(x, degree = degree, Boundary.knots = c(-5, 5)))
  } else {
    fit <- loess(y ~ x, control = loess.control(surface = "direct"))
  }
  
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
#' @param cubic_assoc
#' @return SynSurr results.
SynSurrGWAS <- function(data, cubic_assoc = TRUE) {
  
  if (cubic_assoc) {
    covar <- cbind(1, data$g, poly(data$x, degree = 3))
    colnames(covar) <- c("int", "g", "x1", "x2", "x3")
  } else {
    covar <- cbind(1, data$g, data$x)
    colnames(covar) <- c("int", "g", "x")
  }

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
#' @param miss Target missingness.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set. 
#' @return Data.frame.
SimInstance <- function(
  miss,  
  bg = 0.0,
  n = 1e3
) {
  
  # Training and GWAS data. 
  data_train <- DGP(n = n, bg = bg, miss = 0.0)
  data_gwas <- DGP(n = n, bg = bg, miss = miss)
  
  # Fit SynSurr model. 
  synsurr_train <- FitSynSurr(data_train)
  synsurr_gwas <- FitSynSurr(data_gwas)
  
  # Generate surrogate.
  df_base <- data_gwas %>% dplyr::select(yobs, g, x)
  df_train <- df_gwas <- df_base
  df_train$yhat <- synsurr_train(data_gwas)
  df_gwas$yhat <- synsurr_gwas(data_gwas)
  
  # Run SynSurr.
  gwas_train <- SynSurrGWAS(df_train)
  gwas_train$indep_train <- TRUE
  gwas_gwas <- SynSurrGWAS(df_gwas)
  gwas_gwas$indep_train <- FALSE
  
  out <- rbind(gwas_train, gwas_gwas)
  out$miss <- miss
  return(out)
}


#' Simulation
#' 
#' @param miss Target missingness.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set. 
#' @param reps Simulation replicates.
#' @return Data.frame.
Sim <- function(
    miss,
    bg = 0.0,
    n = 1e3,
    reps = 500
) {
  
  out <- lapply(seq_len(reps), function(i) {
    
    results <- SimInstance(miss = miss, bg = bg, n = n)
    
  })
  out <- do.call(rbind, out)
  return(out)
  
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

# results_file <- "results/subject_reuse/subject_reuse_sim_t1e.tsv"
# results_file <- "results/subject_reuse/subject_reuse_sim_misspec_t1e.tsv"
# results_file <- "results/subject_reuse/subject_reuse_sim_high_corr_t1e.tsv"
# results_file <- "results/subject_reuse/subject_reuse_sim_loess_t1e.tsv"
results_file <- "results/subject_reuse/subject_reuse_sim_cubic_assoc_t1e.tsv"

if (file.exists(results_file)) {
  
  t1e_sim <- data.table::fread(results_file)
  
} else {
  
  file_dir <- "results/subject_reuse/t1e_cubic_assoc/"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(miss, indep_train, bg_p) %>%
    dplyr::mutate(chi2 = qchisq(bg_p, df = 1)) %>%
    dplyr::group_by(miss, indep_train) %>%
    dplyr::summarise(
      mu = mean(chi2),
      se = sqrt(var(chi2) / dplyr::n()),
      lower = mu - 2 * se,
      upper = mu + 2 * se
    )
  data.table::fwrite(
    x = agg,
    file = results_file,
    sep = "\t"
  )
  
}


# Run simulation.
if (FALSE) {
  
  LOOPS <- 5
  REPS <- 1e4
  stem <- "results/subject_reuse/t1e_loess/"
  miss_rates <- c(0.00, 0.25, 0.50, 0.75, 0.90)
  
  sink <- lapply(seq_len(LOOPS), function(i) {
    
    results <- lapply(miss_rates, function(m) {
      
      sim <- Sim(miss = m, reps = REPS)
      
      date <- Sys.Date()
      time <- format(Sys.time(), "%H%M%S")
      fout <- glue::glue("{stem}/m{m}-{date}-{time}.tsv")
      
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


# -----------------------------------------------------------------------------
# Bias simulation.
# -----------------------------------------------------------------------------

# results_file <- "results/subject_reuse/subject_reuse_sim_bias.tsv"
# results_file <- "results/subject_reuse/subject_reuse_misspec_sim_bias.tsv"
# results_file <- "results/subject_reuse/subject_reuse_high_corr_sim_bias.tsv"
# results_file <- "results/subject_reuse/subject_reuse_sim_loess_bias.tsv"
results_file <- "results/subject_reuse/subject_reuse_cubic_assoc_bias.tsv"

if (file.exists(results_file)) {
  
  bias_sim <- data.table::fread(results_file)
  
} else {
  
  file_dir <- "results/subject_reuse/bias_cubic_assoc/"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(miss, indep_train, bg) %>%
    dplyr::group_by(miss, indep_train) %>%
    dplyr::summarise(
      mu = mean(bg),
      se = sqrt(var(bg) / dplyr::n()),
      lower = mu - 2 * se,
      upper = mu + 2 * se
    )
  data.table::fwrite(
    x = agg,
    file = results_file,
    sep = "\t"
  )
  
}


if (FALSE) {
  
  LOOPS <- 2
  REPS <- 1e4
  stem <- "results/subject_reuse/bias_cubic_assoc/"
  miss_rates <- c(0.00, 0.25, 0.50, 0.75, 0.90)
  
  sink <- lapply(seq_len(LOOPS), function(i) {
    
    results <- lapply(miss_rates, function(m) {
      
      sim <- Sim(miss = m, bg = 0.1, reps = REPS)
      
      date <- Sys.Date()
      time <- format(Sys.time(), "%H%M%S")
      fout <- glue::glue("{stem}/m{m}-{date}-{time}.tsv")
      
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


# -----------------------------------------------------------------------------
# Plotting.
# -----------------------------------------------------------------------------

#' Expect Chi2 Plot
#'
#' @param df Data.frame.
#' @param y_lim Y-axis limits.
ExpectedChi2Plot <- function(df, y_lim = c(0, 1.4)) {
  q <- ggplot() + 
    theme_bw() +
    theme(legend.position = "top") +
    geom_hline(
      yintercept = 1.0,
      linetype = "dashed",
      color = "gray",
      linewidth = 1.2
    ) +
    geom_col(
      data = df,
      aes(x = miss_factor, y = mu, fill = indep_train),
      position = position_dodge()
    ) +
    geom_errorbar(
      data = df,
      aes(x = miss_factor, ymin = lower, ymax = upper, group = indep_train),
      position = position_dodge(0.9),
      width = 0.5
    ) + 
    scale_y_continuous(
      name = expression(Expected~X^2),
      limits = y_lim,
      breaks = seq(from = 0.0, to = max(y_lim), by = 0.2)
    ) +
    scale_x_discrete(
      name = "Target Missingness"
    ) + 
    ggsci::scale_fill_nejm(
      name = "Independent\nTraining Data"
    )
  return(q)
}


#' Bias plot
#'
#' @param df Data.frame.
#' @param y_lim Y-axis limits.
BiasPlot <- function(df, y_lim = c(0, 0.14)) {
  q <- ggplot() + 
    theme_bw() +
    theme(legend.position = "top") +
    geom_hline(
      yintercept = 0.1,
      linetype = "dashed",
      color = "gray",
      linewidth = 1.2
    ) +
    geom_col(
      data = df,
      aes(x = miss_factor, y = mu, fill = indep_train),
      position = position_dodge()
    ) +
    geom_errorbar(
      data = df,
      aes(x = miss_factor, ymin = lower, ymax = upper, group = indep_train),
      position = position_dodge(0.9),
      width = 0.5
    ) + 
    scale_y_continuous(
      name = expression(Estimated~beta[G]),
      limits = y_lim,
      breaks = seq(from = 0.0, to = max(y_lim), by = 0.02)
    ) +
    scale_x_discrete(
      name = "Target Missingness"
    ) + 
    ggsci::scale_fill_nejm(
      name = "Independent\nTraining Data"
    )
  return(q)
}


#' Pannel
#' 
#' @param bias_df Data.frame.
#' @param t1e_df Data.frame.
CreatePanel <- function(bias_df, t1e_df) {
  
  # Bias.
  q1 <- BiasPlot(bias_df) + ggtitle("Estimated Genetic Effect")
  
  # Type I error.
  q2 <- ExpectedChi2Plot(t1e_df) + ggtitle(expression(Expected~X^2~under~the~Null))
  
  legend <- ggpubr::get_legend(q1)
  q_legend <- ggpubr::as_ggplot(legend)
  
  q1 <- q1 + theme(legend.position = "none")
  q2 <- q2 + theme(legend.position = "none")
  
  q_panel <- cowplot::plot_grid(
    plotlist = list(q1, q_legend, q2),
    labels = NULL,
    axis = "left",
    ncol = 1,
    rel_heights = c(5, 1, 5)
  )
  return(q_panel)
}

# -----------------------------------------------------------------------------

t1e_sim$miss_factor <- factor(
  t1e_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)

bias_sim$miss_factor <- factor(
  bias_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)



if (FALSE) {
  ggsave(
    plot = q_panel,
    file = "results/data_reuse.pdf",
    width = 9.0,
    height = 6.0
  )
}


# -----------------------------------------------------------------------------
# Supplemental pannel.
# -----------------------------------------------------------------------------

bias_sim_1 <- data.table::fread("results/subject_reuse/subject_reuse_sim_bias.tsv")
bias_sim_0 <- data.table::fread("results/subject_reuse/subject_reuse_misspec_sim_bias.tsv")

t1e_sim_1 <- data.table::fread("results/subject_reuse/subject_reuse_sim_t1e.tsv")
t1e_sim_0 <- data.table::fread("results/subject_reuse/subject_reuse_sim_misspec_t1e.tsv")

# Bias.
bias_sim_0$miss_factor <- factor(
  t1e_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)

bias_sim_1$miss_factor <- factor(
  t1e_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)

# Type I error.
t1e_sim_0$miss_factor <- factor(
  t1e_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)

t1e_sim_1$miss_factor <- factor(
  t1e_sim$miss,
  levels = c(0.00, 0.25, 0.50, 0.75, 0.90),
  labels = c("0%", "25%", "50%", "75%", "90%")
)

# Overall panel.
panel_0 <- CreatePanel(bias_sim_0, t1e_sim_0)
panel_1 <- CreatePanel(bias_sim_1, t1e_sim_1)

overall <- cowplot::plot_grid(
  plotlist = list(panel_0, panel_1),
  ncol = 2,
  labels = c("A", "B")
)

if (FALSE) {
  
  ggsave(
    plot = overall,
    file = "results/data_reuse.pdf",
    width = 12.0,
    height = 6.0
  )
  
}
