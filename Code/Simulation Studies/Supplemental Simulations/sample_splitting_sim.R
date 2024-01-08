#' Purpose: Examine the trade-off between allocating subjects to the 
#' model-building versus GWAS data sets. 
#' Updated: 2023-08-23
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
  rho = 0.5
) {
  
  # Total sample size.
  n <- n0 / (1 - miss)
  
  # Genotype and covariates. G and X have correlation rho.
  g <- stats::rnorm(n)
  x <- sqrt(1 - rho^2) * stats::rnorm(n) + rho * g
  
  # Generate polynomial basis for x. 
  x_basis <- cbind(sin(x), cos(x), sin(x) * cos(x))
  x_basis <- scale(x_basis, center = FALSE)
  x_deg <- ncol(x_basis)
  
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


#' Split Data
#' 
#' @param data Data.
#' @param prop_train Proportion of subjects with observed target outcomes to use
#'   for training.
#' @return List of data.frames.
SplitData <- function(data, prop_train = 0.10) {
  
  # Split into complete and incomplete.
  complete <- data[!is.na(data$yobs), ]
  incomplete <- data[is.na(data$yobs), ]
  
  # Split complete data into model-building and gwas.
  n0 <- nrow(complete)
  n_train <- round(n0 * prop_train)

  draw <- sort(sample(n0, size = n_train))
  key <- (seq_len(n0) %in% draw)
  
  train <- complete[key, ]
  gwas <- rbind(complete[!key, ], incomplete)
  
  # Output.
  out <- list(
    train = train,
    gwas = gwas
  )
  return(out)
}


#' Fit SynSurr Model
#' 
#' @param data Data.frame.
#' @return Prediction function.
FitSynSurr <- function(data) {
  
  df <- data %>% dplyr::filter(complete.cases(.))
  x <- df %>% dplyr::pull(x)
  y <- df %>% dplyr::pull(y)
  
  # Fit.
  fit <- lm(y ~ splines::bs(x, Boundary.knots = c(-5, 5)))
  
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
#' @param train_props Training proportions.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set. 
#' @return Data.frame.
SimInstance <- function(
  train_props,  
  bg = 0.0,
  n = 1e3
) {
  
  # Overall data. 
  data <- DGP(n = n, bg = bg, miss = 0.5)
  
  # Try different sample splits.
  out <- lapply(train_props, function(m) {
    
    split_data <- SplitData(data, prop_train = m)
    data_train <- split_data$train
    data_gwas <- split_data$gwas
    
    # Fit SynSurr model. 
    model <- FitSynSurr(data_train)
    
    # Generate surrogate.
    data_gwas$yhat <- model(data_gwas)
    
    # Run SynSurr.E
    results <- SynSurrGWAS(data_gwas)
    results$prop_train <- m
    results <- results %>% dplyr::relocate(prop_train, .before = "bg")
    results$rho <- cor(data_gwas$yobs, data_gwas$yhat, use = "pairwise.complete")
    return(results)
    
  })
  out <- do.call(rbind, out)
  return(out)
}


#' Simulation
#' 
#' @param train_props Training proportions.
#' @param bg Genetic effect.
#' @param n Size of the GWAS data set.
#' @param reps Simulation replicates.
#' @return Data.frame.
Sim <- function(
    train_props,
    bg = 0.1,
    n = 1e3,
    reps = 500
) {
  
  out <- lapply(seq_len(reps), function(i) {
    
    results <- SimInstance(train_props = train_props, bg = bg, n = n)
    
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
# Run simulation.
# -----------------------------------------------------------------------------

if (FALSE) {
  
  LOOPS <- 10
  REPS <- 1e3
  stem <- "results/sample_splitting/power"
  train_props <- seq(from = 0.05, to = 0.25, by = 0.05)
  
  sink <- lapply(seq_len(LOOPS), function(i) {
    
    sim <- Sim(train_props = train_props, reps = REPS)
    
    date <- Sys.Date()
    time <- format(Sys.time(), "%H%M%S")
    fout <- glue::glue("{stem}/power-{date}-{time}.tsv")
    
    data.table::fwrite(
      x = sim,
      file = fout,
      sep = "\t"
    )
    
    return(NULL)
  })
  
}


# -----------------------------------------------------------------------------
# Power simulation.
# -----------------------------------------------------------------------------

results_file <- "results/sample_splitting/sample_splitting_power.tsv"

if (file.exists(results_file)) {
  
  power_sim <- data.table::fread(results_file)
  power_sim$prop_train <- factor(
    power_sim$prop_train,
    levels = c(0.05, 0.1, 0.15, 0.2, 0.25),
    labels = c("5%", "10%", "15%", "20%", "25%")
  )
  
} else {
  
  file_dir <- "results/sample_splitting/power"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(prop_train, bg_p) %>%
    dplyr::mutate(chi2 = qchisq(bg_p, df = 1, lower.tail = FALSE)) %>%
    dplyr::group_by(prop_train) %>%
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

# Size of GWAS data set.
df <- data.frame(
  prop_train = c(0.05, 0.1, 0.15, 0.2, 0.25),
  n = c(1950, 1900, 1850, 1800, 1750),
  n0 = c(950, 900, 850, 800, 750)
)
df$prop_train <- factor(
  df$prop_train,
  levels = c(0.05, 0.1, 0.15, 0.2, 0.25),
  labels = c("5%", "10%", "15%", "20%", "25%")
)
df <- df %>%
  tidyr::pivot_longer(
    n:n0,
    names_to = "set",
    values_to = "n"
  )
df$set <- factor(
  df$set,
  levels = c("n", "n0"),
  labels = c("Total", "Target\nObserved")
)

# Plot sample sizes.
q_n <- ggplot() + 
  theme_bw() +
  theme(
    legend.position = "top"
  ) +
  geom_col(
    data = df,
    aes(x = prop_train, y = n, fill = set),
    position = position_dodge()
  ) + 
  geom_text(
    data = df,
    aes(x = prop_train, y = n, label = n, group = set),
    position = position_dodge(width = 0.9),
    vjust = -0.2
  ) +
  ggsci::scale_fill_nejm(
    name = NULL
  ) +
  scale_x_discrete(
    name = "Proportion of Data\nAllocated to Model-Building"
  ) +
  scale_y_continuous(
    name = "Number of Subjects",
    limits = c(0, 2100)
  ) + 
  ggtitle("Size of GWAS Data Set")


# Plot expected chi2.
q_pwr <- ggplot() +
  theme_bw() +
  geom_col(
    data = power_sim,
    aes(x = prop_train, y = mu, fill = prop_train),
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = power_sim,
    aes(x = prop_train, ymin = lower, ymax = upper),
    width = 0.5
  ) + 
  scale_fill_brewer(
    palette = "Blues"
  ) +
  scale_x_discrete(
    name = "Proportion of Data\nAllocated to Model-Building"
  ) +
  scale_y_continuous(
    name = expression(Expected~chi^2),
    limits = c(0, 12)
  ) +
  ggtitle("Power")

# -----------------------------------------------------------------------------

# Bias
results_file <- "results/sample_splitting/sample_splitting_bias.tsv"

if (file.exists(results_file)) {
  
  bias_sim <- data.table::fread(results_file)
  bias_sim$prop_train <- factor(
    bias_sim$prop_train,
    levels = c(0.05, 0.1, 0.15, 0.2, 0.25),
    labels = c("5%", "10%", "15%", "20%", "25%")
  )
  
} else {
  
  file_dir <- "results/sample_splitting/power"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(prop_train, bg) %>%
    dplyr::group_by(prop_train) %>%
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
  
q_bias <- ggplot() + 
  theme_bw() +
  geom_hline(
    yintercept = 0.1,
    linetype = "dashed",
    color = "gray",
    linewidth = 1.2
  ) +
  geom_col(
    data = bias_sim,
    aes(x = prop_train, y = mu, fill = prop_train),
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = bias_sim,
    aes(x = prop_train, ymin = lower, ymax = upper),
    width = 0.5
  ) + 
  scale_x_discrete(
    name = "Proportion of Data\nAllocated to Model-Building"
  ) +
  scale_y_continuous(
    name = expression(Estimated~beta[G]),
    limits = c(0.0, 0.14),
    breaks = seq(from = 0.0, to = 0.14, by = 0.02)
  ) +
  scale_fill_brewer(
    palette = "Blues"
  ) +
  ggtitle("Unbiasedness")


# -----------------------------------------------------------------------------

# Bias
results_file <- "results/sample_splitting/target_surrogate_correaltion.tsv"

if (file.exists(results_file)) {
  
  rho_sim <- data.table::fread(results_file)
  rho_sim$prop_train <- factor(
    rho_sim$prop_train,
    levels = c(0.05, 0.1, 0.15, 0.2, 0.25),
    labels = c("5%", "10%", "15%", "20%", "25%")
  )
  
} else {
  
  file_dir <- "results/sample_splitting/power"
  unagg <- Compile(file_dir)
  agg <- unagg %>%
    dplyr::select(prop_train, rho) %>%
    dplyr::group_by(prop_train) %>%
    dplyr::summarise(
      mu = mean(rho),
      se = sqrt(var(rho) / dplyr::n()),
      lower = mu - 2 * se,
      upper = mu + 2 * se
    )
  data.table::fwrite(
    x = agg,
    file = results_file,
    sep = "\t"
  )
  
}

q_rho <- ggplot() + 
  theme_bw() +
  geom_hline(
    yintercept = 0.1,
    linetype = "dashed",
    color = "gray",
    linewidth = 1.2
  ) +
  geom_col(
    data = rho_sim,
    aes(x = prop_train, y = mu, fill = prop_train),
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = rho_sim,
    aes(x = prop_train, ymin = lower, ymax = upper),
    width = 0.5
  ) + 
  scale_x_discrete(
    name = "Proportion of Data\nAllocated to Model-Building"
  ) +
  scale_y_continuous(
    name = "Target-Surrogate Correlation"
  ) +
  scale_fill_brewer(
    palette = "Blues"
  ) +
  ggtitle("Correlation")


# -----------------------------------------------------------------------------

q_panel <- cowplot::plot_grid(
  plotlist = list(q_n, q_bias, q_rho, q_pwr),
  labels = c("A", "B", "C", "D")
)

ggsave(
  plot = q_panel,
  file = "results/sample_splitting.pdf",
  width = 12.0,
  height = 8.0
)
