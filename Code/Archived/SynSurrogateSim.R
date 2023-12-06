# Purpose: Generate synthetic data.
# Updated: 2022-02-16


#' Generate Data
#'
#' Generates a data.frame containing the target and surrogate outcomes
#' accompanied by standard normal covariates. The oracle outcome is
#' identical to the target outcome except for containing no missing values.
#' Beta is the target regression coefficient, while alpha is the surrogate
#' regression coefficient. If omitted, alpha defaults to beta.
#'
#' @param n Sample size.
#' @param beta Target regression coefficients.
#' @param alpha Surrogate regression coefficients.
#' @param include_intercept Should an intercept be included?
#' @param miss Target missingness.
#' @param rho Target-surrogate correlation.
#' @return Data.frame.
#' @export

GenPheno <- function(n,
                     beta,
                     x = NULL,
                     alpha = NULL,
                     include_intercept = FALSE,
                     miss = 0.0,
                     rho = 0.0,
                     beta_g,
                     g,
                     INT = TRUE) {
  if (is.null(alpha)) {
    alpha <- beta
  }

  # Covariates.
  if (include_intercept) {
    n_covar <- length(beta) - 1
  } else {
    n_covar <- length(beta)
  }

  # If x is not provided, generate x from Normal
  if (is.null(x)) {
    x <- rbind(replicate(n_covar, rnorm(n)))
    x <- scale(x)
  }

  if (include_intercept) {
    x <- cbind(1, x)
    colnames(x) <- paste0("x", seq(from = 0, to = n_covar))
  } else {
    colnames(x) <- paste0("x", seq(from = 1, to = n_covar))
  }

  # Means.
  mu_target <- x %*% beta + beta_g * g
  mu_surrogate <- x %*% alpha

  # Outcomes.
  target <- mu_target + rnorm(n)
  surrogate <- mu_surrogate + rho * (target - mu_target) + rnorm(n, sd = sqrt(1 - rho^2))

  # Missingness.
  oracle <- target
  target[sample.int(n, size = round(miss * n))] <- NA

  # INT
  if (INT) {
    oracle <- RNOmni::RankNorm(as.numeric(oracle))
    surrogate <- RNOmni::RankNorm(as.numeric(surrogate))
    target[!is.na(target)] <- RNOmni::RankNorm(target[!is.na(target)])
  }

  # Output.
  data <- data.frame(target, oracle, surrogate, x)
  return(data)
}


#' Generate Genotypes
#'
#' Generates a subject by SNP genotype matrix.
#'
#' @param n Sample size.
#' @param snps Number of SNPs.
#' @param m0 Lower bound on minor allele frequency.
#' @param m1 Upper bound on minor allele frequency.
#' @return Data.frame.
#' @export

GenGeno <- function(n, snps, m0 = 0.05, m1 = 0.50) {
  maf <- runif(n = snps, min = m0, max = m1)
  out <- lapply(maf, function(m) {
    return(rbinom(n = n, size = 2, prob = m))
  })
  out <- do.call(cbind, out)
  storage.mode(out) <- "numeric"
  return(out)
}

SampleSizes <- function(n0, mt = 0, ms = 0) {
  
  # Input check.
  m <- mt + ms
  if (m < 0 || m >= 1) {
    stop("mt+ms must belong to the interval [0,1).")
  }
  
  # Sample size.
  out <- list()
  
  # Surrogate only.
  out$n1 <- ceiling(n0 * mt / (1 - mt - ms))
  
  # Target only.
  out$n2 <- ceiling(n0 * ms / (1 - mt - ms))
  
  # Overall sample size.
  out$n <- (n0 + out$n1 + out$n2)
  
  # output
  return(out)
}

CalcEffectSizes <- function(
  n_covar, 
  pve_covar = 0.2, 
  pve_geno, 
  n_pcs, 
  pve_pcs = 0.05, 
  resid_var = 1.00
) {
  
  # Proportion variation explained by residual.
  pve_resid <- 1 - pve_covar - pve_geno - pve_pcs
  
  if (pve_resid < 0) {
    stop("Total variation explained by covars, geno, and PCs cannot exceed 1.")
  }
  
  # Coefficients. 
  bg <- sqrt(resid_var * pve_geno / pve_resid)
  bx <- sqrt(resid_var * pve_covar / (n_covar * pve_resid))
  bs <- sqrt(resid_var * pve_pcs / (n_pcs * pve_resid))
  
  # Output.
  out <- list(
    "bg" = bg,
    "bx" = rep(bx, times = n_covar),
    "bs" = rep(bs, times = n_pcs)
  )
  return(out)
}