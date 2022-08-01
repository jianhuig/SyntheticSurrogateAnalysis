source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)


#' Simulation
#'
#' Case of correctly specified imputation model.
#'
#' @param n Sample size for the number of complete cases.
#' @param reps Simulation replicates.
#' @param missing_rate Vector of missingness rates to consider.
#' @param rho Vector of synthetic surrogate-target correlations.
#' @return Data.frame.
Sim <- function(n = 1e3, reps = 1e3, missing_rate = c(0, 0.25, 0.5, 0.75),
                rho = c(0, 0.25, 0.5, 0.75)) {
  
  cl <- makeCluster(type = "MPI")
  registerDoParallel(cl)
  
  result <- foreach(m = missing_rate) %:%
    foreach(p = rho) %:%
    foreach(j = 1:reps, .combine = rbind, 
            .packages = c("doParallel", "dplyr", "RNOmni"),
            .errorhandling = "pass") %do% {
              n <- 10^3 # number of complete cases
              beta <- c(1, 0.42, 0.11) # beta0, beta1, beta2
              beta_g <- 0
              X <- cbind(rnorm(n), rbinom(n, 1, 0.5))
              
              pheno <- GenPheno(n = n, beta = beta, x = X, miss = m, 
                                rho = p, include_intercept = TRUE)
              g <- GenGeno(n = n, snps = 1)
              
              oracle <- pheno$oracle + beta_g * g
              target <- pheno$target + beta_g * g
              
              oracle <- RNOmni::RankNorm(as.numeric(oracle))
              target[!is.na(target)] <- RNOmni::RankNorm(target[!is.na(target)])
              
              assoc.oracle <- lm(oracle ~ g + pheno$x1 + pheno$x2)
              
              out <- summary(assoc.oracle)$coefficients["g", c(1, 2, 4)]
              
              assoc.target <- lm(target ~ g + pheno$x1 + pheno$x2)
              
              out <- c(out, summary(assoc.target)$coefficients["g", c(1, 2, 4)])
              
              fit.binormal <- SurrogateRegression::Fit.BNLS(
                t = target,
                s = RNOmni::RankNorm(pheno$surrogate),
                X = cbind(pheno %>% dplyr::select(starts_with("x")), g)
              )
              out <- c(out, fit.binormal@Regression.tab %>%
                         filter(Outcome == "Target" & Coefficient == "g") %>%
                         select(Point, SE, p) %>% as.numeric())
              
              vas <- c("oracle", "target", "bi")
              vis <- c("beta", "se", "p")
              names(out) <- as.vector(t(outer(vas, vis, paste, sep = ".")))
              out
            }
  
  stopCluster(cl)
  
  return(result)
}


# Run the simulation.
sim_result <- Sim(n = 1e3, reps = 1e3, missing_rate = c(0, 0.25, 0.5, 0.75),
                  rho = c(0, 0.25, 0.5, 0.75))

saveRDS(sim_result, file = "t1e.rds")
