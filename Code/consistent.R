source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)

n.sim <- 10^4 # number of replicates
missing_rate <- c(0, 0.25, 0.5, 0.75)
rho <- c(0, 0.25, 0.5, 0.75)

cl <- makeCluster(type = "MPI")
# cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)

result <- foreach(
  j = 1:n.sim,
  .packages = c("doParallel", "dplyr", "RNOmni"),
  .errorhandling = "pass"
) %dopar% {

  # paramter set-up
  n0 <- 10^3 # number of complete cases
  beta <- c(1, 0.42, 0.11) # beta0, beta1, beta2
  X <- cbind(rnorm(n0), rbinom(n0, 1, 0.5))
  g <- rbinom(n = n0, size = 2, prob = 0.25)
  para <- expand.grid(missing_rate, rho)
  beta_g = 0.11576

  # generate phenotype, INT already applied
  pheno <- apply(para, 1, function(m) {
    GenPheno(
      n = n0, beta = beta, x = X, miss = m[1], rho = m[2],
      include_intercept = TRUE, beta_g = beta_g, g = g, INT = TRUE
    )
  })

  # target
  assoc.target <- lapply(pheno, function(x) lm(x$target ~ g + X))
  out <- lapply(assoc.target, function(x) summary(x)$coefficients["g", c(1, 2, 4)])
  out <- matrix(unlist(out), ncol = 3, byrow = TRUE) %>% data.frame()

  # Bivariate
  fit.binormal <- lapply(pheno, function(x) {
    SurrogateRegression::Fit.BNLS(
      t = x$target,
      s = x$surrogate,
      X = cbind(x %>% dplyr::select(starts_with("x")), g)
    )
  })
  out.bi <- lapply(fit.binormal, function(x) {
    x@Regression.tab %>%
      filter(Outcome == "Target" & Coefficient == "g") %>%
      select(Point, SE, p) %>%
      as.numeric()
  })
  out.bi <- matrix(unlist(out.bi), ncol = 3, byrow = TRUE) %>% data.frame()
  out <- cbind(para, out, out.bi)
  vas <- c("target", "bi")
  vis <- c("beta", "se", "p")
  colnames(out) <- c("missing", "rho", as.vector(t(outer(vas, vis, paste, sep = "."))))
  out
}

stopCluster(cl)

saveRDS(result, file = "robustness.rds")