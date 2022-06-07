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
  beta <- c(0.11) # beta0, beta1, beta2
  para <- expand.grid(missing_rate, rho) %>% data.frame() %>%
    mutate(n = SampleSizes(n0 = n0, mt = Var1)$n)
  beta_g = 0.11576
  X <- apply(para, 1, function (x) rnorm(x[3]))
  g <- apply(para, 1, function (x) rbinom(n = x[3], size = 2, prob = 0.25))

  # generate phenotype, INT already applied
  pheno <- lapply(1:nrow(para), function(i) {
    GenPheno(
      n = para[i,3], beta = beta, miss = para[i,1], rho = para[i,2],
      include_intercept = FALSE, beta_g = beta_g, 
      g = g[[i]],
      x = cbind(X[[i]]),
      INT = TRUE
    )
  })

  # target
  assoc.target <- lapply(1:nrow(para), function(i) lm(pheno[[i]]$target ~ g[[i]] + X[[i]]))
  out <- lapply(assoc.target, function(x) summary(x)$coefficients[2, c(1, 2, 4)])
  out <- matrix(unlist(out), ncol = 3, byrow = TRUE) %>% data.frame()

  # Bivariate
  fit.binormal <- lapply(1:nrow(para), function(i) {
    SurrogateRegression::Fit.BNLS(
      t = pheno[[i]]$target,
      s = pheno[[i]]$surrogate,
      X = cbind(pheno[[i]] %>% dplyr::select(starts_with("x")), g[[i]])
    )
  })
  out.bi <- lapply(fit.binormal, function(x) {
    x@Regression.tab %>%
      filter(Outcome == "Target" & Coefficient == "g[[i]]") %>%
      select(Point, SE, p) %>%
      as.numeric()
  })
  out.bi <- matrix(unlist(out.bi), ncol = 3, byrow = TRUE) %>% data.frame()
  out <- cbind(para, out, out.bi)
  vas <- c("target", "bi")
  vis <- c("beta", "se", "p")
  colnames(out) <- c("missing", "rho", 'n.total',as.vector(t(outer(vas, vis, paste, sep = "."))))
  out
}

stopCluster(cl)

saveRDS(result, file = "robustness.rds")