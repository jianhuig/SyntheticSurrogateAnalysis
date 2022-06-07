source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)

n.sim <- 10^4 # number of replicates
missing_rate <- c(0, 0.25, 0.5, 0.75)
rho <- c(0, 0.25, 0.5, 0.75)

# calculate beta_g
cal_betag <- function(h2, maf=0.25){
  var = 2*maf*(1-maf)
  return(sqrt(h2/(var*(1-h2))))
}


cl <- makeCluster(type = "MPI")
#cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)

result <- foreach(
  j = 1:n.sim,
  .packages = c("doParallel", "dplyr", "RNOmni"),
  .errorhandling = "pass"
) %dopar% {
  
  # paramter set-up
  n0 <- 10^3 # number of complete cases
  beta <- c(0.11) # beta0, beta1, beta2
  beta_g = sapply(seq(0.001,0.01, 0.001), function(x) cal_betag(h2 =x))
  para <- expand.grid(missing_rate, rho, beta_g) %>% data.frame() %>%
    mutate(n = SampleSizes(n0 = n0, mt = Var1)$n)
  X <- apply(para, 1, function (x) rnorm(x[4]))
  g <- apply(para, 1, function (x) rbinom(n = x[4], size = 2, prob = 0.25))
  
  
  # generate phenotype, INT already applied
  pheno <- lapply(1:nrow(para), function(i) {
    GenPheno(
      n = para[i,]$n, beta = beta, miss = para[i, 1], rho = para[i, 2],
      include_intercept = FALSE, beta_g = para[i, 3], 
      g = g[[i]],
      x = cbind(X[[i]]),
      INT = TRUE
    )
  })
  
  
  # target
  assoc.target <- lapply(1:nrow(para), function(i) lm(pheno[[i]]$target ~ g[[i]] + X[[i]] - 1))
  out <- lapply(assoc.target, function(x) summary(x)$coefficients["g[[i]]", c(1, 2, 4)])
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
  colnames(out) <- c("missing", "rho","beta_g", "n.total", 
                     as.vector(t(outer(vas, vis, paste, sep = "."))))
  out
}

stopCluster(cl)

#saveRDS(result, file = "../Data/power.rds")
saveRDS(result, file = "power.rds")