source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)

n.sim <- 10^4 # number of replicates
rho <- c(0, 0.25, 0.5, 0.75)

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

result <- foreach(p = rho) %:%
  foreach(j = 1:n.sim, .combine = rbind, .packages = c("doParallel"), .errorhandling = "pass") %do% {
    
    beta_g = 0.1
    n0 = 10^3 # complete case
    g <- GenGeno(n = n0, snps = 1)
    
    # oracle
    target <- rnorm(g * beta_g, 1)
    assoc.oracle <- lm(target ~ g -1)
    out <- summary(assoc.oracle)$coefficients["g", c(1,2,4)]
    
    # correct model + known
    yhat <- rnorm(g * beta_g + rho * (target - g * beta_g), 1-rho^2)
    assoc.surrogate <- lm(yhat ~ g -1)
    out <- c(out, summary(assoc.surrogate)$coefficients["g", c(1,2,4)])
    
    # correct model + unknown
    index <- sample(1:n0, 500)
    beta_hat <- summary(lm(target[index]~g[index]-1))[["coefficients"]][,1]
    yhat <- rnorm(g * beta_hat + rho * (target - g * beta_hat), 1-rho^2)
    assoc.surrogate <- lm(yhat ~ g -1)
    out <- c(out, summary(assoc.surrogate)$coefficients["g", c(1,2,4)])
    # incorrect model + known
    
    
    
    vas <- c("oracle", "target", "surrogate", "bi")
    vis <- c("beta", "se", "p")
    names(out) <- as.vector(t(outer(vas, vis, paste, sep=".")))
    out}

saveRDS(result, file = "t1e.rds")