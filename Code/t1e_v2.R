source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)

n.sim <- 10^4 # number of replicates
missing_rate <- c(0, 0.25, 0.5, 0.75)
rho <- c(0, 0.25, 0.5, 0.75)

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

result <- foreach(m = missing_rate) %:%
  foreach(p = rho) %:%
  foreach(j = 1:n.sim, .combine = rbind, .packages = c("doParallel"), .errorhandling = "pass") %do% {
    
pheno <- GenPheno(n = 10^3, beta = c(1,2), miss = m, rho = p)
g <- GenGeno(n = 10^3, snps = 1)

assoc.oracle <- lm(oracle ~ . -1, data = data.frame(
  cbind(pheno %>% select(oracle, starts_with("x")), g)
))
out <- summary(assoc.oracle)$coefficients["g", c(1,2,4)]

assoc.observed <- lm(target ~ . -1, data = data.frame(
  cbind(pheno %>% select(target, starts_with("x")), g)
))
out <- c(out, summary(assoc.observed)$coefficients["g", c(1,2,4)])

assoc.surrogate <- lm(surrogate ~ . -1, data = data.frame(
  cbind(pheno %>% select(surrogate, starts_with("x")), g)
))
out <- c(out, summary(assoc.surrogate)$coefficients["g", c(1,2,4)])

fit.binormal <- SurrogateRegression::Fit.BNLS(
  t = pheno$target,
  s = pheno$surrogate,
  X = cbind(pheno %>% select(starts_with("x")),g)
)
out <- c(out, fit.binormal@Regression.tab %>%
           filter(Outcome == "Target" & Coefficient == "g") %>%
           select(Point, SE, p) %>% as.numeric())

vas <- c("oracle", "target", "surrogate", "bi")
vis <- c("beta", "se", "p")
names(out) <- as.vector(t(outer(vas, vis, paste, sep=".")))
out}

saveRDS(result, file = "t1e.rds")
