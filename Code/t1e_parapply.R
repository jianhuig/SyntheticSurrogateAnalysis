library(dplyr)
library(doParallel)
source("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Code/SynSurrogateSim.R")
n.sim <- 10^5 # number of replicates


cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)

pheno <- lapply(c(0, 0.25, 0.5, 0.75), function(m) {
  inner <- lapply(c(0, 0.25, 0.5, 0.75), function(p) {
    return(GenPheno(n = 10^3, beta = c(1,1,2), miss = m, rho = p, include_intercept = TRUE))})
  return(inner)}
)
G <- GenGeno(n = 10^3, snps = n.sim)

clusterExport(cl, list("GenGeno", "GenPheno","pheno","G"))
clusterEvalQ(cl, {
  library(dplyr)
})


result <- parLapply(cl, X = 1:n.sim, fun = function(i) {
  outer <- lapply(1:4, function(m) {
    inner <- lapply(1:4, function(p) {
      
      pheno <- pheno[[m]][[p]] %>% data.frame()
      g <- G[,i]

      assoc.oracle <- lm(oracle ~ ., data = data.frame(
        cbind(pheno %>% select(oracle, starts_with("x")), g)
      ))
      out <- summary(assoc.oracle)$coefficients["g", c(1, 2, 4)]

      assoc.observed <- lm(target ~ ., data = data.frame(
        cbind(pheno %>% select(target, starts_with("x")), g)
      ))
      out <- c(out, summary(assoc.observed)$coefficients["g", c(1, 2, 4)])

      assoc.surrogate <- lm(surrogate ~ ., data = data.frame(
        cbind(pheno %>% select(surrogate, starts_with("x")), g)
      ))
      out <- c(out, summary(assoc.surrogate)$coefficients["g", c(1, 2, 4)])

      fit.binormal <- SurrogateRegression::Fit.BNLS(
        t = pheno$target,
        s = pheno$surrogate,
        X = cbind(pheno %>% select(starts_with("x")), g)
      )
      out <- c(out, fit.binormal@Regression.tab %>%
        filter(Outcome == "Target" & Coefficient == "g") %>%
        select(Point, SE, p) %>% as.numeric())

      vas <- c("oracle", "target", "surrogate", "bi")
      vis <- c("beta", "se", "p")
      names(out) <- as.vector(t(outer(vas, vis, paste, sep = ".")))
      return(out)
    })
    return(inner)
  })
  return(outer)
})

stopCluster(cl)

saveRDS(result, file = paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data/", "t1e,nsim=",n.sim,".rds"))

t1e <- readRDS(paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data/", "t1e,nsim=",n.sim,".rds"))

t1e <-result
prop <- c()
for(m in 1:4){
  for(j in 1:4){
    rate <- c(0, 0.25, 0.5, 0.75)[m]
    p <- c(0, 0.25, 0.5, 0.75)[j]
    n <- ncol(G)
    
    prop.oracle <- sum(unlist(lapply(1:n, function(i){
      return(t1e[[i]][[m]][[j]][["oracle.p"]])}))<0.05)/n
    prop.target <- sum(unlist(lapply(1:n, function(i){
      return(t1e[[i]][[m]][[j]][["target.p"]])}))<0.05)/n
    prop.surrogate <- sum(unlist(lapply(1:n, function(i){
      return(t1e[[i]][[m]][[j]][["surrogate.p"]])}))<0.05)/n
    prop.bi <- sum(unlist(lapply(1:n, function(i){
      return(t1e[[i]][[m]][[j]][["bi.p"]])}))<0.05)/n
    
    prop <- rbind(prop, c(rate, p, prop.oracle, prop.target, prop.surrogate,prop.bi))
}}
prop <- prop %>% data.frame()
colnames(prop) = c("mssing","rho","oracle","target","surrogate","bivariate")
prop



#### Version 2
m = 0
p = 0
set.seed(123)

G <- GenGeno(n = 10^3, snps = 10^4)

cl <- makeCluster(parallel::detectCores()-1)
registerDoParallel(cl)

clusterExport(cl, list("GenGeno", "GenPheno", "m", "p","pheno","G"))
clusterEvalQ(cl, {
  library(dplyr)
})


result <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
  outer <- lapply(c(0, 0.25, 0.5, 0.75), function(m) {
    inner <- lapply(c(0, 0.25, 0.5, 0.75), function(p) {
      pheno <- GenPheno(n = 10^3, beta = c(1,1,2), miss = m, rho = p, include_intercept = TRUE)
      g = G[,i]
      
      assoc.oracle <- lm(oracle ~ ., data = data.frame(
        cbind(pheno %>% select(oracle, starts_with("x")), g)
      ))
      out <- summary(assoc.oracle)$coefficients["g", c(1, 2, 4)]
      
      assoc.observed <- lm(target ~ ., data = data.frame(
        cbind(pheno %>% select(target, starts_with("x")), g)
      ))
      out <- c(out, summary(assoc.observed)$coefficients["g", c(1, 2, 4)])
      
      assoc.surrogate <- lm(surrogate ~ ., data = data.frame(
        cbind(pheno %>% select(surrogate, starts_with("x")), g)
      ))
      out <- c(out, summary(assoc.surrogate)$coefficients["g", c(1, 2, 4)])
      
      fit.binormal <- SurrogateRegression::Fit.BNLS(
        t = pheno$target,
        s = pheno$surrogate,
        X = cbind(1, cbind(pheno %>% select(starts_with("x")), g))
      )
      out <- c(out, fit.binormal@Regression.tab %>%
                 filter(Outcome == "Target" & Coefficient == "g") %>%
                 select(Point, SE, p) %>% as.numeric())
      
      vas <- c("oracle", "target", "surrogate", "bi")
      vis <- c("beta", "se", "p")
      names(out) <- as.vector(t(outer(vas, vis, paste, sep = ".")))
      return(out)
    })
    return(inner)
  })
  return(outer)
})