library(data.table)
library(dplyr)
library(bigsnpr)
library(doParallel)
library(BEDMatrix)

create_missing_phenotypes <- function(path, rate, seed = 123) {
  pdata <- fread(path) # phenotype
  if (sum(is.na(pdata)) > 0) {
    stop("cannot contain missing data")
  }
  nmissing <- floor(nrow(pdata) * rate) # number of missing phenotypes
  set.seed(123)
  missing_index <- sample(1:nrow(pdata), size = nmissing)
  pdata[missing_index, ]$V3 <- NA # set phenotype missing value
  pdata <- pdata[, 2:3] # only keep id and phenotype
  colnames(pdata) <- c("id", "target")
  return(data.frame(pdata))
}

merge_data <- function(pdata, sdata, cdata) {
  combined <- cbind(pdata, cdata) # add covariate
  combined <- combined %>% inner_join(sdata, by = "id") # add surrogate data
  return(combined)
}


args=(commandArgs(TRUE))
print(args)
for(k in 1:length(args)){
  eval(parse(text=args[[k]]))
}

load("covariate.RData") # load covariate data
load("ypred.RData") # load rf model
pdata <- create_missing_phenotypes(path = "phenotype.txt", rate = missing_rate)
pheno <- merge_data(pdata = pdata, sdata = y_pred, cdata = covar)
G <- BEDMatrix::BEDMatrix(path = "final.bed", simple_names = T) # read genetic data


cl <- makeCluster(type="MPI")
registerDoParallel(cl)

clusterExport(cl, list("pheno"))
clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "final.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i){
  g <- as.numeric(G[which(rownames(G) %in% pheno$id), i]) # snp i
  # observed phenotype ~ intercept + age + sex + 10 genetic PC + SNP_i
  assoc.observed <- lm(target ~ . - 1, data = data.frame(
    cbind(pheno %>% select(c("target", starts_with("x"))), g)
  ))
  out <- summary(assoc.observed)$coefficients["g", 1:2] # only beta and se to save memory
  
  # surrogate ~ intercept + age + sex + 10 genetic PC + SNP_i
  assoc.imputed <- lm(surrogate ~ . - 1, data = data.frame(
    cbind(pheno %>% select(c("surrogate", starts_with("x"))), g)
  ))
  out <- c(out, summary(assoc.imputed)$coefficients["g", 1:2])
  
  # Bivariate model
  g_complete <- g[!is.na(g)]
  fit.binormal <- SurrogateRegression::Fit.BNLS(
    t = pheno$target[!is.na(g)],
    s = pheno$surrogate[!is.na(g)],
    X = cbind(g_complete, (pheno %>% select(starts_with("x")))[!is.na(g), ])
  )
  out <- c(out, fit.binormal@Regression.tab %>%
             filter(Outcome == "Target" & Coefficient == "g_complete") %>%
             select(Point, SE) %>% as.numeric())
  names(out) <- c("beta_g_obs", "se_obs", "beta_g_surrogate", "se_surrogate", "beta_g_bi", "se_bi")
  return(out)
})

stopCluster(cl)

results <- do.call(rbind, results)
results <- data.frame(results)

save(results, file = paste0("binormal_results",missing_rate,".RData"))