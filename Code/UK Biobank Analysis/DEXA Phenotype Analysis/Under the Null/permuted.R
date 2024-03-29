library(data.table)
library(dplyr)
library(bigsnpr)
library(doParallel)
library(BEDMatrix)

library(ranger)
library(readr)
library(dplyr)
library(data.table)
library(tidyr)

#setwd("~/Desktop/SyntheticSurrogateAnalysis/Data")

cov_id <- c("f.21022.0.0","f.22001.0.0","f.50.0.0","f.23098.0.0","f.23104.0.0") #age,sex,height, weight, BMI
pheno_id <- c("f.23283.2.0")

# split data into train and test
file <- read.table("CombinedTotalMass.tab", header=TRUE, sep="\t")
in_id <- fread("final.king.cutoff.in.id") %>% rename(f.eid = `#FID`) %>% select(-IID)
out_id <- fread("final.king.cutoff.out.id") %>% rename(f.eid = `#FID`) %>% select(-IID)

cov_column <- file %>% select(cov_id)
train <- file %>% inner_join(out_id) %>% tidyr::drop_na(pheno_id) %>% tidyr::drop_na(all_of(cov_id)) %>% select(all_of(pheno_id), all_of(cov_id))
test <- file %>% inner_join(in_id) %>% tidyr::drop_na(all_of(cov_id))%>% select(f.eid, all_of(pheno_id), all_of(cov_id))

# Random Forest model
model.rf <- ranger::ranger(
  data = train,
  as.formula(paste0(pheno_id, "~."))
)

yhat <- predict(model.rf, data = test %>% select(-f.eid,-pheno_id))$predictions
test$yhat <- yhat

# Merge genetic PC
pcs <- fread("UKB_PC.tab")[,1:11]
colnames(pcs) <- c("f.eid", paste0("PC",1:10))
test <- test %>% inner_join(pcs, by = "f.eid")

# Inverse Normal Transformation
INT <- function(data, pheno, k = 0.375){
  missing_index <- which(is.na(data[[pheno]]))
  data.complete <- data[-missing_index, ]
  data.missing <-data[missing_index, ]
  
  n <- nrow(data.complete)
  r <- rank(data.complete %>% pull(pheno))
  data.complete$int <- qnorm((r - k) / (n - 2 * k + 1))
  data.missing$int <- NA
  
  return(rbind(data.complete, data.missing))
}
test <- INT(test, pheno_id) # INT of phenotype
n <- nrow(test)
r <- rank(test %>% pull(yhat))
test$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat

# Randomly permute
index <- sample(nrow(test))
int_permuted <- test$int[index]
yhat_int_permuted <- test$yhat_int[index]
test$int <- int_permuted
test$yhat_int <- yhat_int_permuted

# Genetic Association Test
G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

clusterExport(cl, list("test","cov_id"))
clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
    g <- as.numeric(G[as.character(test$f.eid), i]) # snp i
    # observed phenotype ~ intercept + age + sex + 10 genetic PC + SNP_i
    assoc.observed <- lm(int ~ ., data = data.frame(
      cbind(test %>% select(int, "f.21022.0.0","f.22001.0.0", starts_with("PC")), g)
    ))
    out <- summary(assoc.observed)$coefficients["g", c(1,2,4)] # only beta, se, p-val
  
    # Bivariate model
    g_complete <- g[!is.na(g)]
    X.cov <- cbind(g_complete, (test %>% select("f.21022.0.0","f.22001.0.0",starts_with("PC")))[!is.na(g), ])
    X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept
  
    fit.binormal <- SurrogateRegression::Fit.BNLS(
      t = test$int[!is.na(g)],
      s = test$yhat_int[!is.na(g)],
      X = X.cov
    )
    out <- c(out, fit.binormal@Regression.tab %>%
               filter(Outcome == "Target" & Coefficient == "g_complete") %>%
               select(Point, SE, p) %>% as.numeric())
  
    vas <- c("marginal","bi")
    vis <- c("beta", "se", "p")
    out <- c(colnames(G)[i], out)
    names(out) <- c("rsid", as.vector(t(outer(vas, vis, paste, sep = "."))))
    return(out)
})

stopCluster(cl)

results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, file = paste0("pheno=", pheno_id, "_result_permuted.rds"))