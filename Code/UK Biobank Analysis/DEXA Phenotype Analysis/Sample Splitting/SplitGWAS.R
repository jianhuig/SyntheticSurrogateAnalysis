library(data.table)
library(dplyr)
library(bigsnpr)
library(doParallel)
library(BEDMatrix)
library(ranger)
library(dplyr)
library(data.table)
library(tidyr)

# setwd("~/Desktop/SyntheticSurrogateAnalysis/Data")

cov_id <- c("f.21022.0.0","f.22001.0.0","f.50.0.0","f.23098.0.0","f.23104.0.0",
             "f.23106.0.0", "f.23107.0.0", "f.23108.0.0", "f.23109.0.0",
             "f.23110.0.0") #age,sex,height, weight, BMI, and impedance measures
# pheno_id <- c("f.23248.2.0")
pheno_id <- (commandArgs(TRUE))
print(pheno_id)

# split data into train and test
file <- read.table("DEXA.tab", header=TRUE, sep="\t")
test_id <- readRDS("test_id.rds")
train_id <- readRDS("train_id.rds")

train <- file %>% 
  filter(!is.na(f.22006.0.0)) %>% # Filter UKB caucasian
  filter(f.eid %in% train_id) %>% 
  tidyr::drop_na(pheno_id) %>% 
  tidyr::drop_na(all_of(cov_id)) %>% 
  select(all_of(c(pheno_id,cov_id)))
test <- file %>% 
filter(!is.na(f.22006.0.0)) %>% # Filter UKB caucasian 
filter(f.eid %in% test_id) %>% 
tidyr::drop_na(all_of(cov_id))%>% 
select(all_of(c("f.eid", pheno_id, cov_id)))

# Random Forest model
model.rf <- ranger::ranger(
  data = train,
  as.formula(paste0(pheno_id, "~."))
)

yhat <- predict(model.rf, data = test %>% select(-f.eid, -pheno_id))$predictions
test$yhat <- yhat

# Merge genetic PC
pcs <- fread("UKB_PC.tab")[, 1:11]
colnames(pcs) <- c("f.eid", paste0("PC", 1:10))
test <- test %>% inner_join(pcs, by = "f.eid")

# Inverse Normal Transformation
INT <- function(data, pheno, k = 0.375) {
  missing_index <- which(is.na(data[[pheno]]))
  data.complete <- data[-missing_index, ]
  data.missing <- data[missing_index, ]

  n <- nrow(data.complete)
  r <- rank(data.complete %>% pull(pheno))
  data.complete$int <- qnorm((r - k) / (n - 2 * k + 1))
  data.missing$int <- NA

  return(rbind(data.complete, data.missing))
}
test <- INT(test, pheno_id) # INT of phenotype
n <- nrow(test)
r <- rank(test %>% pull(yhat))
test$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) # INT of yhat

# Indicator of int missingness
test$int_missing <- ifelse(is.na(test$int), 1, 0)

# Split GWAS
# Stratified by target phenotype
set.seed(123)
discovery <- test %>%
  group_by(int_missing) %>%
  sample_frac(size = .8) %>%
  ungroup()
validation <- test %>% anti_join(discovery)

# Genetic Association Test
G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

clusterExport(cl, list("discovery", "validation", "cov_id"))
clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
  # Standard GWAS on discovery
  g <- as.numeric(G[as.character(discovery$f.eid), i]) # snp i in discovery
  standard_discovery <- lm(int ~ ., data = data.frame(
    cbind(discovery %>% select(int, "f.21022.0.0", "f.22001.0.0", starts_with("PC")), g)
  ))
  out <- summary(standard_discovery)$coefficients["g", c(1, 2, 4)]

  # Standard GWAS on validation
  g <- as.numeric(G[as.character(validation$f.eid), i]) # snp i
  standard_validation <- lm(int ~ ., data = data.frame(
    cbind(validation %>% select(int, "f.21022.0.0", "f.22001.0.0", starts_with("PC")), g)
  ))
  out <- c(out, summary(standard_validation)$coefficients["g", c(1, 2, 4)])

  # SynSurr on discovery
  g <- as.numeric(G[as.character(discovery$f.eid), i]) # snp i
  g_complete <- g[!is.na(g)]
  X.cov <- cbind(g_complete, (discovery %>% select("f.21022.0.0", "f.22001.0.0", starts_with("PC")))[!is.na(g), ])
  X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept

  fit.binormal.discovery <- SurrogateRegression::FitBNR(
    t = discovery$int[!is.na(g)],
    s = discovery$yhat_int[!is.na(g)],
    X = X.cov
  )
  out <- c(out, fit.binormal.discovery@Regression.tab %>%
    filter(Outcome == "Target" & Coefficient == "g_complete") %>%
    select(Point, SE, p) %>% as.numeric())

  # SynSurr on validation
  g <- as.numeric(G[as.character(validation$f.eid), i]) # snp i
  g_complete <- g[!is.na(g)]
  X.cov <- cbind(g_complete, (validation %>% select("f.21022.0.0", "f.22001.0.0", starts_with("PC")))[!is.na(g), ])
  X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept

  fit.binormal.validation <- SurrogateRegression::FitBNR(
    t = validation$int[!is.na(g)],
    s = validation$yhat_int[!is.na(g)],
    X = X.cov
  )
  out <- c(out, fit.binormal.validation@Regression.tab %>%
    filter(Outcome == "Target" & Coefficient == "g_complete") %>%
    select(Point, SE, p) %>% as.numeric())
  out <- c(colnames(G)[i], out)

  names <- c("rsid",
    "standard_discovery_beta", "standard_discovery_se", "standard_discovery_p",
    "standard_validation_beta", "standard_validation_se", "standard_validation_p",
    "synsurr_discovery_beta", "synsurr_discovery_se", "synsurr_discovery_p",
    "synsurr_validation_beta", "synsurr_validation_se", "synsurr_validation_p"
  )
  names(out) <- names
  return(out)
})

stopCluster(cl)

results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, file = paste0("pheno=", pheno_id, "_split_result.rds"))
