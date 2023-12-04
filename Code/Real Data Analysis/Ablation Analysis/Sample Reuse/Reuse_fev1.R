library(data.table)
library(dplyr)
library(doParallel)
library(randomForest)

create_missing_phenotypes <- function(field, in_id_file, rate, seed = 123) {
  temp <- readRDS(paste0("field_", field, "_cleaned.rds")) # cleaned pheno + cov
  in_id <- fread(in_id_file) %>%
    rename(f.eid = `#FID`) %>%
    select(-IID)
  temp <- temp %>%
    inner_join(in_id) %>%
    filter(!is.na(int))

  # number of missing phenotypes
  nmissing <- floor(nrow(temp) * rate)

  set.seed(seed)
  missing_index <- sample(which(!is.na(temp$int)), size = nmissing)

  # oracle
  temp$oracle <- temp$int

  # set phenotype missing value
  temp[missing_index, ]$int <- NA
  return(data.frame(temp %>% select(-paste0("f.", field, ".0.0"))))
}

# FEV1 Analysis ============================================================
setwd("Data/Old/")
missing_rate <- as.numeric(commandArgs(trailingOnly=TRUE))
field <- 20150

# White British
ancestry <- fread("Ancestry.tab") %>% filter(!is.na(f.22006.0.0))

# train data
train <- create_missing_phenotypes(
  field = field,
  in_id_file = "final.king.cutoff.in.id", rate = missing_rate
)

# merge with White British
train <- train %>% inner_join(ancestry, by = "f.eid")

# Get original FEV1 data
fev1 <- data.table::fread("/scratch/l/leisun/jianhuig/UKBB/Data/Old/FEV1.tab") %>% select(f.eid, f.20150.0.0)

# Joint with training data
train <- train %>% inner_join(fev1, by = "f.eid")

# remove missing covaraites and phenotypes for training
train <- train %>%
  tidyr::drop_na(all_of(grep("0.0", colnames(train))))

# Use observed data to train
train_sub <- train %>% select(grep("0.0", colnames(train))) %>% select(-f.22006.0.0) 

# Random Forest Model
rf <- ranger::ranger(
    data = train_sub,
    f.20150.0.0 ~ .
  )
# Prediction
train$yhat <- predict(rf, data = train %>% select(grep("0.0", colnames(train))))$predictions

# Inverse normal transformation
n <- nrow(train)
r <- rank(train %>% pull(yhat))
train$yhat <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

G <- BEDMatrix::BEDMatrix(path = "/scratch/l/leisun/jianhuig/UKBB/Data/allchromosome.bed", simple_names = T) # read genetic data

clusterExport(cl, list("train"))

clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "/scratch/l/leisun/jianhuig/UKBB/Data/allchromosome.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
    g <- as.numeric(G[as.character(train$f.eid), i]) # snp i
  
    # Bivariate model
    g_complete <- g[!is.na(g)]
    X.cov <- cbind(g_complete, (train %>% select("f.21022.0.0","f.22001.0.0",starts_with("PC")))[!is.na(g), ])
    X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept
  
    fit.binormal <- SurrogateRegression::FitBNR(
      t = train$int[!is.na(g)],
      s = train$yhat[!is.na(g)],
      X = X.cov
    )
    out <- c(colnames(G)[i], fit.binormal@Regression.tab %>%
               filter(Outcome == "Target" & Coefficient == "g_complete") %>%
               select(Point, SE, p) %>% as.numeric())
  
    vas <- c("bi")
    vis <- c("beta", "se", "p")
    names(out) <- c("rsid", as.vector(t(outer(vas, vis, paste, sep = "."))))
    return(out)
})

stopCluster(cl)

results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, file = paste0("pheno=", field, "_reuse2_missing=", missing_rate,".rds"))