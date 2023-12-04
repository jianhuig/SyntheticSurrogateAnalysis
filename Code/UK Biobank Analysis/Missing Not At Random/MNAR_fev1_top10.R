library(dplyr)
library(data.table)
library(doParallel)

#setwd("/Users/jianhuigao/Desktop/SyntheticSurrogateAnalysis/Data")

# read height data
dat <- readRDS('field_20150_cleaned.rds')

# Filter UKB caucasian 
white <- data.table::fread("Ancestry.tab")
dat <- dat %>% inner_join(white) %>% filter(!is.na(f.22006.0.0))

# remove NA heights
dat <- dat %>% filter(!is.na(f.20150.0.0))
# upper and lower 10%
upper <- quantile(dat$f.20150.0.0, 0.9, na.rm = TRUE)
lower <- quantile(dat$f.20150.0.0, 0, na.rm = TRUE)
dat <- dat %>% mutate(ind = between(f.20150.0.0, lower, upper))

# sample split
train_id <- fread("final.king.cutoff.out.id") %>% rename(f.eid = FID) %>% select(-IID)
test_id <- fread("final.king.cutoff.in.id") %>% rename(f.eid = FID) %>% select(-IID)

train <- dat %>% inner_join(train_id)
test <- dat %>% inner_join(test_id)  %>% mutate(age = f.21022.0.0, 
                                                sex = f.22001.0.0,
                                                waist = f.48.0.0, 
                                                height = f.50.0.0,
                                                hip = f.49.0.0,
                                                smoking = f.20160.0.0,
                                                bmi = f.21001.0.0,
                                                basal = f.23105.0.0,
                                                impedance = f.23106.0.0,
                                                fev1 = f.20150.0.0) %>% select(f.eid, height, age, sex, waist, hip, smoking, bmi, basal, impedance, fev1, int, starts_with("PC"), ind) %>%
  rename(oracle_int = int) %>%
  na.omit()

# Random Forest Model
train <- train %>% mutate(age = f.21022.0.0, 
                          sex = f.22001.0.0,
                          waist = f.48.0.0, 
                          height = f.50.0.0,
                          hip = f.49.0.0,
                          smoking = f.20160.0.0,
                          bmi = f.21001.0.0,
                          basal = f.23105.0.0,
                          impedance = f.23106.0.0,
                          fev1 = f.20150.0.0) %>%
  filter(ind) %>%
  select(fev1, age, sex, waist, height, smoking, bmi, basal, impedance) %>% 
  na.omit() 

model.rf <- ranger::ranger(
  data = train,
  fev1 ~ .
)

yhat <- predict(model.rf, 
                data = test %>%
                  select(age, sex, waist, height, smoking, bmi, basal, impedance))$predictions
test$yhat <- yhat

# Set NAs
test$int <- test$fev1
test$int[which(test$ind == FALSE)] <- NA

# Inverse Normal Transformation
k = 0.375
data.complete <- test[which(test$ind == TRUE), ]
n <- nrow(data.complete)
r <- rank(data.complete %>% pull(int))
test$int[which(test$ind == TRUE)] <- qnorm((r - k) / (n - 2 * k + 1))

# Inverse Normal Transformation
n <- nrow(test)
r <- rank(test %>% pull(yhat))
test$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) #INT of yhat

# Merge genetic PC
#pcs <- fread("UKB_PC.tab")[,1:11]
#colnames(pcs) <- c("f.eid", paste0("PC",1:10))
#test <- test %>% inner_join(pcs, by = "f.eid")

# Genetic Association Test
G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

clusterExport(cl, list("test"))
clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
  g <- as.numeric(G[as.character(test$f.eid), i]) # snp i
  # oracle
  oracle <- lm(oracle_int ~ ., data = data.frame(
    cbind(test %>% select(oracle_int, age, sex, starts_with("PC")), g)
  ))
  out <- summary(oracle)$coefficients["g", c(1,2,4)] # only beta, se, p-val
  
  # observed phenotype ~ intercept + age + sex + 10 genetic PC + SNP_i
  assoc.observed <- lm(int ~ ., data = data.frame(
    cbind(test %>% select(int, age, sex, starts_with("PC")), g)
  ))
  out <- c(out, summary(assoc.observed)$coefficients["g", c(1,2,4)]) # only beta, se, p-val
  
  # Bivariate model
  g_complete <- g[!is.na(g)]
  X.cov <- cbind(g_complete, (test %>% select(age, sex,starts_with("PC")))[!is.na(g), ])
  X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept
  
  fit.binormal <- SurrogateRegression::FitBNR(
    t = test$int[!is.na(g)],
    s = test$yhat_int[!is.na(g)],
    X = X.cov
  )
  out <- c(out, fit.binormal@Regression.tab %>%
             filter(Outcome == "Target" & Coefficient == "g_complete") %>%
             select(Point, SE, p) %>% as.numeric())
  
  vas <- c("oracle","marginal","bi")
  vis <- c("beta", "se", "p")
  out <- c(colnames(G)[i], out)
  names(out) <- c("rsid", as.vector(t(outer(vas, vis, paste, sep = "."))))
  return(out)
})

stopCluster(cl)
results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, file = paste0("fev1_remove_top10.rds"))