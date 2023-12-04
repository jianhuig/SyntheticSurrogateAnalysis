library(CALIBERrfimpute)
library(data.table)
library(dplyr)
library(randomForest)

styler::style_file("Code/Write_Imputation.R")


# Helper Functions ============================================================
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

merge_data <- function(pdata, field = field, int = TRUE) {
  rf_model <- readRDS(paste0("prediected_", field, ".rds"))
  if (int) {
    k <- 0.375
    n <- nrow(rf_model)
    r <- rank(rf_model$yhat)
    rf_model$yhat <- qnorm((r - k) / (n - 2 * k + 1))
  }
  return(rf_model %>% inner_join(pdata, by = "f.eid"))
}

# Height Analysis ============================================================
setwd("Data/Old/")
missing_rate <- 0.5
field <- 50

# White British
ancestry <- fread("Ancestry.tab") %>% filter(!is.na(f.22006.0.0))

# train data
train <- create_missing_phenotypes(
  field = field,
  in_id_file = "final.king.cutoff.out.id", rate = missing_rate
)

# merge with White British
train <- train %>% inner_join(ancestry, by = "f.eid")

# test data
test <- create_missing_phenotypes(
  field = field,
  in_id_file = "final.king.cutoff.in.id", rate = missing_rate
)

test <- merge_data(pdata = test, field = field, int = TRUE)

# merge with White British
test <- test %>% inner_join(ancestry, by = "f.eid")

# add covriates
weight <- fread("height.tab") %>% select(f.eid, f.21002.0.0, f.50.0.0)

# merge with train and test
train <- train %>% inner_join(weight, by = "f.eid")
test <- test %>% inner_join(weight, by = "f.eid")

# covariate columns
cov_column <- c("f.21022.0.0", "f.22001.0.0", "f.21002.0.0", "f.48.0.0")

# remove missing covaraites and phenotypes for training
train <- train %>%
  tidyr::drop_na(paste0("f.", field, ".0.0")) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column))

test_temp <- test %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column)) %>%
  mutate(!!paste0("f.", field, ".0.0") := NA) # Only use train data to train

# Multiple Random Forest Imputation
# impute 5 times
test_imputed_rf <- mice(
  data = rbind(train, test_temp),
  method = "rf",
  m = 5,
  seed = 123
)

# Multiple linear regression imputation
# impute 5 times
test_imputed_linear <- mice(
  data = rbind(
    train %>% select(!!paste0("f.", field, ".0.0"), f.48.0.0),
    test_temp %>% select(!!paste0("f.", field, ".0.0"), f.48.0.0)
  ),
  method = "norm",
  seed = 123
)

# permute random forest imputation
test_imputed_rf_permute <- test_imputed_rf
for (m in 1:5) {
  test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][, m] <- sample(test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][, m])
}

# negate random forest imputation
test_imputed_rf_negate <- test_imputed_rf
for (m in 1:5) {
  test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][, m] <- -test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][, m]
}

# Add mean imputation
test_imputed_mean <- mice(
  data = rbind(
    train %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0),
    test_temp %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0)
  ),
  method = "mean",
  seed = 123,
  m = 1
)

# missing index in test
missing_index <- which(is.na(test$int))

# Merge imputed data
for (m in 1:5) {
  test <- test %>% mutate(!!paste0("imputed_rf_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_linear_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_permute_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_negate_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  # replace NAs with imputed values
  test[[paste0("imputed_rf_", m)]][missing_index] <-
    test_imputed_rf$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with linear imputed values
  test[[paste0("imputed_linear_", m)]][missing_index] <-
    test_imputed_linear$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_linear_", m)]] <-
    qnorm((rank(test[[paste0("imputed_linear_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with permuted random forest imputed values
  test[[paste0("imputed_rf_permute_", m)]][missing_index] <-
    test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_permute_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_permute_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with negated random forest imputed values
  test[[paste0("imputed_rf_negate_", m)]][missing_index] <-
    test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_negate_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_negate_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
}

# Add mean imputation
test <- test %>% mutate(!!paste0("imputed_mean") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_mean")]][missing_index] <- test_imputed_mean$imp[[paste0("f.", field, ".0.0")]][missing_index, 1]
# inverse normal transformation to imputed data
test[[paste0("imputed_mean")]] <-
  qnorm((rank(test[[paste0("imputed_mean")]]) - 0.375) /
    (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Imputation
model <- ranger::ranger(
    data = train ,
    as.formula(paste0("f.", field, ".0.0", "~."))
  )
imputed <- predict(model, test_temp)$predictions
test <- test %>% mutate(!!paste0("imputed_rf") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf"]] <- qnorm((rank(test[["imputed_rf"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Negate Imputation
imputed <- -imputed
test <- test %>% mutate(!!paste0("imputed_rf_negate") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf_negate")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf_negate"]] <- qnorm((rank(test[["imputed_rf_negate"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Permute Imputation
imputed <- sample(-imputed)
test <- test %>% mutate(!!paste0("imputed_rf_permute") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf_permute")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf_permute"]] <- qnorm((rank(test[["imputed_rf_permute"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))


# Add Single Linear Imputation
model <- lm(as.formula(paste0("f.", field, ".0.0", " ~ f.48.0.0")), data = train)
imputed <- predict(model, test_temp)
test <- test %>% mutate(!!paste0("imputed_linear") := imputed)
# inverse normal transformation to imputed data
test[["imputed_linear"]] <- qnorm((rank(test[["imputed_linear"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))


# Write imputed data to txt file
write.table(test %>%
  select(f.eid, starts_with("imputed_"), oracle) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, starts_with("imputed"), oracle) %>%
  rename(`#FID` = f.eid), "height_imputed.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# write covariate data to txt file
write.table(test %>%
  select(f.eid, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  rename(`#FID` = f.eid), "height_covariate.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# write imputed data to Run SynSurr
saveRDS(test %>%
  select(
    f.eid, int, yhat, imputed_linear, imputed_mean,
    f.22001.0.0, f.21022.0.0, starts_with("PC")
  ), "height_imputed.rds")

# FEV1 Analysis ============================================================
setwd("Data/Old/")
missing_rate <- 0.5
field <- 20150


# White British
ancestry <- fread("Ancestry.tab") %>% filter(!is.na(f.22006.0.0))

# train data
train <- create_missing_phenotypes(
  field = field,
  in_id_file = "final.king.cutoff.out.id", rate = missing_rate
)

# merge with White British
train <- train %>% inner_join(ancestry, by = "f.eid")

# test data
test <- create_missing_phenotypes(
  field = field,
  in_id_file = "final.king.cutoff.in.id", rate = missing_rate
)

test <- merge_data(pdata = test, field = 20150)

# merge with White British
test <- test %>% inner_join(ancestry, by = "f.eid")

# load FEV1 tab
fev1 <- fread("FEV1.tab") %>% select(f.eid, f.20150.0.0)
train <- train %>% inner_join(fev1, by = "f.eid")
test <- test %>% inner_join(fev1, by = "f.eid")


# covariate columns
cov_column <- c("f.21022.0.0", "f.22001.0.0", "f.50.0.0", "f.48.0.0", "f.49.0.0", "f.20160.0.0", "f.21001.0.0", "f.23105.0.0", "f.23106.0.0")

# remove missing covaraites and phenotypes for training
train <- train %>%
  tidyr::drop_na(paste0("f.", field, ".0.0")) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column))

test_temp <- test %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column)) %>%
  mutate(!!paste0("f.", field, ".0.0") := NA) # Only use train data to train

# Multiple Random Forest Imputation
# impute 5 times
test_imputed_rf <- mice(
  data = rbind(train, test_temp),
  method = "rf",
  m = 5,
  seed = 123
)

# Multiple linear regression imputation
# impute 5 times
test_imputed_linear <- mice(
  data = rbind(
    train %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0),
    test_temp %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0)
  ),
  method = "norm",
  seed = 123
)

# permute random forest imputation
test_imputed_rf_permute <- test_imputed_rf
for (m in 1:5) {
  test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][, m] <- sample(test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][, m])
}

# negate random forest imputation
test_imputed_rf_negate <- test_imputed_rf
for (m in 1:5) {
  test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][, m] <- -test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][, m]
}

# Add mean imputation
test_imputed_mean <- mice(
  data = rbind(
    train %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0),
    test_temp %>% select(!!paste0("f.", field, ".0.0"), f.21022.0.0, f.22001.0.0)
  ),
  method = "mean",
  seed = 123,
  m = 1
)

# missing index in test
missing_index <- which(is.na(test$int))

# Merge imputed data
for (m in 1:5) {
  test <- test %>% mutate(!!paste0("imputed_rf_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_linear_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_permute_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_negate_", m) :=
    .data[[paste0("f.", field, ".0.0")]])
  # replace NAs with imputed values
  test[[paste0("imputed_rf_", m)]][missing_index] <-
    test_imputed_rf$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with linear imputed values
  test[[paste0("imputed_linear_", m)]][missing_index] <-
    test_imputed_linear$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_linear_", m)]] <-
    qnorm((rank(test[[paste0("imputed_linear_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with permuted random forest imputed values
  test[[paste0("imputed_rf_permute_", m)]][missing_index] <-
    test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_permute_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_permute_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with negated random forest imputed values
  test[[paste0("imputed_rf_negate_", m)]][missing_index] <-
    test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_negate_", m)]] <-
    qnorm((rank(test[[paste0("imputed_rf_negate_", m)]]) - 0.375) /
      (nrow(test) - 2 * 0.375 + 1))
}

# Add mean imputation
test <- test %>% mutate(!!paste0("imputed_mean") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_mean")]][missing_index] <- test_imputed_mean$imp[[paste0("f.", field, ".0.0")]][missing_index, 1]
# inverse normal transformation to imputed data
test[[paste0("imputed_mean")]] <-
  qnorm((rank(test[[paste0("imputed_mean")]]) - 0.375) /
    (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Imputation
model <- ranger::ranger(
    data = train ,
    as.formula(paste0("f.", field, ".0.0", "~."))
  )
imputed <- predict(model, test_temp)$predictions
test <- test %>% mutate(!!paste0("imputed_rf") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf"]] <- qnorm((rank(test[["imputed_rf"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Negate Imputation
imputed <- -imputed
test <- test %>% mutate(!!paste0("imputed_rf_negate") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf_negate")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf_negate"]] <- qnorm((rank(test[["imputed_rf_negate"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))

# Add Single RF Permute Imputation
imputed <- sample(-imputed)
test <- test %>% mutate(!!paste0("imputed_rf_permute") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_rf_permute")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_rf_permute"]] <- qnorm((rank(test[["imputed_rf_permute"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))


# Add Single Linear Imputation
model <- lm(as.formula(paste0("f.", field, ".0.0", " ~ f.21022.0.0 + f.22001.0.0")), data = train)
imputed <- predict(model, test_temp)
test <- test %>% mutate(!!paste0("imputed_linear") :=
  .data[[paste0("f.", field, ".0.0")]])
test[[paste0("imputed_linear")]][missing_index] <- imputed[missing_index]
# inverse normal transformation to imputed data
test[["imputed_linear"]] <- qnorm((rank(test[["imputed_linear"]]) - 0.375) /
  (nrow(test) - 2 * 0.375 + 1))


# Write imputed data to txt file
write.table(test %>%
  select(f.eid, starts_with("imputed_"), oracle) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, starts_with("imputed"), oracle) %>%
  rename(`#FID` = f.eid), "fev1_imputed.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# write covariate data to txt file
write.table(test %>%
  select(f.eid, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  rename(`#FID` = f.eid), "fev1_covariate.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# write imputed data to Run SynSurr
saveRDS(test %>%
  select(
    f.eid, int, yhat, imputed_linear, imputed_mean,
    f.22001.0.0, f.21022.0.0, starts_with("PC")
  ), "fev1_imputed.rds")