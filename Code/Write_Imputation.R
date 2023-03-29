library(CALIBERrfimpute)
library(data.table)
library(dplyr)
library(randomForest)

styler::style_file("Code/Write_Imputation.R")

# Height Analysis ============================================================
# Load data
height <- fread("Data/height.tab")
field <- 50

# Split data into train and test
in_id <- fread("Data/final.king.cutoff.in.id") %>%
  rename(f.eid = `#FID`) %>%
  select(-IID)

out_id <- fread("Data/final.king.cutoff.out.id") %>%
  rename(f.eid = `#FID`) %>%
  select(-IID)

ancestry <- fread("Data/Ancestry.tab") %>% filter(!is.na(f.22006.0.0))

# drop if covariates are missing
cov_column <- colnames(height %>%
  dplyr::select(-paste0("f.", field, ".0.0")))[grep(
  ".0.0$",
  colnames(height %>% dplyr::select(-paste0("f.", field, ".0.0")))
)]
# remove ethnicity
cov_column <- cov_column[!cov_column %in% c("f.22006.0.0", "f.21000.0.0")]

train <- height %>%
  inner_join(out_id) %>%
  inner_join(ancestry) %>%
  tidyr::drop_na(paste0("f.", field, ".0.0")) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column))

test <- height %>%
  inner_join(in_id) %>%
  inner_join(ancestry) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(f.eid, paste0("f.", field, ".0.0"), all_of(cov_column))

# oracle
test$oracle <- test %>% pull(paste0("f.", field, ".0.0"))

# Create 50% missing data
set.seed(123)
nmissing <- floor(nrow(test) * 0.5)
missing_index <- sample(which(!is.na(test %>%
  pull(paste0("f.", field, ".0.0")))), nmissing)
test[missing_index, paste0("f.", field, ".0.0")] <- NA
sum(is.na(test %>% pull(paste0("f.", field, ".0.0"))))
test_temp <- test %>%
  select(-f.eid, -oracle) %>%
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
train_mean <- train %>%
  pull(paste0("f.", field, ".0.0")) %>%
  mean()

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
test[[paste0("imputed_mean")]][missing_index] <- train_mean
# inverse normal transformation to imputed data
test[[paste0("imputed_mean")]] <-
  qnorm((rank(test[[paste0("imputed_mean")]]) - 0.375) /
    (nrow(test) - 2 * 0.375 + 1))

# Inverse normal transformation to oracle data
oracle_missing <- which(is.na(test$oracle))
test_complete <- test[-oracle_missing, ]
test_missing <- test[oracle_missing, ]
test_complete$oracle_int <- qnorm((rank(test_complete %>%
  pull(oracle)) - 0.375) / (nrow(test_complete) - 2 * 0.375 + 1))
test_missing$oracle_int <- NA
test <- rbind(test_complete, test_missing)

# Inverse normal transformation to original data
missing_index <- which(is.na(test %>% pull(paste0("f.", field, ".0.0"))))
test_complete <- test[-missing_index, ]
test_missing <- test[missing_index, ]
test_complete$int <- qnorm((rank(test_complete %>%
  pull(paste0("f.", field, ".0.0"))) - 0.375) /
  (nrow(test_complete) - 2 * 0.375 + 1))
test_missing$int <- NA
test <- rbind(test_complete, test_missing)

# Write imputed data to txt file
write.table(test %>%
  select(f.eid, starts_with("imputed_"), oracle_int) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, starts_with("imputed"), oracle_int) %>%
  rename(`#FID` = f.eid), "Data/height_imputed.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# Merge genetic PC
pcs <- fread("UKB_PC.tab")[, 1:11]
colnames(pcs) <- c("f.eid", paste0("PC", 1:10))
test <- test %>% inner_join(pcs, by = "f.eid")

# writr covariate data to txt file
write.table(test %>%
  select(f.eid, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  rename(`#FID` = f.eid), "Data/height_covariate.txt",
sep = "\t", row.names = FALSE, quote = FALSE
)

# write imputed data to Run SynSurr
saveRDS(test %>%
  select(
    f.eid, oracle_int, int, imputed_rf_1, imputed_linear_1,
    imputed_rf_permute_1, imputed_rf_negate_1, imputed_mean,
    f.22001.0.0, f.21022.0.0, starts_with("PC")
  ), "Data/height_imputed.rds")

# FEV1 Analysis ============================================================
# Load data
fev1 <- fread("Data/FEV1.tab")
field <- 20150

# Split data into train and test
in_id <- fread("Data/final.king.cutoff.in.id") %>%
  rename(f.eid = `#FID`) %>%
  select(-IID)

out_id <- fread("Data/final.king.cutoff.out.id") %>%
  rename(f.eid = `#FID`) %>%
  select(-IID)

ancestry <- fread("Data/Ancestry.tab") %>% filter(!is.na(f.22006.0.0))

# drop if covariates are missing
cov_column <- colnames(fev1 %>% dplyr::select(-paste0("f.", field, ".0.0")))[grep(".0.0$", colnames(fev1 %>% dplyr::select(-paste0("f.", field, ".0.0"))))]
train <- fev1 %>%
  inner_join(out_id) %>%
  inner_join(ancestry) %>%
  tidyr::drop_na(paste0("f.", field, ".0.0")) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(paste0("f.", field, ".0.0"), all_of(cov_column))

test <- fev1 %>%
  inner_join(in_id) %>%
  inner_join(ancestry) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(f.eid, paste0("f.", field, ".0.0"), all_of(cov_column))

# Create 50% missing data
set.seed(123)
nmissing <- floor(nrow(test) * 0.5)
missing_index <- sample(which(!is.na(test %>% pull(paste0("f.", field, ".0.0")))), nmissing)
test[missing_index, paste0("f.", field, ".0.0")] <- NA
test_temp <- test %>%
  select(-f.eid) %>%
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

# Merge imputed data
for (m in 1:5) {
  test <- test %>% mutate(!!paste0("imputed_rf_", m) := .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_linear_", m) := .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_permute_", m) := .data[[paste0("f.", field, ".0.0")]])
  test <- test %>% mutate(!!paste0("imputed_rf_negate_", m) := .data[[paste0("f.", field, ".0.0")]])
  # replace NAs with imputed values
  test[[paste0("imputed_rf_", m)]][missing_index] <- test_imputed_rf$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_", m)]] <- qnorm((rank(test[[paste0("imputed_rf_", m)]]) - 0.375) / (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with linear imputed values
  test[[paste0("imputed_linear_", m)]][missing_index] <- test_imputed_linear$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_linear_", m)]] <- qnorm((rank(test[[paste0("imputed_linear_", m)]]) - 0.375) / (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with permuted random forest imputed values
  test[[paste0("imputed_rf_permute_", m)]][missing_index] <- test_imputed_rf_permute$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_permute_", m)]] <- qnorm((rank(test[[paste0("imputed_rf_permute_", m)]]) - 0.375) / (nrow(test) - 2 * 0.375 + 1))
  # replace NAs with negated random forest imputed values
  test[[paste0("imputed_rf_negate_", m)]][missing_index] <- test_imputed_rf_negate$imp[[paste0("f.", field, ".0.0")]][missing_index, m]
  # inverse normal transformation to imputed data
  test[[paste0("imputed_rf_negate_", m)]] <- qnorm((rank(test[[paste0("imputed_rf_negate_", m)]]) - 0.375) / (nrow(test) - 2 * 0.375 + 1))
}

# Inverse normal transformation to original data
test_complete <- test[-missing_index, ]
test_missing <- test[missing_index, ]
test_complete$int <- qnorm((rank(test_complete %>% pull(paste0("f.", field, ".0.0"))) - 0.375) / (nrow(test_complete) - 2 * 0.375 + 1))
test_missing$int <- NA
test <- rbind(test_complete, test_missing)


# Write imputed data to txt file
write.table(test %>%
  select(f.eid, starts_with("imputed_")) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, starts_with("imputed")) %>%
  rename(`#FID` = f.eid), "Data/fev1_imputed.txt", sep = "\t", row.names = F, quote = F)

# Merge genetic PC
pcs <- fread("UKB_PC.tab")[, 1:11]
colnames(pcs) <- c("f.eid", paste0("PC", 1:10))
test <- test %>% inner_join(pcs, by = "f.eid")

# writr covariate data to txt file
write.table(test %>%
  select(f.eid, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  mutate(IID = f.eid) %>%
  select(f.eid, IID, f.21022.0.0, f.22001.0.0, starts_with("PC")) %>%
  rename(`#FID` = f.eid), "Data/fev1_covariate.txt", sep = "\t", row.names = F, quote = F)

# write imputed data to Run SynSurr
saveRDS(test %>%
  select(int, imputed_rf_1, imputed_linear_1, imputed_rf_permute_1, imputed_rf_negate_1, f.22001.0.0, f.21022.0.0, starts_with("PC")), "Data/fev1_imputed.rds")
