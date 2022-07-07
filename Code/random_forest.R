library(ranger)
library(readr)
library(dplyr)
library(data.table)
library(tidyr)

setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data")

# split data into train and test
split_data <- function(file, field, in_id_file, out_id_file) {
  if (tools::file_ext(file) == "tab") {
    temp <- read.table(file, header = TRUE)
  }
  if (tools::file_ext(file) == "rds") {
    temp <- readRDS(file)
  }
  in_id <- readRDS(in_id_file)
  out_id <- readRDS(out_id_file)

  # drop if covariates are missing
  cov_column <- colnames(temp %>% select(-paste0("f.", field, ".0.0")))[grep("0.0", colnames(temp %>% select(-paste0("f.", field, ".0.0"))))]

  return(list(
    train = temp %>% filter(f.eid %in% out_id) %>% tidyr::drop_na(paste0("f.", field, ".0.0")) %>% tidyr::drop_na(all_of(cov_column)),
    test = temp %>% filter(f.eid %in% in_id) %>% tidyr::drop_na(all_of(cov_column))
  ))
}

# Prediction model
# Random Forest model
rf <- function(data, field) {
  temp <- data %>% select(grep("0.0", colnames(data)))
  return(ranger::ranger(
    data = temp,
    as.formula(paste0("f.", field, ".0.0", "~."))
  ))
}

# Prediction
rf_prediction <- function(data, field, model) {
  temp <- data %>%
    select(grep("0.0", colnames(data))) %>%
    select(-paste0("f.", field, ".0.0"))
  return(predict(model, data = temp)$predictions)
}

# ========================== Run Model =========================================
field_id <- 23283
combined <- split_data(
  file = "TotalMass.tab", field = field_id,
  in_id_file = "in.id.rds", out_id_file = "out.id.rds"
)
model.rf <- rf(data = combined$train, field = field_id)
y_pred <- data.frame(
  f.eid = combined$test$f.eid,
  yhat = rf_prediction(
    data = combined$test,
    field = field_id,
    model = model.rf
  )
)
saveRDS(y_pred, file = paste0("predicted_", field_id, ".rds")) # prediction data

# ==========================Make Plot ==========================================
phenotype <- "FEV1"
setwd(paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Result/", phenotype))

png("rf_prediction.png", width = 1080)
plot(combined$test %>% pull(paste0("f.", field_id, ".0.0")), y_pred$yhat,
  xlab = paste("Observed", phenotype),
  ylab = paste("Predicted", phenotype),
  main = paste0("Random-Forest"),
)
abline(0, 1, col = "red")
dev.off()

sink(file = paste0(phenotype, ".model.summary.txt"))
print(model.rf)
sink()

colnames(train) <- c("y", "sex", "bmi", "age", "smoking", "tester")

spec <- feature_spec(train, y ~ .) %>%
  step_numeric_column(c("bmi", "age", "tester"), normalizer_fn = scaler_standard()) %>%
  step_categorical_column_with_identity(c("sex", "smoking"), num_buckets = 2) %>%
  fit()
