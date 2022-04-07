library(ranger)
library(readr)
library(dplyr)
library(data.table)
library(tidyr)

setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data")

# split data into train and test
split_data <- function(field, in_id_file, out_id_file) {
  temp <- readRDS(paste0("field_", field, "_cleaned.rds"))
  in_id <- readRDS(in_id_file)
  out_id <- readRDS(out_id_file)
  
  # drop if covariates are missing
  y_index <- grep("f.3063.0.0", colnames(data))
  cov_column <- colnames(data)[-y_index][grep("0.0", colnames(data)[-y_index])]
  
  return(list(
    train = temp %>% filter(f.eid %in% out_id) %>% tidyr::drop_na(),
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
field_id <- 3063
combined <- split_data(field = field_id, in_id_file = "in.id.rds", out_id_file = "out.id.rds")
model.rf <- rf(data = combined$train, field = field_id)
y_pred <- data.frame(
  f.eid = combined$test$f.eid,
  yhat = rf_prediction(
    data = combined$test,
    field = field_id,
    model = model.rf
  )
)
saveRDS(y_pred, file = paste0("prediected_", field_id, ".rds")) # prediction data

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

sink(file = paste0(phenotype,".model.summary.txt"))
print(model.rf)
sink()
