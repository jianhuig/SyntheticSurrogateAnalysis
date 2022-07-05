library(keras)
library(tfdatasets)
library(tidyr)
library(dplyr)

setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data")

file <- "BMD.tab"
field <- 3148

temp <- read.table(file, header = TRUE)
in_id <- readRDS("in.id.rds")
out_id <- readRDS("out.id.rds")
# drop if covariates are missing
cov_column <- colnames(temp %>% select(-paste0("f.", field, ".0.0")))[grep("0.0", colnames(temp %>% select(-paste0("f.", field, ".0.0"))))]

train <- temp %>%
  filter(f.eid %in% out_id) %>%
  tidyr::drop_na(paste0("f.", field, ".0.0")) %>%
  tidyr::drop_na(all_of(cov_column)) %>%
  select(grep("0.0", colnames(temp))) %>%
  rename(label = paste0("f.", field, ".0.0"))

# normalize features
spec <- feature_spec(train, label ~ .) %>%
  step_numeric_column(f.20116.0.0, f.22001.0.0) %>%
  step_numeric_column(-f.20116.0.0, -f.22001.0.0, normalizer_fn = scaler_standard()) %>%
  fit()

spec

layer <- layer_dense_features(
  feature_columns = dense_features(spec),
  dtype = tf$float32
)
layer(train)

input <- layer_input_from_dataset(train %>% select(-label))

output <- input %>%
  layer_dense_features(dense_features(spec)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 1)

model <- keras_model(input, output)

summary(model)

model %>% 
  compile(
    loss = "mse",
    optimizer = optimizer_rmsprop(),
    metrics = list("mean_absolute_error")
  )

build_model <- function() {
  input <- layer_input_from_dataset(train_df %>% select(-label))
  
  output <- input %>% 
    layer_dense_features(dense_features(spec)) %>% 
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1) 
  
  model <- keras_model(input, output)
  
  model %>% 
    compile(
      loss = "mse",
      optimizer = optimizer_rmsprop(),
      metrics = list("mean_absolute_error")
    )
  
  model
}

# Display training progress by printing a single dot for each completed epoch.
print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 80 == 0) cat("\n")
    cat(".")
  }
)    

model <- build_model()

history <- model %>% fit(
  x = train %>% select(-label),
  y = train$label,
  epochs = 500,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(print_dot_callback)
)

plot(history)

test <- temp %>% 
  filter(f.eid %in% in_id) %>% 
  tidyr::drop_na(all_of(cov_column)) %>% 
  select(grep("0.0", colnames(temp))) %>%
  rename(label = paste0("f.", field, ".0.0"))


test_predictions <- model %>% predict(test %>% select(-label))
plot(test$label, test_predictions[ , 1])
abline(a = 0, b = 1, col = 'red')

summary(lm(test$label~test_predictions[ , 1]))
