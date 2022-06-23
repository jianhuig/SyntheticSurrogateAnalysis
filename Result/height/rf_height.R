library(ranger)


field_id <- 50

dat <- readRDS(paste0("field_", field_id, "_cleaned.rds"))
in_id <- readRDS("in.id.rds")
out_id <- readRDS("out.id.rds")

train <- dat %>% filter(f.eid %in% out_id)
weight <- read.table("Weight.txt", fill = TRUE, header = TRUE)
train <- train %>%
  left_join(weight %>% select(f.eid, "f.21002.0.0"), by = "f.eid") %>%
  select(
    f.eid, f.21022.0.0, f.22001.0.0,
    f.48.0.0, f.21002.0.0, f.49.0.0, f.50.0.0
  ) %>%
  tidyr::drop_na()


train.X <- train %>% select(
  f.21022.0.0, f.22001.0.0,
  f.48.0.0, f.21002.0.0, f.49.0.0
)

train.Y <- train$'f.50.0.0'

rf.model <- ranger::ranger(
  data = cbind(y = train.Y, x = train.X),
  formula = y ~ .
)

test <- dat %>% filter(f.eid %in% in_id)
test <- test %>%
  left_join(weight %>% select(f.eid, "f.21002.0.0"), by = "f.eid") %>%
  select(
    f.eid, f.21022.0.0, f.22001.0.0,
    f.48.0.0, f.21002.0.0, f.49.0.0, f.50.0.0
  ) %>%
  tidyr::drop_na()
test.X <- test %>% select(
  f.21022.0.0, f.22001.0.0,
  f.48.0.0, f.21002.0.0, f.49.0.0
)

yhat <- predict(rf.model, data = data.frame(x = test.X))$predictions


saveRDS(cbind(f.eid = test$f.eid, yhat = yhat) %>% data.frame(), 
        file = "predicted_50.rds")