library(ranger)
library(readr)
library(dplyr)
library(data.table)
library(tidyr)

setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/PHD_thesis/SyntheticSurrogateAnalysis/Data")

# load Dataset
pheno <- fread("Phenotype.tab")
waist <- fread("waist.tab")
weight <- read_table2("Weight.txt")
pheno <- cbind(pheno, weight[,2], waist[,2])

# Only keep measurement from initial visit
field <- "50" # weight
fieldlist <- c(field, "22001", "21022", "21002", "48") # height, sex, age, weight, waist
pheno <- pheno %>% select(c("f.eid", paste0("f.", fieldlist, ".0.0")))
colnames(pheno) <- c("id", "height", "sex", "age", "weight", "waist")

# Read individuals removed for association analysis
ID <- fread("final.king.cutoff.out.id")[, 1]
colnames(ID) <- c("id")
ID$id <- as.numeric(ID$id)
train <- inner_join(pheno, ID)%>% drop_na()

# Read inviduals in association analysis
load("keep_id.RData")
test <- pheno %>% filter(id %in% keep_id) %>% drop_na()

# Build Prediction model
rf <- function(center = T, ranknormalize = F, train, test){
  # center numerical predictors
  if(center){
    train <- train %>% mutate(age = scale(age, scale = F), 
                            weight = scale(weight, scale = F), 
                            waist = scale(waist, scale = F))
    
    test <- test %>% mutate(age = scale(age, scale = F), 
                              weight = scale(weight, scale = F), 
                              waist = scale(waist, scale = F))
  }
  # rank normal transformation of height
  if(ranknormalize){
    k <- 0.375 # default offset
    n <- nrow(train)
    r <- rank(train$height)
    train <- train %>% mutate(height = qnorm((r - k) / (n - 2 * k + 1)))
    n <- nrow(test)
    r <- rank(test$height)
    test <- test %>% mutate(height = qnorm((r - k) / (n - 2 * k + 1)))
  }
  
  # Random Forest model
  rf <- ranger(
      height ~ sex + age + weight + waist,
      data = train,
    )
  
  return(list(result = c(rf$prediction.error, rf$r.squared),
              height = test$height,
              yhat = predict(rf, data = test[, 2:6])$predictions))
  }


setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/PHD_thesis/SyntheticSurrogateAnalysis/Result/")

results <- rf(train = train, test = test)
y_pred <- data.frame(id = test$id,yhat = results$yhat)
save(y_pred, file = "ypred.RData")


png("rf_prediction.png", width = 1080)

  plot(results[["height"]], results[["yhat"]],
    xlab = "Observed Height",
    ylab = "Predicted Height",
    main = paste0("Random-Forest"),
    xlim = c(120,210),
    ylim = c(120,210)
  )
  abline(0, 1, col = "red")
dev.off()

# rank normalize height
png("rf_prediction_nornamlized.png", width = 1080)
  results <- rf(train = train, test = test, ranknormalize = T)
  plot(results[["height"]], results[["yhat"]],
       xlab = "Observed Height",
       ylab = "Predicted Height",
       main = paste0("Random-Forest")
  )
  abline(0, 1, col = "red")
dev.off()
