library(data.table)
library(dplyr)

cov_id <- c("f.21022.0.0", "f.22001.0.0", "f.50.0.0", "f.23098.0.0", "f.23104.0.0") # age,sex,height, weight, BMI # nolint

yhat_all <- c()

for (pheno_id in c("f.23248.2.0", "f.23260.2.0", "f.23265.2.0", "f.23277.2.0", "f.23283.2.0", "f.23287.2.0")) { # nolint: line_length_linter.

  # split data into train and test
  file <- read.table("CombinedTotalMass.tab", header = TRUE, sep = "\t")
  in_id <- fread("final.king.cutoff.in.id") %>%
    rename(f.eid = `#FID`) %>%
    select(-IID)
  out_id <- fread("final.king.cutoff.out.id") %>%
    rename(f.eid = `#FID`) %>%
    select(-IID)

  cov_column <- file %>% select(cov_id)

  # Filter UKB caucasian
  white <- data.table::fread("Ancestry.tab")
  file <- file %>%
    inner_join(white) %>%
    filter(!is.na(f.22006.0.0))


  train <- file %>%
    inner_join(out_id) %>%
    tidyr::drop_na(pheno_id) %>%
    tidyr::drop_na(all_of(cov_id)) %>%
    select(all_of(pheno_id), all_of(cov_id))
  test <- file %>%
    inner_join(in_id) %>%
    tidyr::drop_na(all_of(cov_id)) %>%
    select(f.eid, all_of(pheno_id), all_of(cov_id))

  # Random Forest model
  model_rf <- ranger::ranger(
    data = train,
    as.formula(paste0(pheno_id, "~."))
  )

  yhat <- predict(model_rf,
    data = test %>% select(-f.eid, -pheno_id)
  )$predictions
  test$yhat <- yhat

  # Merge genetic PC
  pcs <- fread("UKB_PC.tab")[, 1:11]
  colnames(pcs) <- c("f.eid", paste0("PC", 1:10))
  test <- test %>% inner_join(pcs, by = "f.eid")

  # Inverse Normal Transformation
  int_trans <- function(data, pheno, k = 0.375) {
    missing_index <- which(is.na(data[[pheno]]))
    data_complete <- data[-missing_index, ]
    data_missing <- data[missing_index, ]

    n <- nrow(data_complete)
    r <- rank(data_complete %>% pull(pheno))
    data_complete$int <- qnorm((r - k) / (n - 2 * k + 1))
    data_missing$int <- NA

    return(rbind(data_complete, data_missing))
  }
  test <- int_trans(test, pheno_id) # INT of phenotype
  n <- nrow(test)
  r <- rank(test %>% pull(yhat))
  test$yhat_int <- qnorm((r - 0.375) / (n - 2 * 0.375 + 1)) # INT of yhat

  yhat_all <- cbind(yhat_all, test$yhat_int)
}

# prepare for GWAS
pheno_pred <- test %>%
  select(f.eid) %>%
  rename(IID = f.eid) %>%
  mutate(`#FID` = IID) %>%
  select(`#FID`, IID)

pheno_pred <- cbind(pheno_pred, yhat_all)
colnames(pheno_pred) <- c("FID", "IID", c("f.23248.2.0", "f.23260.2.0", "f.23265.2.0", "f.23277.2.0", "f.23283.2.0", "f.23287.2.0"))

# write to txt file
write.table(pheno_pred,
  paste0("pheno_pred.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# prepare covariates file

cov <- test %>%
  select(f.eid, "f.21022.0.0", "f.22001.0.0", starts_with("PC")) %>%
  rename(IID = f.eid) %>%
  mutate(`#FID` = IID) %>%
  select(`#FID`, IID, "f.21022.0.0", "f.22001.0.0", starts_with("PC"))

# write to DEXA covariates txt file
write.table(cov, "DEXA_cov.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# run GWAS association with predicted phenotype
system("sbatch association.sh")

# read GWAS and write to csv file
for(pheno_id in c("f.23248.2.0", "f.23260.2.0", "f.23265.2.0", "f.23277.2.0", "f.23283.2.0", "f.23287.2.0")){
  gwas <- fread(paste0("plink2.", pheno_id, ".glm.linear")) %>%
    select(`ID`,`#CHROM`, `POS`, `REF`, `ALT`, `P`) %>%
    filter(P < 1e-5)
  write.csv(gwas, paste0("gwas_", pheno_id, ".csv"), row.names = FALSE)
}


