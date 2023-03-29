library(dplyr)

# read all files with .linear extension
filenames <- list.files(
    path = "Results/ImputationAnalysis",
    pattern = ".linear",
    full.names = TRUE
)

# get the names of the files
# only keep everything after plink2
# and before .linear

for (filename in filenames) {
    # get the name of the file
    name <- gsub("(?:[^.]+\\.)([^.]+).*", "\\1", filename)
    temp <- data.table::fread(filename)
    assign(paste0(name), temp)
}

# read SynSurr results
synsurr_results <- readRDS("Results/ImputationAnalysis/SynSurr_height.rds")

# SNPs with p < 5e-8
oracle_snps <- oracle_int %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# Single Imputation Mean
single_mean <- imputed_mean %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# recover rate
recover_rate <- c()
recover_rate <- c(recover_rate, mean(oracle_snps$ID %in% single_mean$ID))

# Single Impuation Random Forest
single_rf <- imputed_rf_1 %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# recover rate
recover_rate <- c(recover_rate, mean(oracle_snps$ID %in% single_rf$ID))

# Single Imputation Linear Regression
single_lr <- imputed_linear_1 %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# recover rate
recover_rate <- c(recover_rate, mean(oracle_snps$ID %in% single_lr$ID))

# Single Imputation Permutation
single_perm <- imputed_rf_permute_1 %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# recover rate
recover_rate <- c(recover_rate, mean(oracle_snps$ID %in% single_perm$ID))

# Single Imputation Negation
single_neg <- imputed_rf_negate_1 %>%
    data.frame() %>%
    select(ID, BETA, SE, P) %>%
    filter(P < 5e-8)

# recover rate
recover_rate <- c(recover_rate, mean(oracle_snps$ID %in% single_neg$ID))

names(recover_rate) <- c("SI Mean", "SI RF", "SI LR", "SI Perm", "SI Neg")

# SynSurr Random Forest
synsurr_rf <- synsurr_results %>%
    select(rf.p) %>%
    filter(rf.p < 5e-8)