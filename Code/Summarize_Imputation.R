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
single_lr <- imputed_linear_5 %>%
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

# SynSurr Mean
synsurr_mean <- synsurr_results %>%
select(rsid, mean.beta, mean.se, mean.p) %>%
    rename(ID = rsid, BETA = mean.beta, SE = mean.se, P = mean.p) %>%
    mutate(BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
    filter(P < 5e-8)
# recover rate
recover_rate_synsurr <- mean(oracle_snps$ID %in% synsurr_mean$ID)

# SynSurr Random Forest
synsurr_rf <- synsurr_results %>%
select(rsid, rf.beta, rf.se, rf.p) %>%
    rename(ID = rsid, BETA = rf.beta, SE = rf.se, P = rf.p) %>%
    mutate(BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
    filter(P < 5e-8)
# recover rate
recover_rate_synsurr <- c(recover_rate_synsurr, mean(oracle_snps$ID %in% synsurr_rf$ID))

# SynSurr Linear Regression
synsurr_lr <- synsurr_results %>%
select(rsid, linear.beta, linear.se, linear.p) %>%
    rename(ID = rsid, BETA = linear.beta, SE = linear.se, P = linear.p) %>%
    mutate(BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
    filter(P < 5e-8)
# recover rate
recover_rate_synsurr <- c(recover_rate_synsurr, mean(oracle_snps$ID %in% synsurr_lr$ID))

# SynSurr Permutation
synsurr_perm <- synsurr_results %>%
select(rsid, permute.beta, permute.se, permute.p) %>%
    rename(ID = rsid, BETA = permute.beta, SE = permute.se, P = permute.p) %>%
    mutate(BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
    filter(P < 5e-8)
# recover rate
recover_rate_synsurr <- c(recover_rate_synsurr, mean(oracle_snps$ID %in% synsurr_perm$ID))

# SynSurr Negation
synsurr_neg <- synsurr_results %>%
select(rsid, negate.beta, negate.se, negate.p) %>%
    rename(ID = rsid, BETA = negate.beta, SE = negate.se, P = negate.p) %>%
    mutate(BETA = as.numeric(BETA), SE = as.numeric(SE), P = as.numeric(P)) %>%
    filter(P < 5e-8)
# recover rate
recover_rate_synsurr <- c(recover_rate_synsurr, mean(oracle_snps$ID %in% synsurr_neg$ID))
names(recover_rate_synsurr) <- c("SynSurr Mean", "SynSurr RF", "SynSurr LR", "SynSurr Perm", "SynSurr Neg")
rbind(recover_rate, recover_rate_synsurr)

# Multiple Imputation Mean
# recover rate
recover_rate_mi <- c(NA)

# Multiple Impuation Random Forest
# combine all the results using rubin's rule
mi_rf <- imputed_rf_1 %>% select(ID, BETA, SE, P) %>%
inner_join(imputed_rf_2 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_3 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_4 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_5 %>% select(ID, BETA, SE, P), by = "ID") %>%
mutate(beta_pooled = (BETA.x + BETA.y + BETA.x.x + BETA.y.y + BETA)/5,
       within_variance = (SE.x^2 + SE.y^2 + SE.x.x^2 + SE.y.y^2 + SE^2)/5,
       between_variance = (BETA.x - beta_pooled)^2/4 + (BETA.y - beta_pooled)^2/4 + (BETA.x.x - beta_pooled)^2/4 + (BETA.y.y - beta_pooled)^2/4 + (BETA - beta_pooled)^2/4,
       se_pooled = sqrt(within_variance + between_variance + between_variance/5),
       t_pooled = beta_pooled/se_pooled,
       lambda = (between_variance + between_variance/5)/(within_variance + between_variance + between_variance/5),
       p_pooled = 2*pt(-abs(t_pooled), df = 4/lambda^2)) %>%
select(ID, beta_pooled, se_pooled, p_pooled) %>%
rename(BETA = beta_pooled, SE = se_pooled, P = p_pooled) %>%
filter(P < 5e-8)
# recover rate
recover_rate_mi <- c(recover_rate_mi, mean(oracle_snps$ID %in% mi_rf$ID))

# Multiple Imputation Linear Regression
# combine all the results using rubin's rule
mi_lr <- imputed_linear_1 %>% select(ID, BETA, SE, P) %>%
inner_join(imputed_linear_2 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_linear_3 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_linear_4 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_linear_5 %>% select(ID, BETA, SE, P), by = "ID") %>%
mutate(beta_pooled = (BETA.x + BETA.y + BETA.x.x + BETA.y.y + BETA)/5,
       within_variance = (SE.x^2 + SE.y^2 + SE.x.x^2 + SE.y.y^2 + SE^2)/5,
       between_variance = (BETA.x - beta_pooled)^2/4 + (BETA.y - beta_pooled)^2/4 + (BETA.x.x - beta_pooled)^2/4 + (BETA.y.y - beta_pooled)^2/4 + (BETA - beta_pooled)^2/4,
       se_pooled = sqrt(within_variance + between_variance + between_variance/5),
       t_pooled = beta_pooled/se_pooled,
       lambda = (between_variance + between_variance/5)/(within_variance + between_variance + between_variance/5),
       p_pooled = 2*pt(-abs(t_pooled), df = 4/lambda^2)) %>%
select(ID, beta_pooled, se_pooled, p_pooled) %>%
rename(BETA = beta_pooled, SE = se_pooled, P = p_pooled) %>%
filter(P < 5e-8)
# recover rate
recover_rate_mi <- c(recover_rate_mi, mean(oracle_snps$ID %in% mi_lr$ID))

# Multiple Imputation Permutation
# combine all the results using rubin's rule
mi_perm <- imputed_rf_permute_1 %>% select(ID, BETA, SE, P) %>%
inner_join(imputed_rf_permute_2 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_permute_3 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_permute_4 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_permute_5 %>% select(ID, BETA, SE, P), by = "ID") %>%
mutate(beta_pooled = (BETA.x + BETA.y + BETA.x.x + BETA.y.y + BETA)/5,
       within_variance = (SE.x^2 + SE.y^2 + SE.x.x^2 + SE.y.y^2 + SE^2)/5,
       between_variance = (BETA.x - beta_pooled)^2/4 + (BETA.y - beta_pooled)^2/4 + (BETA.x.x - beta_pooled)^2/4 + (BETA.y.y - beta_pooled)^2/4 + (BETA - beta_pooled)^2/4,
       se_pooled = sqrt(within_variance + between_variance + between_variance/(5 - 1)),
       t_pooled = beta_pooled/se_pooled,
       lambda = (between_variance + between_variance/(5 - 1))/(within_variance + between_variance + between_variance/(5 - 1)),
       p_pooled = 2*pt(-abs(t_pooled), df = 4/lambda^2)) %>%
select(ID, beta_pooled, se_pooled, p_pooled) %>%
rename(BETA = beta_pooled, SE = se_pooled, P = p_pooled) %>%
filter(P < 5e-8)
# recover rate
recover_rate_mi <- c(recover_rate_mi, mean(oracle_snps$ID %in% mi_perm$ID))

# Multiple Imputation Negation
# combine all the results using rubin's rule
mi_neg <- imputed_rf_negate_1 %>% select(ID, BETA, SE, P) %>%
inner_join(imputed_rf_negate_2 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_negate_3 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_negate_4 %>% select(ID, BETA, SE, P), by = "ID") %>%
inner_join(imputed_rf_negate_5 %>% select(ID, BETA, SE, P), by = "ID") %>%
mutate(beta_pooled = (BETA.x + BETA.y + BETA.x.x + BETA.y.y + BETA)/5,
       within_variance = (SE.x^2 + SE.y^2 + SE.x.x^2 + SE.y.y^2 + SE^2)/5,
       between_variance = (BETA.x - beta_pooled)^2/4 + (BETA.y - beta_pooled)^2/4 + (BETA.x.x - beta_pooled)^2/4 + (BETA.y.y - beta_pooled)^2/4 + (BETA - beta_pooled)^2/4,
       se_pooled = sqrt(within_variance + between_variance + between_variance/(5 - 1)),
       t_pooled = beta_pooled/se_pooled,
       lambda = (between_variance + between_variance/(5 - 1))/(within_variance + between_variance + between_variance/(5 - 1)),
       p_pooled = 2*pt(-abs(t_pooled), df = 4/lambda^2)) %>%
select(ID, beta_pooled, se_pooled, p_pooled) %>%
rename(BETA = beta_pooled, SE = se_pooled, P = p_pooled) %>%
filter(P < 5e-8)
# recover rate
recover_rate_mi <- c(recover_rate_mi, mean(oracle_snps$ID %in% mi_neg$ID))
names(recover_rate_mi) <- c("MI mean", "MI RF", "MI LR", "MI Perm", "MI Neg")

# combine all the results
c(recover_rate, recover_rate_synsurr, recover_rate_mi)

est <- cbind(imputed_linear_1$BETA, imputed_linear_2$BETA, imputed_linear_3$BETA, imputed_linear_4$BETA, imputed_linear_5$BETA)
se <- cbind(imputed_linear_1$SE, imputed_linear_2$SE, imputed_linear_3$SE, imputed_linear_4$SE, imputed_linear_5$SE)
m <- rowMeans(est)
var_within <- rowMeans(se^2)
var_between <- rowSums((est - m)^2)/4
var_total <- var_within + (1+1/5)* var_between
se_pooled <- sqrt(var_total)
t_pooled <- m/se_pooled
lambda <- (var_between + var_between/5)/(var_total)
df = 4/lambda^2
p_pooled <- 2*pt(-abs(t_pooled), df = df)


