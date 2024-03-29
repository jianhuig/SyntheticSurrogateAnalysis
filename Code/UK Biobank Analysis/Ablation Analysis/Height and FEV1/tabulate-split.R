library(dplyr)

# read all files with .rds extension
filenames <- list.files(
    path = "Results/SplitAnalysis",
    pattern = ".rds",
    full.names = TRUE
)
output <- c()
for(filename in filenames){
    # get the name of the file
    name <- paste0("pheno_", gsub("(?:[^.]+\\.)([^.]+).*", "\\1", filename))
    temp <- readRDS(filename)

    # assign(paste0(name), temp)

    # results for standard GWAS
    std_discov <- temp %>%
        filter(standard_discovery_p < 5e-8)
    std_valid <- temp %>%
        filter(standard_validation_p < 0.05/nrow(std_discov))
    # recover rate
    std_recov <- nrow(std_valid %>% filter(standard_discovery_p < 5e-8))
    std_recov_rate <- std_recov/nrow(std_discov)
    # results for Synsurr GWAS
    synsurr_discov <- temp %>%
        filter(synsurr_discovery_p < 5e-8)
    synsurr_valid <- temp %>%
        filter(synsurr_validation_p < 0.05/nrow(synsurr_discov))
    # recover rate
    synsurr_recov <- nrow(synsurr_valid %>% filter(synsurr_discovery_p < 5e-8))
    synsurr_recov_rate <- synsurr_recov/nrow(synsurr_discov)

    # append to output
    output <- rbind(output, c(name, nrow(std_discov), std_recov, std_recov_rate, nrow(synsurr_discov), synsurr_recov, synsurr_recov_rate))
}

