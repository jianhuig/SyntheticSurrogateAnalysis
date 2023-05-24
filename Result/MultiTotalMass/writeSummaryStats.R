bim <- read.table("~/Desktop/SyntheticSurrogateAnalysis/Data/final.bim", header = TRUE)
library(dplyr)
library("xlsx")

for(pheno in c("f.23265.2.0", "f.23260.2.0", "f.23283.2.0", "f.23277.2.0",
               "f.23287.2.0", "f.23248.2.0")){
  result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/pheno=", pheno, "_result.rds"))
  result <- result %>% inner_join(bim, by = "rsid") %>% rename(SNP = rsid, POS = BP,
                                                               BETA = bi.beta,
                                                               SE = bi.se,
                                                               P = bi.p
                                                        ) %>%
    mutate(Z = as.numeric(BETA)/as.numeric(SE))
  write.table(result %>% select(SNP, CHR, POS, Ref, Alt, BETA, SE, Z, P),
             file = paste0("SumStats_pheno=",pheno, ".txt"),
             row.names = FALSE,
            quote = FALSE)
}
