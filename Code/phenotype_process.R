library(data.table)
library(dplyr)

### ======================== phenotype Data =================###
pheno <- fread("Phenotype.tab")
# select measurement at baseline
field <- "50" # height
fieldlist <- c(field, "21000", "21022", "22001") # height, ethnicity, age, sex
pheno <- pheno %>% select(c("f.eid", paste0("f.", fieldlist, ".0.0")))
# subset British population
pheno <- pheno %>% filter(f.21000.0.0 == 1001) # 442,551 individuals left
load("keep_id.RData") # IDs that remove closed individuals
pheno <- pheno %>% filter(f.eid %in% keep_id)

# Inverse-normal transformation for continuous phenotype
k <- 0.375 # default offset
n <- nrow(pheno)
r <- rank(pheno %>% pull(paste0("f.", field, ".0.0")))
pheno <- pheno %>% mutate(int = qnorm((r - k) / (n - 2 * k + 1)))

# Save phenotype and covariate data
eigenvec <- fread("plink2.eigenvec") %>% rename(f.eid = `#FID`)
pheno <- pheno %>% left_join(eigenvec, by = "f.eid")
write.table(pheno %>% select(c("f.eid", "IID", "int")), file = "phenotype.txt", sep = " ", row.names = F, quote = F, col.names = F)
write.table(pheno %>% select(c("f.eid", "IID", paste0("f.", c("21022", "22001"), ".0.0"), paste0("PC", 1:10))), file = "covariate.txt", sep = " ", row.names = F, quote = F, col.names = F)
