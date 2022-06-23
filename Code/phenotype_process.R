library(data.table)
library(dplyr)

# obatin phenotype
## cd /home/jianhuig/projects/def-leisun/jianhuig
## system("./ukbconv ukb47570.enc_ukb r -iphenotype_list.txt -oFEV1")

get_baseline_measure <- function(file_name){
	temp <- fread(file_name)
	return(temp%>% select(f.eid, grep(".0.0", colnames(temp)))) # select measurement at baseline
}

INT <- function(data, field, k = 0.375){
  missing_index <- which(is.na(data[[paste0("f.", field, ".0.0")]]))
  data.complete <- data[!missing_index, ]
  data.missing <-data[missing_index, ]
  
	n <- nrow(data.complete)
	r <- rank(data.complete %>% pull(paste0("f.", field, ".0.0")))
	data.complete$int <- qnorm((r - k) / (n - 2 * k + 1))
	data.missing$int <- NA
	
	return(rbind(data.complete, data.missing))
}

merge_pc <- function(data, pc_file, field){
	eigenvec <- fread(pc_file) %>% rename(f.eid = `#FID`) %>% select(-IID)
	temp <- data %>% inner_join(eigenvec, by = "f.eid")
	saveRDS(temp, file = paste0("field_",field,"_cleaned.rds"))
}


pheno <- get_baseline_measure(file_name = "BMD.tab")
pheno <- INT(data = pheno, field = 3148) # Inverse-normal transformation for continuous phenotype
merge_pc(data = pheno, pc_file = "plink2.eigenvec", field = 3148)
