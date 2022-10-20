library(data.table)
library(dplyr)
library(bigsnpr)
library(doParallel)
library(BEDMatrix)

create_missing_phenotypes <- function(field, in_id_file, rate, seed = 123) {
  temp <- readRDS(paste0("field_", field, "_cleaned.rds")) # cleaned pheno + cov
  in_id <- fread("final.king.cutoff.in.id") %>% rename(f.eid = `#FID`) %>% select(-IID)
  temp <- temp %>% inner_join(in_id) %>% filter(!is.na(int))

  # number of missing phenotypes
  nmissing <- floor(nrow(temp) * rate)

  set.seed(seed)
  missing_index <- sample(which(!is.na(temp$int)), size = nmissing)
  
  #oracle
  temp$oracle <- temp$int
  
  # set phenotype missing value
  temp[missing_index, ]$int <- NA
  return(data.frame(temp %>% select(-paste0("f.", field, ".0.0"))))
}

merge_data <- function(pdata, field = field_id, int = TRUE) {
  rf_model <- readRDS(paste0("prediected_", field, ".rds"))
  if (int) {
    k <- 0.375
    n <- nrow(rf_model)
    r <- rank(rf_model$yhat)
    rf_model$yhat <- qnorm((r - k) / (n - 2 * k + 1))
  }
  return(rf_model %>% inner_join(pdata, by = "f.eid"))
}


# args <- (commandArgs(TRUE))
# print(args)
# for (k in 1:length(args)) {
#   eval(parse(text = args[[k]]))
# }
missing_rate = 0.75

field_id <- 50

pdata <- create_missing_phenotypes(field = field_id, 
                                   in_id_file = "in.id.rds", rate = missing_rate)
pheno <- merge_data(pdata = pdata)
G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data

cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

clusterExport(cl, list("pheno"))
clusterEvalQ(cl, {
  library(dplyr)
  G <- BEDMatrix::BEDMatrix(path = "allchromosome.bed", simple_names = T) # read genetic data
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
	tryCatch({
  g <- as.numeric(G[as.character(pheno$f.eid), i]) # snp i
  # observed phenotype ~ intercept + age + sex + 10 genetic PC + SNP_i
  assoc.observed <- lm(int ~ ., data = data.frame(
    cbind(pheno %>% select(int, f.21022.0.0, f.22001.0.0, starts_with("PC")), g)
  ))
  out <- summary(assoc.observed)$coefficients["g", c(1,2,4)] # only beta, se, p-val

  # orcale ~ intercept + age + sex + 10 genetic PC + SNP_i
  assoc.oracle <- lm(oracle ~ ., data = data.frame(
    cbind(pheno %>% select(oracle, f.21022.0.0, f.22001.0.0, starts_with("PC")), g)
  ))
  out <- c(out, summary(assoc.oracle)$coefficients["g", c(1,2,4)])

  # Bivariate model
  g_complete <- g[!is.na(g)]
  X.cov <- cbind(g_complete, (pheno %>% select(f.21022.0.0, f.22001.0.0, starts_with("PC")))[!is.na(g), ])
  X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept
  
  fit.binormal <- SurrogateRegression::Fit.BNLS(
    t = pheno$int[!is.na(g)],
    s = pheno$yhat[!is.na(g)],
    X = X.cov
  )
  out <- c(out, fit.binormal@Regression.tab %>%
    filter(Outcome == "Target" & Coefficient == "g_complete") %>%
    select(Point, SE, p) %>% as.numeric())
	
   vas <- c("marginal","oracle","bi")
    vis <- c("beta", "se", "p")
  names(out) <- as.vector(t(outer(vas, vis, paste, sep = ".")))},
  error=function(e) {
    NULL
  })
  return(out)
})

stopCluster(cl)

results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, file = paste0("binormal_field=", field_id, "_nmissing=",missing_rate, ".rds"))
