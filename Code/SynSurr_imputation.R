library(data.table)
library(dplyr)
library(bigsnpr)
library(doParallel)
library(BEDMatrix)

# read in the data
pheno <- readRDS("Data/FEV1_imputed.rds")

# read in the genetic data
G <- BEDMatrix::BEDMatrix(path = "Data/allchromosome.bed", simple_names = TRUE)

# set up the parallel computing
cl <- makeCluster(type = "MPI")
registerDoParallel(cl)

clusterExport(cl, list("pheno"))
clusterEvalQ(cl, {
    library(dplyr)
    G <- BEDMatrix::BEDMatrix(
        path = "Data/allchromosome.bed",
        simple_names = TRUE
    )
})

results <- parLapply(cl, X = 1:ncol(G), fun = function(i) {
    tryCatch(
        {
            g <- as.numeric(G[as.character(pheno$f.eid), i]) # snp i
            g_complete <- g[!is.na(g)]
            X.cov <- cbind(
                g_complete,
                (pheno %>%
                    select(
                        f.21022.0.0, f.22001.0.0,
                        starts_with("PC")
                    ))[!is.na(g), ]
            )
            X.cov <- cbind(X.cov, rep(1, nrow(X.cov))) # append intercept
            # SynSurr with random forest
            fit.binormal <- SurrogateRegression::FitBNR(
                t = pheno$int[!is.na(g)],
                s = pheno$imputed_rf_1[!is.na(g)],
                X = X.cov
            )
            out <- fit.binormal@Regression.tab %>%
                filter(Outcome == "Target" & Coefficient == "g_complete") %>%
                select(Point, SE, p) %>%
                as.numeric()
            # SynSurr with linear regression
            fit.binormal <- SurrogateRegression::FitBNR(
                t = pheno$int[!is.na(g)],
                s = pheno$imputed_linear_1[!is.na(g)],
                X = X.cov
            )
            out <- c(out, fit.binormal @Regression.tab %>%
                filter(Outcome == "Target" & Coefficient == "g_complete") %>%
                select(Point, SE, p) %>% as.numeric())
            # SynSurr with permuted
            fit.binormal <- SurrogateRegression::FitBNR(
                t = pheno$int[!is.na(g)],
                s = pheno$imputed_rf_permute_1[!is.na(g)],
                X = X.cov
            )
            out <- c(out, fit.binormal @Regression.tab %>%
                filter(Outcome == "Target" & Coefficient == "g_complete") %>%
                select(Point, SE, p) %>% as.numeric())

            # SynSurr with negative random forest
            fit.binormal <- SurrogateRegression::FitBNR(
                t = pheno$int[!is.na(g)],
                s = -pheno$imputed_rf_1[!is.na(g)],
                X = X.cov
            )
            out <- c(out, fit.binormal @Regression.tab %>%
                filter(Outcome == "Target" & Coefficient == "g_complete") %>%
                select(Point, SE, p) %>% as.numeric())

            # SynSurr with mean
            fit.binormal <- SurrogateRegression::FitBNR(
                t = pheno$int[!is.na(g)],
                s = pheno$imputed_mean[!is.na(g)],
                X = X.cov
            )
            out <- c(out, fit.binormal @Regression.tab %>%
                filter(Outcome == "Target" & Coefficient == "g_complete") %>%
                select(Point, SE, p) %>% as.numeric())

            # add oracle
            assoc.oracle <- lm(oracle_int ~ ., data = data.frame(
                cbind(pheno %>% select(oracle_int, f.21022.0.0, f.22001.0.0, starts_with("PC")), g)
            ))
            out <- c(out, summary(assoc.oracle)$coefficients["g", c(1, 2, 4)])


            vas <- c("rf", "linear", "permute", "negate", "mean", "oracle")
            vis <- c("beta", "se", "p")
            out <- c(colnames(G)[i], out)
            names(out) <- c(
                "rsid",
                as.vector(t(outer(vas, vis, paste, sep = ".")))
            )
        },
        error = function(e) {
            NULL
        }
    )

    return(out)
})

stopCluster(cl)
results <- do.call(rbind, results)
results <- data.frame(results)

saveRDS(results, "Data/SynSurr_FEV1.rds")
