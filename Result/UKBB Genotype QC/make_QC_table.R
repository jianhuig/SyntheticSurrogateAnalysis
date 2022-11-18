# Set working directory.
setwd('~/Documents/GitHub/SyntheticSurrogateAnalysis/Result/UKBB Genotype QC')

# Load necessary packages.
library(xtable)

# Read in data.
qc_data <- read.csv('Data QC - Sheet1.csv', header = T)
head(qc_data)

# Clean data.
colnames(qc_data) <- c(" ",
                       "Initial Dataset",
                       "Variants with < 90% genotype rate",
                       "Variants left",
                       "Variants with HWE < 1e-5",
                       "Variants left",
                       "Variants with MAF < 1%", # Needs to be updated.
                       "Analysis Dataset"
                       )
dim(qc_data)

qc_data <- qc_data[, -which(colnames(qc_data) == "Variants left")]
dim(qc_data)

# Make table for the manuscript.
print(xtable(qc_data,
             align=c(
               "|p{2.5cm}|",
               "|p{1cm}|",
               rep("p{2.2cm}|", 5))
             ),
      include.rownames=FALSE)




