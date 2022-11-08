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
                       "Variants with MAF < 5%", # Needs to be updated.
                       "Variants left",
                       "Variants with LD Pruning > 0.9",
                       "Variants left"
                       )

# Make table for the manuscript.
print(xtable(qc_data,
             align=c(
               "|p{2.5cm}|",
               "|p{1cm}|",
               rep("p{1.4cm}|", 9))
             ),
      include.rownames=FALSE)




