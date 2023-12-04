library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scattermore)
library(fastman)

setwd("/Users/jianhuigao/OneDrive - University of Toronto/EHR_research/SyntheticSurrogateAnalysis/Data/")

# read summary statistics from nealab
nealab <- fread("50_irnt.gwas.imputed_v3.both_sexes.tsv")
# filter out maf > 0,05
nealab <- nealab %>%
  filter(minor_allele > 0.05) %>%
  select(c("variant","minor_allele","beta", "pval"))
# seperate variant into chr, BP, allele1, allele2
split <- str_split_fixed(nealab$variant, ":", 4)
nealab$CHR <- as.numeric(split[, 1])
nealab$BP <- as.numeric(split[, 2])
nealab$a1 <- split[, 3]
nealab$a2 <- split[, 4]
nealab <- nealab %>%
  filter(a1 %in% c("A", "T", "C", "G")) %>% # exclude all In-Del variants
  filter(a2 %in% c("A", "T", "C", "G")) %>% # exclude all In-Del variants
  rename(P = pval) 

# If minor allele != alt2, then the sign of beta needs to be changed. As plink forces alt2 to be minor allele.
nealab <- nealab %>% 
  mutate(flip = ifelse(a2 == minor_allele, FALSE, TRUE))%>% 
  mutate(beta = ifelse(flip, -beta, beta)) %>% 
  select(c("CHR", "BP", "beta", "P")) %>%
  filter(P > 0) %>%
  drop_na()
save(nealab, file = "height_nlab.RData")

# Load Data
setwd("/Users/jianhuigao/OneDrive - University of Toronto/EHR_research/SyntheticSurrogateAnalysis/Data/")
gwas <- fread("plink2.PHENO1.glm.linear") %>% filter(TEST == "ADD") # only need genetic coef
gwas <- gwas %>% rename(CHR = `#CHROM`, BP = POS) %>% filter(P > 0)
gwas$P <- as.numeric(gwas$P)
load("height_nlab.Rdata")

# Manhattan plots
setwd("/Users/jianhuigao/OneDrive - University of Toronto/EHR_research/SyntheticSurrogateAnalysis/Result/")
png("height_manhattan.png", width = 10, height = 6, units = "in", res = 300)
fastman(gwas, suggest_line = FALSE) # Fast Manhattan plot
dev.off()

png("height_manhattan_nlab.png", width = 10, height = 6, units = "in", res = 300)
fastman(nealab, suggest_line = FALSE, main = "Nealab Result Manhattan plot") # Fast manhattan plot
dev.off()

# Plot p-value contrast
combined <- inner_join(gwas, nealab, by = c("CHR", "BP")) # only include snps in both results
png("pval_scatter.png")
ggplot(data = combined, aes(x = -log10(P.x), y = -log10(P.y))) +
  geom_scattermore(col = rev(heat.colors(nrow(combined), alpha=0.1))) + # fast scatter plot from scattermore
  geom_abline(intercept = 0, color = "red") + # add 45 degree line
  xlab("pval") +
  ylab("pval from nealab")
dev.off()

# Plot effect size contrast
png("beta_scatter.png")
ggplot(data = combined, aes(x = BETA, y = beta)) +
  geom_scattermore(col = rev(heat.colors(nrow(combined), alpha=0.1))) + # fast scatter plot from scattermore
  geom_abline(intercept = 0, color = "red", linetype =2) + # add 45 degree line
  xlab("Effect Size") +
  ylab("Effect Size from Nealab")
dev.off()
