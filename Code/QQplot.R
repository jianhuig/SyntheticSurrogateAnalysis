library(dplyr)

file <- list.files("/Users/jianhuigao/Desktop/SyntheticSurrogateAnalysis/Data",
                   pattern = "*.permuted", full.names = T)

df <- data.frame()
pheno <- c("Android total mass", "Arms total mass", "Gynoid total mass",
           "Total mass", "Legs total mass", "Trunk total mass")

for(i in 1:length(file)){
  bi.p <- readRDS(file[i]) %>% pull(bi.p) %>% as.numeric()
  df <- rbind(df, data.frame(pval = bi.p, trait = pheno[i]))
}

p <- df %>% 
  group_by(trait) %>%
  mutate(expected = -log10(ppoints(n())), observed = -log10(sort(pval)),
         clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n(), shape2 = n():1)),
         cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n(), shape2 = n():1))) %>%
  ggplot()+
  geom_point(aes(expected, observed, color = factor(trait)))+
  xlab(expression(paste("Expected -log"[10], plain(P))))+
  ylab(expression(paste("Observed -log"[10], plain(P)))) +
  guides(color = guide_legend(title = "Phenotypes")) +
  geom_abline(intercept = 0, slope = 1, color = "red")+
  geom_line(aes(expected, cupper), linetype = 2) +
  geom_line(aes(expected, clower), linetype = 2) +
  theme_bw()+
  theme(legend.position = "bottom")

ggsave(
  plot = p,
  filename = "qq_ukb_permuted.png",
  device = "png",
  height = 6,
  width = 7.5,
  units = "in",
  dpi = 360
)
