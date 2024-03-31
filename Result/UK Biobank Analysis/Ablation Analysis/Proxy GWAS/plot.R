library(data.table)

oracle <- fread("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Result/UKBB Imputation(1)/Height/plink2.oracle.glm.linear") %>%
  select(ID, BETA, SE, P) %>%
  filter(P < 5e-8)
proxy <- fread("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/plink2.height_proxy.glm.linear") %>% filter(ID %in% oracle$ID)


model <- lm(y ~ x, data = data.frame(x = oracle$BETA, 
                                     y = proxy$BETA))

# Construct the label for slope and R^2
slope <- round(coef(model)[2], 2)
r_squared <- round(summary(model)$r.squared, 2)
label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")
# Now use ggplot to make the plot
p1 <- ggplot(data.frame(x = oracle$BETA, y = proxy$BETA),
             aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  labs(x = "Oracle Effect Size", y = "Proxy GWAS Effect Size") +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  annotate("text", x = -0.02, y = 0.07, label = label_text, parse = TRUE,
           hjust = 1, vjust = 1, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_blank(),
        legend.key = element_blank())+
  xlim(-0.12, 0.12) +
  ylim(-0.12, 0.12)+
  ggtitle("Height")+ theme(plot.title = element_text(hjust = 0.5))

oracle <- fread("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Result/UKBB Imputation(1)/FEV1/plink2.oracle.glm.linear") %>%
  select(ID, BETA, SE, P) %>%
  filter(P < 5e-8)
proxy <- fread("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/plink2.fev1_proxy.glm.linear") %>% filter(ID %in% oracle$ID)


model <- lm(y ~ x, data = data.frame(x = oracle$BETA, 
                                     y = proxy$BETA))

# Construct the label for slope and R^2
slope <- round(coef(model)[2], 2)
r_squared <- round(summary(model)$r.squared, 2)
label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")
p2 <- ggplot(data.frame(x = oracle$BETA, y = proxy$BETA),
             aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  labs(x = "Oracle Effect Size", y = "Proxy GWAS Effect Size") +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  annotate("text", x = -0.02, y = 0.07, label = label_text, parse = TRUE,
           hjust = 1, vjust = 1, size = 3.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_blank(),
        legend.key = element_blank())+
  xlim(-0.12, 0.12) +
  ylim(-0.12, 0.12)+
  ggtitle("FEV1")+ theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1),
       filename = 'proxy_height_fev1_effect_size.pdf',
       width = 8, height = 4, dpi = 360,
       units = 'in')

