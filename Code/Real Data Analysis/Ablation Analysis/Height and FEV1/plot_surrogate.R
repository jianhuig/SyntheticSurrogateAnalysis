library(ggplot2)
library(dplyr)

field_id = 3063
phenotype = "FEV1"

data_path <- "/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Data/"
result_path <- paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/SyntheticSurrogateAnalysis/Result/", phenotype,"/")


for(i in 1:4){
  j = c(0.1, 0.25, 0.5, 0.75)[i]
  results <- readRDS(paste0(data_path,"binormal_field=",field_id,"_nmissing=",j,".rds"))
  p1 <- results %>% ggplot(aes(x=beta_g_bi, y = beta_g_obs)) + 
    scattermore::geom_scattermore() + 
    geom_abline(slope = 1, col = 'red') + 
    labs(x = expression(beta[bivariate]), y = expression(beta[observed])) + 
    ggtitle(paste0("Missing rate = ", j))
  p2 <- results %>% ggplot(aes(x=beta_g_bi, y = beta_g_surrogate)) + 
    scattermore::geom_scattermore() + 
    geom_abline(slope = 1, col = 'red') + 
    labs(x = expression(beta[bivariate]), y = expression(beta[predicted])) + 
    ggtitle(paste0("Missing rate = ",j))
  
  ggsave(file = paste0(result_path,"beta_comparision_missing",j,".png"),
         plot = ggpubr::ggarrange(p1,p2))
}

longdata <- data.frame()
for(i in 1:4){
  j = c(0.1, 0.25, 0.5, 0.75)[i]
  results <- readRDS(paste0(data_path,"binormal_field=",field_id,"_nmissing=",j,".rds"))
  longdata <- rbind(longdata,cbind((results$se_obs/results$se_bi)^2,rep("bivariate/observed", nrow(results)), rep(j, nrow(results))))
  longdata <- rbind(longdata,cbind((results$se_surrogate/results$se_bi)^2,rep("bivariate/surrogate", nrow(results)), rep(j, nrow(results))))
}
longdata <- data.frame(longdata)
longdata$V1 <- as.numeric(longdata$V1)
longdata %>% ggplot() + geom_boxplot(aes(y= V1, color = V2, x = factor(V3))) + labs(y = "RE", x = "Percent of missing", color = "Ratio") + ylim(c(0,1.2))
ggsave(filename = paste0(result_path, "relative_efficiency.png"))
