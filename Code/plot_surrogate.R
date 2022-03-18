setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/EHR_research/PHD_thesis/SyntheticSurrogateAnalysis/Data")

library(ggplot2)
library(dplyr)

result <- list()
for(i in 1:4){
  j = c(0.1, 0.25, 0.5, 0.75)[i]
  load(paste0("binormal_results",j,".RData"))
  p1 <- results %>% ggplot(aes(x=beta_g_bi, y = beta_g_obs)) + scattermore::geom_scattermore() + geom_abline(slope = 1, col = 'red') + xlim(-0.1, 0.1) + ylim(-0.1, 0.1) + labs(x = expression(beta[bivariate]), y = expression(beta[observed])) + ggtitle(paste0("Missing rate = ", j))
  p2 <- results %>% ggplot(aes(x=beta_g_bi, y = beta_g_surrogate)) + scattermore::geom_scattermore() + geom_abline(slope = 1, col = 'red') + xlim(-0.1, 0.1) + ylim(-0.1, 0.1)+ labs(x = expression(beta[bivariate]), y = expression(beta[predicted])) + ggtitle(paste0("Missing rate = ",j))
  
  png(paste0("beta_comparision_missing",j,".png"))
  ggpubr::ggarrange(p1,p2)
  dev.off()
}

longdata <- data.frame()
for(i in 1:4){
  j = c(0.1, 0.25, 0.5, 0.75)[i]
  load(paste0("binormal_results",j,".RData"))
  longdata <- rbind(longdata,cbind((results$se_obs/results$se_bi)^2,rep("bivariate/observed", nrow(results)), rep(j, nrow(results))))
  longdata <- rbind(longdata,cbind((results$se_surrogate/results$se_bi)^2,rep("bivariate/surrogate", nrow(results)), rep(j, nrow(results))))
}
longdata <- data.frame(longdata)
longdata$V1 <- as.numeric(longdata$V1)
longdata %>% ggplot() + geom_boxplot(aes(y= V1, color = V2, x = factor(V3))) + labs(y = "RE", x = "Percent of missing", color = "Ratio")
