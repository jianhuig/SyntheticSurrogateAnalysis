library(ranger)
library(dplyr)
library(ggplot2)
# Random Forest Plot ------------
setwd("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data")
## Height -------------------------
field_id <- 50
cov_id <- c("f.21022.0.0","f.22001.0.0","f.48.0.0","f.21002.0.0","f.49.0.0") #age,sex,height, weight, BMI
pheno_id <- "f.50.0.0"

dat <- readRDS(paste0("field_", field_id, "_cleaned.rds"))
in_id <- data.table::fread("final.king.cutoff.in.id")$IID
out_id <- data.table::fread("final.king.cutoff.out.id")$IID
weight <- read.table("Weight.txt", fill = TRUE, header = TRUE)
white <- data.table::fread("Ancestry.tab")

train <- dat %>% left_join(weight %>% select(f.eid, "f.21002.0.0"), by = "f.eid") %>% filter(f.eid %in% out_id) %>% inner_join(white) %>% filter(!is.na(f.22006.0.0))%>% tidyr::drop_na(pheno_id) %>% tidyr::drop_na(all_of(cov_id)) %>% select(pheno_id, cov_id) %>% filter(f.50.0.0 > 10)
test <- dat %>% left_join(weight %>% select(f.eid, "f.21002.0.0"), by = "f.eid")  %>% inner_join(white) %>% filter(!is.na(f.22006.0.0)) %>% filter(f.eid %in% in_id) %>% tidyr::drop_na(all_of(cov_id)) %>% select(f.eid,pheno_id, cov_id)  %>% filter(f.50.0.0 > 10)


model.rf <- ranger::ranger(
  data = train,
  as.formula(paste0(pheno_id, "~."))
)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
yhat <- predict(model.rf, data = train %>% select(-pheno_id))$predictions
train$yhat <- yhat
p1 <- train %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed Height (cm)") + ylab("Predicted Height (cm)") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 180, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA) + ggtitle("Model-building dataset")+xlim(c(120,225))+ylim(c(120,225))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        legend.key = element_blank())

yhat <- predict(model.rf, data = test %>% select(-f.eid,-pheno_id))$predictions
test$yhat <- yhat
p2 <- test %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed Height (cm)") + ylab("Predicted Height (cm)") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 180, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA)+ ggtitle("GWAS dataset")+xlim(c(120,225))+ylim(c(120,225)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        legend.key = element_blank())

ggsave("rf_pred_height_heat.pdf", plot = ggpubr::ggarrange(p1, p2), device = "pdf", width = 8, height = 4)

## FEV1 --------------------------------------------------
field_id <- 20150
cov_id <- c("f.21022.0.0","f.22001.0.0","f.48.0.0","f.21002.0.0","f.49.0.0") #age,sex,height, weight, BMI
pheno_id <- "f.20150.0.0"

dat <- readRDS(paste0("field_", field_id, "_cleaned.rds"))
#in_id <- read.table("final.king.cutoff.in.id", header = TRUE) %>% pull(IID)
#out_id <- read.table("final.king.cutoff.out.id", header = TRUE) %>% pull(IID)
weight <- read.table("Weight.txt", fill = TRUE, header = TRUE)
white <- data.table::fread("Ancestry.tab")

train <- dat %>% filter(f.eid %in% out_id) %>% inner_join(white) %>% filter(!is.na(f.22006.0.0))%>% tidyr::drop_na(pheno_id) %>% tidyr::drop_na(all_of(cov_id)) %>% select(pheno_id, cov_id)
test <- dat %>% inner_join(white) %>% filter(!is.na(f.22006.0.0)) %>% filter(f.eid %in% in_id) %>% tidyr::drop_na(all_of(cov_id)) %>% select(f.eid,pheno_id, cov_id)%>% tidyr::drop_na(pheno_id)


model.rf <- ranger::ranger(
  data = train,
  as.formula(paste0(pheno_id, "~."))
)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
yhat <- predict(model.rf, data = train %>% select(-pheno_id))$predictions
train$yhat <- yhat
p1 <- train %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed FEV1") + ylab("Predicted FEV1") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 8, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA) + ggtitle("Model-building dataset")+xlim(c(0,15))+ylim(c(0,15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        legend.key = element_blank())

yhat <- predict(model.rf, data = test %>% select(-f.eid,-pheno_id))$predictions
test$yhat <- yhat
p2 <- test %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed FEV1") + ylab("Predicted FEV1") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 8, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA) + ggtitle("GWAS dataset")+xlim(c(0,15))+ylim(c(0,15)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave("rf_pred_FEV1_heat.pdf", plot = ggpubr::ggarrange(p1, p2), device = "pdf", width = 8, height = 4)

# Effect Size ------------------------------------------------------------------
## Height -------------------------
id <- 50

result <- readRDS(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", id, "_nmissing=", 0, ".rds"))

true_snp_index <- which(result$oracle.p < 5e-8)
betas <- c()
for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  result <- readRDS(paste0("~/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", id, "_nmissing=", i, ".rds"))
  result <- result[true_snp_index, ]
  betas <- rbind(betas, cbind(i, result$oracle.beta, result$marginal.beta, result$bi.beta))
}
colnames(betas) <- c("missing", "oracle", "marginal", "bivariate")

missing_names <- c(
  `0` = "Missing Rate = 0%",
  `0.25` = "Missing Rate = 25%",
  `0.5` = "Missing Rate = 50%",
  `0.75` = "Missing Rate = 75%",
  `0.9` = "Missing Rate = 90%"
)

betas <- betas %>% data.frame()
plist <- list()
for(i in c(0, 0.25, 0.5, 0.75, 0.9)){
  model <- lm(betas[betas$missing==i, "marginal"] ~ betas[betas$missing==i, "oracle"])
  slope <- round(coef(model)[2], 2)
  r_squared <- round(summary(model)$r.squared, 2)
  label_text <- label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")

  plist[[length(plist) + 1]] <- betas %>%
    filter(missing == i) %>%
    ggplot(aes(x = oracle, y = marginal)) +
    geom_point() +
    xlim(c(-0.15, 0.15)) +
    ylim(c(-0.15, 0.15)) +
    geom_smooth(method = "lm") +
    annotate("text", x = -0.02, y = 0.1, label = label_text, parse = TRUE,
             hjust = 1, vjust = 1, size = 3.5)+
    xlab(expression(beta[oracle])) +
    ylab(expression(beta[standard]))  +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "white", colour = NA),
          legend.background = element_blank(),
          legend.key = element_blank())+
    ggtitle(paste0("Missing Rate = ", i*100, "%"))+
    theme(plot.title = element_text(hjust = 0.5))
  }

ggsave(
  plot = ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 2),
  filename = "effectsize_height_marginal.pdf",
  device = "pdf",
  height = 6,
  width = 7.5,
  units = "in"
)

betas <- betas %>% data.frame()
plist <- list()
for(i in c(0, 0.25, 0.5, 0.75, 0.9)){
  model <- lm(betas[betas$missing==i, "bivariate"] ~ betas[betas$missing==i, "oracle"])
  slope <- round(coef(model)[2], 2)
  r_squared <- round(summary(model)$r.squared, 2)
  label_text <- label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")
  
  plist[[length(plist) + 1]] <- betas %>%
    filter(missing == i) %>%
    ggplot(aes(x = oracle, y = bivariate)) +
    geom_point() +
    xlim(c(-0.15, 0.15)) +
    ylim(c(-0.15, 0.15)) +
    geom_smooth(method = "lm") +
    annotate("text", x = -0.02, y = 0.1, label = label_text, parse = TRUE,
             hjust = 1, vjust = 1, size = 3.5)+
    xlab(expression(beta[oracle])) +
    ylab(expression(beta[SynSurr]))  +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "white", colour = NA),
          legend.background = element_blank(),
          legend.key = element_blank())+
    ggtitle(paste0("Missing Rate = ", i*100, "%"))+
    theme(plot.title = element_text(hjust = 0.5))
}
ggsave(
  plot = ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 2),
  filename = "effectsize_height.pdf",
  device = "pdf",
  height = 6,
  width = 7.5,
  units = "in"
)


## FEV1 -------------------------
id <- 20150

result <- readRDS(paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", id, "_nmissing=", 0, ".rds"))

true_snp_index <- which(result$oracle.p < 5e-8)
betas <- c()
for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  result <- readRDS(paste0("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", id, "_nmissing=", i, ".rds"))
  result <- result[true_snp_index, ]
  betas <- rbind(betas, cbind(i, result$oracle.beta, result$marginal.beta, result$bi.beta))
}
colnames(betas) <- c("missing", "oracle", "marginal", "bivariate")

betas <- betas %>% data.frame()
plist <- list()
for(i in c(0, 0.25, 0.5, 0.75, 0.9)){
  model <- lm(betas[betas$missing==i, "marginal"] ~ betas[betas$missing==i, "oracle"])
  slope <- round(coef(model)[2], 2)
  r_squared <- round(summary(model)$r.squared, 2)
  label_text <- label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")
  
  plist[[length(plist) + 1]] <- betas %>%
    filter(missing == i) %>%
    ggplot(aes(x = oracle, y = marginal)) +
    geom_point() +
    xlim(c(-0.15, 0.15)) +
    ylim(c(-0.15, 0.15)) +
    geom_smooth(method = "lm") +
    annotate("text", x = -0.03, y = 0.1, label = label_text, parse = TRUE,
             hjust = 1, vjust = 1, size = 3.5)+
    xlab(expression(beta[oracle])) +
    ylab(expression(beta[standard]))  +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "white", colour = NA),
          legend.background = element_blank(),
          legend.key = element_blank())+
    ggtitle(paste0("Missing Rate = ", i*100, "%"))+
    theme(plot.title = element_text(hjust = 0.5))
}
ggsave(
  plot = ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 2),
  filename = "effectsize_FEV1_marginal.pdf",
  device = "pdf",
  height = 6,
  width = 7.5,
  units = "in"
)

betas <- betas %>% data.frame()
plist <- list()
for(i in c(0, 0.25, 0.5, 0.75, 0.9)){
  model <- lm(betas[betas$missing==i, "bivariate"] ~ betas[betas$missing==i, "oracle"])
  slope <- round(coef(model)[2], 2)
  r_squared <- round(summary(model)$r.squared, 2)
  label_text <- paste0("atop(italic(R)^2 == ", r_squared, ", ~slope == ", slope, ")")
  
  plist[[length(plist) + 1]] <- betas %>%
    filter(missing == i) %>%
    ggplot(aes(x = oracle, y = bivariate)) +
    geom_point() +
    xlim(c(-0.15, 0.15)) +
    ylim(c(-0.15, 0.15)) +
    geom_smooth(method = "lm") +
    annotate("text", x = -0.02, y = 0.1, label = label_text, parse = TRUE,
             hjust = 1, vjust = 1, size = 3.5)+
    xlab(expression(beta[oracle])) +
    ylab(expression(beta[SynSurr]))  +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "white", colour = NA),
          legend.background = element_blank(),
          legend.key = element_blank())+
    ggtitle(paste0("Missing Rate = ", i*100, "%"))+
    theme(plot.title = element_text(hjust = 0.5))
}
ggsave(
  plot = ggpubr::ggarrange(plotlist = plist, ncol = 3, nrow = 2),
  filename = "effectsize_FEV1.pdf",
  device = "pdf",
  height = 6,
  width = 7.5,
  units = "in"
)


# Tables ---------------------------------------------------------------------
## FPR --------------------------------------------
metric <- c()
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", 0, ".rds"))
true_snp_height <- which(result$oracle.p<5e-8)
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", 0, ".rds"))
true_snp_fev1 <- which(result$oracle.p<5e-8)

for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  height <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", i, ".rds"))
  fev1 <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", i, ".rds"))
  metric <- rbind(metric, c(
    i,
    paste0(sum(height$marginal.p[-true_snp_height]<5e-8),"(",round(sum(height$marginal.p[-true_snp_height]<5e-8)/nrow(height)*100,3), "%)"),
    paste0(sum(height$bi.p[-true_snp_height]<5e-8),"(", round(sum(height$bi.p[-true_snp_height]<5e-8)/nrow(height)*100,3), "%)"),
    paste0(sum(fev1$marginal.p[-true_snp_fev1]<5e-8),"(", round(sum(fev1$marginal.p[-true_snp_fev1]<5e-8)/nrow(fev1)*100,3), "%)"),
    paste0(sum(fev1$marginal.p[-true_snp_fev1]<5e-8),"(", round(sum(fev1$marginal.p[-true_snp_fev1]<5e-8)/nrow(fev1)*100,3), "%)")
  ))
}
colnames(metric) <- c("Missing Rate(%)", "Standard", "SynSurr", "Standard","SynSurr")

metric %>%
  kable(booktab = TRUE, align = 'c', format = 'latex') %>%
  add_header_above(c(" " = 1, "Height" = 2, "FEV1" = 2))

## MSE  -------------------------------------------------------------
metric <- c()
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", 0, ".rds"))
true_snp_height <- which(result$oracle.p<5e-8)
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", 0, ".rds"))
true_snp_fev1 <- which(result$oracle.p<5e-8)

for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  height <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", i, ".rds"))
  fev1 <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", i, ".rds"))
  metric <- rbind(metric, c(
    i,
    signif(mean(height$oracle.se[true_snp_height]^2),2),
    signif(mean(height$marginal.se[true_snp_height]^2),2),
    signif(mean(height$bi.se[true_snp_height]^2),2),
    round(1/mean(height$bi.se[true_snp_height]^2/height$marginal.se[true_snp_height]^2),2),
    signif(mean(fev1$oracle.se[true_snp_fev1]^2),2),
    signif(mean(fev1$marginal.se[true_snp_fev1]^2),2),
    signif(mean(fev1$bi.se[true_snp_fev1]^2),2),
    round(1/mean(fev1$bi.se[true_snp_fev1]^2/fev1$marginal.se[true_snp_fev1]^2),2)
  ))
}
colnames(metric) <- c("Missing Rate(%)","Oracle","Standard", "SynSurr","Relative Efficiency","Oracle","Standard","SynSurr", "Relative Efficiency")

metric %>%
  kable(booktab = TRUE, format = "latex",align = 'c') %>%
  add_header_above(c(" " = 1, "MSE" = 3, " "=1,"MSE" = 3, " "=1))%>%
  add_header_above(c(" " = 1, "Height" = 4, "FEV1"=4))

## Chisq  -------------------------------------------------------------
metric <- c()
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", 0, ".rds"))
true_snp_height <- which(result$oracle.p<5e-8)
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", 0, ".rds"))
true_snp_fev1 <- which(result$oracle.p<5e-8)

for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  height <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", i, ".rds"))
  fev1 <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", i, ".rds"))
  metric <- rbind(metric, c(
    i,
    round(mean((height$oracle.beta[true_snp_height]/height$oracle.se[true_snp_height])^2),2),
    round(mean((height$marginal.beta[true_snp_height]/height$marginal.se[true_snp_height])^2),2),
    round(mean((height$bi.beta[true_snp_height]/height$bi.se[true_snp_height])^2),2),
    round(mean(((height$bi.beta[true_snp_height]/height$bi.se[true_snp_height])^2))/mean((height$marginal.beta[true_snp_height]/height$marginal.se[true_snp_height])^2),2),
    round(mean((fev1$oracle.beta[true_snp_fev1]/fev1$oracle.se[true_snp_fev1])^2),2),
    round(mean((fev1$marginal.beta[true_snp_fev1]/fev1$marginal.se[true_snp_fev1])^2),2),
    round(mean((fev1$bi.beta[true_snp_fev1]/fev1$bi.se[true_snp_fev1])^2),2),
    round(mean(((fev1$bi.beta[true_snp_fev1]/fev1$bi.se[true_snp_fev1])^2))/mean((fev1$marginal.beta[true_snp_fev1]/fev1$marginal.se[true_snp_fev1])^2),2)
  ))
}
colnames(metric) <- c("Missing Rate(%)","Oracle","Standard", "SynSurr","Relative Efficiency","Oracle","Standard","SynSurr", "Relative Efficiency")

metric %>%
  kable(booktab = TRUE, format = "latex",align = 'c') %>%
  add_header_above(c(" " = 1, "MSE" = 3, " "=1,"MSE" = 3, " "=1))%>%
  add_header_above(c(" " = 1, "Height" = 4, "FEV1"=4))

# Number of Significant SNPs

metric <- c()
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", 0, ".rds"))
true_snp_height <- which(result$oracle.p<5e-8)
result <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", 0, ".rds"))
true_snp_fev1 <- which(result$oracle.p<5e-8)

for (i in c(0, 0.25, 0.5, 0.75, 0.9)) {
  height <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 50, "_nmissing=", i, ".rds"))
  fev1 <- readRDS(paste0("~/Desktop/SyntheticSurrogateAnalysis/Data/binormal_field=", 20150, "_nmissing=", i, ".rds"))
  metric <- rbind(metric, c(
    i,
    nrow(height),
    sum(height$oracle.p<5e-8),
    sum(height$marginal.p<5e-8),
    sum()
    round(mean((height$bi.beta[true_snp_height]/height$bi.se[true_snp_height])^2),2),
    round(mean(((height$bi.beta[true_snp_height]/height$bi.se[true_snp_height])^2))/mean((height$marginal.beta[true_snp_height]/height$marginal.se[true_snp_height])^2),2),
    round(mean((fev1$oracle.beta[true_snp_fev1]/fev1$oracle.se[true_snp_fev1])^2),2),
    round(mean((fev1$marginal.beta[true_snp_fev1]/fev1$marginal.se[true_snp_fev1])^2),2),
    round(mean((fev1$bi.beta[true_snp_fev1]/fev1$bi.se[true_snp_fev1])^2),2),
    round(mean(((fev1$bi.beta[true_snp_fev1]/fev1$bi.se[true_snp_fev1])^2))/mean((fev1$marginal.beta[true_snp_fev1]/fev1$marginal.se[true_snp_fev1])^2),2)
  ))
}
colnames(metric) <- c("Missing Rate(%)","Oracle","Standard", "SynSurr","Relative Efficiency","Oracle","Standard","SynSurr", "Relative Efficiency")

metric %>%
  kable(booktab = TRUE, format = "latex",align = 'c') %>%
  add_header_above(c(" " = 1, "MSE" = 3, " "=1,"MSE" = 3, " "=1))%>%
  add_header_above(c(" " = 1, "Height" = 4, "FEV1"=4))
