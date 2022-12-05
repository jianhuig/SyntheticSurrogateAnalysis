library(ranger)

setwd("/Users/jianhuigao/Desktop/SyntheticSurrogateAnalysis/Data/")
field_id <- 50
cov_id <- c("f.21022.0.0","f.22001.0.0","f.48.0.0","f.21002.0.0","f.49.0.0") #age,sex,height, weight, BMI
pheno_id <- "f.50.0.0"

dat <- readRDS(paste0("field_", field_id, "_cleaned.rds"))
in_id <- read.table("final.king.cutoff.in.id", header = TRUE) %>% pull(IID)
out_id <- read.table("final.king.cutoff.out.id", header = TRUE) %>% pull(IID)
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
p1 <- train %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed Height (cm)") + ylab("Predicted Height (cm)") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 180, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA)+ theme_bw() + ggtitle("Model-building dataset")+xlim(c(120,200))+ylim(c(120,220))

yhat <- predict(model.rf, data = test %>% select(-f.eid,-pheno_id))$predictions
test$yhat <- yhat
p2 <- test %>% rename(pheno = pheno_id) %>% mutate(density = get_density(pheno, yhat)) %>% ggplot(aes(x = pheno, y = yhat, color = density)) + geom_point(alpha = 0.2) + xlab("Observed Height (cm)") + ylab("Predicted Height (cm)") +geom_smooth(method = "lm")+ ggpubr::stat_regline_equation(label.y = 180, aes(label = ..rr.label..))+scale_color_gradient(low = "red", high = "yellow", na.value = NA)+ theme_bw() + ggtitle("Model-building dataset")+xlim(c(120,200))+ylim(c(120,220))
ggpubr::ggarrange(p1, p2)