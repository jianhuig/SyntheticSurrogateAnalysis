# Purpose: Plot DEXA distribution predictions.

library(dplyr)
library(ggplot2)

files <- c(
  "results/surrogates/yhat_android_total_mass.rds",
  "results/surrogates/yhat_arms_total_mass.rds",
  "results/surrogates/yhat_gynoid_total_mass.rds",
  "results/surrogates/yhat_legs_total_mass.rds",
  "results/surrogates/yhat_total_mass.rds",
  "results/surrogates/yhat_trunk_total_mass.rds"
)

data <- lapply(files, function(fp) {
  df <- readRDS(fp)
  name <- basename(fp)
  name <- unlist(strsplit(name, split = "\\."))[1]
  name <- gsub("yhat_", "", x = name)
  colnames(df) <- c("y", "is_dexa_obs")
  df$trait <- name
  df <- df %>% dplyr::relocate(trait)
  return(df)
})
data <- do.call(rbind, data)

if (!is.factor(data$trait)) {
  data$trait <- factor(
    x = data$trait,
    levels = sort(unique(data$trait)),
    labels = c(
      "Android\nmass",
      "Arms\nmass",
      "Gynoid\nmass",
      "Legs\nmass",
      "Total\nmass",
      "Trunk\nmass"
    )
  )
}

if (!is.factor(data$is_dexa_obs)) {
  data$is_dexa_obs <- factor(
    x = data$is_dexa_obs,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
}


q <- ggplot(data = data) + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "top",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = 'black')
  ) + 
  geom_violin(
    aes(x = trait, y = y, fill = is_dexa_obs)
  ) + 
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) + 
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Predicted Mass [kg]",
    limits = c(0, 155)
  )

if (FALSE) {
  ggsave(
    plot = q,
    file = "results/surrogate_distributions_by_dexa.pdf",
    width = 9.0,
    height = 4.5
  )
}
