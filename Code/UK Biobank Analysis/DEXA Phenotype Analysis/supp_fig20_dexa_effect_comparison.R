# Purpose: Compare effect sizes.

library(dplyr)
library(ggplot2)

data <- data.table::fread(file = "compiled_results.tsv")
data <- data %>% dplyr::filter(synsurr_p <= 5e-8)

data$trait <- factor(
  data$trait,
  levels = c(
    "android_23248",
    "arms_23260",
    "gynoid_23265",
    "legs_23277",
    "total_23283",
    "trunk_23287"
  ),
  labels = c(
    "Android Mass",
    "Arm Mass",
    "Gynoid Mass",
    "Leg Mass",
    "Total Mass",
    "Trunk Mass"
  )
)
split_data <- split(data, data$trait)

plot_list <- lapply(split_data, function(df) {
  
  fit <- lm(synsurr_beta ~ 0 + standard_beta, data = df)
  slope <- round(coef(fit), digits = 2)
  r2 <- round(cor(df$synsurr_beta, df$standard_beta)^2, digits = 2)
  
  # Create plot.
  q <- ggplot() + 
    theme_bw() + 
    theme(
      strip.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) + 
    theme(
      axis.line = element_line(color = 'black')
    ) +
    geom_abline(
      linetype = "dashed",
      intercept = 0,
      slope = 1
    ) + 
    geom_point(
      data = df,
      aes(x = standard_beta, y = synsurr_beta),
      alpha = 0.8
    ) + 
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.2,
      vjust = 1.5,
      label = bquote(atop(
        R^2 == .(r2),
        "Slope" == .(slope)
      ))
    ) + 
    stat_smooth(
      data = data,
      aes(x = standard_beta, y = synsurr_beta),
      method = "lm",
      formula = formula(y ~ 0 + x),
      color = "#0072B5FF",
      se = FALSE
    ) +
    scale_x_continuous(
      name = "Standard Effect Size"
    ) + 
    scale_y_continuous(
      name = "SynSurr Effect Size"
    ) + ggtitle(
      unique(df$trait)
    )
  return(q)
  
})


# Overall.
q <- cowplot::plot_grid(
  plotlist = plot_list,
  rows = 2
)

if (FALSE) {
  ggsave(
    plot = q,
    file = "results/effect_size_comparison.pdf",
    width = 12.0,
    height = 8.0
  )
}
