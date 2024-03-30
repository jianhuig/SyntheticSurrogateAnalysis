# Purpose: Intersect SynSurr results with GWAS catalog.
# Updated: 2023-07-25
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------

data <- data.table::fread("results/ukbb/v2/ldsc_overall_results.tsv")
data$pheno <- factor(
  data$pheno,
  levels = c(23248, 23260, 23265, 23277, 23283, 23287),
  labels = c("Android\nmass", "Arms\nmass", "Gynoid\nmass", "Legs\nmass", "Total\nmass", "Trunk\nmass")
)
data$method <- factor(
  data$method,
  levels = c("standard", "synsurr"),
  labels = c("Standard", "SynSurr")
)
z <- stats::qnorm(0.975)
data$lower <- data$rg - z * data$se
data$upper <- data$rg + z * data$se

q <- ggplot(data = data) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = 'black')
  ) + 
  geom_col(
    aes(x = pheno, y = rg, fill = method),
    position = position_dodge()
  ) + 
  geom_errorbar(
    aes(x = pheno, ymin = lower, ymax = upper, group = method),
    position = position_dodge(width = 0.9)
  ) +
  ggsci::scale_fill_nejm(
    name = "Method"
  ) + 
  scale_x_discrete(
    name = "Trait"
  ) +
  scale_y_continuous(
    name = "Genetic Correlation"
  ) + 
  geom_hline(
    linetype = "dashed",
    color = "gray",
    yintercept = 1.0
  )

if (FALSE) {
  ggsave(
    plot = q,
    file = "results/genetic_correlations.pdf",
    width = 9.0,
    height = 6.0
  )
}

