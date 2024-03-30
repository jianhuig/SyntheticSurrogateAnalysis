# Purpose: Number of significant associations and average chi2.
# Updated: 2023-07-23

library(dplyr)
library(ggplot2)


# -----------------------------------------------------------------------------
# Data preparation.
# -----------------------------------------------------------------------------

data <- data.table::fread(file = "results/ukbb/v2/compiled_results.tsv")
split_data <- split(data, data$trait)

df <- lapply(split_data, function(df) {
  
  sig_standard <- (df$standard_p <= 5e-8)
  sig_synsurr <- (df$synsurr_p <= 5e-8)
  sig_any <- (sig_standard | sig_synsurr)
  
  n_standard <- sum(sig_standard)
  n_synsurr <- sum(sig_synsurr)
  
  chi_standard <- qchisq(df$standard_p[sig_any], df = 1, lower.tail = FALSE)
  chi_synsurr <- qchisq(df$synsurr_p[sig_any], df = 1, lower.tail = FALSE)
  
  out <- data.frame(
    trait = unique(df$trait),
    n_standard = n_standard,
    n_synsurr = n_synsurr,
    mu_standard = mean(chi_standard),
    se_standard = sd(chi_standard) / sqrt(n_standard),
    mu_synsurr = mean(chi_synsurr),
    se_synsurr = sd(chi_synsurr) / sqrt(n_synsurr)
  )
  return(out)
})
df <- do.call(rbind, df)
rownames(df) <- NULL

df$trait <- gsub("_\\d+", "", x = df$trait)
df$trait <- factor(
  df$trait,
  levels = c("android", "arms", "gynoid", "legs", "total", "trunk"),
  labels = c("Android\nmass", "Arms\nmass", "Gynoid\nmass", "Legs\nmass", "Total\nmass", "Trunk\nmass")
)

df1 <- df %>%
  dplyr::select(trait, n_standard, n_synsurr) %>%
  tidyr::pivot_longer(
    cols = c("n_standard", "n_synsurr"),
    names_to = "method",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    method = gsub("n_", "", method)
  )
df1$method <- factor(
  df1$method,
  levels = c("standard", "synsurr"),
  labels = c("Standard", "SynSurr")
)

df2a <- df %>%
  dplyr::select(trait, mu_standard, mu_synsurr) %>%
  tidyr::pivot_longer(
    cols = c("mu_standard", "mu_synsurr"),
    names_to = "method",
    values_to = "chi2"
  ) %>%
  dplyr::mutate(
    method = gsub("mu_", "", method)
  )

df2b <- df %>%
  dplyr::select(trait, se_standard, se_synsurr) %>%
  tidyr::pivot_longer(
    cols = c("se_standard", "se_synsurr"),
    names_to = "method",
    values_to = "se"
  ) %>%
  dplyr::mutate(
    method = gsub("se_", "", method)
  )

z <- qnorm(0.975)
df2 <- merge(df2a, df2b, by = c("trait", "method"))
df2$lower <- df2$chi2 - z * df2$se
df2$upper <- df2$chi2 + z * df2$se

df2$method <- factor(
  df2$method,
  levels = c("standard", "synsurr"),
  labels = c("Standard", "SynSurr")
)

# -----------------------------------------------------------------------------

q1 <- ggplot(data = df1) +
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
    aes(x = trait, y = count, fill = method),
    position = position_dodge()
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Genome-wide\nSignificant Associations"
  ) +
  ggsci::scale_fill_nejm(
    name = "Method"
  ) + 
  geom_text(
    aes(x = trait, y = count, label = count, group = method),
    position = position_dodge(width = 0.9),
    vjust = -0.2
  )

legend <- ggpubr::get_legend(q1)
legend <- ggpubr::as_ggplot(legend)
q1 <- q1 + theme(legend.position = "none")

q2 <- ggplot(data = df2) +
  theme_bw() + 
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = 'black')
  ) + 
  geom_col(
    aes(x = trait, y = chi2, fill = method),
    position = position_dodge()
  ) +
  geom_errorbar(
    aes(x = trait, ymin = lower, ymax = upper, group = method),
    position = position_dodge(width = 0.9),
    width = 0.5
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = bquote(atop("Average"~chi^2~"at",
                       "Genome-wide Significant Variants"))
  ) +
  ggsci::scale_fill_nejm(
    name = "Method"
  )

q <- cowplot::plot_grid(
  plotlist = list(legend, q1, q2),
  labels = c("", "A", "B"),
  ncol = 1,
  rel_heights = c(2, 10, 10)
)

if (FALSE) {
  ggsave(
    plot = q,
    file = "results/ukbb/v2_body_comp_results.pdf",
    width = 12.0,
    height = 9.0,
    units = "in",
    dpi = 480
  )
}
