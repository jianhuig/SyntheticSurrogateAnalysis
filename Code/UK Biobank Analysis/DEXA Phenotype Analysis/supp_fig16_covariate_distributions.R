# Purpose: Plot covariate distributions of patients with and without DEXA.
# Updated: 2023-11-20

library(dplyr)
library(ggplot2)

covar <- data.table::fread(file = "data/simple_covar.tsv")
dexa <- data.table::fread(file = "data/dexa.tsv")
covar <- covar %>%
  dplyr::inner_join(dexa, by = "eid")
covar$is_dexa_obs <- !is.na(covar$total_arm_mass_23260)

age <- data.table::fread(file = "data/age.tsv")
covar <- covar %>%
  dplyr::inner_join(age, by = "eid")

height <- data.table::fread(file = "data/height_50.tsv")
weight <- data.table::fread(file = "data/weight_21002.tsv")
bmi <- data.table::fread(file = "data/bmi_21001.tsv")

covar <- covar %>%
  dplyr::left_join(height, by = "eid") %>%
  dplyr::left_join(weight, by = "eid") %>%
  dplyr::left_join(bmi, by = "eid")

# ------------------------------------------------------------------------------

impedance <- data.table::fread(file = "data/impedance.tsv")
covar <- covar %>%
  dplyr::left_join(impedance, by = "eid")

# ------------------------------------------------------------------------------
# Genetic sex
# ------------------------------------------------------------------------------

z <- qnorm(0.975)
df <- covar %>%
  dplyr::group_by(is_dexa_obs) %>%
  dplyr::summarise(
    prop_female = mean(1 - gsex)
  )

q <- ggplot() + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = "black")
  ) +
  geom_hline(
    yintercept = 0.50,
    linetype = "dashed",
    color = "gray"
  ) + 
  geom_col(
    data = df,
    aes(x = is_dexa_obs, y = prop_female, fill = is_dexa_obs),
  ) + 
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) + 
  scale_x_discrete(
    name = "DEXA Measured",
    labels = c("No", "Yes")
  ) + 
  scale_y_continuous(
    name = "Proportion Female",
    limits = c(0, 1.0)
  ) 
  # ggtitle("Genetic Sex by DEXA Measurement")

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/covariates/gsex_by_dexa.pdf",
    width = 6.0,
    height = 3.0
  )
}

# ------------------------------------------------------------------------------
# Age.
# ------------------------------------------------------------------------------

q <- ggplot(data = covar) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = "black")
  ) + 
  geom_violin(
    aes(x = is_dexa_obs, y = age_21022, fill = is_dexa_obs),
    alpha = 0.8
  ) +
  geom_boxplot(
    aes(x = is_dexa_obs, y = age_21022, fill = is_dexa_obs),
    width = 0.1,
  ) + 
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) + 
  scale_x_discrete(
    name = "DEXA Measured",
    labels = c("No", "Yes")
  ) + 
  scale_y_continuous(
    name = "Age",
    limits = c(38, 72)
  )
  # ggtitle("Age by DEXA Measurement")

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/covariates/age_by_dexa.pdf",
    width = 6.0,
    height = 3.0
  )
}


# ------------------------------------------------------------------------------
# Body measures.
# ------------------------------------------------------------------------------

df <- covar %>%
  dplyr::select(is_dexa_obs, height, weight, bmi) %>%
  dplyr::filter(!is.na(bmi)) %>%
  tidyr::pivot_longer(
    height:bmi,
    names_to = "pheno",
    values_to = "y"
  )

if (!is.factor(df$pheno)) {
  df$pheno <- factor(
    x = df$pheno,
    levels = c("height", "weight", "bmi"),
    labels = c("Height", "Weight", "BMI")
  )
}

q <- ggplot(data = df) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = "black")
  ) +
  geom_violin(
    aes(x = pheno, y = y, fill = is_dexa_obs),
    position = position_dodge(width = 0.9),
    alpha = 0.8
  ) +
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) +
  scale_x_discrete(
    name = NULL,
  ) + 
  scale_y_continuous(
    name = "Phenotypic Value",
  )
# ggtitle("Age by DEXA Measurement")

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/covariates/height_weight_by_dexa.pdf",
    width = 9.0,
    height = 4.5
  )
}


# Impedance.
df <- covar %>%
  dplyr::select(
    is_dexa_obs, 
    left_arm,
    left_leg,
    right_arm,
    right_leg,
    whole_body
  ) %>%
  dplyr::filter(!is.na(whole_body)) %>%
  tidyr::pivot_longer(
    left_arm:whole_body,
    names_to = "pheno",
    values_to = "y"
  ) %>%
  dplyr::filter(!is.na(y))

if (!is.factor(df$pheno)) {
  df$pheno <- factor(
    x = df$pheno,
    levels = c("left_arm", "left_leg", "right_arm", "right_leg", "whole_body"),
    labels = c("Left\nArm", "Left\nLeg", "Right\nArm", "Right\nLeg", "Whole\nBody")
  )
}

q <- ggplot(data = df) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = "black")
  ) + 
  geom_violin(
    aes(x = pheno, y = y, fill = is_dexa_obs),
    position = position_dodge(width = 0.9),
    alpha = 0.8
  ) +
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) +
  scale_x_discrete(
    name = NULL,
  ) + 
  scale_y_continuous(
    name = expression(Impedance~(Omega)),
  )
# ggtitle("Age by DEXA Measurement")

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/covariates/impedance_by_dexa.pdf",
    width = 9.0,
    height = 4.5
  )
}


# ------------------------------------------------------------------------------
# Principal Components
# ------------------------------------------------------------------------------

df <- covar %>%
  dplyr::select(is_dexa_obs, gpc1:gpc20) %>%
  tidyr::pivot_longer(
    gpc1:gpc20,
    names_to = "which_pc",
    values_to = "y"
  )

if (!is.factor(df$which_pc)) {
  df$which_pc <- factor(
    df$which_pc,
    levels = paste0("gpc", seq(1:20)),
    ordered = TRUE
  )
}

if (!is.factor(df$is_dexa_obs)) {
  df$is_dexa_obs <- factor(
    df$is_dexa_obs,
    levels = c(FALSE, TRUE),
    label = c("No", "Yes")
  )
}

q <- ggplot(data = df) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "top",
    panel.border = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = "black")
  ) + 
  geom_violin(
    aes(x = which_pc, y = y, fill = is_dexa_obs),
    position = position_dodge(width = 0.9),
    alpha = 0.8
  ) +
  ggsci::scale_fill_nejm(
    name = "DEXA Measured?"
  ) +
  scale_x_discrete(
    name = "Genetic PC",
    labels = seq(1:20)
  ) + 
  scale_y_continuous(
    name = "PC Value",
    limits = c(-14, 14)
  )  
  # ggtitle("Genetic PCs by DEXA Measurement")

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/covariates/pcs_by_dexa.pdf",
    width = 9.0,
    height = 4.5
  )
}

