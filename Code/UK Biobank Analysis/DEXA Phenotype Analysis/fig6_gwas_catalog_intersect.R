# Purpose: Intersect SynSurr results with GWAS catalog.
# Updated: 2022-11-11
BASE_DIR <- NULL

# devtools::install_github(repo = "zrmacc/CompareLoci")
library(CompareLoci)
library(dplyr)
library(ggplot2)
FLANK <- 250e3

# -----------------------------------------------------------------------------

#' Load Results
#' @param file File name.
#' @return Data.frame
LoadResults <- function(file) {
  # df <- data.table::fread(
  #   file,
  #   select = c("rsid", "CHR", "BP", "marginal.p", "bi.p"),
  #   col.names = c("rsid", "chr", "bp", "standard_p", "synsurr_p")
  # )
  df <- data.table::fread(
    file,
    select = c(
      "SNP", "CHR", "POS",
      "BETA_Standard", "SE_Standard", "P_Standard",
      "BETA_SynSurr", "SE_SynSurr", "P_SynSurr"),
    col.names = c(
      "rsid", "chr", "bp", 
      "standard_beta", "standard_se", "standard_p", 
      "synsurr_beta", "synsurr_se", "synsurr_p"
    )
  )
  return(df)  
}

# -----------------------------------------------------------------------------
# Prepare data.
# -----------------------------------------------------------------------------

if (FALSE) {
  
  setwd(file.path(BASE_DIR, "results/ukbb/v2"))
  files <- dir()
  files <- files[grepl("^[a-z]+_[0-9]+\\.tsv", files)]
  
  data <- lapply(files, function(fp) {
    df <- LoadResults(fp)
    trait <- unlist(strsplit(fp, split = ".", fixed = TRUE))[1]
    df$trait <- trait
    df <- df %>%
      dplyr::relocate(trait, .before = rsid)
    return(df)
  })
  
  data <- do.call(rbind, data)
  
  data.table::fwrite(
    x = data,
    file = file.path(BASE_DIR, "results/ukbb/v2/compiled_results.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# Intersections.
# -----------------------------------------------------------------------------

if (FALSE) {
  
  setwd(BASE_DIR)
  
  data <- data.table::fread(file = "results/ukbb/v2/compiled_results.tsv")
  split_data <- split(data, data$trait)
  
  gwas <- data.table::fread(file = "gwas_catalog/compiled.tsv") %>%
    dplyr::select(chr, bp) %>%
    dplyr::mutate(
      start = bp - FLANK,
      end = bp + FLANK
    )
  
  overlaps <- lapply(split_data, function(df) {
    
    standard <- df %>%
      dplyr::filter(standard_p <= 5e-8) %>%
      dplyr::select(chr, bp) %>%
      dplyr::mutate(
        start = bp - FLANK,
        end = bp + FLANK
      ) 
    
    synsurr <- df %>%
      dplyr::filter(synsurr_p <= 5e-8) %>%
      dplyr::select(chr, bp) %>%
      dplyr::mutate(
        start = bp - FLANK,
        end = bp + FLANK
      ) 
    
    synsurr_only <- df %>%
      dplyr::filter(synsurr_p <= 5e-8) %>%
      dplyr::filter(standard_p > 5e-8) %>%
      dplyr::select(chr, bp) %>%
      dplyr::mutate(
        start = bp - FLANK,
        end = bp + FLANK
      ) 
    
    over_synsurr <- CalcOverlapStats(gwas, synsurr)
    over_synsurr_only <- CalcOverlapStats(gwas, synsurr_only)
    over_standard <- CalcOverlapStats(gwas, standard)
    
    out <- data.frame(
      trait = unique(df$trait),
      standard_total = nrow(standard),
      standard_rep = length(over_standard@right_overlapped),
      synsurr_total = nrow(synsurr),
      synsurr_rep = length(over_synsurr@right_overlapped),
      synsurr_only_total = nrow(synsurr_only),
      synsurr_only_rep = length(over_synsurr_only@right_overlapped)
    )
    return(out)
  })
  
  overlaps <- do.call(rbind, overlaps)
  row.names(overlaps) <- NULL
  
  data.table::fwrite(
    x = overlaps,
    file = "gwas_catalog/v2_ukbb_gwas_catalog_overlaps.tsv",
    sep = "\t"
  )
  
}

# -----------------------------------------------------------------------------
# Proportion of SNPs replicated by GWAS catalog.
# -----------------------------------------------------------------------------

setwd(BASE_DIR)

overlaps <- data.table::fread(file = "gwas_catalog/v2_ukbb_gwas_catalog_overlaps.tsv")

df <- overlaps %>%
  dplyr::mutate(
    trait = gsub("_\\d+", "", x = overlaps$trait),
    standard_prop_rep = standard_rep / standard_total,
    synsurr_prop_rep = synsurr_rep / synsurr_total,
    synsurr_only_prop_rep = synsurr_only_rep / synsurr_only_total
  )
df$trait <- factor(
  df$trait,
  levels = c("android", "arms", "gynoid", "legs", "total", "trunk"),
  labels = c("Android\nmass", "Arms\nmass", "Gynoid\nmass", "Legs\nmass", "Total\nmass", "Trunk\nmass")
)

df1 <- df %>%
  dplyr::select(trait, standard_rep, synsurr_rep, synsurr_only_rep) %>%
  tidyr::pivot_longer(
    cols = c("standard_rep", "synsurr_rep", "synsurr_only_rep"),
    names_to = "method",
    values_to = "count"
  )
df1$method <- factor(
  x = df1$method,
  levels = c("standard_rep", "synsurr_rep", "synsurr_only_rep"),
  labels = c("Standard", "SynSurr", "SynSurr\nOnly")
)

df2 <- df %>%
  dplyr::select(trait, standard_prop_rep, synsurr_prop_rep, synsurr_only_prop_rep) %>%
  tidyr::pivot_longer(
    cols = c("standard_prop_rep", "synsurr_prop_rep", "synsurr_only_prop_rep"),
    names_to = "method",
    values_to = "count"
  )
df2$method <- factor(
  x = df2$method,
  levels = c("standard_prop_rep", "synsurr_prop_rep", "synsurr_only_prop_rep"),
  labels = c("Standard", "SynSurr", "SynSurr\nOnly")
)

q1 <- ggplot(data = df1) +
  theme_bw() + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) + 
  geom_col(
    aes(x = trait, y = count, fill = method),
    position = position_dodge()
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Number\nOverlapped"
  ) +
  ggsci::scale_fill_nejm(
    name = "Method"
  )

legend <- ggpubr::get_legend(q1)
legend <- ggpubr::as_ggplot(legend)
q1 <- q1 + theme(legend.position = "none")

q2 <- ggplot(data = df2) +
  theme_bw() + 
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) + 
  geom_col(
    aes(x = trait, y = count, fill = method),
    position = position_dodge()
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Proportion\nOverlapped"
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
    file = "results/gwas_catalog/v2_overlaps.pdf",
    width = 8.0,
    height = 7.0,
    units = "in",
    dpi = 480
  )
}

# -----------------------------------------------------------------------------
# Figure without SynSurr Only.
# -----------------------------------------------------------------------------

df1 <- df %>%
  dplyr::select(trait, standard_rep, synsurr_rep) %>%
  tidyr::pivot_longer(
    cols = c("standard_rep", "synsurr_rep"),
    names_to = "method",
    values_to = "count"
  )
df1$method <- factor(
  x = df1$method,
  levels = c("standard_rep", "synsurr_rep"),
  labels = c("Standard", "SynSurr")
)

df2 <- df %>%
  dplyr::select(trait, standard_prop_rep, synsurr_prop_rep) %>%
  tidyr::pivot_longer(
    cols = c("standard_prop_rep", "synsurr_prop_rep"),
    names_to = "method",
    values_to = "count"
  )
df2$method <- factor(
  x = df2$method,
  levels = c("standard_prop_rep", "synsurr_prop_rep"),
  labels = c("Standard", "SynSurr")
)

q1 <- ggplot(data = df1) +
  theme_bw() + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) + 
  theme(axis.line = element_line(color = 'black')) +
  geom_col(
    aes(x = trait, y = count, fill = method),
    position = position_dodge()
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Number\nOverlapped"
  ) +
  ggsci::scale_fill_nejm(
    name = "Method"
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
  theme(axis.line = element_line(color = 'black')) +
  geom_col(
    aes(x = trait, y = count, fill = method),
    position = position_dodge()
  ) +
  scale_x_discrete(
    name = "Trait"
  ) + 
  scale_y_continuous(
    name = "Proportion\nOverlapped"
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
    file = "gwas_catalog/v2_overlaps.pdf",
    width = 8.0,
    height = 7.0,
    units = "in",
    dpi = 480
  )
}
