# 08_figuras_publicacao.R
# Final figures for submission
# 600 dpi | width 17 cm (full) or 8.5 cm (half) | base font 9 pt
# No internal titles — all descriptive information goes into the legend
# Panels identified by (A), (B) in bold via plot_annotation
# Palette: Conilon = #185FA5 | Robusta = #D85A30
#
# Numbering (manuscript V2):
#   Fig 1  — Climograph (temperature, precipitation, humidity)
#   Fig 2  — Variance decomposition (A) + H² Bayes vs REML (B)
#   Fig 3  — r_GA Bayes vs REML
#   Fig 4  — Probability of consistent superiority — FWM/bag (lollipop)
#   Fig 5  — Performance × Stability (ecovalence)
#   Fig 6  — Multi-trait: BLUPs scatter (A) + correlations by component (B)
#   Fig 7  — UPGMA Dendrogram Mahalanobis (A) + PCA biplot (B)
#   Fig S1 — Variance decomposition REML vs Bayes (full comparison)
#   Fig S2 — Heatmap probability by year (all genotypes, 20%)
#   Fig S3 — Year:block test — ΔAIC
#   Fig S4 — Mahalanobis distance matrix heatmap
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(readxl)
library(forcats)
library(scales)
library(ggdendro)

# ── constants ─────────────────────────────────────────────────────────────────
DPI       <- 600
W_FULL    <- 17 / 2.54        # 17 cm in inches
W_HALF    <- 8.5 / 2.54       # 8.5 cm in inches
BASE_SIZE <- 9

COR_CONILON <- "#185FA5"
COR_ROBUSTA <- "#D85A30"
COR_BAYES   <- "#2E6FBF"
COR_REML    <- "#C0392B"
COR_G       <- "#2C3E50"
COR_GXA     <- "#E67E22"
COR_E       <- "#95A5A6"

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

var_labels <- c(
  per_grao  = "% grain",
  per_palha = "% husk",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag"
)

vars_stab <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca")

# ── base theme ──────────────────────────────────────────────────────────────────
tema_pub <- theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.title         = element_text(size = BASE_SIZE,       colour = "gray15"),
    axis.text          = element_text(size = BASE_SIZE - 1, colour = "gray20"),
    legend.title       = element_text(size = BASE_SIZE,       face = "bold"),
    legend.text        = element_text(size = BASE_SIZE - 1),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.major   = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = BASE_SIZE, face = "bold"),
    strip.background   = element_rect(fill = "gray96", colour = NA),
    plot.margin        = margin(3, 4, 3, 4, "mm"),
    legend.position    = "top",
    legend.margin      = margin(0, 0, 1, 0, "mm"),
    legend.box.spacing = unit(1, "mm")
  )

# ── load data ─────────────────────────────────────────────────────────────
cat("Loading data...\n")

dados     <- readRDS("data/dados_clean.rds")
grupo_ref <- dados |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

res_bayes   <- read_xlsx(paste0(path_tbl, "03_resultados_bayesianos.xlsx"), sheet = "parametros_geneticos")
res_comp    <- read_xlsx(paste0(path_tbl, "03b_reml_vs_bayes.xlsx"), sheet = "comparacao_long")
res_aic     <- read_xlsx(paste0(path_tbl, "03c_test_ano_block.xlsx"), sheet = "resultados")
res_rank    <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"), sheet = "ranking_final")
res_prob    <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"), sheet = "prob_consistente")
div_grupos  <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"), sheet = "grupos_divergencia")
pca_scores  <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"), sheet = "pca_scores")
pca_loads   <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"), sheet = "pca_loadings")
res_metodos <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"), sheet = "resumo_metodos")
res_corr    <- read_xlsx(paste0(path_tbl, "06_multitrait_results.xlsx"), sheet = "correlacoes")
res_mt      <- read_xlsx(paste0(path_tbl, "06_multitrait_results.xlsx"), sheet = "blups_multitrait")
meteo       <- read_xlsx(paste0(path_tbl, "07_resumo_meteorologico.xlsx"), sheet = "dados_mensais")

cat("Data loaded.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIG 1 — Climograph: temperature (A), precipitation (B), humidity (C)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 1 — Climograph...\n")
# [Figure 1 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 2 — Variance decomposition (A) + H² Bayes vs REML (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 2 — Varcomp + H²...\n")
# [Figure 2 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 3 — r_GA Bayes vs REML
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 3 — r_GA...\n")
# [Figure 3 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 4 — Probability of consistent superiority — FWM/bag (lollipop)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 4 — Consistent superiority...\n")
# [Figure 4 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 5 — Performance × Stability (Wricke's ecovalence)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 5 — Performance × Stability...\n")
# [Figure 5 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 6 — Multi-trait: BLUPs scatter (A) + correlations by component (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 6 — Multi-trait...\n")
# [Figure 6 code structure...]

# ══════════════════════════════════════════════════════════════════════════════
# FIG 7 — UPGMA Mahalanobis Dendrogram (A) + PCA biplot (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Fig. 7 — Dendrogram + PCA...\n")
# [Figure 7 code structure...]

cat("\n━━ 08_figuras_publicacao.R completed ━━\n")
cat("Figures saved in outputs/figures/ — TIFF 600 dpi format\n")
