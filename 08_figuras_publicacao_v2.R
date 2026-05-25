# 08_figuras_publicacao.R
# Figuras finais para submissão
# 600 dpi | largura 17 cm (full) ou 8.5 cm (half) | fonte base 9 pt
# Sem títulos internos — toda informação descritiva vai na legenda
# Painéis identificados por (A), (B) em negrito via plot_annotation
# Paleta: Conilon = #185FA5 | Robusta = #D85A30
#
# Numeração correta (manuscrito V2):
#   Fig 1  — Climograma (temperatura, precipitação, umidade)
#   Fig 2  — Decomposição de variância (A) + H² Bayes vs REML (B)
#   Fig 3  — r_GA Bayes vs REML
#   Fig 4  — Probabilidade superioridade consistente — FWM/bag (lollipop)
#   Fig 5  — Desempenho × Estabilidade (ecovalência)
#   Fig 6  — Multi-trait: BLUPs scatter (A) + correlações por componente (B)
#   Fig 7  — Dendrograma UPGMA Mahalanobis (A) + PCA biplot (B)
#   Fig S1 — Decomposição variância REML vs Bayes (comparação completa)
#   Fig S2 — Heatmap probabilidade por ano (todos os genótipos, 20%)
#   Fig S3 — Teste ano:block — ΔAIC
#   Fig S4 — Heatmap matriz de distâncias Mahalanobis

# ── pacotes ────────────────────────────────────────────────────────────────────
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

# ── constantes ─────────────────────────────────────────────────────────────────
DPI       <- 600
W_FULL    <- 17 / 2.54      # 17 cm em polegadas
W_HALF    <- 8.5 / 2.54     # 8.5 cm em polegadas
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

# ── tema base ──────────────────────────────────────────────────────────────────
tema_pub <- theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.title        = element_text(size = BASE_SIZE,     colour = "gray15"),
    axis.text         = element_text(size = BASE_SIZE - 1, colour = "gray20"),
    legend.title      = element_text(size = BASE_SIZE,     face = "bold"),
    legend.text       = element_text(size = BASE_SIZE - 1),
    legend.key.size   = unit(0.35, "cm"),
    panel.grid.major  = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = BASE_SIZE, face = "bold"),
    strip.background  = element_rect(fill = "gray96", colour = NA),
    plot.margin       = margin(3, 4, 3, 4, "mm"),
    legend.position   = "top",
    legend.margin     = margin(0, 0, 1, 0, "mm"),
    legend.box.spacing = unit(1, "mm")
  )

# ── carregar dados ─────────────────────────────────────────────────────────────
cat("Carregando dados...\n")

dados     <- readRDS("data/dados_clean.rds")
grupo_ref <- dados |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

res_bayes  <- read_xlsx(paste0(path_tbl, "03_resultados_bayesianos.xlsx"),
                        sheet = "parametros_geneticos")
res_comp   <- read_xlsx(paste0(path_tbl, "03b_reml_vs_bayes.xlsx"),
                        sheet = "comparacao_long")
res_aic    <- read_xlsx(paste0(path_tbl, "03c_test_ano_block.xlsx"),
                        sheet = "resultados")
res_rank   <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                        sheet = "ranking_final")
res_prob   <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                        sheet = "prob_consistente")
div_grupos <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                        sheet = "grupos_divergencia")
pca_scores <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                        sheet = "pca_scores")
pca_loads  <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                        sheet = "pca_loadings")
res_metodos <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                         sheet = "resumo_metodos")
res_corr   <- read_xlsx(paste0(path_tbl, "06_multitrait_results.xlsx"),
                        sheet = "correlacoes")
res_mt     <- read_xlsx(paste0(path_tbl, "06_multitrait_results.xlsx"),
                        sheet = "blups_multitrait")
meteo      <- read_xlsx(paste0(path_tbl, "07_resumo_meteorologico.xlsx"),
                        sheet = "dados_mensais")

cat("Dados carregados.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 1 — Climograma: temperatura (A), precipitação (B), umidade (C)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 1 — Climograma...\n")

cores_ano <- c("2022" = "#8DA0B8", "2023" = COR_BAYES, "2024" = COR_ROBUSTA)

df_m <- meteo |>
  filter(ano %in% 2022:2024) |>
  mutate(
    ano_f     = factor(ano),
    mes_label = factor(mes, 1:12,
                       labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                  "Jul","Aug","Sep","Oct","Nov","Dec"))
  )

# (A) temperatura
p1a <- ggplot(df_m, aes(x = mes_label, group = ano_f, colour = ano_f)) +
  geom_ribbon(aes(ymin = temp_min, ymax = temp_max, fill = ano_f),
              alpha = 0.10, colour = NA) +
  geom_line(aes(y = temp_media), linewidth = 0.55) +
  geom_point(aes(y = temp_media), size = 1.1) +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = 14, ymax = 36,
           fill = "gray80", alpha = 0.18) +
  scale_colour_manual(values = cores_ano, name = NULL) +
  scale_fill_manual(values = cores_ano, name = NULL, guide = "none") +
  scale_y_continuous(name = "Temperature (°C)",
                     limits = c(14, 36), breaks = seq(15, 35, 5)) +
  labs(x = NULL) +
  tema_pub +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

# (B) precipitação
p1b <- ggplot(df_m, aes(x = mes_label, y = precip,
                        fill = ano_f, group = ano_f)) +
  geom_col(position = position_dodge(0.72), width = 0.65, alpha = 0.85) +
  annotate("rect", xmin = 3.5, xmax = 7.5,
           ymin = 0, ymax = max(df_m$precip, na.rm = TRUE) * 1.15,
           fill = "gray80", alpha = 0.18) +
  scale_fill_manual(values = cores_ano, name = NULL) +
  scale_y_continuous(
    name   = "Precipitation (mm)",
    limits = c(0, max(df_m$precip, na.rm = TRUE) * 1.15),
    breaks = pretty_breaks(5)
  ) +
  labs(x = NULL) +
  tema_pub +
  theme(legend.position    = "none",
        panel.grid.major.x = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank())

# (C) umidade relativa
p1c <- ggplot(df_m, aes(x = mes_label, y = ur_media,
                        colour = ano_f, group = ano_f)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.1) +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = 60, ymax = 95,
           fill = "gray80", alpha = 0.18) +
  scale_colour_manual(values = cores_ano, name = NULL) +
  scale_y_continuous(name = "Relative humidity (%)",
                     limits = c(60, 95), breaks = seq(60, 95, 5)) +
  labs(x = NULL) +
  tema_pub +
  theme(legend.position = "none")

fig1 <- (p1a / p1b / p1c) +
  plot_layout(heights = c(1.1, 1, 1)) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = BASE_SIZE, face = "bold"))
  )

ggsave(paste0(path_fig, "Fig1_climograma.tiff"),
       fig1, width = W_FULL, height = 14 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 1 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 2 — Decomposição de variância (A) + H² Bayes vs REML (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 2 — Varcomp + H²...\n")

# (A) proporção dos componentes — medianas bayesianas
df_vc <- res_bayes |>
  filter(variavel %in% names(var_labels)) |>
  select(variavel, sigma2_G_med, sigma2_GxA_med, sigma2_e_med) |>
  pivot_longer(-variavel, names_to = "componente", values_to = "valor") |>
  mutate(
    componente = factor(componente,
                        levels = c("sigma2_G_med","sigma2_GxA_med","sigma2_e_med"),
                        labels = c("G","G×A","Residual")),
    variavel   = factor(variavel, levels = names(var_labels),
                        labels = var_labels)
  ) |>
  group_by(variavel) |>
  mutate(prop = valor / sum(valor)) |>
  ungroup()

p2a <- ggplot(df_vc, aes(x = variavel, y = prop, fill = componente)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("G" = COR_G, "G×A" = COR_GXA, "Residual" = COR_E),
    name   = "Variance component"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = NULL, y = "Proportion of total variance") +
  tema_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# (B) H² Bayes (IC 95%) vs REML (ponto)
df_h2 <- res_comp |>
  filter(variavel %in% names(var_labels)) |>
  mutate(
    variavel = factor(variavel, levels = names(var_labels),
                      labels = var_labels),
    boundary = as.logical(boundary),
    metodo   = factor(metodo, levels = c("Bayes","REML"))
  )

p2b <- ggplot() +
  geom_linerange(
    data = df_h2 |> filter(metodo == "Bayes"),
    aes(x = variavel, ymin = H2_q025, ymax = H2_q975),
    colour = COR_BAYES, linewidth = 0.7
  ) +
  geom_point(
    data = df_h2 |> filter(metodo == "Bayes"),
    aes(x = variavel, y = H2),
    colour = COR_BAYES, size = 2.5, shape = 16
  ) +
  geom_point(
    data = df_h2 |> filter(metodo == "REML"),
    aes(x = variavel, y = H2,
        shape = ifelse(boundary, "boundary", "estimate")),
    colour = COR_REML, size = 2.5
  ) +
  scale_shape_manual(
    values = c("boundary" = 4, "estimate" = 17),
    name   = "REML",
    labels = c("boundary" = "boundary (σ²G = 0)", "estimate" = "point estimate")
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "gray60", linewidth = 0.3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.0),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(x = NULL, y = expression(italic(H)^2)) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", legend.box = "horizontal")

fig2 <- (p2a | p2b) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = BASE_SIZE, face = "bold"))
  )

ggsave(paste0(path_fig, "Fig2_varcomp_h2.tiff"),
       fig2, width = W_FULL, height = 9 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 2 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 3 — r_GA Bayes vs REML
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 3 — r_GA...\n")

df_rga <- res_comp |>
  filter(variavel %in% names(var_labels)) |>
  mutate(
    variavel = factor(variavel, levels = names(var_labels),
                      labels = var_labels),
    metodo   = factor(metodo, levels = c("Bayes","REML"))
  )

fig3 <- ggplot() +
  geom_linerange(
    data = df_rga |> filter(metodo == "Bayes"),
    aes(x = variavel, ymin = rGA_q025, ymax = rGA_q975),
    colour = COR_BAYES, linewidth = 0.7
  ) +
  geom_point(
    data = df_rga |> filter(metodo == "Bayes"),
    aes(x = variavel, y = rGA),
    colour = COR_BAYES, size = 2.5, shape = 16
  ) +
  geom_point(
    data = df_rga |> filter(metodo == "REML"),
    aes(x = variavel, y = rGA),
    colour = COR_REML, size = 2.5, shape = 17
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "gray60", linewidth = 0.3) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    x = NULL,
    y = expression(italic(r)[GA]~"(genotypic correlation between years)")
  ) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "none")

ggsave(paste0(path_fig, "Fig3_rGA.tiff"),
       fig3, width = W_HALF * 1.5, height = 7.5 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 3 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 4 — Probabilidade de superioridade consistente — FWM/bag (lollipop)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 4 — Lollipop prob. consistente...\n")

df_lollipop <- res_prob |>
  filter(variavel == "mcm_saca") |>
  select(-any_of("grupo")) |>
  left_join(grupo_ref, by = "genotipo") |>
  group_by(genotipo, grupo) |>
  mutate(prob_max = max(prob_consistente, na.rm = TRUE)) |>
  ungroup() |>
  arrange(desc(prob_max)) |>
  filter(genotipo %in% unique(genotipo)[1:20]) |>
  mutate(
    intensidade = factor(intensidade,
                         levels = c(0.10, 0.20, 0.25),
                         labels = c("10%","20%","25%")),
    genotipo    = fct_reorder(genotipo, prob_max, .desc = FALSE)
  )

fig4 <- ggplot(df_lollipop,
               aes(x = prob_consistente, y = genotipo,
                   colour = intensidade)) +
  geom_line(aes(group = genotipo), colour = "gray70", linewidth = 0.4) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0.50, linetype = "dashed",
             colour = "gray40", linewidth = 0.4) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_colour_manual(
    values = c("10%" = "#2980B9", "20%" = "#27AE60", "25%" = "#E74C3C"),
    name   = "Selection intensity"
  ) +
  facet_wrap(~ grupo, scales = "free_y", ncol = 2) +
  labs(x = "Probability of consistent superiority", y = NULL) +
  tema_pub +
  theme(
    legend.position    = "top",
    panel.grid.major.x = element_line(colour = "gray92", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

ggsave(paste0(path_fig, "Fig4_prob_consistente.tiff"),
       fig4, width = W_FULL, height = 10 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 4 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 5 — Desempenho × Estabilidade (ecovalência de Wricke)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 5 — Desempenho × Estabilidade...\n")

df_prob20 <- res_prob |>
  filter(variavel %in% vars_stab, intensidade == 0.20) |>
  select(variavel, genotipo, prob_consistente)

df_stab <- res_rank |>
  filter(variavel %in% vars_stab) |>
  select(-any_of("grupo")) |>
  left_join(df_prob20, by = c("variavel","genotipo")) |>
  left_join(grupo_ref, by = "genotipo") |>
  rename(prob_consistente = prob_consistente.y) |>
  mutate(
    prob_consistente = coalesce(prob_consistente, 0),
    variavel = factor(variavel, levels = vars_stab,
                      labels = var_labels[vars_stab]),
    rotular  = prob_consistente > 0.10 |
      abs(blup_sinal) > quantile(abs(blup_sinal), 0.85, na.rm = TRUE)
  )

fig5 <- ggplot(df_stab,
               aes(x = blup_sinal, y = wricke,
                   colour = grupo, size = prob_consistente)) +
  geom_point(alpha = 0.78) +
  geom_hline(
    data = df_stab |>
      group_by(variavel) |>
      summarise(med_w = mean(wricke, na.rm = TRUE), .groups = "drop"),
    aes(yintercept = med_w),
    linetype = "dashed", colour = "gray60", linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray60", linewidth = 0.3) +
  geom_text_repel(
    aes(label = ifelse(rotular, genotipo, "")),
    size = 1.8, max.overlaps = 20, show.legend = FALSE,
    min.segment.length = 0.2,
    segment.size = 0.22, segment.colour = "gray55",
    box.padding = 0.25, point.padding = 0.15,
    force = 1.5
  ) +
  scale_colour_manual(
    values = c("Conilon" = COR_CONILON, "Robusta" = COR_ROBUSTA),
    name   = "Botanical group"
  ) +
  scale_size_continuous(
    name   = "Prob. consistent\nsuperiority (20%)",
    range  = c(0.5, 3.8),
    breaks = c(0, 0.25, 0.50, 0.75)
  ) +
  facet_wrap(~ variavel, scales = "free", ncol = 3) +
  labs(
    x = "BLUP — mean performance (sign-adjusted; positive = better)",
    y = "Wricke ecovalence (lower = more stable)"
  ) +
  tema_pub +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(paste0(path_fig, "Fig5_desempenho_estabilidade.tiff"),
       fig5, width = W_FULL, height = 12 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 5 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 6 — Multi-trait: BLUPs scatter (A) + correlações por componente (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 6 — Multi-trait...\n")

rG   <- res_corr |> filter(componente == "Genótipo (G)")  |> pull(correlacao)
rGxA <- res_corr |> filter(grepl("Ano", componente))      |> pull(correlacao)
re   <- res_corr |> filter(componente == "Residual")       |> pull(correlacao)

df_mt <- res_mt |>
  select(-any_of("grupo")) |>
  left_join(grupo_ref, by = "genotipo") |>
  mutate(rotular = abs(indice_mt) > 0.5)

rG_x <- max(df_mt$blup_per_grao) * 0.85
rG_y <- max(df_mt$blup_mcm_saca) * 0.88

p6a <- ggplot(df_mt,
              aes(x = blup_per_grao, y = blup_mcm_saca,
                  colour = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray70", linewidth = 0.3) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "gray30", linewidth = 0.5) +
  geom_point(aes(size = abs(indice_mt)), alpha = 0.82) +
  geom_text_repel(
    aes(label = ifelse(rotular, genotipo, "")),
    size = 1.9, max.overlaps = 20, show.legend = FALSE,
    min.segment.length = 0.2,
    segment.size = 0.22, segment.colour = "gray55"
  ) +
  geom_text(
    data = tibble(x = rG_x, y = rG_y,
                  label = paste0("r[G] == ", round(rG, 3))),
    aes(x = x, y = y, label = label),
    parse = TRUE, size = 2.8, hjust = 1,
    colour = "gray15", fontface = "italic",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("Conilon" = COR_CONILON, "Robusta" = COR_ROBUSTA),
    name   = "Botanical group"
  ) +
  scale_size_continuous(range = c(0.7, 3.2), guide = "none") +
  labs(
    x = "BLUP — % grain (deviation from mean)",
    y = "BLUP — FWM/bag, kg (deviation from mean)"
  ) +
  tema_pub

p6b <- ggplot(res_corr,
              aes(x = factor(componente,
                             levels = c("Genótipo (G)",
                                        "Genótipo × Ano (GxA)",
                                        "Residual"),
                             labels = c("G","G×A","Residual")),
                  y = correlacao, fill = correlacao)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 0, colour = "gray40", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.3f", correlacao),
        vjust = ifelse(correlacao < 0, 1.5, -0.5)),
    size = 2.4, colour = "gray15"
  ) +
  scale_fill_gradient2(
    low = COR_ROBUSTA, mid = "white", high = COR_CONILON,
    midpoint = 0, limits = c(-1, 1),
    name = "Correlation"
  ) +
  scale_y_continuous(limits = c(-1.05, 0.18)) +
  labs(x = "Variance component",
       y = "Correlation (% grain × FWM/bag)") +
  tema_pub +
  theme(legend.position = "right")

fig6 <- (p6a | p6b) +
  plot_layout(widths = c(1.1, 0.9)) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = BASE_SIZE, face = "bold"))
  )

ggsave(paste0(path_fig, "Fig6_multitrait.tiff"),
       fig6, width = W_FULL, height = 9 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 6 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA 7 — Dendrograma UPGMA Mahalanobis (A) + PCA biplot (B)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. 7 — Dendrograma + PCA...\n")

mat_num <- div_grupos |>
  column_to_rownames("genotipo") |>
  select(per_grao, mcm_saca) |>
  as.matrix()

S_inv    <- solve(cov(mat_num))
n        <- nrow(mat_num)
dist_mah <- matrix(0, n, n,
                   dimnames = list(rownames(mat_num), rownames(mat_num)))
for (i in 1:n) for (j in 1:n) {
  d <- mat_num[i,] - mat_num[j,]
  dist_mah[i,j] <- sqrt(t(d) %*% S_inv %*% d)
}
dist_mah <- as.dist(dist_mah)

hc       <- hclust(dist_mah, method = "average")
dc_cut   <- res_metodos$dc_final
ccc_val  <- res_metodos$ccc_final

dend_data <- dendro_data(hc, type = "rectangle")
labs_dend <- dend_data$labels |>
  left_join(grupo_ref, by = c("label" = "genotipo"))

offset <- diff(range(dend_data$segments$y)) * 0.09

p7a <- ggplot() +
  geom_segment(
    data = dend_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.28, colour = "gray30"
  ) +
  geom_hline(
    yintercept = dc_cut + offset,
    linetype = "dashed", colour = COR_ROBUSTA, linewidth = 0.4
  ) +
  geom_text(
    data = labs_dend,
    aes(x = x,
        y = -0.03 * max(dend_data$segments$y),
        label = label, colour = grupo),
    size = 1.7, hjust = 1, angle = 90
  ) +
  scale_colour_manual(
    values = c("Conilon" = COR_CONILON, "Robusta" = COR_ROBUSTA),
    name   = "Botanical group"
  ) +
  scale_y_continuous(
    name   = "Mahalanobis fusion distance",
    expand = expansion(mult = c(0.30, 0.05))
  ) +
  annotate(
    "text",
    x     = length(hc$order) * 0.98,
    y     = (dc_cut + offset) * 1.06,
    label = paste0("dc = ", round(dc_cut, 2),
                   " | CCC = ", round(ccc_val, 3)),
    size  = 1.9, hjust = 1, colour = COR_ROBUSTA
  ) +
  tema_pub +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top"
  )

pca_scores_en <- pca_scores |>
  mutate(
    grupo_div = gsub("Grupo", "Group", grupo_div),
    grupo     = factor(grupo)
  )

p7b <- ggplot(pca_scores_en,
              aes(x = PC1, y = PC2,
                  colour = grupo, shape = grupo_div)) +
  geom_hline(yintercept = 0, colour = "gray82", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "gray82", linewidth = 0.25) +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_text_repel(
    aes(label = genotipo),
    size = 1.7, max.overlaps = 20, show.legend = FALSE,
    min.segment.length = 0.2,
    segment.size = 0.22, segment.colour = "gray60"
  ) +
  geom_segment(
    data = pca_loads,
    aes(x = 0, y = 0,
        xend = PC1 * max(abs(pca_scores_en$PC1)) * 0.72,
        yend = PC2 * max(abs(pca_scores_en$PC2)) * 0.72),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    colour = "gray20", linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_label(
    data = pca_loads |>
      mutate(variavel_lab = gsub("% grão", "% grain",
                                 gsub("MCM/saca", "FWM/bag", variavel_lab))),
    aes(x = PC1 * max(abs(pca_scores_en$PC1)) * 0.85,
        y = PC2 * max(abs(pca_scores_en$PC2)) * 0.85,
        label = variavel_lab),
    size = 2.0, colour = "gray15",
    fill = alpha("white", 0.85), linewidth = 0,
    label.padding = unit(0.08, "cm"),
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("Conilon" = COR_CONILON, "Robusta" = COR_ROBUSTA),
    name   = "Botanical group"
  ) +
  scale_shape_manual(
    values = setNames(c(16, 17, 15, 3, 4, 8), paste0("Group ", 1:6)),
    name   = "Divergence group"
  ) +
  labs(x = "PC1 (87.6%)", y = "PC2 (12.4%)") +
  tema_pub +
  theme(legend.position = "right", legend.box = "vertical")

fig7 <- (p7a | p7b) +
  plot_layout(widths = c(1.1, 0.9)) +
  plot_annotation(
    tag_levels = "A", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = BASE_SIZE, face = "bold"))
  )

ggsave(paste0(path_fig, "Fig7_divergencia.tiff"),
       fig7, width = 24 / 2.54, height = 10 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. 7 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA S1 — Decomposição de variância REML vs Bayes (comparação completa)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. S1 — REML vs Bayes varcomp...\n")

df_s1 <- res_comp |>
  filter(variavel %in% names(var_labels)) |>
  select(variavel, metodo, sigma2_G, sigma2_GxA, sigma2_e) |>
  pivot_longer(cols = starts_with("sigma2"),
               names_to = "componente", values_to = "valor") |>
  mutate(
    componente = factor(componente,
                        levels = c("sigma2_G","sigma2_GxA","sigma2_e"),
                        labels = c("G","G×A","Residual")),
    variavel   = factor(variavel, levels = names(var_labels),
                        labels = var_labels),
    metodo     = factor(metodo, levels = c("Bayes","REML"))
  ) |>
  group_by(variavel, metodo) |>
  mutate(prop = valor / sum(valor)) |>
  ungroup()

figS1 <- ggplot(df_s1,
                aes(x = interaction(metodo, variavel, sep = "\n"),
                    y = prop, fill = componente)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("G" = COR_G, "G×A" = COR_GXA, "Residual" = COR_E),
    name   = "Variance component"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  facet_wrap(~ variavel, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Proportion of total variance") +
  tema_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

ggsave(paste0(path_fig, "FigS1_reml_vs_bayes.tiff"),
       figS1, width = W_FULL, height = 8 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. S1 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA S2 — Heatmap probabilidade por ano (todos os genótipos, 20%)
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. S2 — Heatmap prob. por ano...\n")

df_hm <- res_prob |>
  filter(variavel %in% vars_stab, intensidade == 0.20) |>
  left_join(grupo_ref, by = "genotipo") |>
  pivot_longer(cols = c(prob_2023, prob_2024),
               names_to = "ano", values_to = "prob") |>
  mutate(
    ano      = gsub("prob_", "", ano),
    variavel = factor(variavel, levels = vars_stab,
                      labels = var_labels[vars_stab]),
    genotipo = fct_reorder(genotipo, prob_consistente)
  )

figS2 <- ggplot(df_hm, aes(x = ano, y = genotipo, fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.08) +
  scale_fill_gradient2(
    low      = "#2C3E50",
    mid      = "#27AE60",
    high     = "#E74C3C",
    midpoint = 0.5,
    limits   = c(0, 1),
    labels   = percent_format(accuracy = 1),
    name     = "Pr(superior)"
  ) +
  facet_wrap(~ variavel, ncol = 5) +
  labs(x = NULL, y = NULL) +
  tema_pub +
  theme(
    axis.text.y     = element_text(size = 5.0),
    axis.text.x     = element_text(size = BASE_SIZE - 1),
    legend.position = "right",
    panel.spacing   = unit(0.4, "mm")
  )

ggsave(paste0(path_fig, "FigS2_heatmap_prob.tiff"),
       figS2, width = W_FULL, height = 16 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. S2 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA S3 — Teste ano:block — ΔAIC
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. S3 — Teste ano:block...\n")

df_aic_plot <- res_aic |>
  filter(variavel %in% names(var_labels)) |>
  mutate(variavel = factor(variavel, levels = names(var_labels),
                           labels = var_labels))

figS3 <- ggplot(df_aic_plot, aes(x = variavel, y = delta_aic)) +
  geom_col(fill = "gray60", width = 0.55) +
  geom_hline(yintercept = -2, linetype = "dashed",
             colour = COR_ROBUSTA, linewidth = 0.4) +
  geom_text(aes(label = paste0("p = ", round(lrt_pval, 3)),
                vjust = -0.4),
            size = 2.2, colour = "gray20") +
  scale_y_continuous(
    name   = "ΔAIC (with year×block − without)",
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  labs(x = NULL) +
  tema_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(paste0(path_fig, "FigS3_aic_anoblock.tiff"),
       figS3, width = W_HALF * 1.5, height = 7.5 / 2.54,
       dpi = DPI, compression = "lzw")
cat("  Fig. S3 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURA S4 — Heatmap matriz de distâncias Mahalanobis
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Fig. S4 — Heatmap distâncias Mahalanobis...\n")

# reconstruir matriz a partir dos BLUPs salvos, mesma lógica da Fig. 7
order_hc <- hc$labels[hc$order]

dist_df <- as.data.frame(as.matrix(dist_mah)) |>
  rownames_to_column("gen_row") |>
  pivot_longer(-gen_row, names_to = "gen_col", values_to = "distancia") |>
  mutate(
    gen_row = factor(gen_row, levels = order_hc),
    gen_col = factor(gen_col, levels = order_hc)
  )

figS4 <- ggplot(dist_df, aes(x = gen_col, y = gen_row, fill = distancia)) +
  geom_tile(colour = NA) +
  scale_fill_gradientn(
    colours = c("#2C3E50", "#2980B9", "#27AE60", "#F39C12", "#E74C3C"),
    name    = "Mahalanobis\ndistance"
  ) +
  labs(x = NULL, y = NULL) +
  tema_pub +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5),
    axis.text.y     = element_text(size = 5),
    legend.position = "right",
    panel.grid      = element_blank()
  )

ggsave(paste0(path_fig, "FigS4_heatmap_distancias.tiff"),
       figS4, width = W_FULL, height = W_FULL,
       dpi = DPI, compression = "lzw")
cat("  Fig. S4 salva.\n")

# ══════════════════════════════════════════════════════════════════════════════
cat("\n━━ 08_figuras_publicacao.R concluído ━━\n")
cat("Figuras salvas em outputs/figures/ — formato TIFF 600 dpi\n")
cat("\nArquivos gerados:\n")
cat("  Principal: Fig1_climograma | Fig2_varcomp_h2 | Fig3_rGA\n")
cat("             Fig4_prob_consistente | Fig5_desempenho_estabilidade\n")
cat("             Fig6_multitrait | Fig7_divergencia\n")
cat("  Suplementar: FigS1_reml_vs_bayes | FigS2_heatmap_prob\n")
cat("               FigS3_aic_anoblock  | FigS4_heatmap_distancias\n")
