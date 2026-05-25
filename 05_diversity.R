# 05_diversity.R
# Divergência genética entre genótipos de C. canephora
# baseada nos BLUPs bayesianos das variáveis de rendimento
#
# Métodos:
#   - Distância de Mahalanobis (com fallback para euclidiana padronizada)
#   - Agrupamento hierárquico: UPGMA e Ward.D2
#   - Critério de Mojena (k = 1.25) para ponto de corte
#   - Coeficiente de correlação cofenética (CCC)
#   - PCA complementar
#   - Contribuição de variáveis para a divergência
#   - Análise intra-grupo (Conilon e Robusta separados)
#
# Saídas:
#   outputs/tables/05_diversity_results.xlsx
#   outputs/figures/05a_dendrograma_upgma.png
#   outputs/figures/05b_dendrograma_ward.png
#   outputs/figures/05c_pca_biplot.png
#   outputs/figures/05d_contribuicao_variaveis.png
#   outputs/figures/05e_heatmap_distancias.png

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readxl)
library(ggplot2)
library(patchwork)
library(writexl)
library(ggdendro)
library(ggrepel)
library(factoextra)
library(cluster)

# ── dados ──────────────────────────────────────────────────────────────────────
path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# carregar BLUPs
blups_G <- read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "blups_genotipo"
) |>
  mutate(
    genotipo = if_else(genotipo == "Pé.de.Ouro", "Pé de Ouro", genotipo),
    grupo    = if_else(genotipo == "Pé de Ouro" & is.na(grupo),
                       "Robusta", grupo)
  )

# variáveis para divergência — excluir vcm_mcm e vcm_saca
# vcm_mcm: sem estrutura genética
# vcm_saca: H² muito baixo, maioria ambiental
# foco nas variáveis com base genética mais sólida
vars_div <- c("per_grao", "mcm_mgb", "mcm_saca")

cat("\n── Variáveis usadas na divergência:", vars_div, "──\n")

# ── construir matriz de BLUPs ──────────────────────────────────────────────────
mat_blups <- blups_G |>
  filter(variavel %in% vars_div) |>
  select(genotipo, variavel, blup_orig_med) |>
  pivot_wider(
    names_from  = variavel,
    values_from = blup_orig_med
  ) |>
  arrange(genotipo)

# guardar grupo separado
grupo_ref <- blups_G |>
  filter(variavel == vars_div[1]) |>
  select(genotipo, grupo) |>
  mutate(genotipo = if_else(genotipo == "Pé.de.Ouro",
                            "Pé de Ouro", genotipo))

# matriz numérica com rownames
gen_names <- mat_blups$genotipo
mat_num   <- mat_blups |>
  select(all_of(vars_div)) |>
  as.matrix()
rownames(mat_num) <- gen_names

cat("Dimensões da matriz de BLUPs:", nrow(mat_num), "×", ncol(mat_num), "\n")
cat("Genótipos:", nrow(mat_num), "\n\n")

# ── 1. teste de adequação — verificar invertibilidade ─────────────────────────
S     <- cov(mat_num)
det_S <- det(S)

cat("── Matriz de covariância dos BLUPs ──\n")
print(round(S, 4))
cat("Determinante:", round(det_S, 6), "\n")
cat("Condição (kappa):", round(kappa(S), 2), "\n\n")

# definir método de distância com base na invertibilidade
usar_mahalanobis <- det_S > 1e-10 && kappa(S) < 1000

if (usar_mahalanobis) {
  cat("Matriz invertível — usando distância de Mahalanobis\n\n")
  S_inv  <- solve(S)
  dist_mat <- as.dist(
    outer(seq_len(nrow(mat_num)), seq_len(nrow(mat_num)),
          FUN = Vectorize(function(i, j) {
            diff_ij <- mat_num[i, ] - mat_num[j, ]
            sqrt(t(diff_ij) %*% S_inv %*% diff_ij)
          }))
  )
  metodo_dist <- "Mahalanobis"
} else {
  cat("Matriz singular ou mal condicionada — usando euclidiana padronizada\n\n")
  mat_scaled <- scale(mat_num)
  dist_mat   <- dist(mat_scaled, method = "euclidean")
  metodo_dist <- "Euclidiana padronizada"
}

# ── 2. agrupamento hierárquico ─────────────────────────────────────────────────

# UPGMA
hc_upgma <- hclust(dist_mat, method = "average")

# Ward.D2
hc_ward  <- hclust(dist_mat, method = "ward.D2")

# ── 3. coeficiente cofenético ──────────────────────────────────────────────────
ccc_upgma <- cor(dist_mat, cophenetic(hc_upgma))
ccc_ward  <- cor(dist_mat, cophenetic(hc_ward))

cat("── Coeficiente de correlação cofenética (CCC) ──\n")
cat("  UPGMA:", round(ccc_upgma, 4), "\n")
cat("  Ward.D2:", round(ccc_ward, 4), "\n")
cat("  Método preferido:", if_else(ccc_upgma >= ccc_ward,
                                   "UPGMA", "Ward.D2"), "\n\n")

# método com melhor CCC
hc_melhor  <- if (ccc_upgma >= ccc_ward) hc_upgma else hc_ward
nome_melhor <- if (ccc_upgma >= ccc_ward) "UPGMA" else "Ward.D2"

# ── 4. ponto de corte — critério de Mojena ────────────────────────────────────
mojena_k <- 1.25

fusao      <- hc_melhor$height
dc_mojena  <- mean(fusao) + mojena_k * sd(fusao)
n_grupos   <- length(unique(cutree(hc_melhor, h = dc_mojena)))

cat("── Critério de Mojena (k =", mojena_k, ") ──\n")
cat("  Distância de corte:", round(dc_mojena, 4), "\n")
cat("  Número de grupos:", n_grupos, "\n\n")

# grupos por genótipo
grupos_gen <- cutree(hc_melhor, h = dc_mojena)

df_grupos <- tibble(
  genotipo     = names(grupos_gen),
  grupo_div    = paste0("Grupo ", grupos_gen)
) |>
  left_join(grupo_ref, by = "genotipo")

cat("── Composição dos grupos de divergência ──\n")
df_grupos |>
  count(grupo_div, grupo) |>
  arrange(grupo_div, grupo) |>
  print()

# ── 5. contribuição de cada variável para a divergência ───────────────────────
# variância explicada de cada variável na separação dos grupos
contrib_var <- tibble(
  variavel    = vars_div,
  variancia   = apply(mat_num, 2, var),
  cv          = apply(mat_num, 2, sd) / abs(apply(mat_num, 2, mean)) * 100
) |>
  mutate(
    prop_var = variancia / sum(variancia) * 100
  )

cat("\n── Contribuição das variáveis ──\n")
print(contrib_var)

# ── 6. PCA ────────────────────────────────────────────────────────────────────
pca_res <- prcomp(mat_num, scale. = TRUE, center = TRUE)

# variância explicada
var_exp <- summary(pca_res)$importance[2, ] * 100
cat("\n── PCA — variância explicada ──\n")
print(round(var_exp, 2))

# scores com grupos
df_pca <- as_tibble(pca_res$x[, 1:3]) |>
  mutate(
    genotipo  = gen_names,
    grupo_div = paste0("Grupo ", grupos_gen)
  ) |>
  left_join(grupo_ref, by = "genotipo") |>
  left_join(
    df_grupos |> select(genotipo, grupo_div),
    by = "genotipo"
  )

# loadings
df_loadings <- as_tibble(pca_res$rotation[, 1:2],
                         rownames = "variavel") |>
  mutate(
    variavel_lab = c(
      per_grao = "% grão",
      mcm_mgb  = "MCM/MGB",
      mcm_saca = "MCM/saca"
    )[variavel]
  )

# ── 7. visualizações ───────────────────────────────────────────────────────────

cores_grupo  <- c("Conilon" = "#185FA5", "Robusta" = "#D85A30")
escala_grupo <- scale_colour_manual(values = cores_grupo, name = "grupo botânico")

# helper — dendrograma com ggplot
plot_dendro <- function(hc, titulo, ccc, dc, n_grp) {
  
  dend_data <- ggdendro::dendro_data(hc, type = "rectangle")
  
  # adicionar grupo botânico aos labels
  labels_df <- dend_data$labels |>
    left_join(grupo_ref, by = c("label" = "genotipo"))
  
  ggplot() +
    geom_segment(
      data = dend_data$segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.35, colour = "gray30"
    ) +
    geom_hline(
      yintercept = dc,
      linetype   = "dashed",
      colour     = "#D85A30",
      linewidth  = 0.5
    ) +
    geom_text(
      data = labels_df,
      aes(x = x, y = -0.02 * max(dend_data$segments$y),
          label = label, colour = grupo),
      size  = 2.2, hjust = 1, angle = 90
    ) +
    escala_grupo +
    scale_y_continuous(
      expand = expansion(mult = c(0.25, 0.05))
    ) +
    annotate(
      "text",
      x    = length(hc$order) * 0.98,
      y    = dc * 1.05,
      label = paste0("Mojena k=", mojena_k,
                     " | dc=", round(dc, 2)),
      size  = 3, hjust = 1, colour = "#D85A30"
    ) +
    labs(
      title    = titulo,
      subtitle = paste0("CCC = ", round(ccc, 3),
                        " | ", n_grp, " grupos | ",
                        metodo_dist),
      x = NULL, y = "distância de fusão"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "top"
    )
}

p_upgma <- plot_dendro(hc_upgma,
                       "Dendrograma — UPGMA",
                       ccc_upgma, dc_mojena, n_grupos)

p_ward  <- plot_dendro(hc_ward,
                       "Dendrograma — Ward.D2",
                       ccc_ward, dc_mojena, n_grupos)

ggsave(paste0(path_fig, "05a_dendrograma_upgma.png"),
       p_upgma, width = 14, height = 7, dpi = 300)
ggsave(paste0(path_fig, "05b_dendrograma_ward.png"),
       p_ward,  width = 14, height = 7, dpi = 300)

# PCA biplot
escala_div <- scale_shape_manual(
  values = setNames(1:n_grupos,
                    paste0("Grupo ", 1:n_grupos)),
  name   = "grupo divergência"
)

p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, colour = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "gray70", linewidth = 0.3) +
  geom_point(aes(colour = grupo, shape = grupo_div.y),
             size = 2.5, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = genotipo, colour = grupo),
    size = 2.5, max.overlaps = 20, show.legend = FALSE
  ) +
  # vetores dos loadings
  geom_segment(
    data = df_loadings,
    aes(x = 0, y = 0,
        xend = PC1 * max(abs(df_pca$PC1)) * 0.8,
        yend = PC2 * max(abs(df_pca$PC2)) * 0.8),
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    colour    = "gray20",
    linewidth = 0.5
  ) +
  geom_label(
    data = df_loadings,
    aes(x = PC1 * max(abs(df_pca$PC1)) * 0.85,
        y = PC2 * max(abs(df_pca$PC2)) * 0.85,
        label = variavel_lab),
    size = 3, colour = "gray20",
    fill = alpha("white", 0.8),
    label.size = 0
  ) +
  escala_grupo +
  escala_div +
  labs(
    x        = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    title    = "PCA — divergência genética para rendimento",
    subtitle = paste0("BLUPs bayesianos | variáveis: ",
                      paste(vars_div, collapse = ", "))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave(paste0(path_fig, "05c_pca_biplot.png"),
       p_pca, width = 10, height = 8, dpi = 300)

# contribuição das variáveis
p_contrib <- ggplot(contrib_var,
                    aes(x = reorder(variavel, prop_var),
                        y = prop_var,
                        fill = variavel)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(prop_var, 1), "%")),
            hjust = -0.2, size = 3.5) +
  scale_fill_manual(
    values = c(per_grao  = "#185FA5",
               mcm_mgb   = "#1D9E75",
               mcm_saca  = "#D85A30")
  ) +
  scale_y_continuous(
    limits = c(0, max(contrib_var$prop_var) * 1.2),
    labels = scales::percent_format(scale = 1)
  ) +
  coord_flip() +
  labs(
    x     = NULL,
    y     = "contribuição para variância total (%)",
    title = "Contribuição das variáveis para a divergência genética"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(paste0(path_fig, "05d_contribuicao_variaveis.png"),
       p_contrib, width = 7, height = 4, dpi = 300)

# heatmap de distâncias
dist_df <- as.matrix(dist_mat) |>
  as_tibble(rownames = "gen_i") |>
  pivot_longer(-gen_i, names_to = "gen_j",
               values_to = "distancia") |>
  left_join(grupo_ref |> rename(gen_i = genotipo,
                                grupo_i = grupo),
            by = "gen_i") |>
  left_join(grupo_ref |> rename(gen_j = genotipo,
                                grupo_j = grupo),
            by = "gen_j") |>
  mutate(
    gen_i = factor(gen_i, levels = gen_names[hc_melhor$order]),
    gen_j = factor(gen_j, levels = gen_names[hc_melhor$order])
  )

p_heat_dist <- ggplot(dist_df,
                      aes(x = gen_j, y = gen_i,
                          fill = distancia)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("#042C53", "#185FA5", "#1D9E75",
                "#D85A30", "#993C1D"),
    name    = metodo_dist
  ) +
  labs(
    x = NULL, y = NULL,
    title    = paste0("Matriz de distâncias — ", metodo_dist),
    subtitle = "Genótipos ordenados pelo dendrograma"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1,
                                    vjust = 0.5, size = 6),
    axis.text.y      = element_text(size = 6),
    legend.position  = "right",
    panel.grid       = element_blank()
  )

ggsave(paste0(path_fig, "05e_heatmap_distancias.png"),
       p_heat_dist, width = 12, height = 11, dpi = 300)

# ── 8. análise intra-grupo ─────────────────────────────────────────────────────
cat("\n── Divergência intra-grupo ──\n")

for (grp in c("Conilon", "Robusta")) {
  
  gen_grp <- grupo_ref |>
    filter(grupo == grp) |>
    pull(genotipo)
  
  mat_grp <- mat_num[rownames(mat_num) %in% gen_grp, ]
  
  cat("\n", grp, "—", nrow(mat_grp), "genótipos\n")
  
  if (nrow(mat_grp) < 4) {
    cat("  Poucos genótipos para análise separada\n")
    next
  }
  
  S_grp   <- cov(mat_grp)
  det_grp <- det(S_grp)
  
  if (det_grp > 1e-10 && kappa(S_grp) < 1000) {
    S_inv_grp <- solve(S_grp)
    dist_grp  <- as.dist(
      outer(seq_len(nrow(mat_grp)), seq_len(nrow(mat_grp)),
            FUN = Vectorize(function(i, j) {
              d <- mat_grp[i, ] - mat_grp[j, ]
              sqrt(t(d) %*% S_inv_grp %*% d)
            }))
    )
    cat("  Distância: Mahalanobis\n")
  } else {
    dist_grp <- dist(scale(mat_grp), method = "euclidean")
    cat("  Distância: Euclidiana padronizada (matriz singular)\n")
  }
  
  hc_grp  <- hclust(dist_grp, method = "average")
  ccc_grp <- cor(dist_grp, cophenetic(hc_grp))
  dc_grp  <- mean(hc_grp$height) + mojena_k * sd(hc_grp$height)
  ng_grp  <- length(unique(cutree(hc_grp, h = dc_grp)))
  
  cat("  CCC:", round(ccc_grp, 3), "\n")
  cat("  Grupos (Mojena):", ng_grp, "\n")
}

# ── 9. salvar tabelas ──────────────────────────────────────────────────────────
write_xlsx(
  list(
    grupos_divergencia  = df_grupos |>
      left_join(
        as_tibble(mat_num, rownames = "genotipo"),
        by = "genotipo"
      ),
    contribuicao_vars   = contrib_var,
    pca_scores          = df_pca,
    pca_loadings        = df_loadings,
    resumo_metodos      = tibble(
      metodo_distancia  = metodo_dist,
      ccc_upgma         = ccc_upgma,
      ccc_ward          = ccc_ward,
      metodo_preferido  = nome_melhor,
      mojena_k          = mojena_k,
      dc_mojena         = dc_mojena,
      n_grupos          = n_grupos
    )
  ),
  paste0(path_tbl, "05_diversity_results.xlsx")
)

cat("\n05_diversity.R concluído.\n")
cat("Figuras salvas em outputs/figures/\n")
cat("Tabelas salvas em outputs/tables/05_diversity_results.xlsx\n")








# divergência completa com per_grao e mcm_saca apenas ---------------------
# ── refazer divergência com per_grao + mcm_saca (Mahalanobis) ─────────────────
vars_div2 <- c("per_grao", "mcm_saca")

mat_num2   <- mat_num[, vars_div2]
S2         <- cov(mat_num2)
S2_inv     <- solve(S2)

cat("Matriz de covariância:\n")
print(round(S2, 4))
cat("Kappa:", round(kappa(S2), 2), "\n\n")

# ── distância de Mahalanobis ───────────────────────────────────────────────────
dist_mah <- as.dist(
  outer(seq_len(nrow(mat_num2)), seq_len(nrow(mat_num2)),
        FUN = Vectorize(function(i, j) {
          d <- mat_num2[i, ] - mat_num2[j, ]
          sqrt(t(d) %*% S2_inv %*% d)
        }))
)
attr(dist_mah, "Labels") <- rownames(mat_num2)

# ── agrupamento ────────────────────────────────────────────────────────────────
hc_upgma2 <- hclust(dist_mah, method = "average")
hc_ward2  <- hclust(dist_mah, method = "ward.D2")

# ── CCC ───────────────────────────────────────────────────────────────────────
ccc_upgma2 <- cor(dist_mah, cophenetic(hc_upgma2))
ccc_ward2  <- cor(dist_mah, cophenetic(hc_ward2))

cat("CCC UPGMA:", round(ccc_upgma2, 4), "\n")
cat("CCC Ward.D2:", round(ccc_ward2, 4), "\n")
cat("Método preferido:", if_else(ccc_upgma2 >= ccc_ward2,
                                 "UPGMA", "Ward.D2"), "\n\n")

# ── Mojena — dc calculado separado por método ─────────────────────────────────
mojena_k <- 1.25

dc_upgma2 <- mean(hc_upgma2$height) + mojena_k * sd(hc_upgma2$height)
dc_ward2  <- mean(hc_ward2$height)  + mojena_k * sd(hc_ward2$height)

n_grupos_upgma2 <- length(unique(cutree(hc_upgma2, h = dc_upgma2)))
n_grupos_ward2  <- length(unique(cutree(hc_ward2,  h = dc_ward2)))

cat("── Mojena UPGMA — dc:", round(dc_upgma2, 4),
    "| grupos:", n_grupos_upgma2, "\n")
cat("── Mojena Ward   — dc:", round(dc_ward2, 4),
    "| grupos:", n_grupos_ward2, "\n\n")

# usar método com melhor CCC
if (ccc_upgma2 >= ccc_ward2) {
  hc_final2   <- hc_upgma2
  dc_final2   <- dc_upgma2
  n_final2    <- n_grupos_upgma2
  nome_final2 <- "UPGMA"
  ccc_final2  <- ccc_upgma2
} else {
  hc_final2   <- hc_ward2
  dc_final2   <- dc_ward2
  n_final2    <- n_grupos_ward2
  nome_final2 <- "Ward.D2"
  ccc_final2  <- ccc_ward2
}

cat("Método final:", nome_final2, "\n")
cat("CCC final:", round(ccc_final2, 4), "\n")
cat("Grupos:", n_final2, "\n\n")

# ── grupos por genótipo ────────────────────────────────────────────────────────
grupos_gen2 <- cutree(hc_final2, h = dc_final2)

df_grupos2 <- tibble(
  genotipo  = names(grupos_gen2),
  grupo_div = paste0("Grupo ", grupos_gen2)
) |>
  left_join(grupo_ref, by = "genotipo")

cat("── Composição dos grupos ──\n")
df_grupos2 |>
  count(grupo_div, grupo) |>
  arrange(grupo_div, grupo) |>
  print()

# ── PCA corrigida — contribuição pelos loadings ────────────────────────────────
pca2 <- prcomp(mat_num2, scale. = TRUE, center = TRUE)
var_exp2 <- summary(pca2)$importance[2, ] * 100

cat("\n── PCA — variância explicada ──\n")
print(round(var_exp2, 2))

cat("\n── Loadings PC1 e PC2 ──\n")
print(round(pca2$rotation, 4))

# contribuição correta — via loadings^2 na PC1
contrib_corr <- tibble(
  variavel    = vars_div2,
  loading_pc1 = pca2$rotation[, 1],
  contrib_pc1 = pca2$rotation[, 1]^2 / sum(pca2$rotation[, 1]^2) * 100,
  loading_pc2 = pca2$rotation[, 2],
  contrib_pc2 = pca2$rotation[, 2]^2 / sum(pca2$rotation[, 2]^2) * 100
)

cat("\n── Contribuição das variáveis (via loadings) ──\n")
print(contrib_corr)

# ── dendrograma final ──────────────────────────────────────────────────────────
dend_data2 <- ggdendro::dendro_data(hc_final2, type = "rectangle")

labels_df2 <- dend_data2$labels |>
  left_join(grupo_ref, by = c("label" = "genotipo"))

p_dend_final <- ggplot() +
  geom_segment(
    data = dend_data2$segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.35, colour = "gray30"
  ) +
  geom_hline(
    yintercept = dc_final2,
    linetype   = "dashed",
    colour     = "#D85A30",
    linewidth  = 0.5
  ) +
  geom_text(
    data = labels_df2,
    aes(x = x,
        y = -0.02 * max(dend_data2$segments$y),
        label = label, colour = grupo),
    size = 2.2, hjust = 1, angle = 90
  ) +
  scale_colour_manual(
    values = c("Conilon" = "#185FA5", "Robusta" = "#D85A30"),
    name   = "grupo botânico"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
  annotate(
    "text",
    x     = length(hc_final2$order) * 0.98,
    y     = dc_final2 * 1.05,
    label = paste0("Mojena k=", mojena_k,
                   " | dc=", round(dc_final2, 2)),
    size  = 3, hjust = 1, colour = "#D85A30"
  ) +
  labs(
    title    = paste0("Dendrograma — ", nome_final2),
    subtitle = paste0("CCC = ", round(ccc_final2, 3),
                      " | ", n_final2, " grupos",
                      " | Distância de Mahalanobis",
                      " | Variáveis: per_grao + mcm_saca"),
    x = NULL, y = "distância de fusão"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top"
  )

ggsave(paste0(path_fig, "05f_dendrograma_final_mahalanobis.png"),
       p_dend_final, width = 14, height = 7, dpi = 300)

# ── PCA biplot final ───────────────────────────────────────────────────────────
df_pca2 <- as_tibble(pca2$x[, 1:2]) |>
  mutate(genotipo = rownames(mat_num2)) |>
  left_join(grupo_ref, by = "genotipo") |>
  left_join(df_grupos2 |> select(genotipo, grupo_div),
            by = "genotipo")

df_load2 <- as_tibble(pca2$rotation[, 1:2],
                      rownames = "variavel") |>
  mutate(
    variavel_lab = c(per_grao  = "% grão",
                     mcm_saca  = "MCM/saca")[variavel]
  )

escala_div2 <- scale_shape_manual(
  values = setNames(c(16, 17, 15, 3, 4, 8)[seq_len(n_final2)],
                    paste0("Grupo ", seq_len(n_final2))),
  name   = "grupo divergência"
)

p_pca2 <- ggplot(df_pca2, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, colour = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "gray70", linewidth = 0.3) +
  geom_point(aes(colour = grupo, shape = grupo_div),
             size = 2.8, alpha = 0.85) +
  ggrepel::geom_text_repel(
    aes(label = genotipo, colour = grupo),
    size = 2.5, max.overlaps = 25, show.legend = FALSE
  ) +
  geom_segment(
    data = df_load2,
    aes(x = 0, y = 0,
        xend = PC1 * max(abs(df_pca2$PC1)) * 0.75,
        yend = PC2 * max(abs(df_pca2$PC2)) * 0.75),
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    colour    = "gray20", linewidth = 0.5
  ) +
  geom_label(
    data = df_load2,
    aes(x = PC1 * max(abs(df_pca2$PC1)) * 0.82,
        y = PC2 * max(abs(df_pca2$PC2)) * 0.82,
        label = variavel_lab),
    size = 3.2, colour = "gray20",
    fill = alpha("white", 0.8), linewidth = 0
  ) +
  scale_colour_manual(
    values = c("Conilon" = "#185FA5", "Robusta" = "#D85A30"),
    name   = "grupo botânico"
  ) +
  escala_div2 +
  labs(
    x        = paste0("PC1 (", round(var_exp2[1], 1), "%)"),
    y        = paste0("PC2 (", round(var_exp2[2], 1), "%)"),
    title    = "PCA — divergência genética para rendimento",
    subtitle = "BLUPs bayesianos | Mahalanobis | per_grao + mcm_saca"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave(paste0(path_fig, "05g_pca_biplot_final.png"),
       p_pca2, width = 10, height = 8, dpi = 300)

# ── atualizar xlsx ─────────────────────────────────────────────────────────────
write_xlsx(
  list(
    grupos_divergencia = df_grupos2 |>
      left_join(as_tibble(mat_num2, rownames = "genotipo"),
                by = "genotipo"),
    contribuicao_vars  = contrib_corr,
    pca_scores         = df_pca2,
    pca_loadings       = df_load2,
    resumo_metodos     = tibble(
      metodo_distancia = "Mahalanobis",
      variaveis        = paste(vars_div2, collapse = " + "),
      kappa_S          = kappa(S2),
      ccc_upgma        = ccc_upgma2,
      ccc_ward         = ccc_ward2,
      metodo_preferido = nome_final2,
      ccc_final        = ccc_final2,
      mojena_k         = mojena_k,
      dc_final         = dc_final2,
      n_grupos         = n_final2
    )
  ),
  paste0(path_tbl, "05_diversity_results_final.xlsx")
)

cat("\nConcluído. Figuras 05f e 05g salvas.\n")
cat("Tabela atualizada em 05_diversity_results_final.xlsx\n")


