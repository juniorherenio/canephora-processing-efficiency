# 09_tabela3_publicacao.R
# Tabela 3 — Desempenho, estabilidade e recomendação de genótipos
#
# Conteúdo: para cada genótipo × variável (cinco caracteres retidos):
#   - BLUP médio na escala original (desvio da média geral)
#   - Ecovalência de Wricke relativa (%)
#   - Probabilidade de superioridade consistente a 20%
#   - Grupo de divergência genética (UPGMA Mahalanobis)
#   - Grupo botânico
#
# Estrutura de saída:
#   Tab 1 — "Table3_FWMbag"  : tabela completa para FWM/bag (caráter principal)
#   Tab 2 — "Table3_all"     : todos os cinco caracteres em formato longo
#   Tab 3 — "TableS1_complete": todos os genótipos × todos os caracteres
#                               para material suplementar
#
# Inserção no manuscrito:
#   Tabela 3 (principal) — após o parágrafo de resultados de estabilidade
#   (seção 3.3), antes da seção 3.4 (correlações genéticas).
#   TableS1 — material suplementar, citada no mesmo parágrafo.
#
# Formatação:
#   - BLUP: duas casas decimais
#   - Wricke (%): uma casa decimal
#   - Probabilidade: três casas decimais
#   - Ordenação: decrescente por prob_consistente dentro de cada variável

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(writexl)
library(openxlsx)

path_tbl <- "outputs/tables/"

# ── labels de variáveis (inglês, para o artigo) ────────────────────────────────
var_labels <- c(
  per_grao  = "% grain",
  per_palha = "% husk",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag"
)

vars_resp <- names(var_labels)

# ── carregar dados ─────────────────────────────────────────────────────────────
cat("Carregando dados...\n")

# estabilidade — ranking final com todas as colunas necessárias
df_rank <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                     sheet = "ranking_final") |>
  filter(variavel %in% vars_resp)

# probabilidade consistente a 20%
df_prob <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                     sheet = "prob_consistente") |>
  filter(variavel %in% vars_resp, intensidade == 0.20) |>
  dplyr::select(variavel, genotipo, prob_consistente)

# grupos de divergência — da análise de diversidade
df_grupos <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                       sheet = "grupos_divergencia") |>
  dplyr::select(genotipo, grupo_div) |>
  mutate(genotipo = as.character(genotipo))

cat("  Dados carregados.\n\n")

# ── montar tabela base ─────────────────────────────────────────────────────────
cat("Montando tabela base...\n")

df_base <- df_rank |>
  dplyr::select(variavel, genotipo, grupo,
         blup_orig_med, wricke_rel,
         rank_desempenho, rank_estab, rank_composto) |>
  left_join(df_prob, by = c("variavel", "genotipo")) |>
  left_join(df_grupos, by = "genotipo") |>
  mutate(
    # média geral para cada variável (para calcular valor absoluto do BLUP)
    variavel_label = var_labels[variavel],
    # flag de recomendação: prob > 0.25 a 20% de intensidade
    recomendado    = prob_consistente > 0.25
  ) |>
  arrange(variavel, desc(prob_consistente))

cat("  Dimensões:", nrow(df_base), "×", ncol(df_base), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# TABELA 3 — FWM/bag (caráter principal) — formato publicação
# ══════════════════════════════════════════════════════════════════════════════
cat("Gerando Table 3 — FWM/bag...\n")

tab3_fwmbag <- df_base |>
  filter(variavel == "mcm_saca") |>
  transmute(
    Genotype             = genotipo,
    `Botanical group`    = grupo,
    `Divergence group`   = grupo_div,
    `BLUP (kg saca⁻¹)`  = round(blup_orig_med, 2),
    `Wricke (% Wi)`      = round(wricke_rel,    1),
    `Pr(consistent superiority)` = round(prob_consistente, 3),
    `Performance rank`   = rank_desempenho,
    `Stability rank`     = rank_estab,
    `Composite rank`     = rank_composto
  ) |>
  arrange(`Composite rank`)

cat("  Genótipos em Table 3:", nrow(tab3_fwmbag), "\n")
cat("  Top 10:\n")
print(tab3_fwmbag |> head(10))

# ══════════════════════════════════════════════════════════════════════════════
# TABELA SUPLEMENTAR S5 — todos os genótipos × todos os caracteres
# ══════════════════════════════════════════════════════════════════════════════
cat("\nGerando Table S5 — todos os caracteres...\n")

tabS5 <- df_base |>
  transmute(
    Trait                        = var_labels[variavel],
    Genotype                     = genotipo,
    `Botanical group`            = grupo,
    `Divergence group`           = grupo_div,
    `BLUP (original scale)`      = round(blup_orig_med, 3),
    `Wricke ecovalence (% Wi)`   = round(wricke_rel,    1),
    `Pr(consistent superiority)` = round(prob_consistente, 3),
    `Performance rank`           = rank_desempenho,
    `Stability rank`             = rank_estab,
    `Composite rank`             = rank_composto
  ) |>
  arrange(Trait, `Composite rank`)

cat("  Linhas em Table S5:", nrow(tabS5), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# SALVAR — xlsx com formatação
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSalvando xlsx formatado...\n")

wb <- createWorkbook()

# ── estilos ────────────────────────────────────────────────────────────────────
st_header <- createStyle(
  fontName   = "Arial", fontSize = 10, fontColour = "white",
  fgFill     = "#2C3E50", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold", wrapText = TRUE,
  border = "TopBottomLeftRight", borderColour = "#FFFFFF"
)

st_body <- createStyle(
  fontName = "Arial", fontSize = 9,
  halign   = "CENTER", valign = "CENTER",
  border   = "TopBottomLeftRight", borderColour = "#CCCCCC"
)

st_body_left <- createStyle(
  fontName = "Arial", fontSize = 9,
  halign   = "LEFT", valign = "CENTER",
  border   = "TopBottomLeftRight", borderColour = "#CCCCCC"
)

st_highlight <- createStyle(
  fontName = "Arial", fontSize = 9, fontColour = "#FFFFFF",
  fgFill   = "#185FA5", halign = "CENTER", valign = "CENTER",
  textDecoration = "bold",
  border   = "TopBottomLeftRight", borderColour = "#CCCCCC"
)

st_zebra <- createStyle(
  fontName = "Arial", fontSize = 9,
  fgFill   = "#F5F7FA",
  halign   = "CENTER", valign = "CENTER",
  border   = "TopBottomLeftRight", borderColour = "#CCCCCC"
)

st_zebra_left <- createStyle(
  fontName = "Arial", fontSize = 9,
  fgFill   = "#F5F7FA",
  halign   = "LEFT", valign = "CENTER",
  border   = "TopBottomLeftRight", borderColour = "#CCCCCC"
)

# ── aba Table 3 — FWM/bag ──────────────────────────────────────────────────────
addWorksheet(wb, "Table3_FWMbag")
writeData(wb, "Table3_FWMbag", tab3_fwmbag, startRow = 3, startCol = 1,
          headerStyle = st_header)

# título e legenda
writeData(wb, "Table3_FWMbag",
          "Table 3. Mean performance, temporal stability and probability of consistent superiority for FWM/bag in 48 Coffea canephora genotypes evaluated over two years (2023–2024) in Jaguaré, Espírito Santo, Brazil.",
          startRow = 1, startCol = 1)
mergeCells(wb, "Table3_FWMbag", cols = 1:9, rows = 1)

writeData(wb, "Table3_FWMbag",
          "BLUP: best linear unbiased prediction (deviation from overall mean, kg bag⁻¹). Wricke (%Wi): relative Wricke ecovalence (lower = more stable). Pr(consistent superiority): Bayesian probability of maintaining superior rank in both years at 20% dplyr::selection intensity. Divergence group: UPGMA cluster based on Mahalanobis distance (% grain + FWM/bag BLUPs). Ranks: 1 = best.",
          startRow = 2, startCol = 1)
mergeCells(wb, "Table3_FWMbag", cols = 1:9, rows = 2)

# estilo cabeçalho
addStyle(wb, "Table3_FWMbag", st_header,
         rows = 3, cols = 1:ncol(tab3_fwmbag), gridExpand = TRUE)

# corpo — zebra striping + destacar prob > 0.50
n_rows <- nrow(tab3_fwmbag)
for (i in seq_len(n_rows)) {
  row_i   <- i + 3   # offset pelo título
  prob_i  <- tab3_fwmbag$`Pr(consistent superiority)`[i]
  
  if (prob_i >= 0.50) {
    # destaque azul para genótipos com prob >= 50%
    addStyle(wb, "Table3_FWMbag", st_highlight,
             rows = row_i, cols = 1:ncol(tab3_fwmbag), gridExpand = TRUE)
  } else if (i %% 2 == 0) {
    addStyle(wb, "Table3_FWMbag", st_zebra,
             rows = row_i, cols = 3:ncol(tab3_fwmbag), gridExpand = TRUE)
    addStyle(wb, "Table3_FWMbag", st_zebra_left,
             rows = row_i, cols = 1:2, gridExpand = TRUE)
  } else {
    addStyle(wb, "Table3_FWMbag", st_body,
             rows = row_i, cols = 3:ncol(tab3_fwmbag), gridExpand = TRUE)
    addStyle(wb, "Table3_FWMbag", st_body_left,
             rows = row_i, cols = 1:2, gridExpand = TRUE)
  }
}

# larguras de coluna
setColWidths(wb, "Table3_FWMbag",
             cols = 1:ncol(tab3_fwmbag),
             widths = c(14, 14, 16, 13, 11, 20, 15, 13, 13))
setRowHeights(wb, "Table3_FWMbag", rows = 1, heights = 40)
setRowHeights(wb, "Table3_FWMbag", rows = 2, heights = 55)
setRowHeights(wb, "Table3_FWMbag", rows = 3, heights = 30)

# congelar painel no cabeçalho
freezePane(wb, "Table3_FWMbag", firstActiveRow = 4)

# ── aba Table S5 — todos os caracteres ────────────────────────────────────────
addWorksheet(wb, "TableS5_all_traits")
writeData(wb, "TableS5_all_traits", tabS5, startRow = 3, startCol = 1,
          headerStyle = st_header)

writeData(wb, "TableS5_all_traits",
          "Table S5. Mean performance, temporal stability and probability of consistent superiority for five processing efficiency traits in 48 Coffea canephora genotypes.",
          startRow = 1, startCol = 1)
mergeCells(wb, "TableS5_all_traits", cols = 1:10, rows = 1)

writeData(wb, "TableS5_all_traits",
          "Trait abbreviations: % grain = ratio of grain mass to total dry fruit mass; % husk = ratio of husk mass to total dry fruit mass; FWM/GW = fresh fruit mass per unit of processed grain mass; FWM/bag = fresh fruit mass required per 60-kg processed bag; FVol/bag = fresh fruit volume required per 60-kg processed bag. Genotypes ordered by composite rank within each trait.",
          startRow = 2, startCol = 1)
mergeCells(wb, "TableS5_all_traits", cols = 1:10, rows = 2)

addStyle(wb, "TableS5_all_traits", st_header,
         rows = 3, cols = 1:ncol(tabS5), gridExpand = TRUE)

# corpo com zebra por bloco de trait
traits_order <- unique(tabS5$Trait)
row_offset   <- 4
for (trait in traits_order) {
  rows_trait <- which(tabS5$Trait == trait)
  for (k in seq_along(rows_trait)) {
    row_i <- row_offset + rows_trait[k] - 1
    if (k %% 2 == 0) {
      addStyle(wb, "TableS5_all_traits", st_zebra,
               rows = row_i, cols = 3:ncol(tabS5), gridExpand = TRUE)
      addStyle(wb, "TableS5_all_traits", st_zebra_left,
               rows = row_i, cols = 1:2, gridExpand = TRUE)
    } else {
      addStyle(wb, "TableS5_all_traits", st_body,
               rows = row_i, cols = 3:ncol(tabS5), gridExpand = TRUE)
      addStyle(wb, "TableS5_all_traits", st_body_left,
               rows = row_i, cols = 1:2, gridExpand = TRUE)
    }
  }
}

setColWidths(wb, "TableS5_all_traits",
             cols = 1:ncol(tabS5),
             widths = c(12, 14, 14, 16, 16, 16, 20, 15, 13, 13))
setRowHeights(wb, "TableS5_all_traits", rows = 1, heights = 40)
setRowHeights(wb, "TableS5_all_traits", rows = 2, heights = 55)
setRowHeights(wb, "TableS5_all_traits", rows = 3, heights = 30)
freezePane(wb, "TableS5_all_traits", firstActiveRow = 4)

# ── salvar ─────────────────────────────────────────────────────────────────────
out_path <- paste0(path_tbl, "09_Table3_publicacao.xlsx")
saveWorkbook(wb, out_path, overwrite = TRUE)

cat("\nArquivo salvo:", out_path, "\n")
cat("\nResumo Table 3 (FWM/bag):\n")
cat("  Total de genótipos:", nrow(tab3_fwmbag), "\n")
cat("  Com prob >= 0.50:  ",
    sum(tab3_fwmbag$`Pr(consistent superiority)` >= 0.50), "\n")
cat("  Com prob >= 0.25:  ",
    sum(tab3_fwmbag$`Pr(consistent superiority)` >= 0.25), "\n")

cat("\nResumo Table S5 (todos os caracteres):\n")
cat("  Total de linhas:", nrow(tabS5), "\n")
cat("  Traits:", paste(unique(tabS5$Trait), collapse = ", "), "\n")

cat("\n━━ 09_tabela3_publicacao.R concluído ━━\n")