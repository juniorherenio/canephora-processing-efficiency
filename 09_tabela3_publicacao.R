# 09_table3_publication.R
# Table 3 — Performance, stability, and genotype recommendation
#
# Content: for each genotype × trait (five retained traits):
#   - Mean BLUP on the original scale (deviation from overall mean)
#   - Relative Wricke's ecovalence (%)
#   - Probability of consistent superiority at 20%
#   - Genetic divergence group (UPGMA Mahalanobis)
#   - Botanical group
#
# Output structure:
#   Tab 1 — "Table3_FWMbag"  : full table for FWM/bag (main trait)
#   Tab 2 — "Table3_all"     : all five traits in long format
#   Tab 3 — "TableS5_all_traits": all genotypes × all traits
#                               for supplementary material
#
# Insertion in the manuscript:
#   Table 3 (main) — after the stability results paragraph
#   (section 3.3), before section 3.4 (genetic correlations).
#   Table S5 — supplementary material, cited in the same paragraph.
#
# Formatting:
#   - BLUP: two decimal places
#   - Wricke (%): one decimal place
#   - Probability: three decimal places
#   - Sorting: descending by prob_consistente within each trait
#
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(writexl)
library(openxlsx)

path_tbl <- "outputs/tables/"

# ── variable labels (English, for the article) ────────────────────────────────
var_labels <- c(
  per_grao  = "% grain",
  per_palha = "% husk",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag"
)

vars_resp <- names(var_labels)

# ── load data ─────────────────────────────────────────────────────────────
cat("Loading data...\n")

# stability — final ranking with all necessary columns
df_rank <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                     sheet = "ranking_final") |>
  filter(variavel %in% vars_resp)

# consistent probability at 20%
df_prob <- read_xlsx(paste0(path_tbl, "04_stability_results.xlsx"),
                     sheet = "prob_consistente") |>
  filter(variavel %in% vars_resp, intensidade == 0.20) |>
  dplyr::select(variavel, genotipo, prob_consistente)

# divergence groups — from diversity analysis
df_groups <- read_xlsx(paste0(path_tbl, "05_diversity_results_final.xlsx"),
                       sheet = "grupos_divergencia") |>
  dplyr::select(genotipo, grupo_div) |>
  mutate(genotipo = as.character(genotipo))

cat("  Data loaded.\n\n")

# ── assemble base table ─────────────────────────────────────────────────────────
cat("Assembling base table...\n")

df_base <- df_rank |>
  dplyr::select(variavel, genotipo, grupo,
         blup_orig_med, wricke_rel,
         rank_desempenho, rank_estab, rank_composto) |>
  left_join(df_prob, by = c("variavel", "genotipo")) |>
  left_join(df_groups, by = "genotipo") |>
  mutate(
    # overall mean for each trait (to calculate absolute BLUP value)
    variavel_label = var_labels[variavel],
    # recommendation flag: prob > 0.25 at 20% intensity
    recommended    = prob_consistente > 0.25
  ) |>
  arrange(variavel, desc(prob_consistente))

cat("  Dimensions:", nrow(df_base), "×", ncol(df_base), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE 3 — FWM/bag (main trait) — publication format
# ══════════════════════════════════════════════════════════════════════════════
cat("Generating Table 3 — FWM/bag...\n")

tab3_fwmbag <- df_base |>
  filter(variavel == "mcm_saca") |>
  transmute(
    Genotype                   = genotipo,
    `Botanical group`          = grupo,
    `Divergence group`         = grupo_div,
    `BLUP (kg bag⁻¹)`          = round(blup_orig_med, 2),
    `Wricke (% Wi)`            = round(wricke_rel,    1),
    `Pr(consistent superiority)` = round(prob_consistente, 3),
    `Performance rank`         = rank_desempenho,
    `Stability rank`           = rank_estab,
    `Composite rank`           = rank_composto
  ) |>
  arrange(`Composite rank`)

cat("  Genotypes in Table 3:", nrow(tab3_fwmbag), "\n")
cat("  Top 10:\n")
print(tab3_fwmbag |> head(10))

# ══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY TABLE S5 — all genotypes × all traits
# ══════════════════════════════════════════════════════════════════════════════
cat("\nGenerating Table S5 — all traits...\n")

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

cat("  Rows in Table S5:", nrow(tabS5), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# SAVE — xlsx with formatting
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSaving formatted xlsx...\n")

wb <- createWorkbook()

# ── styles ────────────────────────────────────────────────────────────────────
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

# ── Table 3 tab — FWM/bag ──────────────────────────────────────────────────────
addWorksheet(wb, "Table3_FWMbag")
writeData(wb, "Table3_FWMbag", tab3_fwmbag, startRow = 3, startCol = 1,
          headerStyle = st_header)

# title and legend
writeData(wb, "Table3_FWMbag",
          "Table 3. Mean performance, temporal stability and probability of consistent superiority for FWM/bag in 48 Coffea canephora genotypes evaluated over two years (2023–2024) in Jaguaré, Espírito Santo, Brazil.",
          startRow = 1, startCol = 1)
mergeCells(wb, "Table3_FWMbag", cols = 1:9, rows = 1)

writeData(wb, "Table3_FWMbag",
          "BLUP: best linear unbiased prediction (deviation from overall mean, kg bag⁻¹). Wricke (%Wi): relative Wricke ecovalence (lower = more stable). Pr(consistent superiority): Bayesian probability of maintaining superior rank in both years at 20% selection intensity. Divergence group: UPGMA cluster based on Mahalanobis distance (% grain + FWM/bag BLUPs). Ranks: 1 = best.",
          startRow = 2, startCol = 1)
mergeCells(wb, "Table3_FWMbag", cols = 1:9, rows = 2)

# header style
addStyle(wb, "Table3_FWMbag", st_header,
         rows = 3, cols = 1:ncol(tab3_fwmbag), gridExpand = TRUE)

# body — zebra striping + highlight prob > 0.50
n_rows <- nrow(tab3_fwmbag)
for (i in seq_len(n_rows)) {
  row_i  <- i + 3   # offset due to the title
  prob_i <- tab3_fwmbag$`Pr(consistent superiority)`[i]
  
  if (prob_i >= 0.50) {
    # blue highlight for genotypes with prob >= 50%
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

# column widths
setColWidths(wb, "Table3_FWMbag",
             cols = 1:ncol(tab3_fwmbag),
             widths = c(14, 14, 16, 13, 11, 20, 15, 13, 13))
setRowHeights(wb, "Table3_FWMbag", rows = 1, heights = 40)
setRowHeights(wb, "Table3_FWMbag", rows = 2, heights = 55)
setRowHeights(wb, "Table3_FWMbag", rows = 3, heights = 30)

# freeze pane at the header
freezePane(wb, "Table3_FWMbag", firstActiveRow = 4)

# ── Table S5 tab — all traits ────────────────────────────────────────
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

# body with zebra striping by trait block
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

# ── save ─────────────────────────────────────────────────────────────────────
out_path <- paste0(path_tbl, "09_Table3_publication.xlsx")
saveWorkbook(wb, out_path, overwrite = TRUE)

cat("\nFile saved:", out_path, "\n")
cat("\nSummary Table 3 (FWM/bag):\n")
cat("  Total genotypes:", nrow(tab3_fwmbag), "\n")
cat("  With prob >= 0.50:  ",
    sum(tab3_fwmbag$`Pr(consistent superiority)` >= 0.50), "\n")
cat("  With prob >= 0.25:  ",
    sum(tab3_fwmbag$`Pr(consistent superiority)` >= 0.25), "\n")

cat("\nSummary Table S5 (all traits):\n")
cat("  Total rows:", nrow(tabS5), "\n")
cat("  Traits:", paste(unique(tabS5$Trait), collapse = ", "), "\n")

cat("\n━━ 09_table3_publication.R completed ━━\n")
