# 01_load_explore.R
# Carregamento, verificação de estrutura e exploração descritiva
# Saídas: data/dados_clean.rds, outputs/tables/descritivas.csv

# ── pacotes ────────────────────────────────────────────────────────────────────
library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(skimr)
library(scales)

# ── paths ──────────────────────────────────────────────────────────────────────
path_data     <- "data/dados.xlsx"
path_out_rds  <- "data/dados_clean.rds"
path_out_fig  <- "outputs/figures/"
path_out_tbl  <- "outputs/tables/"

# ── 1. leitura ─────────────────────────────────────────────────────────────────
raw <- read_excel(
  path_data,
  col_types = c("text", "text", "text", "text",
                "numeric", "numeric", "numeric",
                "numeric", "numeric", "numeric")
)

# ── 2. limpeza e tipagem ───────────────────────────────────────────────────────
dados <- raw |>
  janitor::clean_names() |>
  mutate(
    genotipo = as.factor(genotipo),
    ano      = as.factor(ano),
    block    = as.factor(block),
    grupo    = if_else(
      genotipo %in% c("Peneirão", "Clementino", "Pirata", "AD1",
                      "K61", "CM1", "Z21", "L80"),
      "Conilon", "Robusta"
    ) |> as.factor()
  )

# ── 3. verificação de balanceamento ───────────────────────────────────────────
cat("\n── Balanceamento genótipo × ano ──\n")
print(table(dados$genotipo, dados$ano))

cat("\n── Observações por grupo ──\n")
dados |>
  count(grupo, ano) |>
  print()

# ── 4. estatísticas descritivas ───────────────────────────────────────────────
vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca", "vcm_mcm")

cat("\n── Descritivas globais ──\n")
dados |>
  select(all_of(vars_resp)) |>
  skim() |>
  print()

cat("\n── Descritivas por ano ──\n")
desc_ano <- dados |>
  group_by(ano) |>
  summarise(
    across(
      all_of(vars_resp),
      list(
        media = \(x) mean(x, na.rm = TRUE),
        dp    = \(x) sd(x, na.rm = TRUE),
        min   = \(x) min(x, na.rm = TRUE),
        max   = \(x) max(x, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )
  )
print(desc_ano)

# ── 5. correlação entre anos por variável ─────────────────────────────────────
# Média por genótipo × ano → wide → cor entre anos
# Isso indica quão consistente é a classificação dos genótipos entre os dois anos
cor_anos <- dados |>
  group_by(genotipo, ano) |>
  summarise(
    across(all_of(vars_resp), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = ano,
    values_from = all_of(vars_resp),
    names_glue  = "{.value}__{ano}"
  )

cat("\n── Correlação de Pearson entre anos (média por genótipo) ──\n")
cor_tab <- sapply(vars_resp, function(v) {
  x <- cor_anos[[paste0(v, "__2023")]]
  y <- cor_anos[[paste0(v, "__2024")]]
  cor(x, y, use = "complete.obs")
})
print(round(cor_tab, 3))

# ── 6. plots exploratórios ────────────────────────────────────────────────────

# 6a. Distribuição por ano (densidade + boxplot)
p_dist <- dados |>
  pivot_longer(all_of(vars_resp), names_to = "variavel", values_to = "valor") |>
  ggplot(aes(x = valor, fill = ano, colour = ano)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~variavel, scales = "free", ncol = 3) +
  scale_fill_manual(values  = c("2023" = "#378ADD", "2024" = "#D85A30")) +
  scale_colour_manual(values = c("2023" = "#185FA5", "2024" = "#993C1D")) +
  labs(x = NULL, y = "densidade", fill = "ano", colour = "ano") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text  = element_text(size = 9, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(paste0(path_out_fig, "01a_distribuicao_por_ano.png"),
       p_dist, width = 10, height = 7, dpi = 300)

# 6b. Média por genótipo em 2023 vs. 2024 (scatter de consistência)
p_scatter <- cor_anos |>
  select(genotipo, starts_with("mcm_saca")) |>
  ggplot(aes(x = mcm_saca__2023, y = mcm_saca__2024,
             label = genotipo)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "gray60", linewidth = 0.4) +
  geom_point(size = 2.5, colour = "#378ADD", alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.8, colour = "gray40",
                           max.overlaps = 15) +
  labs(
    x = "MCM/saca 2023 (kg)",
    y = "MCM/saca 2024 (kg)",
    title = "Consistência do ranking entre anos — MCM/saca"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(paste0(path_out_fig, "01b_scatter_2023_2024_mcm_saca.png"),
       p_scatter, width = 7, height = 6, dpi = 300)

# ── 7. salvar objeto limpo ────────────────────────────────────────────────────
saveRDS(dados, path_out_rds)
cat("\nDados salvos em", path_out_rds, "\n")
cat("Dimensões:", nrow(dados), "×", ncol(dados), "\n")


# ── parcela permanente (genotipo × bloco — mesmo indivíduo nos dois anos) ──────
dados <- dados |>
  mutate(parcela = interaction(genotipo, block, drop = TRUE))

# verificação: cada parcela deve ter exatamente 2 observações (uma por ano)
check_parcela <- dados |>
  count(parcela) |>
  filter(n != 2)

if (nrow(check_parcela) == 0) {
  cat("\nParcela permanente OK — 144 parcelas × 2 anos\n")
} else {
  cat("\nATENÇÃO — parcelas com contagem irregular:\n")
  print(check_parcela)
}

saveRDS(dados, path_out_rds)
