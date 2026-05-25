# 01_load_explore.R
# Loading, structure verification, and descriptive exploration
# Outputs: data/dados_clean.rds, outputs/tables/descritivas.csv
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
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

# ── 1. reading ─────────────────────────────────────────────────────────────────
raw <- read_excel(
  path_data,
  col_types = c("text", "text", "text", "text",
                "numeric", "numeric", "numeric",
                "numeric", "numeric", "numeric")
)

# ── 2. cleaning and typing ───────────────────────────────────────────────────────
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

# ── 3. balancing verification ───────────────────────────────────────────
cat("\n── Genotype × Year balance ──\n")
print(table(dados$genotipo, dados$ano))

cat("\n── Observations by group ──\n")
dados |>
  count(grupo, ano) |>
  print()

# ── 4. descriptive statistics ───────────────────────────────────────────────
vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca", "vcm_mcm")

cat("\n── Global descriptives ──\n")
dados |>
  select(all_of(vars_resp)) |>
  skim() |>
  print()

cat("\n── Descriptives by year ──\n")
desc_ano <- dados |>
  group_by(ano) |>
  summarise(
    across(
      all_of(vars_resp),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE),
        min  = \(x) min(x, na.rm = TRUE),
        max  = \(x) max(x, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    )
  )
print(desc_ano)

# ── 5. correlation between years by variable ─────────────────────────────────────
# Mean per genotype × year → wide → correlation between years
# This indicates how consistent the genotype ranking is between the two years
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

cat("\n── Pearson correlation between years (mean per genotype) ──\n")
cor_tab <- sapply(vars_resp, function(v) {
  x <- cor_anos[[paste0(v, "__2023")]]
  y <- cor_anos[[paste0(v, "__2024")]]
  cor(x, y, use = "complete.obs")
})
print(round(cor_tab, 3))

# ── 6. exploratory plots ────────────────────────────────────────────────────

# 6a. Distribution by year (density + boxplot)
# Defining display labels for variables
var_labels <- c(
  per_grao  = "Grain proportion (% grain)",
  per_palha = "Husk proportion (% husk)",
  mcm_mgb   = "Fruit fresh mass per grain mass (FWM/GW)",
  mcm_saca  = "Fruit fresh mass per bag (FWM/bag)",
  vcm_saca  = "Fruit volume per bag (FVol/bag)",
  vcm_mcm   = "Fruit volume-to-fresh mass ratio (FVol/FWM)"
)

p_dist <- dados |>
  pivot_longer(all_of(vars_resp), names_to = "variavel", values_to = "valor") |>
  mutate(variavel = recode(variavel, !!!var_labels)) |>
  ggplot(aes(x = valor, fill = ano, colour = ano)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  facet_wrap(~variavel, scales = "free", ncol = 3) +
  scale_fill_manual(values  = c("2023" = "#378ADD", "2024" = "#D85A30")) +
  scale_colour_manual(values = c("2023" = "#185FA5", "2024" = "#993C1D")) +
  labs(x = NULL, y = "Density", fill = "Year", colour = "Year") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text  = element_text(size = 9, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(paste0(path_out_fig, "01a_distribuicao_por_ano.png"),
       p_dist, width = 10, height = 7, dpi = 300)

# 6b. Mean per genotype in 2023 vs. 2024 (consistency scatter plot)
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
    x = "FWM/bag 2023 (kg)",
    y = "FWM/bag 2024 (kg)",
    title = "Genotype ranking consistency between years — FWM/bag"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(paste0(path_out_fig, "01b_scatter_2023_2024_mcm_saca.png"),
       p_scatter, width = 7, height = 6, dpi = 300)

# ── 7. saving cleaned object ────────────────────────────────────────────────────
saveRDS(dados, path_out_rds)
cat("\nData saved in", path_out_rds, "\n")
cat("Dimensions:", nrow(dados), "×", ncol(dados), "\n")


# ── permanent plot (genotype × block — same individual in both years) ──────
dados <- dados |>
  mutate(parcela = interaction(genotipo, block, drop = TRUE))

# verification: each plot must have exactly 2 observations (one per year)
check_parcela <- dados |>
  count(parcela) |>
  filter(n != 2)

if (nrow(check_parcela) == 0) {
  cat("\nPermanent plot OK — 144 plots × 2 years\n")
} else {
  cat("\nWARNING — plots with irregular count:\n")
  print(check_parcela)
}

saveRDS(dados, path_out_rds)
