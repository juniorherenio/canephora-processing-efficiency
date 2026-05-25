# 04_stability.R
# Stability analysis: Wricke's ecovalence and Bayesian probability of 
# consistent superiority
#
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readxl)
library(ggplot2)
library(patchwork)
library(writexl)
library(brms)      
library(posterior) 

# ── data and paths ──────────────────────────────────────────────────────────────
dados   <- readRDS("data/dados_clean.rds")
modelos <- readRDS("data/modelos_finais.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# target variables — excluding vcm_mcm (lacks genetic structure)
vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca")

# display labels
var_labels <- c(
  per_grao  = "% grain",
  per_palha = "% husk",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag"
)

# desirable direction
# "alto" (high) = higher is better | "baixo" (low) = lower is better
var_direcao <- c(
  per_grao  = "alto",  
  per_palha = "baixo", 
  mcm_mgb   = "baixo", 
  mcm_saca  = "baixo", 
  vcm_saca  = "baixo"  
)

# ── load pre-calculated BLUPs ──────────────────────────────────────────────────
blups_G <- read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "blups_genotipo"
) |>
  filter(variavel %in% vars_resp)

blups_GxA <- read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "blups_gxa"
) |>
  filter(variavel %in% vars_resp)

# ── reconstruct genotypic values by year ───────────────────────────────────────
media_geral <- dados |>
  summarise(across(all_of(vars_resp), \(x) mean(x, na.rm = TRUE))) |>
  pivot_longer(everything(),
               names_to  = "variavel",
               values_to = "media_geral")

df_valores <- blups_GxA |>
  dplyr::select(variavel, genotipo, ano, blup_gxa_med, grupo) |>
  left_join(
    blups_G |> dplyr::select(variavel, genotipo, blup_orig_med),
    by = c("variavel", "genotipo")
  ) |>
  left_join(media_geral, by = "variavel") |>
  mutate(
    valor_ano = media_geral + blup_orig_med + blup_gxa_med,
    ano       = as.character(ano)
  )

# ── Wricke's ecovalence ────────────────────────────────────────────────────────
df_blups_G <- blups_G |>
  dplyr::select(variavel, genotipo, grupo,
                blup_orig_med, blup_orig_q025, blup_orig_q975,
                prob_positivo) |>
  left_join(media_geral, by = "variavel") |>
  mutate(
    valor_medio = media_geral + blup_orig_med,
    sinal = if_else(var_direcao[variavel] == "alto", 1, -1),
    blup_sinal = blup_orig_med * sinal
  )

df_wricke <- df_valores |>
  group_by(variavel, genotipo, grupo) |>
  summarise(
    blup_gxa_2023  = blup_gxa_med[ano == "2023"],
    blup_gxa_2024  = blup_gxa_med[ano == "2024"],
    wricke         = blup_gxa_2023^2 + blup_gxa_2024^2,
    amplitude_gxa  = abs(blup_gxa_2023 - blup_gxa_2024),
    inversao       = sign(blup_gxa_2023) != sign(blup_gxa_2024),
    .groups        = "drop"
  ) |>
  left_join(
    df_blups_G |>
      dplyr::select(variavel, genotipo, blup_orig_med, blup_sinal, valor_medio),
    by = c("variavel", "genotipo")
  ) |>
  group_by(variavel) |>
  mutate(
    wricke_total    = sum(wricke),
    wricke_rel      = wricke / wricke_total * 100,
    rank_estab      = rank(wricke),
    rank_desempenho = rank(-blup_sinal)
  ) |>
  ungroup()

# ── Bayesian probability of consistent superiority ────────────────────────────
intensidades <- c(0.10, 0.20, 0.25)
n_gen        <- 48L

calcular_prob_consistente <- function(var, intensidades, modelos, media_orig) {
  m   <- modelos[[var]]$modelo
  dp  <- modelos[[var]]$dp_orig
  mu  <- media_orig |> filter(variavel == var) |> pull(media_geral)
  post <- as_draws_df(m)
  
  cols_G   <- grep("^r_genotipo\\[", names(post), value = TRUE)
  cols_G   <- cols_G[!grepl("ano", cols_G)]
  cols_GxA <- grep("^r_genotipo:ano\\[", names(post), value = TRUE)
  gen_names <- gsub("r_genotipo\\[(.+),Intercept\\]", "\\1", cols_G)
  
  mat_G <- as.matrix(post[, cols_G]) * dp
  cols_2023 <- grep("_2023", cols_GxA, value = TRUE)
  cols_2024 <- grep("_2024", cols_GxA, value = TRUE)
  
  gen_order_2023 <- gsub("r_genotipo:ano\\[(.+)_2023,Intercept\\]", "\\1", cols_2023)
  gen_order_2024 <- gsub("r_genotipo:ano\\[(.+)_2024,Intercept\\]", "\\1", cols_2024)
  
  mat_GxA_2023 <- as.matrix(post[, cols_2023]) * dp
  mat_GxA_2024 <- as.matrix(post[, cols_2024]) * dp
  
  mat_GxA_2023 <- mat_GxA_2023[, match(gen_names, gen_order_2023)]
  mat_GxA_2024 <- mat_GxA_2024[, match(gen_names, gen_order_2024)]
  
  mat_2023 <- mu + mat_G + mat_GxA_2023
  mat_2024 <- mu + mat_G + mat_GxA_2024
  
  map_dfr(intensidades, function(intens) {
    n_top <- ceiling(n_gen * intens)
    direcao <- var_direcao[var]
    
    sup_2023 <- apply(mat_2023, 1, function(vals) if (direcao == "baixo") rank(vals) <= n_top else rank(-vals) <= n_top)
    sup_2024 <- apply(mat_2024, 1, function(vals) if (direcao == "baixo") rank(vals) <= n_top else rank(-vals) <= n_top)
    
    tibble(
      variavel         = var,
      genotipo         = gen_names,
      intensidade      = intens,
      prob_consistente = rowMeans(sup_2023 & sup_2024),
      prob_qualquer    = rowMeans(sup_2023 | sup_2024),
      prob_2023        = rowMeans(sup_2023),
      prob_2024        = rowMeans(sup_2024)
    )
  })
}

df_prob <- map_dfr(vars_resp, calcular_prob_consistente, 
                   intensidades = intensidades, modelos = modelos, media_orig = media_geral) |>
  mutate(genotipo = gsub("\\.", " ", genotipo)) |>
  left_join(df_wricke |> dplyr::select(variavel, genotipo, grupo, wricke, wricke_rel, rank_estab, rank_desempenho, inversao),
            by = c("variavel", "genotipo"))

# ── final ranking and export ──────────────────────────────────────────────────
df_ranking <- df_wricke |>
  left_join(df_prob |> filter(intensidade == 0.20) |> select(variavel, genotipo, prob_consistente, prob_2023, prob_2024),
            by = c("variavel", "genotipo")) |>
  group_by(variavel) |>
  mutate(
    rank_composto   = (rank_desempenho + rank_estab) / 2,
    quadrante       = paste(if_else(rank_desempenho <= 24, "superior", "inferior"), 
                            if_else(rank_estab <= 24, "estável", "instável"), sep = " + ")
  ) |>
  ungroup()

write_xlsx(list(wricke_ecovalencia = df_wricke, prob_consistente = df_prob, ranking_final = df_ranking),
           "outputs/tables/04_stability_results.xlsx")
