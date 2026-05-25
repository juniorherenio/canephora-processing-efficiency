# 03b_reml_comparison.R
# Comparação entre estimativas de componentes de variância e parâmetros
# genéticos via REML (lme4) e via Inferência Bayesiana (brms)
#
# Objetivo: demonstrar que REML produz estimativas de fronteira (boundary)
# para várias variáveis e subestima a incerteza dos parâmetros, justificando
# o uso da abordagem bayesiana
#
# Modelo: y ~ ano + block + (1|genotipo) + (1|genotipo:ano)
#
# Saídas:
#   outputs/tables/03b_reml_vs_bayes.xlsx
#   outputs/figures/03b_comparacao_h2.png
#   outputs/figures/03b_comparacao_rga.png
#   outputs/figures/03b_varcomp_reml_vs_bayes.png

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lme4)
library(ggplot2)
library(patchwork)
library(writexl)

# ── dados e paths ──────────────────────────────────────────────────────────────
dados        <- readRDS("data/dados_clean.rds")
bayes_params <- readRDS("data/modelos_finais.rds")

# carregar tabela de parâmetros bayesianos já calculada
params_bayes <- readxl::read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "parametros_geneticos"
)

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

n_anos   <- 2
n_blocos <- 3

# ── controle do otimizador ─────────────────────────────────────────────────────
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ── funções auxiliares ─────────────────────────────────────────────────────────

# extrai componentes de variância do lmer e calcula parâmetros genéticos
extrair_reml <- function(var, data) {
  
  f <- as.formula(
    paste(var, "~ ano + block + (1|genotipo) + (1|genotipo:ano)")
  )
  
  m <- tryCatch(
    lmer(f, data = data, REML = TRUE, control = ctrl),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lmer(f, data = data, REML = TRUE, control = ctrl)
    )
  )
  
  if (is.null(m)) {
    return(tibble(variavel = var, metodo = "REML",
                  sigma2_G = NA, sigma2_GxA = NA, sigma2_e = NA,
                  H2 = NA, acuracia = NA, rGA = NA,
                  CVg = NA, CVe = NA, razao_CVg_CVe = NA,
                  singular = NA, boundary = NA))
  }
  
  vc       <- as.data.frame(VarCorr(m))
  vc_named <- setNames(vc$vcov, vc$grp)
  
  sigma2_G   <- vc_named["genotipo"]
  sigma2_GxA <- vc_named["genotipo:ano"]
  sigma2_e   <- vc_named["Residual"]
  
  # herdabilidade na base da média
  sigma2_P <- sigma2_G + sigma2_GxA / n_anos + sigma2_e / (n_anos * n_blocos)
  H2       <- sigma2_G / sigma2_P
  
  # acurácia
  acuracia <- sqrt(H2)
  
  # correlação genotípica entre anos
  rGA <- sigma2_G / (sigma2_G + sigma2_GxA)
  
  # CVg e CVe
  media_var <- mean(data[[var]], na.rm = TRUE)
  CVg       <- sqrt(sigma2_G)  / media_var * 100
  CVe       <- sqrt(sigma2_e)  / media_var * 100
  
  # diagnóstico
  sing <- isSingular(m)
  bnd  <- any(vc$vcov <= 1e-6)
  
  tibble(
    variavel      = var,
    metodo        = "REML",
    sigma2_G      = sigma2_G,
    sigma2_GxA    = sigma2_GxA,
    sigma2_e      = sigma2_e,
    H2            = H2,
    acuracia      = acuracia,
    rGA           = rGA,
    CVg           = CVg,
    CVe           = CVe,
    razao_CVg_CVe = CVg / CVe,
    singular      = sing,
    boundary      = bnd
  )
}

# ── 1. REML para todas as variáveis ───────────────────────────────────────────
cat("\n── Estimativas REML ──\n\n")

reml_results <- map_dfr(vars_resp, function(v) {
  cat("Variável:", v, "\n")
  res <- extrair_reml(v, dados)
  cat("  Singular:", res$singular,
      "| Boundary:", res$boundary,
      "| H²:", round(res$H2, 3),
      "| rGA:", round(res$rGA, 3), "\n")
  res
})

print(reml_results |> select(variavel, sigma2_G, sigma2_GxA, sigma2_e,
                             H2, rGA, singular, boundary))

# ── 2. Tabela comparativa REML vs Bayes ───────────────────────────────────────

# formatar bayes para comparação
bayes_comp <- params_bayes |>
  select(variavel,
         sigma2_G   = sigma2_G_med,
         sigma2_GxA = sigma2_GxA_med,
         sigma2_e   = sigma2_e_med,
         H2         = H2_med,
         H2_q025,
         H2_q975,
         acuracia   = acuracia_med,
         rGA        = rGA_med,
         rGA_q025,
         rGA_q975,
         CVg        = CVg_med,
         CVe        = CVe_med,
         razao_CVg_CVe) |>
  mutate(metodo = "Bayes", singular = FALSE, boundary = FALSE)

# juntar
comp_full <- bind_rows(
  reml_results |>
    mutate(H2_q025 = NA, H2_q975 = NA,
           rGA_q025 = NA, rGA_q975 = NA),
  bayes_comp
) |>
  arrange(variavel, metodo)

# tabela resumo lado a lado
comp_wide <- comp_full |>
  select(variavel, metodo, sigma2_G, sigma2_GxA, sigma2_e,
         H2, H2_q025, H2_q975, rGA, rGA_q025, rGA_q975,
         singular, boundary) |>
  pivot_wider(
    names_from  = metodo,
    values_from = c(sigma2_G, sigma2_GxA, sigma2_e,
                    H2, H2_q025, H2_q975,
                    rGA, rGA_q025, rGA_q975,
                    singular, boundary),
    names_glue  = "{.value}__{metodo}"
  )

cat("\n── Comparação REML vs Bayes ──\n")
print(comp_wide |> select(variavel,
                          H2__REML, H2__Bayes, H2_q025__Bayes, H2_q975__Bayes,
                          rGA__REML, rGA__Bayes, rGA_q025__Bayes, rGA_q975__Bayes,
                          boundary__REML))

# ── 3. visualizações ───────────────────────────────────────────────────────────

# labels para o eixo x
var_labels <- c(
  per_grao  = "% grão",
  per_palha = "% palha",
  mcm_mgb   = "MCM/MGB",
  mcm_saca  = "MCM/saca",
  vcm_saca  = "VCM/saca",
  vcm_mcm   = "VCM/MCM"
)

# paleta
cores <- c("REML" = "#D85A30", "Bayes" = "#185FA5")

# 3a. comparação H²
df_h2 <- comp_full |>
  filter(!is.na(H2)) |>
  mutate(
    variavel_lab = var_labels[variavel],
    variavel_lab = factor(variavel_lab, levels = var_labels)
  )

p_h2 <- ggplot(df_h2, aes(x = variavel_lab, y = H2, colour = metodo,
                          shape = metodo)) +
  geom_errorbar(
    data = df_h2 |> filter(metodo == "Bayes"),
    aes(ymin = H2_q025, ymax = H2_q975),
    width = 0.15, linewidth = 0.5
  ) +
  geom_point(size = 3.5) +
  geom_point(
    data = df_h2 |> filter(metodo == "REML" & boundary == TRUE),
    aes(x = variavel_lab, y = H2),
    shape = 4, size = 5, colour = "#D85A30", stroke = 1.2
  ) +
  scale_colour_manual(values = cores) +
  scale_shape_manual(values = c("REML" = 17, "Bayes" = 16)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    x = NULL,
    y = "herdabilidade (H²)",
    colour = "método",
    shape  = "método",
    title  = "Herdabilidade — REML vs Bayes",
    subtitle = "× indica estimativa de fronteira (boundary) no REML"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    axis.text.x        = element_text(angle = 30, hjust = 1)
  )

# 3b. comparação rGA
p_rga <- ggplot(df_h2, aes(x = variavel_lab, y = rGA, colour = metodo,
                           shape = metodo)) +
  geom_errorbar(
    data = df_h2 |> filter(metodo == "Bayes"),
    aes(ymin = rGA_q025, ymax = rGA_q975),
    width = 0.15, linewidth = 0.5
  ) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "gray50", linewidth = 0.4) +
  scale_colour_manual(values = cores) +
  scale_shape_manual(values = c("REML" = 17, "Bayes" = 16)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "correlação genotípica entre anos (rGA)",
    colour = "método",
    shape  = "método",
    title  = "Correlação genotípica entre anos — REML vs Bayes",
    subtitle = "linha tracejada = rGA = 0.5"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    axis.text.x        = element_text(angle = 30, hjust = 1)
  )

# 3c. decomposição da variância — REML vs Bayes lado a lado
df_vc <- comp_full |>
  select(variavel, metodo, sigma2_G, sigma2_GxA, sigma2_e) |>
  pivot_longer(
    cols      = c(sigma2_G, sigma2_GxA, sigma2_e),
    names_to  = "componente",
    values_to = "variancia"
  ) |>
  group_by(variavel, metodo) |>
  mutate(
    prop = variancia / sum(variancia, na.rm = TRUE),
    componente = factor(componente,
                        levels = c("sigma2_G", "sigma2_GxA", "sigma2_e"),
                        labels = c("Genótipo", "G×A", "Resíduo")),
    variavel_lab = var_labels[variavel],
    metodo_lab   = paste0(var_labels[variavel], "\n(", metodo, ")")
  ) |>
  ungroup()

p_vc <- ggplot(df_vc |> filter(!is.na(prop)),
               aes(x = metodo_lab, y = prop, fill = componente)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Genótipo" = "#185FA5",
      "G×A"      = "#D85A30",
      "Resíduo"  = "#B4B2A9"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~variavel_lab, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    y = "proporção da variância total",
    fill = "componente",
    title = "Decomposição da variância — REML vs Bayes"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    axis.text.x        = element_text(angle = 40, hjust = 1, size = 7),
    strip.text         = element_text(size = 8)
  )

# salvar figuras
ggsave(paste0(path_fig, "03b_comparacao_h2.png"),
       p_h2,  width = 8, height = 5, dpi = 300)

ggsave(paste0(path_fig, "03b_comparacao_rga.png"),
       p_rga, width = 8, height = 5, dpi = 300)

ggsave(paste0(path_fig, "03b_varcomp_reml_vs_bayes.png"),
       p_vc,  width = 14, height = 5, dpi = 300)

# ── 4. salvar tabelas ──────────────────────────────────────────────────────────
write_xlsx(
  list(
    reml_parametros    = reml_results,
    comparacao_wide    = comp_wide,
    comparacao_long    = comp_full |>
      select(variavel, metodo, sigma2_G, sigma2_GxA, sigma2_e,
             H2, H2_q025, H2_q975, acuracia,
             rGA, rGA_q025, rGA_q975,
             CVg, CVe, razao_CVg_CVe, singular, boundary)
  ),
  paste0(path_tbl, "03b_reml_vs_bayes.xlsx")
)

cat("\nConcluído. Resultados salvos em outputs/tables/03b_reml_vs_bayes.xlsx\n")
