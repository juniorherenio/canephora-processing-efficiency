# 02_model_selection_v2.R
# Seleção de modelo com estrutura correta para melhoramento:
#   Fixos:     Ano, Bloco
#   Aleatórios: Genótipo, G×A (Genótipo:Ano), Parcela (Genótipo×Bloco)
#
# Testes:
#   LRT1 — necessidade de Parcela permanente
#   LRT2 — necessidade de G×A
#   LRT3 — necessidade de Genótipo
#   Diagnóstico de boundary por variável
#
# Saídas:
#   outputs/tables/02v2_model_selection.csv
#   outputs/figures/02v2_boundary_summary.png
#   outputs/figures/02v2_varcomp_profile.png

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(tibble)
library(ggplot2)
library(patchwork)

# ── dados e paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_out_fig <- "outputs/figures/"
path_out_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca", "vcm_mcm")

# ── controle do otimizador ─────────────────────────────────────────────────────
ctrl_reml <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ctrl_ml   <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ── funções auxiliares ─────────────────────────────────────────────────────────

fit_lmer_safe <- function(formula, data, reml) {
  tryCatch(
    lmer(formula, data = data, REML = reml,
         control = if (reml) ctrl_reml else ctrl_ml),
    error   = function(e) NULL,
    warning = function(w) {
      suppressWarnings(
        lmer(formula, data = data, REML = reml,
             control = if (reml) ctrl_reml else ctrl_ml)
      )
    }
  )
}

extract_vc_df <- function(model) {
  if (is.null(model)) return(tibble(grp = NA_character_, vcov = NA_real_))
  as.data.frame(VarCorr(model)) |>
    select(grp, vcov) |>
    as_tibble()
}

is_boundary <- function(model) {
  if (is.null(model)) return(NA)
  vc <- extract_vc_df(model)
  any(vc$vcov <= 1e-6, na.rm = TRUE)
}

lrt_ml <- function(m_full, m_red) {
  # ambos ajustados com ML (REML = FALSE) — correto para LRT de efeitos fixos
  # e para comparação de estruturas aleatórias aninhadas
  if (is.null(m_full) || is.null(m_red)) {
    return(tibble(delta_aic = NA_real_, delta_bic = NA_real_, pval = NA_real_))
  }
  lt <- anova(m_full, m_red, refit = FALSE)
  tibble(
    delta_aic = lt$AIC[1]  - lt$AIC[2],
    delta_bic = lt$BIC[1]  - lt$BIC[2],
    pval      = lt[["Pr(>Chisq)"]][2]
  )
}

# ── modelos candidatos ─────────────────────────────────────────────────────────
# M_full : modelo completo — todos os componentes aleatórios
# M_noPP : sem parcela permanente
# M_noGxA: sem interação G×A
# M_noG  : sem genótipo (apenas para referência — não faz sentido biológico)

build_formulas <- function(var) {
  list(
    full  = as.formula(paste(var,
                             "~ ano + block + (1|genotipo) + (1|genotipo:ano) + (1|parcela)")),
    noPP  = as.formula(paste(var,
                             "~ ano + block + (1|genotipo) + (1|genotipo:ano)")),
    noGxA = as.formula(paste(var,
                             "~ ano + block + (1|genotipo) + (1|parcela)")),
    noG   = as.formula(paste(var,
                             "~ ano + block + (1|genotipo:ano) + (1|parcela)"))
  )
}

# ── loop principal ─────────────────────────────────────────────────────────────
cat("\n── Seleção de modelo v2 (genótipo aleatório) ──\n\n")

results <- map(vars_resp, function(v) {
  
  cat("────────────────────────────\n")
  cat("Variável:", v, "\n")
  
  fmls <- build_formulas(v)
  
  # ajuste com REML = TRUE — para estimativa final dos componentes
  m_full_reml <- fit_lmer_safe(fmls$full,  dados, reml = TRUE)
  
  # ajuste com REML = FALSE — para LRT
  m_full_ml  <- fit_lmer_safe(fmls$full,  dados, reml = FALSE)
  m_noPP_ml  <- fit_lmer_safe(fmls$noPP,  dados, reml = FALSE)
  m_noGxA_ml <- fit_lmer_safe(fmls$noGxA, dados, reml = FALSE)
  m_noG_ml   <- fit_lmer_safe(fmls$noG,   dados, reml = FALSE)
  
  # ── singularidade e boundary
  sing    <- isTRUE(isSingular(m_full_reml))
  bnd     <- is_boundary(m_full_reml)
  cat("  Singular:", sing, "| Boundary:", bnd, "\n")
  
  # ── componentes de variância (REML)
  vc <- extract_vc_df(m_full_reml)
  cat("  Componentes de variância (REML):\n")
  print(vc)
  
  # ── LRTs (ML)
  lrt_pp  <- lrt_ml(m_full_ml, m_noPP_ml)   # testa parcela permanente
  lrt_gxa <- lrt_ml(m_full_ml, m_noGxA_ml)  # testa G×A
  lrt_g   <- lrt_ml(m_full_ml, m_noG_ml)    # testa genótipo
  
  cat("  LRT Parcela permanente — p:", round(lrt_pp$pval,  6),
      "| ΔAIC:", round(lrt_pp$delta_aic,  2), "\n")
  cat("  LRT G×A               — p:", round(lrt_gxa$pval, 6),
      "| ΔAIC:", round(lrt_gxa$delta_aic, 2), "\n")
  cat("  LRT Genótipo          — p:", round(lrt_g$pval,   6),
      "| ΔAIC:", round(lrt_g$delta_aic,   2), "\n\n")
  
  # ── proporção de cada componente na variância total
  vc_named <- vc |> tibble::deframe()
  v_total  <- sum(vc_named, na.rm = TRUE)
  
  tibble(
    variavel          = v,
    singular          = sing,
    boundary          = bnd,
    # componentes (REML)
    sigma2_G          = vc_named["genotipo"],
    sigma2_GxA        = vc_named["genotipo:ano"],
    sigma2_PP         = vc_named["parcela"],
    sigma2_e          = vc_named["Residual"],
    # proporções
    prop_G            = sigma2_G   / v_total,
    prop_GxA          = sigma2_GxA / v_total,
    prop_PP           = sigma2_PP  / v_total,
    prop_e            = sigma2_e   / v_total,
    # LRTs
    lrt_pp_p          = lrt_pp$pval,
    lrt_pp_daic       = lrt_pp$delta_aic,
    lrt_gxa_p         = lrt_gxa$pval,
    lrt_gxa_daic      = lrt_gxa$delta_aic,
    lrt_g_p           = lrt_g$pval,
    lrt_g_daic        = lrt_g$delta_aic
  )
  
}) |> list_rbind()

print(results |> select(variavel, singular, boundary,
                        sigma2_G, sigma2_GxA, sigma2_PP, sigma2_e,
                        lrt_pp_p, lrt_gxa_p, lrt_g_p))

write.csv(results, paste0(path_out_tbl, "02v2_model_selection.csv"),
          row.names = FALSE)

# ── visualização dos componentes de variância ──────────────────────────────────
vc_long <- results |>
  select(variavel, prop_G, prop_GxA, prop_PP, prop_e) |>
  pivot_longer(-variavel, names_to = "componente", values_to = "proporcao") |>
  mutate(
    componente = factor(componente,
                        levels = c("prop_G", "prop_GxA", "prop_PP", "prop_e"),
                        labels = c("Genótipo", "G×A", "Parcela", "Resíduo"))
  )

p_vc <- ggplot(vc_long, aes(x = variavel, y = proporcao, fill = componente)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Genótipo" = "#185FA5",
      "G×A"      = "#D85A30",
      "Parcela"  = "#1D9E75",
      "Resíduo"  = "#B4B2A9"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = NULL, y = "proporção da variância total",
    fill = "componente",
    title = "Decomposição da variância — modelo completo (REML)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position   = "top",
    axis.text.x       = element_text(angle = 30, hjust = 1)
  )

ggsave(paste0(path_out_fig, "02v2_varcomp_stacked.png"),
       p_vc, width = 9, height = 6, dpi = 300)

# ── tabela resumo de boundary para decisão frequentista vs bayes ───────────────
cat("\n── Resumo para decisão frequentista vs. bayesiano ──\n")
results |>
  select(variavel, singular, boundary, lrt_pp_p, lrt_gxa_p, lrt_g_p) |>
  mutate(
    recomendacao = case_when(
      boundary ~ "Bayesiano",
      !boundary & !singular ~ "LMM-REML",
      TRUE ~ "Verificar"
    )
  ) |>
  print()
