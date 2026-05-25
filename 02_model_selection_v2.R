# 02_model_selection_v2.R
# Model selection with correct structure for breeding:
#   Fixed:   Year, Block
#   Random: Genotype, G×A (Genotype:Year), Plot (Genotype×Block)
#
# Tests:
#   LRT1 — need for permanent Plot
#   LRT2 — need for G×A
#   LRT3 — need for Genotype
#   Boundary diagnostic by variable
#
# Outputs:
#   outputs/tables/02v2_model_selection.csv
#   outputs/figures/02v2_boundary_summary.png
#   outputs/figures/02v2_varcomp_profile.png
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(lme4)
library(tibble)
library(ggplot2)
library(patchwork)

# ── data and paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_out_fig <- "outputs/figures/"
path_out_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca", "vcm_mcm")

# ── optimizer control ─────────────────────────────────────────────────────
ctrl_reml <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ctrl_ml   <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ── auxiliary functions ─────────────────────────────────────────────────────────

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
  # both fitted with ML (REML = FALSE) — correct for fixed effects LRT
  # and for comparison of nested random structures
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

# ── candidate models ─────────────────────────────────────────────────────────
# M_full : full model — all random components
# M_noPP : without permanent plot
# M_noGxA: without G×A interaction
# M_noG  : without genotype (reference only — biologically meaningless)

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

# ── main loop ─────────────────────────────────────────────────────────────
cat("\n── Model selection v2 (random genotype) ──\n\n")

results <- map(vars_resp, function(v) {
  
  cat("────────────────────────────\n")
  cat("Variable:", v, "\n")
  
  fmls <- build_formulas(v)
  
  # fit with REML = TRUE — for final component estimation
  m_full_reml <- fit_lmer_safe(fmls$full,  dados, reml = TRUE)
  
  # fit with REML = FALSE — for LRT
  m_full_ml  <- fit_lmer_safe(fmls$full,  dados, reml = FALSE)
  m_noPP_ml  <- fit_lmer_safe(fmls$noPP,  dados, reml = FALSE)
  m_noGxA_ml <- fit_lmer_safe(fmls$noGxA, dados, reml = FALSE)
  m_noG_ml   <- fit_lmer_safe(fmls$noG,   dados, reml = FALSE)
  
  # ── singularity and boundary
  sing    <- isTRUE(isSingular(m_full_reml))
  bnd     <- is_boundary(m_full_reml)
  cat("  Singular:", sing, "| Boundary:", bnd, "\n")
  
  # ── variance components (REML)
  vc <- extract_vc_df(m_full_reml)
  cat("  Variance components (REML):\n")
  print(vc)
  
  # ── LRTs (ML)
  lrt_pp  <- lrt_ml(m_full_ml, m_noPP_ml)   # tests permanent plot
  lrt_gxa <- lrt_ml(m_full_ml, m_noGxA_ml)  # tests G×A
  lrt_g   <- lrt_ml(m_full_ml, m_noG_ml)    # tests genotype
  
  cat("  LRT permanent Plot — p:", round(lrt_pp$pval,  6),
      "| ΔAIC:", round(lrt_pp$delta_aic,  2), "\n")
  cat("  LRT G×A            — p:", round(lrt_gxa$pval, 6),
      "| ΔAIC:", round(lrt_gxa$delta_aic, 2), "\n")
  cat("  LRT Genotype       — p:", round(lrt_g$pval,   6),
      "| ΔAIC:", round(lrt_g$delta_aic,   2), "\n\n")
  
  # ── proportion of each component in total variance
  vc_named <- vc |> tibble::deframe()
  v_total  <- sum(vc_named, na.rm = TRUE)
  
  tibble(
    variavel        = v,
    singular        = sing,
    boundary        = bnd,
    # components (REML)
    sigma2_G        = vc_named["genotipo"],
    sigma2_GxA      = vc_named["genotipo:ano"],
    sigma2_PP       = vc_named["parcela"],
    sigma2_e        = vc_named["Residual"],
    # proportions
    prop_G          = sigma2_G   / v_total,
    prop_GxA        = sigma2_GxA / v_total,
    prop_PP         = sigma2_PP  / v_total,
    prop_e          = sigma2_e   / v_total,
    # LRTs
    lrt_pp_p        = lrt_pp$pval,
    lrt_pp_daic     = lrt_pp$delta_aic,
    lrt_gxa_p       = lrt_gxa$pval,
    lrt_gxa_daic    = lrt_gxa$delta_aic,
    lrt_g_p         = lrt_g$pval,
    lrt_g_daic      = lrt_g$delta_aic
  )
  
}) |> list_rbind()

print(results |> select(variavel, singular, boundary,
                        sigma2_G, sigma2_GxA, sigma2_PP, sigma2_e,
                        lrt_pp_p, lrt_gxa_p, lrt_g_p))

write.csv(results, paste0(path_out_tbl, "02v2_model_selection.csv"),
          row.names = FALSE)

# ── variance components visualization ──────────────────────────────────
# Variable display labels for plotting
var_labels <- c(
  per_grao  = "Grain proportion (% grain)",
  per_palha = "Husk proportion (% husk)",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag",
  vcm_mcm   = "FVol/FWM"
)

vc_long <- results |>
  select(variavel, prop_G, prop_GxA, prop_PP, prop_e) |>
  pivot_longer(-variavel, names_to = "componente", values_to = "proporcao") |>
  mutate(
    componente = factor(componente,
                        levels = c("prop_G", "prop_GxA", "prop_PP", "prop_e"),
                        labels = c("Genotype", "G×A", "Plot", "Residual")),
    variavel = recode(variavel, !!!var_labels)
  )

p_vc <- ggplot(vc_long, aes(x = variavel, y = proporcao, fill = componente)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Genotype" = "#185FA5",
      "G×A"      = "#D85A30",
      "Plot"     = "#1D9E75",
      "Residual" = "#B4B2A9"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = NULL, y = "Proportion of total variance",
    fill = "Component",
    title = "Variance decomposition — full model (REML)"
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

# ── boundary summary table for frequentist vs. bayesian decision ───────────────
cat("\n── Summary for frequentist vs. Bayesian decision ──\n")
results |>
  select(variavel, singular, boundary, lrt_pp_p, lrt_gxa_p, lrt_g_p) |>
  mutate(
    recommendation = case_when(
      boundary ~ "Bayesian",
      !boundary & !singular ~ "LMM-REML",
      TRUE ~ "Verify"
    )
  ) |>
  print()
