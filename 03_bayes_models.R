# 03_bayes_models.R
# Bayesian modeling via brms for all yield traits
#
# Strategy:
#   - Variables standardized (z-score) before fitting
#   - Prior student_t(3, 0, 2.5) for all variance components
#   - Full model: year + block + (1|genotype) + (1|genotype:year) + (1|plot)
#   - Validation: Prior Predictive Check + LOO with/without permanent plot
#   - Back-transformation of BLUPs to original scale
#
# Outputs:
#   outputs/figures/03_ppc_*.png         — prior predictive checks
#   outputs/figures/03_posterior_*.png   — posterior diagnostics
#   outputs/tables/03_varcomp_bayes.csv  — variance components
#   outputs/tables/03_blups_bayes.csv    — BLUPs on original scale
#   outputs/tables/03_loo_parcela.csv    — LOO comparison with/without plot
#   data/modelos_bayes.rds               — list with all fitted models
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)

# ── global settings ──────────────────────────────────────────────────────
options(mc.cores = 12)
bayesplot::color_scheme_set("blue")

# ── data and paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

# ── MCMC parameters ────────────────────────────────────────────────────────────
n_iter   <- 4000
n_warmup <- 2000
n_chains <- 4
n_cores  <- 12
adapt_d  <- 0.95

# ── priors (standardized scale) ────────────────────────────────────────────────
priors_base <- c(
  prior(student_t(3, 0, 2.5), class = sd),
  prior(student_t(3, 0, 2.5), class = sigma),
  prior(normal(0, 5),         class = b)
)

# ── full model formula ─────────────────────────────────────────────────
formula_completa <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano) + (1 | parcela)
)
formula_sem_pp <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano)
)

# ── auxiliary functions ─────────────────────────────────────────────────────────

# standardizes variable and keeps parameters for back-transformation
padronizar <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sg <- sd(x,   na.rm = TRUE)
  list(z = (x - mu) / sg, media = mu, dp = sg)
}

# extracts posterior variance components
extrair_varcomp <- function(modelo, var_nome, media_orig, dp_orig) {
  vc <- as.data.frame(VarCorr(modelo, summary = FALSE))
  
  # column names depend on the model — identify generically
  cols_sd <- grep("^sd_", names(vc), value = TRUE)
  
  map_dfr(cols_sd, function(col) {
    sd_post  <- vc[[col]]
    var_post <- sd_post^2
    
    tibble(
      variavel   = var_nome,
      componente = col,
      # standardized scale
      sd_median_z  = median(sd_post),
      sd_q025_z    = quantile(sd_post, 0.025),
      sd_q975_z    = quantile(sd_post, 0.975),
      var_median_z = median(var_post),
      # back-transformed to original scale
      var_median_orig = median(var_post) * dp_orig^2
    )
  })
}

# extracts BLUPs (random genotype effects) on original scale
extrair_blups <- function(modelo, var_nome, media_orig, dp_orig) {
  re <- ranef(modelo, summary = FALSE)
  
  if (!"genotipo" %in% names(re)) return(NULL)
  
  g_post <- re$genotipo[, , "Intercept"]  # matrix: iter × genotype
  
  tibble(
    variavel  = var_nome,
    genotipo  = colnames(g_post),
    blup_z_median = apply(g_post, 2, median),
    blup_z_q025   = apply(g_post, 2, \(x) quantile(x, 0.025)),
    blup_z_q975   = apply(g_post, 2, \(x) quantile(x, 0.975)),
    # back-transform: BLUP_orig = BLUP_z * dp_orig (deviation from mean)
    blup_orig_median = blup_z_median * dp_orig,
    blup_orig_q025   = blup_z_q025   * dp_orig,
    blup_orig_q975   = blup_z_q975   * dp_orig
  )
}

# ── main loop ─────────────────────────────────────────────────────────────
cat("\n── Bayesian modeling — start ──\n\n")

lista_modelos  <- list()
lista_varcomp  <- list()
lista_blups    <- list()
lista_loo      <- list()

for (v in vars_resp) {
  
  cat("══════════════════════════════════\n")
  cat("Variable:", v, "\n")
  
  # ── standardization
  pad       <- padronizar(dados[[v]])
  dados$y_z <- pad$z
  
  # ── 1. Prior Predictive Check
  cat("  [1/5] Prior predictive check...\n")
  m_ppc <- brm(
    formula_completa,
    data         = dados,
    prior        = priors_base,
    sample_prior = "only",
    chains       = 2,
    iter         = 1000,
    warmup       = 500,
    cores        = n_cores,
    seed         = 42,
    silent       = 2,
    refresh      = 0
  )
  
  ppc_draws <- posterior_predict(m_ppc, ndraws = 200)
  p_ppc <- ppc_dens_overlay(pad$z, ppc_draws[1:50, ]) +
    labs(
      title    = paste0("Prior predictive check — ", v),
      subtitle = "colored lines = prior samples | black line = actual data"
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(paste0(path_fig, "03_ppc_", v, ".png"),
         p_ppc, width = 8, height = 5, dpi = 300)
  
  # ── 2. Full model (with permanent plot)
  cat("  [2/5] Fitting full model...\n")
  m_full <- brm(
    formula_completa,
    data    = dados,
    prior   = priors_base,
    chains  = n_chains,
    iter    = n_iter,
    warmup  = n_warmup,
    cores   = n_cores,
    seed    = 42,
    save_pars = save_pars(all = TRUE),          # needed for moment_match
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 100
  )
  
  # check for divergences and refit if necessary
  nuts   <- nuts_params(m_full)
  n_div  <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  cat("  Divergences:", n_div, "\n")
  
  if (n_div > 10) {
    cat("  Refitting with adapt_delta = 0.99...\n")
    m_full <- update(
      m_full,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.99, max_treedepth = 15),
      silent    = 2,
      refresh   = 100
    )
    nuts  <- nuts_params(m_full)
    n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
    cat("  Divergences after refit:", n_div, "\n")
  }
  
  # ── 3. Model without permanent plot (for LOO)
  cat("  [3/5] Fitting model without permanent plot...\n")
  m_nopp <- brm(
    formula_sem_pp,
    data    = dados,
    prior   = priors_base,
    chains  = n_chains,
    iter    = n_iter,
    warmup  = n_warmup,
    cores   = n_cores,
    seed    = 42,
    save_pars = save_pars(all = TRUE),          # needed for moment_match
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 0
  )
  
  # ── 4. LOO — comparison with/without permanent plot ─────────────────────────
  cat("  [4/5] LOO-CV...\n")
  
  # Standard LOO — without moment_match to avoid memory issues
  loo_full <- loo(m_full)
  loo_nopp <- loo(m_nopp)
  
  # check for problematic points (k > 0.7)
  n_bad_full <- sum(loo_full$diagnostics$pareto_k > 0.7)
  n_bad_nopp <- sum(loo_nopp$diagnostics$pareto_k > 0.7)
  cat("  Points k > 0.7 — full:", n_bad_full,
      "| no plot:", n_bad_nopp, "\n")
  
  # if many problematic points, use WAIC as alternative
  if (n_bad_full > 10) {
    cat("  Many k > 0.7 points — using WAIC as supplement\n")
    waic_full <- waic(m_full)
    waic_nopp <- waic(m_nopp)
    cat("  WAIC full:", round(waic_full$estimates["waic", "Estimate"], 2),
        "| no plot:", round(waic_nopp$estimates["waic", "Estimate"], 2), "\n")
  }
  
  loo_comp <- loo_compare(loo_full, loo_nopp)
  
  lista_loo[[v]] <- tibble(
    variavel    = v,
    modelo      = c("full", "no_plot"),
    elpd_loo    = c(loo_full$estimates["elpd_loo", "Estimate"],
                    loo_nopp$estimates["elpd_loo", "Estimate"]),
    se_elpd_loo = c(loo_full$estimates["elpd_loo", "SE"],
                    loo_nopp$estimates["elpd_loo", "SE"]),
    n_bad_k     = c(n_bad_full, n_bad_nopp)
  )
  
  cat("  LOO comparison:\n")
  print(loo_comp)
  
  # ── 5. Posterior diagnostic and extraction
  cat("  [5/5] Diagnostic and extraction...\n")
  
  p_trace <- mcmc_trace(
    as.array(m_full),
    pars         = grep("^sd_", variables(m_full), value = TRUE),
    facet_args = list(ncol = 2)
  ) +
    labs(title = paste0("Trace plots — standard deviations — ", v)) +
    theme_minimal(base_size = 10)
  
  ggsave(paste0(path_fig, "03_trace_", v, ".png"),
         p_trace, width = 9, height = 6, dpi = 300)
  
  p_post <- pp_check(m_full, ndraws = 100, type = "dens_overlay") +
    labs(title = paste0("Posterior predictive check — ", v)) +
    theme_minimal(base_size = 11)
  
  ggsave(paste0(path_fig, "03_ppcheck_", v, ".png"),
         p_post, width = 8, height = 5, dpi = 300)
  
  lista_varcomp[[v]] <- extrair_varcomp(m_full, v, pad$media, pad$dp)
  lista_blups[[v]]   <- extrair_blups(m_full,   v, pad$media, pad$dp)
  
  lista_modelos[[v]] <- list(
    modelo       = m_full,
    modelo_npp = m_nopp,
    media_orig = pad$media,
    dp_orig    = pad$dp
  )
  
  cat("  Finished:", v, "\n\n")
}

# ── consolidate and save ────────────────────────────────────────────────────────
varcomp_final <- list_rbind(lista_varcomp)
blups_final   <- list_rbind(lista_blups)
loo_final     <- list_rbind(lista_loo)

write.csv(varcomp_final, paste0(path_tbl, "03_varcomp_bayes.csv"),  row.names = FALSE)
write.csv(blups_final,   paste0(path_tbl, "03_blups_bayes.csv"),    row.names = FALSE)
write.csv(loo_final,     paste0(path_tbl, "03_loo_parcela.csv"),    row.names = FALSE)

saveRDS(lista_modelos, "data/modelos_bayes.rds")

cat("\nFinished. Models saved in data/modelos_bayes.rds\n")
cat("Variables processed:", length(lista_modelos), "\n")
