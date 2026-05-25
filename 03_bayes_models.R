# 03_bayes_models.R
# Modelagem bayesiana via brms para todos os caracteres de rendimento
#
# Estratégia:
#   - Variáveis padronizadas (z-score) antes do ajuste
#   - Prior student_t(3, 0, 2.5) para todos os componentes de variância
#   - Modelo completo: ano + bloco + (1|genotipo) + (1|genotipo:ano) + (1|parcela)
#   - Validação: Prior Predictive Check + LOO com/sem parcela permanente
#   - Back-transformação dos BLUPs para escala original
#
# Saídas:
#   outputs/figures/03_ppc_*.png         — prior predictive checks
#   outputs/figures/03_posterior_*.png   — diagnóstico posterior
#   outputs/tables/03_varcomp_bayes.csv  — componentes de variância
#   outputs/tables/03_blups_bayes.csv    — BLUPs na escala original
#   outputs/tables/03_loo_parcela.csv    — comparação LOO com/sem parcela
#   data/modelos_bayes.rds               — lista com todos os modelos ajustados

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)

# ── configurações globais ──────────────────────────────────────────────────────
options(mc.cores = 12)
bayesplot::color_scheme_set("blue")

# ── dados e paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

# ── parâmetros MCMC ────────────────────────────────────────────────────────────
n_iter   <- 4000
n_warmup <- 2000
n_chains <- 4
n_cores  <- 12
adapt_d  <- 0.95

# ── priors (escala padronizada) ────────────────────────────────────────────────
priors_base <- c(
  prior(student_t(3, 0, 2.5), class = sd),
  prior(student_t(3, 0, 2.5), class = sigma),
  prior(normal(0, 5),         class = b)
)

# ── fórmula do modelo completo ─────────────────────────────────────────────────
formula_completa <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano) + (1 | parcela)
)
formula_sem_pp <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano)
)

# ── funções auxiliares ─────────────────────────────────────────────────────────

# padroniza variável e guarda parâmetros para back-transformação
padronizar <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sg <- sd(x,   na.rm = TRUE)
  list(z = (x - mu) / sg, media = mu, dp = sg)
}

# extrai componentes de variância posteriores
extrair_varcomp <- function(modelo, var_nome, media_orig, dp_orig) {
  vc <- as.data.frame(VarCorr(modelo, summary = FALSE))
  
  # nomes das colunas dependem do modelo — identificar genericamente
  cols_sd <- grep("^sd_", names(vc), value = TRUE)
  
  map_dfr(cols_sd, function(col) {
    sd_post  <- vc[[col]]
    var_post <- sd_post^2
    
    tibble(
      variavel   = var_nome,
      componente = col,
      # escala padronizada
      sd_median_z  = median(sd_post),
      sd_q025_z    = quantile(sd_post, 0.025),
      sd_q975_z    = quantile(sd_post, 0.975),
      var_median_z = median(var_post),
      # back-transformado para escala original
      var_median_orig = median(var_post) * dp_orig^2
    )
  })
}

# extrai BLUPs (efeitos aleatórios de genótipo) na escala original
extrair_blups <- function(modelo, var_nome, media_orig, dp_orig) {
  re <- ranef(modelo, summary = FALSE)
  
  if (!"genotipo" %in% names(re)) return(NULL)
  
  g_post <- re$genotipo[, , "Intercept"]  # matriz: iter × genótipo
  
  tibble(
    variavel  = var_nome,
    genotipo  = colnames(g_post),
    blup_z_median = apply(g_post, 2, median),
    blup_z_q025   = apply(g_post, 2, \(x) quantile(x, 0.025)),
    blup_z_q975   = apply(g_post, 2, \(x) quantile(x, 0.975)),
    # back-transform: BLUP_orig = BLUP_z * dp_orig  (desvio em relação à média)
    blup_orig_median = blup_z_median * dp_orig,
    blup_orig_q025   = blup_z_q025   * dp_orig,
    blup_orig_q975   = blup_z_q975   * dp_orig
  )
}

# ── loop principal ─────────────────────────────────────────────────────────────
cat("\n── Modelagem bayesiana — início ──\n\n")

lista_modelos  <- list()
lista_varcomp  <- list()
lista_blups    <- list()
lista_loo      <- list()

for (v in vars_resp) {
  
  cat("══════════════════════════════════\n")
  cat("Variável:", v, "\n")
  
  # ── padronização
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
      subtitle = "linhas coloridas = amostras do prior | linha preta = dados reais"
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(paste0(path_fig, "03_ppc_", v, ".png"),
         p_ppc, width = 8, height = 5, dpi = 300)
  
  # ── 2. Modelo completo (com parcela permanente)
  cat("  [2/5] Ajustando modelo completo...\n")
  m_full <- brm(
    formula_completa,
    data      = dados,
    prior     = priors_base,
    chains    = n_chains,
    iter      = n_iter,
    warmup    = n_warmup,
    cores     = n_cores,
    seed      = 42,
    save_pars = save_pars(all = TRUE),          # necessário para moment_match
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 100
  )
  
  # verificar divergências e reajustar se necessário
  nuts   <- nuts_params(m_full)
  n_div  <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  cat("  Divergências:", n_div, "\n")
  
  if (n_div > 10) {
    cat("  Reajustando com adapt_delta = 0.99...\n")
    m_full <- update(
      m_full,
      save_pars = save_pars(all = TRUE),
      control   = list(adapt_delta = 0.99, max_treedepth = 15),
      silent    = 2,
      refresh   = 100
    )
    nuts  <- nuts_params(m_full)
    n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
    cat("  Divergências após reajuste:", n_div, "\n")
  }
  
  # ── 3. Modelo sem parcela permanente (para LOO)
  cat("  [3/5] Ajustando modelo sem parcela permanente...\n")
  m_nopp <- brm(
    formula_sem_pp,
    data      = dados,
    prior     = priors_base,
    chains    = n_chains,
    iter      = n_iter,
    warmup    = n_warmup,
    cores     = n_cores,
    seed      = 42,
    save_pars = save_pars(all = TRUE),          # necessário para moment_match
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 0
  )
  
  # ── 4. LOO — comparação com/sem parcela permanente ─────────────────────────
  cat("  [4/5] LOO-CV...\n")
  
  # LOO padrão — sem moment_match para evitar problema de memória
  loo_full <- loo(m_full)
  loo_nopp <- loo(m_nopp)
  
  # verificar se há pontos problemáticos (k > 0.7)
  n_bad_full <- sum(loo_full$diagnostics$pareto_k > 0.7)
  n_bad_nopp <- sum(loo_nopp$diagnostics$pareto_k > 0.7)
  cat("  Pontos k > 0.7 — completo:", n_bad_full,
      "| sem parcela:", n_bad_nopp, "\n")
  
  # se muitos pontos problemáticos, usar WAIC como alternativa
  if (n_bad_full > 10) {
    cat("  Muitos pontos k > 0.7 — usando WAIC como complemento\n")
    waic_full <- waic(m_full)
    waic_nopp <- waic(m_nopp)
    cat("  WAIC completo:", round(waic_full$estimates["waic", "Estimate"], 2),
        "| sem parcela:", round(waic_nopp$estimates["waic", "Estimate"], 2), "\n")
  }
  
  loo_comp <- loo_compare(loo_full, loo_nopp)
  
  lista_loo[[v]] <- tibble(
    variavel    = v,
    modelo      = c("completo", "sem_parcela"),
    elpd_loo    = c(loo_full$estimates["elpd_loo", "Estimate"],
                    loo_nopp$estimates["elpd_loo", "Estimate"]),
    se_elpd_loo = c(loo_full$estimates["elpd_loo", "SE"],
                    loo_nopp$estimates["elpd_loo", "SE"]),
    n_bad_k     = c(n_bad_full, n_bad_nopp)
  )
  
  cat("  LOO comparação:\n")
  print(loo_comp)
  
  # ── 5. Diagnóstico posterior e extração
  cat("  [5/5] Diagnóstico e extração...\n")
  
  p_trace <- mcmc_trace(
    as.array(m_full),
    pars       = grep("^sd_", variables(m_full), value = TRUE),
    facet_args = list(ncol = 2)
  ) +
    labs(title = paste0("Trace plots — desvios padrão — ", v)) +
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
    modelo     = m_full,
    modelo_npp = m_nopp,
    media_orig = pad$media,
    dp_orig    = pad$dp
  )
  
  cat("  Concluído:", v, "\n\n")
}

# ── consolidar e salvar ────────────────────────────────────────────────────────
varcomp_final <- list_rbind(lista_varcomp)
blups_final   <- list_rbind(lista_blups)
loo_final     <- list_rbind(lista_loo)

write.csv(varcomp_final, paste0(path_tbl, "03_varcomp_bayes.csv"),  row.names = FALSE)
write.csv(blups_final,   paste0(path_tbl, "03_blups_bayes.csv"),    row.names = FALSE)
write.csv(loo_final,     paste0(path_tbl, "03_loo_parcela.csv"),    row.names = FALSE)

saveRDS(lista_modelos, "data/modelos_bayes.rds")

cat("\nConcluído. Modelos salvos em data/modelos_bayes.rds\n")
cat("Variáveis processadas:", length(lista_modelos), "\n")
