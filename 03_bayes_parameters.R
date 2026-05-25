# 03_bayes_parameters.R
# Modelo bayesiano final — sem parcela permanente
# Extrai: componentes de variância, H², acurácia, correlação genotípica entre
# anos, e BLUPs por genótipo e por genótipo:ano
#
# Modelo: y ~ ano + block + (1|genotipo) + (1|genotipo:ano)
#
# Saídas:
#   outputs/tables/03_varcomp.xlsx
#   outputs/tables/03_parametros_geneticos.xlsx
#   outputs/tables/03_blups_genotipo.xlsx
#   outputs/tables/03_blups_gxa.xlsx
#   outputs/figures/03_trace_*.png
#   outputs/figures/03_ppcheck_*.png
#   data/modelos_finais.rds

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(brms)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(writexl)

# ── dados e paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

# vars_resp <- "per_grao"

# ── configurações MCMC ─────────────────────────────────────────────────────────
options(mc.cores = 12)
bayesplot::color_scheme_set("blue")

n_iter   <- 4000
n_warmup <- 2000
n_chains <- 4
n_cores  <- 12
adapt_d  <- 0.95

# ── modelo final — sem parcela permanente ──────────────────────────────────────
formula_final <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano)
)

# ── priors ─────────────────────────────────────────────────────────────────────
priors_final <- c(
  prior(student_t(3, 0, 2.5), class = sd),
  prior(student_t(3, 0, 2.5), class = sigma),
  prior(normal(0, 5),         class = b)
)

# ── funções auxiliares ─────────────────────────────────────────────────────────

padronizar <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sg <- sd(x,   na.rm = TRUE)
  list(z = (x - mu) / sg, media = mu, dp = sg)
}

# extrai componentes de variância e parâmetros genéticos da posterior
extrair_parametros <- function(modelo, var_nome, media_orig, dp_orig, n_anos = 2, n_blocos = 3) {
  
  # amostras posteriores dos desvios padrão
  post <- as_draws_df(modelo)
  
  # identificar colunas sd
  sd_G   <- post[[ grep("^sd_genotipo__Intercept$",
                        names(post), value = TRUE) ]]
  sd_GxA <- post[[ grep("sd_genotipo:ano__Intercept",
                        names(post), value = TRUE) ]]
  sd_e   <- post[["sigma"]]
  
  # variâncias
  var_G   <- sd_G^2
  var_GxA <- sd_GxA^2
  var_e   <- sd_e^2
  var_P   <- var_G + var_GxA / n_anos + var_e / (n_anos * n_blocos)
  
  # herdabilidade no sentido amplo (base da média)
  H2 <- var_G / var_P
  
  # acurácia seletiva
  acuracia <- sqrt(H2)
  
  # correlação genotípica entre anos
  # rGA = σ²G / (σ²G + σ²GxA)
  # mede quanto do ranking genotípico se preserva entre anos
  rGA <- var_G / (var_G + var_GxA)
  
  # razão CVg/CVe (escala original)
  CVg <- (sqrt(var_G)  * dp_orig) / media_orig * 100
  CVe <- (sd_e         * dp_orig) / media_orig * 100
  razao_CVg_CVe <- CVg / CVe
  
  tibble(
    variavel        = var_nome,
    # componentes — mediana e IC 95% — escala original
    sigma2_G_med    = median(var_G)   * dp_orig^2,
    sigma2_G_q025   = quantile(var_G, 0.025) * dp_orig^2,
    sigma2_G_q975   = quantile(var_G, 0.975) * dp_orig^2,
    sigma2_GxA_med  = median(var_GxA) * dp_orig^2,
    sigma2_GxA_q025 = quantile(var_GxA, 0.025) * dp_orig^2,
    sigma2_GxA_q975 = quantile(var_GxA, 0.975) * dp_orig^2,
    sigma2_e_med    = median(var_e)   * dp_orig^2,
    sigma2_e_q025   = quantile(var_e, 0.025) * dp_orig^2,
    sigma2_e_q975   = quantile(var_e, 0.975) * dp_orig^2,
    # parâmetros genéticos — adimensionais
    H2_med          = median(H2),
    H2_q025         = quantile(H2, 0.025),
    H2_q975         = quantile(H2, 0.975),
    acuracia_med    = median(acuracia),
    acuracia_q025   = quantile(acuracia, 0.025),
    acuracia_q975   = quantile(acuracia, 0.975),
    rGA_med         = median(rGA),
    rGA_q025        = quantile(rGA, 0.025),
    rGA_q975        = quantile(rGA, 0.975),
    CVg_med         = median(CVg),
    CVe_med         = median(CVe),
    razao_CVg_CVe   = median(razao_CVg_CVe)
  )
}

# extrai BLUPs de genótipo (efeito médio entre anos)
extrair_blups_G <- function(modelo, var_nome, media_orig, dp_orig) {
  re   <- as_draws_df(modelo)
  cols <- grep("^r_genotipo\\[", names(re), value = TRUE)
  
  map_dfr(cols, function(col) {
    gen_name <- gsub("r_genotipo\\[(.+),Intercept\\]", "\\1", col)
    samples  <- re[[col]]
    tibble(
      variavel        = var_nome,
      genotipo        = gen_name,
      blup_z_med      = median(samples),
      blup_z_q025     = quantile(samples, 0.025),
      blup_z_q975     = quantile(samples, 0.975),
      blup_orig_med   = median(samples)  * dp_orig,
      blup_orig_q025  = quantile(samples, 0.025) * dp_orig,
      blup_orig_q975  = quantile(samples, 0.975) * dp_orig,
      prob_positivo   = mean(samples > 0)
    )
  })
}

# extrai BLUPs de genótipo:ano (desvio de cada genótipo em cada ano)
extrair_blups_GxA <- function(modelo, var_nome) {
  re   <- as_draws_df(modelo)
  cols <- grep("^r_genotipo:ano\\[", names(re), value = TRUE)
  
  map_dfr(cols, function(col) {
    nm       <- gsub("r_genotipo:ano\\[(.+),Intercept\\]", "\\1", col)
    parts    <- strsplit(nm, "_")[[1]]
    ano_val  <- tail(parts, 1)
    gen_name <- paste(head(parts, -1), collapse = "_")
    samples  <- re[[col]]
    tibble(
      variavel       = var_nome,
      genotipo       = gen_name,
      ano            = ano_val,
      blup_gxa_med   = median(samples),
      blup_gxa_q025  = quantile(samples, 0.025),
      blup_gxa_q975  = quantile(samples, 0.975),
      var_gxa        = var(samples)   # variância posterior — medida de instabilidade
    )
  })
}

# ── loop principal ─────────────────────────────────────────────────────────────
cat("\n── Modelos finais — início ──\n\n")

lista_modelos    <- list()
lista_parametros <- list()
lista_blups_G    <- list()
lista_blups_GxA  <- list()

for (v in vars_resp) {
  
  cat("══════════════════════════════════\n")
  cat("Variável:", v, "\n")
  
  pad       <- padronizar(dados[[v]])
  dados$y_z <- pad$z
  
  # ── ajuste
  cat("  Ajustando modelo final...\n")
  m <- brm(
    formula_final,
    data      = dados,
    prior     = priors_final,
    chains    = n_chains,
    iter      = n_iter,
    warmup    = n_warmup,
    cores     = n_cores,
    seed      = 42,
    save_pars = save_pars(all = TRUE),
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 100
  )
  
  # ── divergências
  nuts  <- nuts_params(m)
  n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  cat("  Divergências:", n_div, "\n")
  
  if (n_div > 10) {
    cat("  Reajustando com adapt_delta = 0.99...\n")
    m <- update(m,
                control   = list(adapt_delta = 0.99, max_treedepth = 15),
                save_pars = save_pars(all = TRUE),
                silent    = 2, refresh = 100)
    nuts  <- nuts_params(m)
    n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
    cat("  Divergências após reajuste:", n_div, "\n")
  }
  
  # ── diagnóstico visual
  p_trace <- mcmc_trace(
    as.array(m),
    pars       = grep("^sd_", variables(m), value = TRUE),
    facet_args = list(ncol = 2)
  ) +
    labs(title = paste0("Trace plots — ", v)) +
    theme_minimal(base_size = 10)
  
  ggsave(paste0(path_fig, "03_trace_", v, ".png"),
         p_trace, width = 9, height = 5, dpi = 300)
  
  p_ppc <- pp_check(m, ndraws = 100, type = "dens_overlay") +
    labs(title = paste0("Posterior predictive check — ", v)) +
    theme_minimal(base_size = 11)
  
  ggsave(paste0(path_fig, "03_ppcheck_", v, ".png"),
         p_ppc, width = 8, height = 5, dpi = 300)
  
  # ── extração
  cat("  Extraindo parâmetros e BLUPs...\n")
  lista_parametros[[v]] <- extrair_parametros(m, v, pad$media, pad$dp)
  lista_blups_G[[v]]    <- extrair_blups_G(m, v, pad$media, pad$dp)
  lista_blups_GxA[[v]]  <- extrair_blups_GxA(m, v)
  
  lista_modelos[[v]] <- list(modelo = m, media_orig = pad$media, dp_orig = pad$dp)
  
  cat("  Concluído:", v, "\n\n")
}

# ── consolidar e salvar ────────────────────────────────────────────────────────
cat("── Salvando resultados...\n")

params_final  <- list_rbind(lista_parametros)
blups_G_final <- list_rbind(lista_blups_G)
blups_GxA_final <- list_rbind(lista_blups_GxA)

# adicionar grupo ao blups
grupo_ref <- dados |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

blups_G_final <- blups_G_final |>
  left_join(grupo_ref, by = "genotipo")

blups_GxA_final <- blups_GxA_final |>
  left_join(grupo_ref, by = "genotipo")

# salvar em xlsx
write_xlsx(
  list(
    parametros_geneticos = params_final,
    blups_genotipo       = blups_G_final,
    blups_gxa            = blups_GxA_final
  ),
  paste0(path_tbl, "03_resultados_bayesianos.xlsx")
)

saveRDS(lista_modelos, "data/modelos_finais.rds")

cat("\nConcluído. Modelos salvos em data/modelos_finais.rds\n")
cat("Variáveis processadas:", length(lista_modelos), "\n")
