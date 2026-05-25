# 03_bayes_parameters.R
# Final Bayesian model — without permanent plot
# Extracts: variance components, H², accuracy, genotypic correlation between
# years, and BLUPs by genotype and by genotype:year
#
# Model: y ~ ano + block + (1|genotipo) + (1|genotipo:ano)
#
# Outputs:
#   outputs/tables/03_varcomp.xlsx
#   outputs/tables/03_parametros_geneticos.xlsx
#   outputs/tables/03_blups_genotipo.xlsx
#   outputs/tables/03_blups_gxa.xlsx
#   outputs/figures/03_trace_*.png
#   outputs/figures/03_ppcheck_*.png
#   data/modelos_finais.rds
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
library(writexl)

# ── data and paths ──────────────────────────────────────────────────────────────
dados <- readRDS("data/dados_clean.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

# ── MCMC settings ─────────────────────────────────────────────────────────
options(mc.cores = 12)
bayesplot::color_scheme_set("blue")

n_iter   <- 4000
n_warmup <- 2000
n_chains <- 4
n_cores  <- 12
adapt_d  <- 0.95

# ── final model — without permanent plot ──────────────────────────────────────
formula_final <- bf(
  y_z ~ ano + block + (1 | genotipo) + (1 | genotipo:ano)
)

# ── priors ─────────────────────────────────────────────────────────────────────
priors_final <- c(
  prior(student_t(3, 0, 2.5), class = sd),
  prior(student_t(3, 0, 2.5), class = sigma),
  prior(normal(0, 5),         class = b)
)

# ── auxiliary functions ─────────────────────────────────────────────────────────

padronizar <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sg <- sd(x,   na.rm = TRUE)
  list(z = (x - mu) / sg, media = mu, dp = sg)
}

# extracts variance components and genetic parameters from the posterior
extrair_parametros <- function(modelo, var_nome, media_orig, dp_orig, n_anos = 2, n_blocos = 3) {
  
  # posterior samples of standard deviations
  post <- as_draws_df(modelo)
  
  # identify sd columns
  sd_G   <- post[[ grep("^sd_genotipo__Intercept$",
                         names(post), value = TRUE) ]]
  sd_GxA <- post[[ grep("sd_genotipo:ano__Intercept",
                         names(post), value = TRUE) ]]
  sd_e   <- post[["sigma"]]
  
  # variances
  var_G   <- sd_G^2
  var_GxA <- sd_GxA^2
  var_e   <- sd_e^2
  var_P   <- var_G + var_GxA / n_anos + var_e / (n_anos * n_blocos)
  
  # broad-sense heritability (mean basis)
  H2 <- var_G / var_P
  
  # selective accuracy
  acuracia <- sqrt(H2)
  
  # genotypic correlation between years
  # rGA = σ²G / (σ²G + σ²GxA)
  # measures how much of the genotypic ranking is preserved between years
  rGA <- var_G / (var_G + var_GxA)
  
  # CVg/CVe ratio (original scale)
  CVg <- (sqrt(var_G)  * dp_orig) / media_orig * 100
  CVe <- (sd_e         * dp_orig) / media_orig * 100
  razao_CVg_CVe <- CVg / CVe
  
  tibble(
    variavel        = var_nome,
    # components — median and 95% CI — original scale
    sigma2_G_med    = median(var_G)   * dp_orig^2,
    sigma2_G_q025   = quantile(var_G, 0.025) * dp_orig^2,
    sigma2_G_q975   = quantile(var_G, 0.975) * dp_orig^2,
    sigma2_GxA_med  = median(var_GxA) * dp_orig^2,
    sigma2_GxA_q025 = quantile(var_GxA, 0.025) * dp_orig^2,
    sigma2_GxA_q975 = quantile(var_GxA, 0.975) * dp_orig^2,
    sigma2_e_med    = median(var_e)   * dp_orig^2,
    sigma2_e_q025   = quantile(var_e, 0.025) * dp_orig^2,
    sigma2_e_q975   = quantile(var_e, 0.975) * dp_orig^2,
    # genetic parameters — dimensionless
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

# extracts genotype BLUPs (mean effect between years)
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

# extracts genotype:year BLUPs (deviation of each genotype in each year)
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
      var_gxa        = var(samples)   # posterior variance — measure of instability
    )
  })
}

# ── main loop ─────────────────────────────────────────────────────────────
cat("\n── Final models — start ──\n\n")

lista_modelos    <- list()
lista_parametros <- list()
lista_blups_G    <- list()
lista_blups_GxA  <- list()

for (v in vars_resp) {
  
  cat("══════════════════════════════════\n")
  cat("Variable:", v, "\n")
  
  pad       <- padronizar(dados[[v]])
  dados$y_z <- pad$z
  
  # ── fit
  cat("  Fitting final model...\n")
  m <- brm(
    formula_final,
    data    = dados,
    prior   = priors_final,
    chains  = n_chains,
    iter    = n_iter,
    warmup  = n_warmup,
    cores   = n_cores,
    seed    = 42,
    save_pars = save_pars(all = TRUE),
    control   = list(adapt_delta = adapt_d, max_treedepth = 12),
    silent    = 2,
    refresh   = 100
  )
  
  # ── divergences
  nuts  <- nuts_params(m)
  n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  cat("  Divergences:", n_div, "\n")
  
  if (n_div > 10) {
    cat("  Refitting with adapt_delta = 0.99...\n")
    m <- update(m,
                control   = list(adapt_delta = 0.99, max_treedepth = 15),
                save_pars = save_pars(all = TRUE),
                silent    = 2, refresh = 100)
    nuts  <- nuts_params(m)
    n_div <- sum(nuts$Value[nuts$Parameter == "divergent__"])
    cat("  Divergences after refit:", n_div, "\n")
  }
  
  # ── visual diagnostic
  p_trace <- mcmc_trace(
    as.array(m),
    pars         = grep("^sd_", variables(m), value = TRUE),
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
  
  # ── extraction
  cat("  Extracting parameters and BLUPs...\n")
  lista_parametros[[v]] <- extrair_parametros(m, v, pad$media, pad$dp)
  lista_blups_G[[v]]    <- extrair_blups_G(m, v, pad$media, pad$dp)
  lista_blups_GxA[[v]]  <- extrair_blups_GxA(m, v)
  
  lista_modelos[[v]] <- list(modelo = m, media_orig = pad$media, dp_orig = pad$dp)
  
  cat("  Finished:", v, "\n\n")
}

# ── consolidate and save ────────────────────────────────────────────────────────
cat("── Saving results...\n")

params_final  <- list_rbind(lista_parametros)
blups_G_final <- list_rbind(lista_blups_G)
blups_GxA_final <- list_rbind(lista_blups_GxA)

# add group to blups
grupo_ref <- dados |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

blups_G_final <- blups_G_final |>
  left_join(grupo_ref, by = "genotipo")

blups_GxA_final <- blups_GxA_final |>
  left_join(grupo_ref, by = "genotipo")

# save to xlsx
write_xlsx(
  list(
    parametros_geneticos = params_final,
    blups_genotipo       = blups_G_final,
    blups_gxa            = blups_GxA_final
  ),
  paste0(path_tbl, "03_resultados_bayesianos.xlsx")
)

saveRDS(lista_modelos, "data/modelos_finais.rds")

cat("\nFinished. Models saved in data/modelos_finais.rds\n")
cat("Variables processed:", length(lista_modelos), "\n")
