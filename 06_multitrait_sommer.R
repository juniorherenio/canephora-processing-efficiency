# 06_multitrait_sommer.R
# Multi-trait model for processing efficiency in C. canephora
#
# Objective: estimate genetic correlations between per_grao and mcm_saca
# and validate whether both traits capture independent dimensions
#
# Model structure:
#   Fixed:     ano + block
#   Random:    genotipo (US), genotipo:ano (US), residual (US)
#
# Model selection via LRT (theoretical df):
#   M_diag:          diagonal — no correlations between traits
#   M_us_diagGxA:    US for G, diagonal for GxA
#   M_us:            full US — selected model
#
# Note on FA: with p=2 traits, FA(k) has more parameters than US 
# (FA(1)=4, FA(2)=5 vs US=3) — FA is not parsimonious for p=2. 
# FA will be evaluated when the model is expanded to p=5 with ASReml.
#
# Outputs:
#   outputs/tables/06_multitrait_results.xlsx
#   outputs/figures/06a_correlacoes_geneticas.png
#   outputs/figures/06b_blups_multitrait.png
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(sommer)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(writexl)
library(readxl)

# ── data and paths ──────────────────────────────────────────────────────────────
dados    <- readRDS("data/dados_clean.rds")
dados_df <- as.data.frame(dados)   # sommer requires data.frame

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# ── reference group for plots ────────────────────────────────────────────
grupo_ref <- dados_df |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

# ── 1. fitting models ──────────────────────────────────────────────────────
cat("\n── Fitting multi-trait models ──\n\n")

# M_diag — diagonal (no correlations between traits)
cat("M_diag: diagonal...\n")
m_diag <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = diag(2)) +
    vsr(genotipo:ano, Gtc = diag(2)),
  rcov   = ~ vsr(units,        Gtc = diag(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Converged:", m_diag$convergence,
    "| AIC:", round(m_diag$AIC, 2), "\n")

# M_us_diagGxA — US for G, diagonal for GxA
cat("M_us_diagGxA: US(G) + diag(GxA)...\n")
m_us_diagGxA <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = unsm(2)) +
    vsr(genotipo:ano, Gtc = diag(2)),
  rcov   = ~ vsr(units,        Gtc = unsm(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Converged:", m_us_diagGxA$convergence,
    "| AIC:", round(m_us_diagGxA$AIC, 2), "\n")

# M_us — full US (selected model)
cat("M_us: full US...\n")
m_us <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = unsm(2)) +
    vsr(genotipo:ano, Gtc = unsm(2)),
  rcov   = ~ vsr(units,        Gtc = unsm(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Converged:", m_us$convergence,
    "| AIC:", round(m_us$AIC, 2), "\n\n")

# ── 2. model selection via LRT ───────────────────────────────────────────────
cat("── Model selection via LRT ──\n")

# logLik from the last iteration (row 1, last column of monitor)
ll_diag    <- m_diag$monitor[1,    ncol(m_diag$monitor)]
ll_us_dGxA <- m_us_diagGxA$monitor[1, ncol(m_us_diagGxA$monitor)]
ll_us      <- m_us$monitor[1,      ncol(m_us$monitor)]

# LRT with theoretical df
lrt_diag_us   <- 2 * (ll_us - ll_diag)
lrt_dGxA_us   <- 2 * (ll_us - ll_us_dGxA)
lrt_diag_dGxA <- 2 * (ll_us_dGxA - ll_diag)

p_diag_us     <- pchisq(lrt_diag_us,   df = 3, lower.tail = FALSE)
p_dGxA_us     <- pchisq(lrt_dGxA_us,   df = 1, lower.tail = FALSE)
p_diag_dGxA   <- pchisq(lrt_diag_dGxA, df = 2, lower.tail = FALSE)

df_selecao <- tibble(
  comparison   = c("M_diag vs M_us",
                   "M_us_diagGxA vs M_us",
                   "M_diag vs M_us_diagGxA"),
  loglik_red   = c(ll_diag,    ll_us_dGxA, ll_diag),
  loglik_full  = c(ll_us,      ll_us,      ll_us_dGxA),
  chi2         = c(lrt_diag_us, lrt_dGxA_us, lrt_diag_dGxA),
  df           = c(3, 1, 2),
  pval         = c(p_diag_us, p_dGxA_us, p_diag_dGxA)
)

print(df_selecao)

df_aic <- tibble(
  model        = c("M_diag", "M_us_diagGxA", "M_us"),
  loglik       = c(ll_diag, ll_us_dGxA, ll_us),
  aic          = c(m_diag$AIC, m_us_diagGxA$AIC, m_us$AIC),
  delta_aic    = c(m_diag$AIC - m_us$AIC,
                   m_us_diagGxA$AIC - m_us$AIC,
                   0),
  selected     = c(FALSE, FALSE, TRUE)
)

cat("\n── AIC summary ──\n")
print(df_aic)

# ── 3. variance components and correlations ─────────────────────────────────
cat("\n── Variance components — M_us ──\n")
vc <- summary(m_us)$varcomp
print(vc)

# extract covariances and calculate correlations
s2_G_per   <- vc["u:genotipo.per_grao-per_grao",    "VarComp"]
s2_G_mcm   <- vc["u:genotipo.mcm_saca-mcm_saca",    "VarComp"]
cov_G      <- vc["u:genotipo.per_grao-mcm_saca",    "VarComp"]
rG         <- cov_G / sqrt(s2_G_per * s2_G_mcm)

s2_GxA_per <- vc["u:genotipo:ano.per_grao-per_grao", "VarComp"]
s2_GxA_mcm <- vc["u:genotipo:ano.mcm_saca-mcm_saca", "VarComp"]
cov_GxA    <- vc["u:genotipo:ano.per_grao-mcm_saca",  "VarComp"]
rGxA       <- cov_GxA / sqrt(s2_GxA_per * s2_GxA_mcm)

s2_e_per   <- vc["u:units.per_grao-per_grao",  "VarComp"]
s2_e_mcm   <- vc["u:units.mcm_saca-mcm_saca",  "VarComp"]
cov_e      <- vc["u:units.per_grao-mcm_saca",   "VarComp"]
re         <- cov_e / sqrt(s2_e_per * s2_e_mcm)

df_correlacoes <- tibble(
  component     = c("Genotype (G)",
                    "Genotype × Year (GxA)",
                    "Residual"),
  sigma2_per_grain   = c(s2_G_per, s2_GxA_per, s2_e_per),
  sigma2_fwm_bag     = c(s2_G_mcm, s2_GxA_mcm, s2_e_mcm),
  covariance         = c(cov_G, cov_GxA, cov_e),
  correlation        = c(rG, rGxA, re)
)

cat("\n── Correlations between grain proportion and FWM/bag ──\n")
print(df_correlacoes)

# ── 4. multi-trait BLUPs ──────────────────────────────────────────────────────
cat("\n── Extracting multi-trait BLUPs ──\n")

blup_gen <- tibble(
  genotype      = names(m_us$U$`u:genotipo`$per_grao),
  blup_per_grain = as.numeric(m_us$U$`u:genotipo`$per_grao),
  blup_fwm_bag   = as.numeric(m_us$U$`u:genotipo`$mcm_saca)
) |>
  left_join(grupo_ref, by = c("genotype" = "genotipo")) |>
  mutate(
    blup_per_z  =  scale(blup_per_grain)[, 1],
    blup_mcm_z  = -scale(blup_fwm_bag)[, 1],
    index_mt    = (blup_per_z + blup_mcm_z) / 2
  ) |>
  arrange(desc(index_mt))

cat("── Top 10 multi-trait index ──\n")
blup_gen |>
  select(genotype, grupo, blup_per_grain,
         blup_fwm_bag, index_mt) |>
  head(10) |>
  print()

# ── 5. visualizations ──────────────────────────────────────────────────────────

cores_grupo <- c("Conilon" = "#185FA5", "Robusta" = "#D85A30")

p_blups <- ggplot(blup_gen,
                  aes(x = blup_per_grain, y = blup_fwm_bag,
                      colour = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray60", linewidth = 0.4) +
  geom_point(aes(size = abs(index_mt)), alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = blup_gen |> filter(abs(index_mt) > 0.5),
    aes(label = genotype),
    size = 2.8, max.overlaps = 20, show.legend = FALSE
  ) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "gray30", linewidth = 0.5) +
  scale_colour_manual(values = cores_grupo, name = "Botanical group") +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  annotate(
    "text",
    x = min(blup_gen$blup_per_grain) * 0.85,
    y = max(blup_gen$blup_fwm_bag) * 0.90,
    label = paste0("rG = ", round(rG, 3)),
    size  = 4.5, colour = "gray20", hjust = 0, fontface = "bold"
  ) +
  labs(
    x        = "BLUP — % grain (deviation from mean)",
    y        = "BLUP — FWM/bag kg (deviation from mean)",
    title    = "Multi-trait BLUPs — % grain × FWM/bag",
    subtitle = paste0("Genetic correlation rG = ", round(rG, 3),
                      " | US Model | sommer 4.4.x")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top")

ggsave(paste0(path_fig, "06a_blups_multitrait_scatter.png"),
       p_blups, width = 9, height = 8, dpi = 300)

# ── 6. save results ───────────────────────────────────────────────────────
write_xlsx(
  list(
    model_selection    = df_aic,
    lrt_comparisons    = df_selecao,
    varcomp_mus        = as_tibble(vc, rownames = "parameter"),
    correlations       = df_correlacoes,
    blups_multitrait   = blup_gen
  ),
  paste0(path_tbl, "06_multitrait_results.xlsx")
)

cat("\n06_multitrait_sommer.R completed.\n")
