# 06_multitrait_sommer.R
# Modelo multi-trait para eficiência de beneficiamento em C. canephora
#
# Objetivo: estimar correlações genéticas entre per_grao e mcm_saca
# e validar se as duas variáveis capturam dimensões independentes
#
# Estrutura do modelo:
#   Fixos:     ano + block
#   Aleatórios: genotipo (US), genotipo:ano (US), residual (US)
#
# Seleção de modelo via LRT (df teóricos):
#   M_diag:      diagonal — sem covariâncias entre traits
#   M_us_diagGxA: US para G, diagonal para GxA
#   M_us:         US completo — modelo selecionado
#
# Nota sobre FA: com p=2 traits, FA(k) tem mais parâmetros que US
# (FA(1)=4, FA(2)=5 vs US=3) — FA não é parcimonioso para p=2
# FA será avaliado quando o modelo for expandido para p=5 com ASReml
#
# Saídas:
#   outputs/tables/06_multitrait_results.xlsx
#   outputs/figures/06a_correlacoes_geneticas.png
#   outputs/figures/06b_blups_multitrait.png

# ── pacotes ────────────────────────────────────────────────────────────────────
library(sommer)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(writexl)
library(readxl)

# ── dados e paths ──────────────────────────────────────────────────────────────
dados    <- readRDS("data/dados_clean.rds")
dados_df <- as.data.frame(dados)   # sommer requer data.frame

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# ── grupo de referência para plots ────────────────────────────────────────────
grupo_ref <- dados_df |>
  distinct(genotipo, grupo) |>
  mutate(genotipo = as.character(genotipo))

# ── 1. ajuste dos modelos ──────────────────────────────────────────────────────
cat("\n── Ajustando modelos multi-trait ──\n\n")

# M_diag — diagonal (sem covariâncias entre traits)
cat("M_diag: diagonal...\n")
m_diag <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = diag(2)) +
    vsr(genotipo:ano, Gtc = diag(2)),
  rcov   = ~ vsr(units,        Gtc = diag(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Convergiu:", m_diag$convergence,
    "| AIC:", round(m_diag$AIC, 2), "\n")

# M_us_diagGxA — US para G, diagonal para GxA
cat("M_us_diagGxA: US(G) + diag(GxA)...\n")
m_us_diagGxA <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = unsm(2)) +
    vsr(genotipo:ano, Gtc = diag(2)),
  rcov   = ~ vsr(units,        Gtc = unsm(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Convergiu:", m_us_diagGxA$convergence,
    "| AIC:", round(m_us_diagGxA$AIC, 2), "\n")

# M_us — US completo (modelo selecionado)
cat("M_us: US completo...\n")
m_us <- mmer(
  cbind(per_grao, mcm_saca) ~ ano + block,
  random = ~ vsr(genotipo,     Gtc = unsm(2)) +
    vsr(genotipo:ano, Gtc = unsm(2)),
  rcov   = ~ vsr(units,        Gtc = unsm(2)),
  data   = dados_df,
  verbose = FALSE
)
cat("  Convergiu:", m_us$convergence,
    "| AIC:", round(m_us$AIC, 2), "\n\n")

# ── 2. seleção de modelo via LRT ───────────────────────────────────────────────
cat("── Seleção de modelo via LRT ──\n")

# logLik da última iteração (linha 1, última coluna do monitor)
ll_diag    <- m_diag$monitor[1,    ncol(m_diag$monitor)]
ll_us_dGxA <- m_us_diagGxA$monitor[1, ncol(m_us_diagGxA$monitor)]
ll_us      <- m_us$monitor[1,      ncol(m_us$monitor)]

# LRT com df teóricos
# M_diag vs M_us: +3 parâmetros (covG, covGxA, covResidual)
# M_us_diagGxA vs M_us: +1 parâmetro (covGxA)
# M_diag vs M_us_diagGxA: +2 parâmetros (covG, covResidual)
lrt_diag_us   <- 2 * (ll_us - ll_diag)
lrt_dGxA_us   <- 2 * (ll_us - ll_us_dGxA)
lrt_diag_dGxA <- 2 * (ll_us_dGxA - ll_diag)

p_diag_us     <- pchisq(lrt_diag_us,   df = 3, lower.tail = FALSE)
p_dGxA_us     <- pchisq(lrt_dGxA_us,   df = 1, lower.tail = FALSE)
p_diag_dGxA   <- pchisq(lrt_diag_dGxA, df = 2, lower.tail = FALSE)

df_selecao <- tibble(
  comparacao   = c("M_diag vs M_us",
                   "M_us_diagGxA vs M_us",
                   "M_diag vs M_us_diagGxA"),
  loglik_red   = c(ll_diag,    ll_us_dGxA, ll_diag),
  loglik_full  = c(ll_us,      ll_us,      ll_us_dGxA),
  chi2         = c(lrt_diag_us, lrt_dGxA_us, lrt_diag_dGxA),
  df           = c(3, 1, 2),
  pval         = c(p_diag_us, p_dGxA_us, p_diag_dGxA)
)

cat("\n")
print(df_selecao)

df_aic <- tibble(
  modelo       = c("M_diag", "M_us_diagGxA", "M_us"),
  loglik       = c(ll_diag, ll_us_dGxA, ll_us),
  aic          = c(m_diag$AIC, m_us_diagGxA$AIC, m_us$AIC),
  delta_aic    = c(m_diag$AIC - m_us$AIC,
                   m_us_diagGxA$AIC - m_us$AIC,
                   0),
  selecionado  = c(FALSE, FALSE, TRUE)
)

cat("\n── Resumo AIC ──\n")
print(df_aic)

# ── 3. componentes de variância e correlações ─────────────────────────────────
cat("\n── Componentes de variância — M_us ──\n")
vc <- summary(m_us)$varcomp
print(vc)

# extrair covariâncias e calcular correlações
s2_G_per   <- vc["u:genotipo.per_grao-per_grao",    "VarComp"]
s2_G_mcm   <- vc["u:genotipo.mcm_saca-mcm_saca",    "VarComp"]
cov_G      <- vc["u:genotipo.per_grao-mcm_saca",     "VarComp"]
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
  componente          = c("Genótipo (G)",
                          "Genótipo × Ano (GxA)",
                          "Residual"),
  sigma2_per_grao     = c(s2_G_per, s2_GxA_per, s2_e_per),
  sigma2_mcm_saca     = c(s2_G_mcm, s2_GxA_mcm, s2_e_mcm),
  covariancia         = c(cov_G, cov_GxA, cov_e),
  correlacao          = c(rG, rGxA, re)
)

cat("\n── Correlações entre per_grao e mcm_saca ──\n")
print(df_correlacoes)

# ── 4. BLUPs multi-trait ──────────────────────────────────────────────────────
cat("\n── Extraindo BLUPs multi-trait ──\n")

blup_gen <- tibble(
  genotipo      = names(m_us$U$`u:genotipo`$per_grao),
  blup_per_grao = as.numeric(m_us$U$`u:genotipo`$per_grao),
  blup_mcm_saca = as.numeric(m_us$U$`u:genotipo`$mcm_saca)
) |>
  left_join(grupo_ref, by = "genotipo") |>
  mutate(
    blup_per_z  =  scale(blup_per_grao)[, 1],
    blup_mcm_z  = -scale(blup_mcm_saca)[, 1],
    indice_mt   = (blup_per_z + blup_mcm_z) / 2
  ) |>
  arrange(desc(indice_mt))

cat("── Top 10 índice multi-trait ──\n")
blup_gen |>
  select(genotipo, grupo, blup_per_grao,
         blup_mcm_saca, indice_mt) |>
  head(10) |>
  print()

cat("Dimensões BLUPs genótipo:", nrow(blup_gen), "×", ncol(blup_gen), "\n")

# verificar top genótipos por eficiência combinada
# índice simples: padronizar e combinar (per_grao alto + mcm_saca baixo)
blup_gen <- blup_gen |>
  mutate(
    blup_per_z  =  scale(blup_per_grao)[, 1],
    blup_mcm_z  = -scale(blup_mcm_saca)[, 1],  # negativo: menor = melhor
    indice_mt   = (blup_per_z + blup_mcm_z) / 2
  ) |>
  arrange(desc(indice_mt))

cat("\n── Top 10 — índice multi-trait ──\n")
blup_gen |>
  select(genotipo, grupo, blup_per_grao,
         blup_mcm_saca, indice_mt) |>
  head(10) |>
  print()

# ── 5. visualizações ──────────────────────────────────────────────────────────

cores_grupo <- c("Conilon" = "#185FA5", "Robusta" = "#D85A30")

# 5a. scatter BLUPs multi-trait
p_blups <- ggplot(blup_gen,
                  aes(x = blup_per_grao, y = blup_mcm_saca,
                      colour = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "gray60", linewidth = 0.4) +
  geom_point(aes(size = abs(indice_mt)), alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = blup_gen |> filter(abs(indice_mt) > 0.5),
    aes(label = genotipo),
    size = 2.8, max.overlaps = 20, show.legend = FALSE
  ) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "gray30", linewidth = 0.5) +
  scale_colour_manual(values = cores_grupo, name = "grupo botânico") +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  annotate(
    "text",
    x = min(blup_gen$blup_per_grao) * 0.85,
    y = max(blup_gen$blup_mcm_saca) * 0.90,
    label = paste0("rG = ", round(rG, 3)),
    size  = 4.5, colour = "gray20", hjust = 0, fontface = "bold"
  ) +
  labs(
    x        = "BLUP — % grão (desvio da média)",
    y        = "BLUP — MCM/saca kg (desvio da média)",
    title    = "BLUPs multi-trait — per_grao × mcm_saca",
    subtitle = paste0("Correlação genética rG = ", round(rG, 3),
                      " | Modelo US | sommer 4.4.x")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top")

ggsave(paste0(path_fig, "06a_blups_multitrait_scatter.png"),
       p_blups, width = 9, height = 8, dpi = 300)
cat("06a salvo\n")

# 5b. correlações por componente
df_cor_plot <- df_correlacoes |>
  mutate(
    componente = factor(componente,
                        levels = c("Genótipo (G)",
                                   "Genótipo × Ano (GxA)",
                                   "Residual")),
    cor_label  = round(correlacao, 3),
    cor_sig    = case_when(
      componente == "Genótipo (G)"         ~ "p < 0.001",
      componente == "Genótipo × Ano (GxA)" ~ "p < 0.001",
      TRUE                                 ~ ""
    )
  )

p_cor <- ggplot(df_cor_plot,
                aes(x = componente, y = correlacao,
                    fill = correlacao)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0, colour = "gray40", linewidth = 0.4) +
  geom_text(aes(label = paste0(cor_label, "\n", cor_sig),
                vjust = ifelse(correlacao < 0, 1.3, -0.3)),
            size = 3.5, colour = "gray20") +
  scale_fill_gradient2(
    low      = "#D85A30",
    mid      = "white",
    high     = "#185FA5",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "correlação"
  ) +
  scale_y_continuous(limits = c(-1.1, 0.2)) +
  labs(
    x        = NULL,
    y        = "correlação entre per_grao e mcm_saca",
    title    = "Correlações entre caracteres por componente de variância",
    subtitle = "Modelo multi-trait US — sommer 4.4.x"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

ggsave(paste0(path_fig, "06b_correlacoes_por_componente.png"),
       p_cor, width = 8, height = 6, dpi = 300)

# ── 6. salvar resultados ───────────────────────────────────────────────────────
write_xlsx(
  list(
    selecao_modelos  = df_aic,
    lrt_comparacoes  = df_selecao,
    varcomp_mus      = as_tibble(vc, rownames = "parametro"),
    correlacoes      = df_correlacoes,
    blups_multitrait = blup_gen
  ),
  paste0(path_tbl, "06_multitrait_results.xlsx")
)

cat("\n06_multitrait_sommer.R concluído.\n")
cat("Figuras salvas em outputs/figures/\n")
cat("Tabelas salvas em outputs/tables/06_multitrait_results.xlsx\n")
