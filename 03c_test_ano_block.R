# 03c_test_ano_block.R
# Tests whether the year:block interaction (years nested within block) 
# improves the model fit for each response trait
#
# Reference: model from the Theobroma grandiflorum paper with repeated 
# measures (Years:Replicates as a fixed effect)
#
# Comparison via LRT with ML:
#   M1: y ~ ano + block + (1|genotipo) + (1|genotipo:ano)
#   M2: y ~ ano + block + ano:block + (1|genotipo) + (1|genotipo:ano)
#
# Outputs:
#   outputs/tables/03c_test_ano_block.xlsx
#   outputs/figures/03c_aic_comparacao.png
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lme4)
library(lmerTest)
library(ggplot2)
library(writexl)

# ── data and paths ──────────────────────────────────────────────────────────────
dados    <- readRDS("data/dados_clean.rds")
path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

vars_resp <- c("per_grao", "per_palha", "mcm_mgb",
               "mcm_saca", "vcm_saca", "vcm_mcm")

ctrl_ml <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# ── labels for plots ───────────────────────────────────────────────────────
var_labels <- c(
  per_grao  = "% grain",
  per_palha = "% husk",
  mcm_mgb   = "FWM/GW",
  mcm_saca  = "FWM/bag",
  vcm_saca  = "FVol/bag",
  vcm_mcm   = "FVol/FWM"
)

# ── test function ────────────────────────────────────────────────────────────
testar_ano_block <- function(var, data) {
  
  cat("────────────────────────────\n")
  cat("Variable:", var, "\n")
  
  f_m1 <- as.formula(paste(
    var, "~ ano + block + (1|genotipo) + (1|genotipo:ano)"
  ))
  f_m2 <- as.formula(paste(
    var, "~ ano + block + ano:block + (1|genotipo) + (1|genotipo:ano)"
  ))
  
  # fitting with ML for correct fixed effects LRT
  m1 <- tryCatch(
    lmer(f_m1, data = data, REML = FALSE, control = ctrl_ml),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lmer(f_m1, data = data, REML = FALSE, control = ctrl_ml)
    )
  )
  
  m2 <- tryCatch(
    lmer(f_m2, data = data, REML = FALSE, control = ctrl_ml),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      lmer(f_m2, data = data, REML = FALSE, control = ctrl_ml)
    )
  )
  
  if (is.null(m1) || is.null(m2)) {
    cat("  ERROR — model failed to converge\n")
    return(tibble(
      variavel       = var,
      aic_m1         = NA, bic_m1 = NA, loglik_m1 = NA,
      aic_m2         = NA, bic_m2 = NA, loglik_m2 = NA,
      delta_aic      = NA, delta_bic = NA,
      lrt_stat       = NA, lrt_df = NA, lrt_pval = NA,
      ano_block_sig  = NA, modelo_preferido = NA
    ))
  }
  
  # LRT
  lrt <- anova(m2, m1, refit = FALSE)
  
  delta_aic <- AIC(m2) - AIC(m1)
  delta_bic <- BIC(m2) - BIC(m1)
  lrt_stat  <- lrt[["Chisq"]][2]
  lrt_df    <- lrt[["Df"]][2]
  lrt_pval  <- lrt[["Pr(>Chisq)"]][2]
  
  sig <- lrt_pval < 0.05
  
  cat("  ΔAIC (M2 - M1):", round(delta_aic, 2),
      "| ΔBIC:", round(delta_bic, 2), "\n")
  cat("  LRT χ²:", round(lrt_stat, 3),
      "| df:", lrt_df,
      "| p:", round(lrt_pval, 4), "\n")
  cat("  year:block significant:", sig, "\n\n")
  
  # coefficients of the year:block term in M2
  coef_m2 <- summary(m2)$coefficients
  coef_ab <- coef_m2[grep("ano.*block|block.*ano", rownames(coef_m2)), ,
                     drop = FALSE]
  
  if (nrow(coef_ab) > 0) {
    cat("  Coefficients year:block:\n")
    print(round(coef_ab, 4))
    cat("\n")
  }
  
  tibble(
    variavel         = var,
    aic_m1           = AIC(m1),
    bic_m1           = BIC(m1),
    loglik_m1        = logLik(m1)[1],
    aic_m2           = AIC(m2),
    bic_m2           = BIC(m2),
    loglik_m2        = logLik(m2)[1],
    delta_aic        = delta_aic,
    delta_bic        = delta_bic,
    lrt_stat         = lrt_stat,
    lrt_df           = lrt_df,
    lrt_pval         = lrt_pval,
    ano_block_sig    = sig,
    modelo_preferido = if_else(sig & delta_aic < -2,
                               "M2 (with year:block)",
                               "M1 (without year:block)")
  )
}

# ── main loop ─────────────────────────────────────────────────────────────
cat("\n── Year:block test ──\n\n")

resultados <- map_dfr(vars_resp, testar_ano_block, data = dados)

cat("── Summary ──\n")
print(resultados |>
        select(variavel, delta_aic, delta_bic,
               lrt_pval, ano_block_sig, modelo_preferido))

# ── visualization ───────────────────────────────────────────────────────────────

df_plot <- resultados |>
  mutate(
    variavel_lab = factor(var_labels[variavel], levels = var_labels),
    sig_label    = if_else(ano_block_sig,
                           paste0("p = ", round(lrt_pval, 3), " *"),
                           paste0("p = ", round(lrt_pval, 3)))
  )

p_aic <- ggplot(df_plot, aes(x = variavel_lab, y = delta_aic,
                             fill = ano_block_sig)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = -2, linetype = "dashed",
             colour = "gray40", linewidth = 0.4) +
  geom_hline(yintercept =  0, colour = "gray60", linewidth = 0.3) +
  geom_text(aes(label = sig_label,
                vjust = if_else(delta_aic < 0, 1.4, -0.5)),
            size = 3.2, colour = "gray20") +
  scale_fill_manual(
    values = c("TRUE" = "#185FA5", "FALSE" = "#B4B2A9"),
    labels = c("TRUE" = "significant (p < 0.05)",
               "FALSE" = "not significant"),
    name   = "year:block"
  ) +
  scale_y_continuous(
    name = "ΔAIC (M2 with year:block − M1 without year:block)",
    breaks = scales::pretty_breaks(6)
  ) +
  labs(
    x        = NULL,
    title    = "Test of year:block effect on model fixed structure",
    subtitle = "Negative ΔAIC indicates improvement with year:block inclusion\ndashed line = ΔAIC = −2 (conventional threshold)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    axis.text.x        = element_text(angle = 30, hjust = 1)
  )

ggsave(paste0(path_fig, "03c_aic_comparacao.png"),
       p_aic, width = 8, height = 5, dpi = 300)

# ── save ─────────────────────────────────────────────────────────────────────
write_xlsx(
  list(resultados = resultados),
  paste0(path_tbl, "03c_test_ano_block.xlsx")
)

cat("\nFinished. Results in outputs/tables/03c_test_ano_block.xlsx\n")
