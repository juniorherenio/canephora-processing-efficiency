# 04_stability.R — Bloco 1
# Carregamento e reconstrução dos valores genotípicos por ano
# a partir dos BLUPs bayesianos
#
# Estrutura:
#   valor_genotipo_ano = media_geral + BLUP_G + BLUP_GxA
#
# Saída deste bloco:
#   df_valores — tibble com valor genotípico por genótipo × ano × variável
#   df_blups_G — tibble com desempenho médio por genótipo × variável

# ── pacotes ────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readxl)
library(ggplot2)
library(patchwork)
library(writexl)
library(brms)      # Necessário se os modelos foram ajustados com brms
library(posterior) # Pacote que contém a função as_draws_df()) 

# ── dados e paths ──────────────────────────────────────────────────────────────
dados   <- readRDS("data/dados_clean.rds")
modelos <- readRDS("data/modelos_finais.rds")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# variáveis de interesse — excluindo vcm_mcm (sem estrutura genética)
vars_resp <- c("per_grao", "per_palha", "mcm_mgb", "mcm_saca", "vcm_saca")

# labels para plots
var_labels <- c(
  per_grao  = "% grão",
  per_palha = "% palha",
  mcm_mgb   = "MCM/MGB",
  mcm_saca  = "MCM/saca",
  vcm_saca  = "VCM/saca"
)

# direção desejável por variável
# "alto" = maior é melhor | "baixo" = menor é melhor
var_direcao <- c(
  per_grao  = "alto",    # mais grão = melhor
  per_palha = "baixo",   # menos palha = melhor
  mcm_mgb   = "baixo",   # menos massa por grão = melhor
  mcm_saca  = "baixo",   # menos massa por saca = melhor
  vcm_saca  = "baixo"    # menos volume por saca = melhor
)

# ── carregar BLUPs já calculados ───────────────────────────────────────────────
blups_G <- read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "blups_genotipo"
) |>
  filter(variavel %in% vars_resp)

blups_GxA <- read_xlsx(
  "outputs/tables/03_resultados_bayesianos.xlsx",
  sheet = "blups_gxa"
) |>
  filter(variavel %in% vars_resp)

# ── reconstruir valores genotípicos por ano ────────────────────────────────────
# media_geral por variável (escala original)
media_geral <- dados |>
  summarise(across(all_of(vars_resp), \(x) mean(x, na.rm = TRUE))) |>
  pivot_longer(everything(),
               names_to  = "variavel",
               values_to = "media_geral")

# juntar BLUP_G + BLUP_GxA + media_geral
df_valores <- blups_GxA |>
  dplyr::select(variavel, genotipo, ano, blup_gxa_med, grupo) |>  # Correção aqui
  left_join(
    blups_G |> dplyr::select(variavel, genotipo, blup_orig_med),  # E correção aqui
    by = c("variavel", "genotipo")
  ) |>
  left_join(media_geral, by = "variavel") |>
  mutate(
    # valor genotípico reconstruído na escala original
    valor_ano = media_geral + blup_orig_med + blup_gxa_med,
    ano       = as.character(ano)
  )

# verificação: estrutura esperada
cat("\n── Estrutura df_valores ──\n")
cat("Dimensões:", nrow(df_valores), "×", ncol(df_valores), "\n")
cat("Variáveis:", unique(df_valores$variavel), "\n")
cat("Anos:", unique(df_valores$ano), "\n")
cat("Genótipos:", n_distinct(df_valores$genotipo), "\n\n")

# verificação: médias reconstruídas vs médias observadas
cat("── Verificação: médias reconstruídas vs observadas ──\n")
check <- df_valores |>
  group_by(variavel) |>
  summarise(
    media_reconstruida = mean(valor_ano),
    .groups = "drop"
  ) |>
  left_join(media_geral, by = "variavel") |>
  mutate(diferenca = media_reconstruida - media_geral)

print(check)

# ── df_blups_G — desempenho médio por genótipo ─────────────────────────────────
df_blups_G <- blups_G |>
  dplyr::select(variavel, genotipo, grupo,
         blup_orig_med, blup_orig_q025, blup_orig_q975,
         prob_positivo) |>
  left_join(media_geral, by = "variavel") |>
  mutate(
    valor_medio = media_geral + blup_orig_med,
    # sinal de seleção: +1 se direção "alto", -1 se direção "baixo"
    # usado para padronizar ranking (maior sempre = melhor)
    sinal = if_else(var_direcao[variavel] == "alto", 1, -1),
    # blup_sinal: positivo sempre indica direção desejável
    blup_sinal = blup_orig_med * sinal
  )

cat("\n── Primeiras linhas df_blups_G ──\n")
print(df_blups_G |>
        filter(variavel == "mcm_saca") |>
        arrange(blup_sinal) |>
        dplyr::select(genotipo, grupo, blup_orig_med, blup_sinal,
               valor_medio, prob_positivo) |>
        head(10))

cat("\nBloco 1 concluído — df_valores e df_blups_G prontos\n")

# 04_stability.R — Bloco 2
# Ecovalência de Wricke adaptada aos BLUPs bayesianos
#
# W_i = soma dos quadrados dos BLUP_GxA por genótipo
# Com dois anos: W_i = BLUP_GxA_2023² + BLUP_GxA_2024²
#
# Interpretação: menor W_i = maior estabilidade temporal
# Genótipos com W_i próximo de zero mantêm desvio consistente entre anos

# ── calcular ecovalência ───────────────────────────────────────────────────────
df_wricke <- df_valores |>
  group_by(variavel, genotipo, grupo) |>
  summarise(
    blup_gxa_2023  = blup_gxa_med[ano == "2023"],
    blup_gxa_2024  = blup_gxa_med[ano == "2024"],
    # ecovalência — soma dos quadrados dos desvios GxA
    wricke         = blup_gxa_2023^2 + blup_gxa_2024^2,
    # amplitude — diferença absoluta entre anos (complementar)
    amplitude_gxa  = abs(blup_gxa_2023 - blup_gxa_2024),
    # inversão de sinal entre anos
    inversao       = sign(blup_gxa_2023) != sign(blup_gxa_2024),
    .groups        = "drop"
  ) |>
  # adicionar desempenho médio
  left_join(
    df_blups_G |>
      dplyr::select(variavel, genotipo, blup_orig_med, blup_sinal, valor_medio),
    by = c("variavel", "genotipo")
  ) |>
  # ecovalência relativa — proporção da ecovalência total por variável
  group_by(variavel) |>
  mutate(
    wricke_total    = sum(wricke),
    wricke_rel      = wricke / wricke_total * 100,
    # ranking de estabilidade: 1 = mais estável
    rank_estab      = rank(wricke),
    # ranking de desempenho: 1 = melhor (considera direção)
    rank_desempenho = rank(-blup_sinal)
  ) |>
  ungroup()

# ── resumo por variável ────────────────────────────────────────────────────────
cat("\n── Ecovalência — resumo por variável ──\n")
df_wricke |>
  group_by(variavel) |>
  summarise(
    wricke_min    = min(wricke),
    wricke_max    = max(wricke),
    wricke_media  = mean(wricke),
    n_inversao    = sum(inversao),
    prop_inversao = mean(inversao),
    .groups       = "drop"
  ) |>
  print()

# ── top 10 mais estáveis por variável ─────────────────────────────────────────
cat("\n── Top 10 mais estáveis — mcm_saca ──\n")
df_wricke |>
  filter(variavel == "mcm_saca") |>
  arrange(wricke) |>
  dplyr::select(genotipo, grupo, wricke, wricke_rel,
         amplitude_gxa, inversao,
         blup_orig_med, rank_estab, rank_desempenho) |>
  head(10) |>
  print()

cat("\n── Top 10 mais instáveis — mcm_saca ──\n")
df_wricke |>
  filter(variavel == "mcm_saca") |>
  arrange(desc(wricke)) |>
  dplyr::select(genotipo, grupo, wricke, wricke_rel,
         amplitude_gxa, inversao,
         blup_orig_med, rank_estab, rank_desempenho) |>
  head(10) |>
  print()

# ── genótipos que invertem sinal entre anos ────────────────────────────────────
cat("\n── Genótipos com inversão de sinal entre anos — mcm_saca ──\n")
df_wricke |>
  filter(variavel == "mcm_saca", inversao == TRUE) |>
  arrange(wricke) |>
  dplyr::select(genotipo, grupo, blup_gxa_2023, blup_gxa_2024,
         wricke, amplitude_gxa) |>
  print()

cat("\nBloco 2 concluído — df_wricke pronto\n")

# 04_stability.R — Bloco 3
# Probabilidade bayesiana de superioridade consistente
#
# Para cada genótipo e cada intensidade de seleção (10%, 20%, 25%),
# calcula a probabilidade de estar no grupo superior
# simultaneamente nos dois anos
#
# Derivado diretamente das amostras MCMC — sem suposições adicionais

# ── recarregar modelos para acessar amostras MCMC ─────────────────────────────
cat("\n── Calculando probabilidades bayesianas de superioridade ──\n\n")

intensidades <- c(0.10, 0.20, 0.25)
n_gen        <- 48L

# ── função principal ───────────────────────────────────────────────────────────
calcular_prob_consistente <- function(var, intensidades, modelos, media_orig) {
  
  cat("  Variável:", var, "\n")
  
  m    <- modelos[[var]]$modelo
  dp   <- modelos[[var]]$dp_orig
  mu   <- media_orig |> filter(variavel == var) |> pull(media_geral)
  
  post <- as_draws_df(m)
  
  # colunas de efeito genotípico principal
  cols_G <- grep("^r_genotipo\\[", names(post), value = TRUE)
  cols_G <- cols_G[!grepl("ano", cols_G)]
  
  # colunas de interação genótipo:ano
  cols_GxA <- grep("^r_genotipo:ano\\[", names(post), value = TRUE)
  
  # nomes dos genótipos
  gen_names <- gsub("r_genotipo\\[(.+),Intercept\\]", "\\1", cols_G)
  
  # efeito fixo de ano (intercepto = 2023, b_ano2024 = desvio de 2024)
  b_intercept <- post[["b_Intercept"]]
  b_ano2024   <- if ("b_ano2024" %in% names(post)) post[["b_ano2024"]] else 0
  
  n_samples <- nrow(post)
  
  # matrizes de valor genotípico por amostra: linhas = amostras, colunas = genótipos
  # escala padronizada → back-transform para escala original
  mat_G <- as.matrix(post[, cols_G]) * dp
  
  # separar GxA por ano
  cols_2023 <- grep("_2023", cols_GxA, value = TRUE)
  cols_2024 <- grep("_2024", cols_GxA, value = TRUE)
  
  # ordenar para garantir mesma ordem dos genótipos
  gen_order_2023 <- gsub("r_genotipo:ano\\[(.+)_2023,Intercept\\]",
                         "\\1", cols_2023)
  gen_order_2024 <- gsub("r_genotipo:ano\\[(.+)_2024,Intercept\\]",
                         "\\1", cols_2024)
  
  mat_GxA_2023 <- as.matrix(post[, cols_2023]) * dp
  mat_GxA_2024 <- as.matrix(post[, cols_2024]) * dp
  
  # reordenar colunas para garantir consistência com gen_names
  idx_2023 <- match(gen_names, gen_order_2023)
  idx_2024 <- match(gen_names, gen_order_2024)
  
  mat_GxA_2023 <- mat_GxA_2023[, idx_2023]
  mat_GxA_2024 <- mat_GxA_2024[, idx_2024]
  
  # valor genotípico total por ano e por amostra
  # valor_ij = mu + BLUP_G + BLUP_GxA_j (na escala original)
  mat_2023 <- mu + mat_G + mat_GxA_2023
  mat_2024 <- mu + mat_G + mat_GxA_2024
  
  # calcular probabilidades para cada intensidade
  map_dfr(intensidades, function(intens) {
    
    n_top <- ceiling(n_gen * intens)
    
    # para cada amostra, quais genótipos estão no top?
    # considera direção: para variáveis "baixo", menor = melhor
    direcao <- var_direcao[var]
    
    sup_2023 <- apply(mat_2023, 1, function(vals) {
      if (direcao == "baixo") {
        rank(vals) <= n_top
      } else {
        rank(-vals) <= n_top
      }
    })  # resultado: matriz n_gen × n_samples
    
    sup_2024 <- apply(mat_2024, 1, function(vals) {
      if (direcao == "baixo") {
        rank(vals) <= n_top
      } else {
        rank(-vals) <= n_top
      }
    })
    
    # probabilidade de ser superior nos DOIS anos simultaneamente
    # sup_2023 e sup_2024 são matrizes n_gen × n_samples
    prob_consistente <- rowMeans(sup_2023 & sup_2024)
    
    # probabilidade de ser superior em pelo menos um ano
    prob_qualquer    <- rowMeans(sup_2023 | sup_2024)
    
    # probabilidade de ser superior em 2023
    prob_2023        <- rowMeans(sup_2023)
    
    # probabilidade de ser superior em 2024
    prob_2024        <- rowMeans(sup_2024)
    
    tibble(
      variavel         = var,
      genotipo         = gen_names,
      intensidade      = intens,
      prob_consistente = prob_consistente,
      prob_qualquer    = prob_qualquer,
      prob_2023        = prob_2023,
      prob_2024        = prob_2024
    )
  })
}

# ── loop por variável ──────────────────────────────────────────────────────────
df_prob <- map_dfr(vars_resp, function(v) {
  calcular_prob_consistente(
    var         = v,
    intensidades = intensidades,
    modelos     = modelos,
    media_orig  = media_geral
  )
})

# adicionar grupo e métricas de wricke
df_prob <- df_prob |>
  mutate(genotipo = gsub("\\.", " ", genotipo)) |>
  left_join(
    df_wricke |>
      dplyr::select(variavel, genotipo, grupo, wricke, wricke_rel,
             rank_estab, rank_desempenho, inversao),
    by = c("variavel", "genotipo")
  )

# ── verificação ────────────────────────────────────────────────────────────────
cat("\n── Verificação — prob_consistente 20% — mcm_saca ──\n")
df_prob |>
  filter(variavel == "mcm_saca", intensidade == 0.20) |>
  arrange(desc(prob_consistente)) |>
  dplyr::select(genotipo, grupo, prob_consistente,
         prob_2023, prob_2024, wricke, rank_desempenho) |>
  head(15) |>
  print()

cat("\n── Quantos genótipos com prob_consistente > 0.50 (20%) ──\n")
df_prob |>
  filter(intensidade == 0.20) |>
  group_by(variavel) |>
  summarise(
    n_prob50  = sum(prob_consistente > 0.50),
    n_prob25  = sum(prob_consistente > 0.25),
    max_prob  = max(prob_consistente),
    .groups   = "drop"
  ) |>
  print()

cat("\nBloco 3 concluído — df_prob pronto\n")



# 04_stability.R — Bloco 4
# Visualizações — figura central do artigo
#
# 4a. Gráfico desempenho × estabilidade (wricke) com prob_consistente
# 4b. Lollipop de prob_consistente por variável e intensidade
# 4c. Heatmap de prob_consistente × ano por genótipo (estilo Chagas)

# ── paleta e temas ─────────────────────────────────────────────────────────────
cores_grupo <- c("Conilon" = "#185FA5", "Robusta" = "#D85A30")

tema_base <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    strip.text         = element_text(size = 9, face = "bold")
  )

# ── correção 04b ───────────────────────────────────────────────────────────────
df_lollipop <- df_prob |>
  filter(variavel == "mcm_saca", genotipo %in% top_gen_mcm) |>
  select(-grupo) |>                          # remover grupo duplicado antes do join
  left_join(
    df_wricke |>
      filter(variavel == "mcm_saca") |>
      select(genotipo, grupo) |>
      distinct(),
    by = "genotipo"
  ) |>
  mutate(
    genotipo    = factor(genotipo, levels = rev(top_gen_mcm)),
    intensidade = factor(paste0(intensidade * 100, "%"),
                         levels = c("10%", "20%", "25%"))
  )

# confirmar que grupo está presente e sem NA
cat("grupo NA:", sum(is.na(df_lollipop$grupo)), "\n")
cat("grupos:", unique(df_lollipop$grupo), "\n")

p_lollipop <- ggplot(df_lollipop,
                     aes(y = genotipo, x = prob_consistente,
                         colour = intensidade)) +
  geom_segment(aes(xend = 0, yend = genotipo),
               linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             colour = "gray40", linewidth = 0.4) +
  scale_colour_manual(
    values = c("10%" = "#185FA5", "20%" = "#1D9E75", "25%" = "#D85A30")
  ) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::percent_format()) +
  facet_wrap(~grupo, scales = "free_y", ncol = 2) +
  labs(
    x        = "probabilidade de superioridade consistente",
    y        = NULL,
    colour   = "intensidade",
    title    = "Probabilidade de superioridade consistente — MCM/saca",
    subtitle = "Top 15 genótipos | linha tracejada = 50%"
  ) +
  tema_base

ggsave(paste0(path_fig, "04b_lollipop_prob_mcm_saca.png"),
       p_lollipop, width = 10, height = 7, dpi = 300)
cat("04b salvo\n")

# ── correção 04a — mesmo problema de grupo duplicado ──────────────────────────
df_plot_ds <- df_wricke |>
  left_join(
    df_prob |>
      filter(intensidade == 0.20) |>
      select(variavel, genotipo, prob_consistente),  # sem grupo
    by = c("variavel", "genotipo")
  ) |>
  mutate(variavel_lab = factor(var_labels[variavel], levels = var_labels))

ref_ds <- df_plot_ds |>
  group_by(variavel) |>
  summarise(
    blup_sinal_med = mean(blup_sinal),
    wricke_med     = mean(wricke),
    .groups        = "drop"
  )

df_plot_ds <- df_plot_ds |>
  left_join(ref_ds, by = "variavel")

cat("grupo NA em df_plot_ds:", sum(is.na(df_plot_ds$grupo)), "\n")

p_desempenho_estab <- ggplot(df_plot_ds,
                             aes(x = blup_sinal, y = wricke)) +
  geom_hline(aes(yintercept = wricke_med),
             linetype = "dashed", colour = "gray60", linewidth = 0.4) +
  geom_vline(aes(xintercept = blup_sinal_med),
             linetype = "dashed", colour = "gray60", linewidth = 0.4) +
  geom_point(aes(colour = grupo, size = prob_consistente,
                 alpha = prob_consistente)) +
  ggrepel::geom_text_repel(
    data = df_plot_ds |> filter(prob_consistente > 0.15),
    aes(label = genotipo, colour = grupo),
    size = 2.6, max.overlaps = 20, show.legend = FALSE
  ) +
  scale_colour_manual(values = cores_grupo) +
  scale_size_continuous(range = c(1, 5), name = "prob. consistente") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  facet_wrap(~variavel_lab, scales = "free", ncol = 3) +
  labs(
    x        = "BLUP — desempenho médio (sinal ajustado: positivo = melhor)",
    y        = "Ecovalência de Wricke (menor = mais estável)",
    colour   = "grupo",
    title    = "Desempenho × Estabilidade temporal",
    subtitle = "Intensidade de seleção: 20% | tamanho = probabilidade de superioridade consistente"
  ) +
  tema_base +
  theme(axis.text.x = element_text(size = 8))

ggsave(paste0(path_fig, "04a_desempenho_estabilidade.png"),
       p_desempenho_estab, width = 13, height = 9, dpi = 300)
cat("04a salvo\n")

# ── 4c. Heatmap prob dentro de cada ano — todas as variáveis ──────────────────
# estilo Chagas et al. (2025, 2026)
# eixo x = ano, eixo y = genótipo, cor = prob_2023 / prob_2024

df_heat <- df_prob |>
  filter(intensidade == 0.20) |>
  select(variavel, genotipo, grupo, prob_2023, prob_2024,
         prob_consistente) |>
  pivot_longer(
    cols      = c(prob_2023, prob_2024),
    names_to  = "ano",
    values_to = "prob"
  ) |>
  mutate(
    ano          = gsub("prob_", "", ano),
    variavel_lab = factor(var_labels[variavel], levels = var_labels)
  ) |>
  # ordenar genótipos por prob_consistente decrescente dentro de cada variável
  group_by(variavel) |>
  mutate(
    genotipo = factor(genotipo,
                      levels = df_prob |>
                        filter(variavel == first(variavel),
                               intensidade == 0.20) |>
                        arrange(prob_consistente) |>
                        pull(genotipo) |>
                        unique())
  ) |>
  ungroup()

p_heat <- ggplot(df_heat,
                 aes(x = ano, y = genotipo, fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours  = c("#042C53", "#185FA5", "#1D9E75", "#D85A30", "#993C1D"),
    values   = scales::rescale(c(0, 0.25, 0.50, 0.75, 1)),
    limits   = c(0, 1),
    name     = "Pr(superior)",
    labels   = scales::percent_format()
  ) +
  facet_wrap(~variavel_lab, nrow = 1) +
  labs(
    x        = NULL,
    y        = NULL,
    title    = "Probabilidade de superioridade por ano",
    subtitle = "Intensidade de seleção: 20%"
  ) +
  tema_base +
  theme(
    axis.text.y     = element_text(size = 6.5),
    axis.text.x     = element_text(size = 9),
    panel.spacing   = unit(0.5, "lines"),
    legend.key.height = unit(1.2, "cm")
  )

ggsave(paste0(path_fig, "04c_heatmap_prob_por_ano.png"),
       p_heat, width = 14, height = 10, dpi = 300)

cat("\nBloco 4 concluído — figuras salvas em outputs/figures/\n")
cat("  04a_desempenho_estabilidade.png\n")
cat("  04b_lollipop_prob_mcm_saca.png\n")
cat("  04c_heatmap_prob_por_ano.png\n")


# 04_stability.R — Bloco 5
# Consolidação, tabelas finais e salvamento
#
# Saídas:
#   outputs/tables/04_stability_results.xlsx
#     aba 1: wricke_ecovalencia     — ecovalência por genótipo × variável
#     aba 2: prob_consistente       — probabilidades por genótipo × variável × intensidade
#     aba 3: ranking_final          — ranking integrado desempenho + estabilidade
#     aba 4: resumo_por_variavel    — estatísticas resumo por variável
#     aba 5: top_genotipos          — genótipos recomendados por variável

# ── ranking integrado ──────────────────────────────────────────────────────────
# combina rank_desempenho + rank_estabilidade em um índice composto
# menor = melhor em ambas as dimensões

df_ranking <- df_wricke |>
  left_join(
    df_prob |>
      filter(intensidade == 0.20) |>
      select(variavel, genotipo, prob_consistente,
             prob_2023, prob_2024),
    by = c("variavel", "genotipo")
  ) |>
  group_by(variavel) |>
  mutate(
    # índice composto: média dos ranks de desempenho e estabilidade
    # rank_desempenho: 1 = melhor desempenho
    # rank_estab: 1 = mais estável
    rank_composto   = (rank_desempenho + rank_estab) / 2,
    rank_prob       = rank(-prob_consistente),
    # classificação por quadrante
    desempenho_cat  = if_else(rank_desempenho <= 24,
                              "superior", "inferior"),
    estabilidade_cat = if_else(rank_estab <= 24,
                               "estável", "instável"),
    quadrante       = paste(desempenho_cat, estabilidade_cat, sep = " + ")
  ) |>
  ungroup() |>
  arrange(variavel, rank_composto)

cat("\n── Distribuição por quadrante — mcm_saca ──\n")
df_ranking |>
  filter(variavel == "mcm_saca") |>
  count(quadrante) |>
  print()

cat("\n── Top 10 ranking composto — mcm_saca ──\n")
df_ranking |>
  filter(variavel == "mcm_saca") |>
  select(genotipo, grupo, rank_desempenho, rank_estab,
         rank_composto, prob_consistente, quadrante) |>
  head(10) |>
  print()

# ── top genótipos recomendados ─────────────────────────────────────────────────
# critério: quadrante superior + estável E prob_consistente > 0.10
df_top <- df_ranking |>
  filter(
    quadrante == "superior + estável",
    prob_consistente > 0.10
  ) |>
  arrange(variavel, rank_composto) |>
  select(variavel, genotipo, grupo,
         blup_orig_med, blup_sinal,
         wricke, wricke_rel,
         rank_desempenho, rank_estab, rank_composto,
         prob_consistente, prob_2023, prob_2024,
         inversao, quadrante)

cat("\n── Genótipos recomendados (superior + estável, prob > 0.10) ──\n")
df_top |>
  select(variavel, genotipo, grupo, rank_composto,
         prob_consistente, wricke_rel) |>
  print(n = 50)

# ── resumo por variável ────────────────────────────────────────────────────────
df_resumo <- df_prob |>
  filter(intensidade == 0.20) |>
  group_by(variavel) |>
  summarise(
    n_total             = n_distinct(genotipo),
    n_prob_gt50         = sum(prob_consistente > 0.50),
    n_prob_gt25         = sum(prob_consistente > 0.25),
    n_prob_gt10         = sum(prob_consistente > 0.10),
    max_prob            = max(prob_consistente),
    gen_max_prob        = genotipo[which.max(prob_consistente)],
    n_inversao          = sum(inversao, na.rm = TRUE),
    prop_inversao       = mean(inversao, na.rm = TRUE),
    .groups             = "drop"
  ) |>
  left_join(
    # adicionar rGA do modelo bayesiano
    readxl::read_xlsx(
      "outputs/tables/03_parametros_completos.xlsx"
    ) |>
      select(variavel, rGA_med, rGA_q025, rGA_q975,
             H2_med, H2_q025, H2_q975,
             prop_GxA_maior_G),
    by = "variavel"
  )

cat("\n── Resumo por variável ──\n")
print(df_resumo)

# ── salvar tudo em xlsx ────────────────────────────────────────────────────────
write_xlsx(
  list(
    wricke_ecovalencia  = df_wricke |>
      select(variavel, genotipo, grupo,
             blup_gxa_2023, blup_gxa_2024,
             wricke, wricke_rel, amplitude_gxa,
             inversao, blup_orig_med, blup_sinal,
             rank_estab, rank_desempenho),
    prob_consistente    = df_prob |>
      select(variavel, genotipo, grupo,
             intensidade, prob_consistente,
             prob_2023, prob_2024, prob_qualquer,
             wricke, rank_estab, rank_desempenho),
    ranking_final       = df_ranking |>
      select(variavel, genotipo, grupo,
             blup_orig_med, blup_sinal,
             wricke, wricke_rel,
             rank_desempenho, rank_estab,
             rank_composto, rank_prob,
             prob_consistente, prob_2023, prob_2024,
             inversao, quadrante,
             desempenho_cat, estabilidade_cat),
    top_genotipos       = df_top,
    resumo_por_variavel = df_resumo
  ),
  "outputs/tables/04_stability_results.xlsx"
)

cat("\nBloco 5 concluído.\n")
cat("Resultados salvos em outputs/tables/04_stability_results.xlsx\n")
cat("\nResumo final do 04_stability.R:\n")
cat("  Genótipos avaliados:", n_distinct(df_ranking$genotipo), "\n")
cat("  Variáveis analisadas:", n_distinct(df_ranking$variavel), "\n")
cat("  Intensidades testadas: 10%, 20%, 25%\n")
cat("  Genótipos recomendados (superior+estável, prob>0.10):",
    n_distinct(df_top$genotipo), "\n")

