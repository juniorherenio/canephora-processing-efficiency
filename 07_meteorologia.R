# 07_meteorologia.R
# Dados climáticos via NASA POWER para Jaguaré-ES
# Coordenadas: lat = -18.99, lon = -40.08
# Período: 2022-2024
# Variáveis: T2M (temp média), T2M_MAX, T2M_MIN, PRECTOTCORR (precip), RH2M (UR)
#
# Saídas:
#   outputs/figures/07_climograma.png
#   outputs/tables/07_resumo_meteorologico.xlsx

# ── pacotes ────────────────────────────────────────────────────────────────────
library(conflicted)

library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(writexl)
library(purrr)

# declarar preferências
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# ── 1. download NASA POWER ────────────────────────────────────────────────────
cat("── Baixando dados NASA POWER — Jaguaré-ES ──\n")

url <- paste0(
  "https://power.larc.nasa.gov/api/temporal/monthly/point?",
  "parameters=T2M,T2M_MAX,T2M_MIN,PRECTOTCORR,RH2M&",
  "community=AG&",
  "longitude=-40.08&latitude=-18.99&",
  "start=2022&end=2024&",
  "format=JSON"
)

resp <- GET(url, timeout(60))
cat("Status:", status_code(resp), "\n")

raw  <- content(resp, as = "text", encoding = "UTF-8")
data <- fromJSON(raw)

# extrair parâmetros
params <- data$properties$parameter

# converter para tibble longo
df_meteo <- map_dfr(names(params), function(p) {
  vals <- params[[p]]
  tibble(
    parametro = p,
    ano_mes   = names(vals),
    valor     = as.numeric(unlist(vals))
  )
}) |>
  filter(grepl("^\\d{6}$", ano_mes)) |>
  mutate(
    ano = as.integer(substr(ano_mes, 1, 4)),
    mes = as.integer(substr(ano_mes, 5, 6))
  ) |>
  filter(between(mes, 1, 12)) |>   # 🔥 remove mês 13
  mutate(
    data = make_date(ano, mes, 1)
  ) |>
  filter(between(ano, 2022, 2024)) |>
  filter(valor != -999)

cat("Dimensões:", nrow(df_meteo), "×", ncol(df_meteo), "\n")
cat("Período:", min(df_meteo$data), "a", max(df_meteo$data), "\n")

# ── correção da precipitação: mm/dia → mm/mês ─────────────────────────────────
df_wide <- df_wide |>
  mutate(
    precip = precip * days_in_month(data)
  )

# verificar totais corrigidos
cat("── Totais anuais de precipitação corrigidos ──\n")
df_wide |>
  group_by(ano) |>
  summarise(precip_anual = sum(precip), .groups = "drop") |>
  print()

# ── recriar resumos com precipitação corrigida ────────────────────────────────
resumo_anual <- df_wide |>
  group_by(ano) |>
  summarise(
    temp_media_anual = mean(temp_media,  na.rm = TRUE),
    temp_max_media   = mean(temp_max,    na.rm = TRUE),
    temp_min_media   = mean(temp_min,    na.rm = TRUE),
    precip_total     = sum(precip,       na.rm = TRUE),
    ur_media_anual   = mean(ur_media,    na.rm = TRUE),
    amplitude_media  = mean(amplitude,   na.rm = TRUE),
    .groups = "drop"
  )

resumo_colheita <- df_wide |>
  filter(mes %in% 4:7) |>
  group_by(ano) |>
  summarise(
    periodo        = "abril-julho",
    temp_media     = mean(temp_media,  na.rm = TRUE),
    temp_max_media = mean(temp_max,    na.rm = TRUE),
    temp_min_media = mean(temp_min,    na.rm = TRUE),
    precip_total   = sum(precip,       na.rm = TRUE),
    ur_media       = mean(ur_media,    na.rm = TRUE),
    .groups = "drop"
  )

cat("\n── Resumo anual corrigido ──\n")
print(resumo_anual)

cat("\n── Período de colheita corrigido ──\n")
print(resumo_colheita)

# ── recriar df_plot com precipitação corrigida ────────────────────────────────
df_plot <- df_wide |>
  filter(ano %in% 2022:2024) |>
  mutate(
    mes_label = factor(mes,
                       levels = 1:12,
                       labels = c("Jan","Fev","Mar","Abr","Mai","Jun",
                                  "Jul","Ago","Set","Out","Nov","Dez"))
  )

# ── recriar p_prec com escala correta ─────────────────────────────────────────
p_prec <- ggplot(df_plot,
                 aes(x = mes_label, y = precip,
                     fill = ano_f, group = ano_f)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, alpha = 0.85) +
  scale_fill_manual(values = cores_ano, name = "ano") +
  scale_y_continuous(
    name   = "Precipitação (mm)",
    limits = c(0, max(df_plot$precip, na.rm = TRUE) * 1.15),
    breaks = scales::pretty_breaks(6)
  ) +
  annotate("rect",
           xmin = 3.5, xmax = 7.5,
           ymin = 0,
           ymax = max(df_plot$precip, na.rm = TRUE) * 1.15,
           fill = "gray80", alpha = 0.25) +
  labs(x = NULL, title = "Precipitação mensal acumulada") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    axis.text.x        = element_text(angle = 0)
  )

# ── remontar e salvar climograma ───────────────────────────────────────────────
p_clima <- p_temp / p_prec / p_ur +
  plot_annotation(
    title    = "Condições climáticas em Jaguaré-ES (2022-2024)",
    subtitle = paste0(
      "NASA POWER — lat: -18.99°, lon: -40.08° | ",
      "retângulo cinza = período de colheita (abr-jul)"
    ),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  colour = "gray40")
    )
  )

ggsave(paste0(path_fig, "07_climograma.png"),
       p_clima, width = 10, height = 11, dpi = 300)
cat("Climograma corrigido salvo.\n")

# ── atualizar xlsx ─────────────────────────────────────────────────────────────
write_xlsx(
  list(
    dados_mensais   = df_wide |>
      select(ano, mes, mes_label, temp_media, temp_max,
             temp_min, amplitude, precip, ur_media),
    resumo_anual    = resumo_anual,
    resumo_colheita = resumo_colheita
  ),
  paste0(path_tbl, "07_resumo_meteorologico.xlsx")
)
cat("Excel atualizado.\n")
