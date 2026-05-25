# 07_meteorologia.R
# Climate data via NASA POWER for Jaguaré-ES
# Coordinates: lat = -18.99, lon = -40.08
# Period: 2022-2024
# Variables: T2M (mean temp), T2M_MAX, T2M_MIN, PRECTOTCORR (precip), RH2M (RH)
#
# Outputs:
#   outputs/figures/07_climograma.png
#   outputs/tables/07_resumo_meteorologico.xlsx
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

# ── packages ────────────────────────────────────────────────────────────────────
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

# declare preferences
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")

path_fig <- "outputs/figures/"
path_tbl <- "outputs/tables/"

# ── 1. download NASA POWER ────────────────────────────────────────────────────
cat("── Downloading NASA POWER data — Jaguaré-ES ──\n")

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

# extract parameters
params <- data$properties$parameter

# convert to long tibble
df_meteo <- map_dfr(names(params), function(p) {
  vals <- params[[p]]
  tibble(
    parameter = p,
    year_month   = names(vals),
    value        = as.numeric(unlist(vals))
  )
}) |>
  filter(grepl("^\\d{6}$", year_month)) |>
  mutate(
    year = as.integer(substr(year_month, 1, 4)),
    month = as.integer(substr(year_month, 5, 6))
  ) |>
  filter(between(month, 1, 12)) |>
  mutate(
    date = make_date(year, month, 1)
  ) |>
  filter(between(year, 2022, 2024)) |>
  filter(value != -999)

cat("Dimensions:", nrow(df_meteo), "×", ncol(df_meteo), "\n")
cat("Period:", min(df_meteo$date), "to", max(df_meteo$date), "\n")

# ── 2. data reshaping and precipitation correction (mm/day → mm/month) ──────────
df_wide <- df_meteo |>
  pivot_wider(names_from = parameter, values_from = value) |>
  rename(
    temp_media = T2M,
    temp_max   = T2M_MAX,
    temp_min   = T2M_MIN,
    precip     = PRECTOTCORR,
    ur_media   = RH2M
  ) |>
  mutate(
    amplitude = temp_max - temp_min,
    precip    = precip * days_in_month(date), # Correction to monthly total
    ano_f     = as.factor(year)
  )

# verify corrected totals
cat("── Corrected annual precipitation totals ──\n")
df_wide |>
  group_by(year) |>
  summarise(precip_annual = sum(precip), .groups = "drop") |>
  print()

# ── 3. summaries ──────────────────────────────────────────────────────────────
resumo_anual <- df_wide |>
  group_by(year) |>
  summarise(
    avg_temp      = mean(temp_media,  na.rm = TRUE),
    avg_max_temp  = mean(temp_max,    na.rm = TRUE),
    avg_min_temp  = mean(temp_min,    na.rm = TRUE),
    total_precip  = sum(precip,        na.rm = TRUE),
    avg_rh        = mean(ur_media,    na.rm = TRUE),
    avg_amplitude = mean(amplitude,   na.rm = TRUE),
    .groups = "drop"
  )

resumo_colheita <- df_wide |>
  filter(month %in% 4:7) |>
  group_by(year) |>
  summarise(
    period         = "April-July",
    avg_temp       = mean(temp_media,  na.rm = TRUE),
    avg_max_temp   = mean(temp_max,    na.rm = TRUE),
    avg_min_temp   = mean(temp_min,    na.rm = TRUE),
    total_precip   = sum(precip,        na.rm = TRUE),
    avg_rh         = mean(ur_media,    na.rm = TRUE),
    .groups = "drop"
  )

# ── 4. visualization ──────────────────────────────────────────────────────────
cores_ano <- c("2022" = "#999999", "2023" = "#56B4E9", "2024" = "#E69F00")

df_plot <- df_wide |>
  mutate(
    month_label = factor(month,
                        levels = 1:12,
                        labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                   "Jul","Aug","Sep","Oct","Nov","Dec"))
  )

p_temp <- ggplot(df_plot, aes(x = month_label, y = temp_media, colour = ano_f, group = ano_f)) +
  geom_line(linewidth = 0.8) + geom_point() +
  scale_colour_manual(values = cores_ano) +
  labs(y = "Temp (°C)", x = NULL, title = "Mean Temperature") + theme_minimal()

p_prec <- ggplot(df_plot, aes(x = month_label, y = precip, fill = ano_f, group = ano_f)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
  scale_fill_manual(values = cores_ano) +
  scale_y_continuous(name = "Precipitation (mm)") +
  annotate("rect", xmin = 3.5, xmax = 7.5, ymin = 0, ymax = Inf, fill = "gray80", alpha = 0.25) +
  labs(x = NULL, title = "Monthly Precipitation") + theme_minimal()

p_ur <- ggplot(df_plot, aes(x = month_label, y = ur_media, colour = ano_f, group = ano_f)) +
  geom_line(linewidth = 0.8) + geom_point() +
  scale_colour_manual(values = cores_ano) +
  labs(y = "RH (%)", x = NULL, title = "Relative Humidity") + theme_minimal()

p_clima <- p_temp / p_prec / p_ur +
  plot_annotation(
    title    = "Climatic conditions in Jaguaré-ES (2022-2024)",
    subtitle = "NASA POWER | shaded area = harvest period (Apr-Jul)",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

ggsave(paste0(path_fig, "07_climograma.png"), p_clima, width = 10, height = 11, dpi = 300)

# ── 5. save files ─────────────────────────────────────────────────────────────
write_xlsx(
  list(
    monthly_data  = df_wide,
    annual_summary = resumo_anual,
    harvest_summary = resumo_colheita
  ),
  paste0(path_tbl, "07_resumo_meteorologico.xlsx")
)
cat("Data saved successfully.\n")
