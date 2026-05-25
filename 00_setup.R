# 00_setup.R
# Configuração global do projeto
# Rodar uma vez para garantir que todos os pacotes estão instalados

pkgs <- c(
  # dados
  "readxl", "dplyr", "tidyr", "forcats", "stringr",
  # visualização
  "ggplot2", "patchwork", "ggdist", "scales", "ggrepel",
  # modelos mistos (frequentista)
  "lme4", "lmerTest", "emmeans", "pbkrtest",
  # bayesiano (reserva — instala mas não carrega ainda)
  "brms", "bayesplot", "tidybayes",
  # diagnóstico e parâmetros
  "performance", "parameters", "see",
  # multivariate / diversidade
  "factoextra", "cluster",
  # misc
  "skimr", "janitor", "purrr", "tibble", "broom.mixed"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

# opções globais
options(
  scipen = 999,
  digits = 4,
  OutDec = "."
)

message("Setup concluído. Pacotes disponíveis: ", length(pkgs) - length(to_install), "/", length(pkgs))


