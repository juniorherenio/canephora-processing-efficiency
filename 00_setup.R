# 00_setup.R
# Global project configuration
# Run once to ensure all packages are installed
# Part of: Gonçalves Júnior et al. (2026), Biology (MDPI)
# Repository: https://github.com/juniorherenio/canephora-processing-efficiency

pkgs <- c(
  # data
  "readxl", "dplyr", "tidyr", "forcats", "stringr",
  # visualization
  "ggplot2", "patchwork", "ggdist", "scales", "ggrepel",
  # mixed models (frequentist)
  "lme4", "lmerTest", "emmeans", "pbkrtest",
  # bayesian (reserve — install but do not load yet)
  "brms", "bayesplot", "tidybayes",
  # diagnostics and parameters
  "performance", "parameters", "see",
  # multivariate / diversity
  "factoextra", "cluster",
  # misc
  "skimr", "janitor", "purrr", "tibble", "broom.mixed"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

# global options
options(
  scipen = 999,
  digits = 4,
  OutDec = "."
)

message("Setup completed. Packages available: ", length(pkgs) - length(to_install), "/", length(pkgs))
