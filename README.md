# canephora-processing-efficiency

R analysis scripts for the manuscript:

**Genetic variation in fruit-to-grain conversion efficiency in *Coffea canephora*:
heritability, temporal instability, and divergence in Robusta hybrids and Conilon**

Deurimar Herênio Gonçalves Júnior, Jéssica Almeida Jorge, Júlio César Pereira Machado,
Danillo Lima Pereira, Weverton Pereira Rodrigues, Fábio Luiz Partelli

*Biology* (MDPI), 2026. https://doi.org/XXXXXXXX

---

## Description

This repository contains the R scripts used for all statistical analyses reported
in the manuscript. The study estimated genetic parameters for five processing
efficiency traits in 48 *Coffea canephora* genotypes (40 Robusta, 8 Conilon)
evaluated over two harvest years (2023–2024) in Jaguaré, Espírito Santo, Brazil.

Raw data are not publicly available due to their origin in an active breeding
program. Requests for data access may be directed to the corresponding author
(partelli@yahoo.com.br).

---

## Script pipeline

| Script | Description |
|--------|-------------|
| `00_setup.R` | Package installation and environment setup |
| `01_load_explore.R` | Data loading, cleaning, and descriptive exploration |
| `02_model_selection_v2.R` | Mixed model structure selection (LRT, boundary diagnostics) |
| `03_bayes_models.R` | Bayesian model fitting via brms/Stan (MCMC) |
| `03_bayes_parameters.R` | Genetic parameter estimation from posterior distributions |
| `03b_reml_comparison.R` | Comparison between REML and Bayesian estimates |
| `03c_test_ano_block.R` | Test of year × block interaction (AIC-based) |
| `04_stability.R` | Wricke ecovalence and Bayesian probability of consistent superiority |
| `05_diversity.R` | Genetic divergence: Mahalanobis distance, UPGMA clustering, PCA |
| `06_multitrait_sommer.R` | Multi-trait model with unstructured covariance (AI-REML, sommer) |
| `07_meteorologia.R` | Climatic data retrieval via NASA POWER API |
| `08_figuras_publicacao_v2.R` | Publication-ready figures (600 dpi TIFF) |
| `09_tabela3_publicacao.R` | Formatted Table 3 for manuscript |

Scripts should be run sequentially (00 → 09). Each script reads inputs from
`data/` and writes outputs to `outputs/tables/` and `outputs/figures/`.

---

## R dependencies

```r
# Core statistical modeling
brms (>= 2.20.0)
lme4
sommer
posterior
loo

# Data manipulation
dplyr
tidyr
purrr
tibble
readxl
writexl

# Visualization
ggplot2
patchwork
ggrepel
ggdendro
factoextra

# Utilities
nasapower
janitor
skimr
scales
cluster
```

---

## Session info

Analysis conducted in R version 4.5.3 (2026-03-11).
Stan backend via cmdstanr or rstan (brms 2.23.0).

---

## Citation

If you use these scripts, please cite:

Gonçalves Júnior, D.H., Jorge, J.A., Machado, J.C.P., Pereira, D.L.,
Rodrigues, W.P., Partelli, F.L. (2026). Genetic variation in fruit-to-grain
conversion efficiency in *Coffea canephora*: heritability, temporal instability,
and divergence in Robusta hybrids and Conilon. *Biology*, XX, XXXX.
https://doi.org/XXXXXXXX

---

## License

Scripts are released under the [MIT License](LICENSE).