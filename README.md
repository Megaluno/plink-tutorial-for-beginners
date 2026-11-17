# PLINK for Livestock Genomics

A step-by-step tutorial on PLINK 1.9 (with PLINK 2.0 and 1.07 notes) using cattle and sheep data, by **Ismaila Emmanuel Oluwasegun**.

**Read it online:** https://megaluno.github.io/plink-tutorial-for-beginners/

## What is covered

Installing PLINK and running it from R/RStudio; file formats; converting, subsetting and merging data; quality control (missingness, allele frequency, Hardy-Weinberg, sex checks, Mendel errors, heterozygosity); linkage disequilibrium; relatedness, the genomic relationship matrix and runs of homozygosity; PCA (eigenvectors, eigenvalues, scree plots, projection and plotting in R); case/control and quantitative-trait association; family tests; permutation; dosages; haplotypes and imputation; CNVs; set-based tests, epistasis, G×E and meta-analysis; how each setting changes the results; a real-data chapter; and a complete QC pipeline.

## Repository layout

| Folder | Contents |
|---|---|
| `*.qmd`, `_quarto.yml`, `custom.scss` | Source of the Quarto book |
| `images/` | All figures, drawn from real PLINK output |
| `data/` | Simulated cattle data (with the answer key in `TRUTH.json`), `data/real/` (public cattle and sheep data) `data/phasing/` (simulated reference panels and chip data for phasing and imputation, with `TRUTH_phasing.json`) and `data/ascertainment/` (sequence and chip SNP sets for the ascertainment-bias example) |
| `qc_lists/` | QC lists used to build the cleaned file |
| `scripts/` | `simulate.py` (makes the simulated data), `chapters/` (the commands of each chapter), `make_figures.py`, `convert_real_data.py`; for Chapter 14 `simulate_phasing.py`, `ls_impute.py` (teaching Li and Stephens phasing and imputation), `eval_imputation.py`, `eval_window_phase.py`, `fig_imputation.py` and `phasing_imputation_pipeline.sh`; for Chapter 13 `lmm_gwas.py` (teaching mixed model); for Chapter 18 `simulate_ascertainment.py` and `ascertainment_pipeline.sh`; `fig_flowcharts.py`, `fig_imputation.py` and `fig_ascertainment.py` draw the diagrams and figures |
| `downloads/` | Everything in one zip |
| `docs/` | The rendered website (served by GitHub Pages) |

## Build the site yourself

Install [Quarto](https://quarto.org), then run `quarto render` in this folder. The website is written to `docs/`.

## Licence

Text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The real datasets keep their original licences: `synbreedData` (GPL-2) and `detectRUNS` (GPL-3).

## Citation

Emmanuel, I. (2026). *PLINK for Livestock Genomics: a step-by-step tutorial with cattle and sheep data.* https://megaluno.github.io/plink-tutorial-for-beginners/
