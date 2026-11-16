# Replication materials: exponential random graph models for exclusive multi-affiliation networks

Replication archive for "Exponential random graph models for exclusive multi-affiliation networks: agencies, players and clubs in European football" (Grund and Dittrich).

## Contents
- `software/ergm.exclusive/`: R package with the exclusive-affiliation terms, the typed and grouped swap proposals and the compiled sampler `ex_fast_simulate()`. Install with `R CMD INSTALL software/ergm.exclusive` (needs `ergm` >= 4.12, `network`, a C compiler).
- `enrichment/analysis/`: R analysis scripts; `sim/` holds the simulation-study world definitions (`sim_world*.R`), the R fits and permutation tests, the summaries and the result summaries (`sim/results/`).
- `enrichment/analysis/output/tp/`: the football network and covariates (compressed, see Data), starting values and result files (`rrefine/`, `fastval/`, `gof/`, `proj/`, `recover*/`).
- `paper/figures/`, `paper/tables/`: scripts for the figures and the generated tables of the paper.
- `enrichment/*.py`: the scraper and table builders for the raw data (see Data).

## Data
The estimation inputs are the derived network and covariates (the football snapshot of 24 August 2026: 3,363 players, 950 agencies, 132 clubs): `nodes`, `edges` and `cov_*` as `.csv.gz` in `enrichment/analysis/output/tp/`. Unpack with `gunzip -k enrichment/analysis/output/tp/*.csv.gz`. Players and agencies are identified by numeric ids only. The raw scraped pages and tables from the data source (Transfermarkt) are not redistributed; `scrape_tm.py`, `parse_agency_clients.py` and `build_tables.py` document how they were obtained, and the scripts that start from them (`ergm_tp_data.R`, `ergm_tp_mkt.R`, `ergm_tp_mkt_n.R`) are included for transparency but need those files. Users of the derived data should respect the source's terms of use.

## Order of scripts (run from the repository root)
1. Estimation (compiled sampler): `enrichment/analysis/ergm_excl_refine_fast.R` (env `MODEL`, `START`, `OUT`, `ITER`, `NS`, `IV`, `DAMP`, `CAP`, `SEED0`, `BURN`); model definitions in `ergm_excl_terms.R` (models `1`, `3m`, `3ms2`, `3s3`, and `3s3c` with env `CUT1`/`CUT2` for the cut-off check). Starting values and refined estimates: `output/tp/rrefine/`. In the paper, Models 1 to 3 correspond to tags `3m`, `3ms2`, `3s3`, and the core model of the appendix to tag `1`.
2. Diagnostics and fit by agency size: `fast_validate.R` (four chains of 200 draws) -> `output/tp/fastval/`; tables via `ergm_excl_diag_fast_tex.R`.
3. Projection analyses of the football network: `football_projection.R` -> `output/tp/proj/`.
4. Recovery at the scale of the application: `ergm_excl_recover_fast.R` (env `REP` = 1..24) -> `output/tp/recover_fast*/`; table via `ergm_excl_recover_tex.R`.
5. Simulation study: summaries in `enrichment/analysis/sim/results/` are processed by `sim_summarise.R`; the permutation tests and R fits are `sim_perm_country.R` and `sim_fit_r.R` (generating worlds: `sim_world*.R`). The joint-model fits of the simulation study were produced with an earlier Stata implementation of the same model, which is not distributed here; their result files are included.
6. Figures: `Rscript paper/figures/make_figures.R`, `make_fig_patterns.R`.

Seeds are fixed in each script. All scripts run from the repository root.

## Licenses
Scripts: MIT (`LICENSE`). R package in `software/`: GPL-3 (see `software/ergm.exclusive`).
