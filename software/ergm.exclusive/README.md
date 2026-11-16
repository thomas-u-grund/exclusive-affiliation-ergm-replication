# ergm.exclusive

Exponential random graph models for **exclusive multi-affiliation networks**: every actor (player) holds exactly one
affiliation of each type (agency, club), and the size of every affiliation is fixed. Development version.

Status (development version; the sampler is verified against exact enumeration, see below)

* **`ergm()` integration**: `constraints = ~exclusive(types)` selects the typed-swap proposal (C, `src/exclusive_swap.c`), so
  that `simulate()` and `ergm()` work on a bipartite network (actors vs. affiliations) with any existing bipartite term, e.g.
  stacking = `b1dsp(2)`. Verified: on the 8-actor toy network the simulated distribution of `b1dsp(2)` matches exact
  enumeration at theta = 0 and theta = 0.8 (60 chains x 10,000 draws, all |z| < 2; `inst/examples/ergm_hook_check.R`).
* `ex_simulate()`: reference typed-swap Metropolis–Hastings sampler (pure R; small networks).
* `ex_stats()`: stacking, nationality homophily within agencies and within clubs, co-location by club country.
* `ex_enumerate()`, `ex_exact()`: exact distribution of any statistic on small networks (560 × 560 arrangements for the toy).
* `ex_reachable()`: breadth-first check that the swaps reach every arrangement.
* `inst/examples/toy_benchmark.R`: sampler against exact enumeration at theta = 0 and theta != 0.
* **Cross-layer terms** (C, `src/terms.c`, tested against direct computation on random arrangements and against exact
  enumeration at nonzero theta, `tests/testthat/test-terms.R`):
  * `coloc(centre, key)`: for every centre actor (e.g. agency) pairs of its players whose other affiliation has the same key
    (club country: the co-location term L; every club its own key with only small agencies as centres: size-dependent stacking);
  * `stack()`: pairs of players sharing their agency and their club (same value as `b1dsp(2)` on this space, about 4x faster in the football application because it scans the smaller partner set);
  * `homcentre(centre, attr)`: pairs of players of the same nationality (`attr`) within the groups of the centre actors (H_A, H_C);
  * `acov(W)`: sum over players of W[agency, club], an agency-club covariate (corridor, value gap, market presence).
* **Football example** (the paper's models (1), (3m), (3s)): `enrichment/analysis/ergm_excl_data.R` (loader),
  `ergm_excl_terms.R` (formulas), `ergm_excl_fit.R` and `ergm_excl_check2.R` (fit and moment check) in the project
  directory (the data are not distributed). On the 4,445-node network: simulation at the fitted estimates reproduces all
  observed statistics for model (3m); `ergm()` started from the fitted estimates reproduces the size model within standard errors.
  The sampler mixes slowly (about 10^6 swaps per independent draw): use `MCMC.interval` of 1e6 and expect hours per fit.

Install from source: `R CMD INSTALL .` (needs `ergm` >= 4.12 and a C compiler). Maintainer email in `DESCRIPTION` is a placeholder.

See the paper in this repository's parent project for the model.

## Fast sampler

`ex_fast_simulate()` is a dedicated compiled sampler for two-affiliation networks (one agency and one club per player). It keeps
count tables so that every supported statistic (`stack`, `coloc`, `homcentre`, `acov`, actor-by-affiliation covariates as in `edgecov`) is updated in O(1)
per proposed swap, and it mixes the plain typed swap with a swap between players that share a partner of the other type (both symmetric). On the
football network it runs 3 million swaps in about 1 s, against about 8 s through `ergm`'s generic simulation; the statistics agree with the `ergm`
terms to rounding error, and the sampler is checked against exact enumeration. The analysis scripts use it through `fast_terms()` / `fast_sim()`
(`enrichment/analysis/ergm_excl_fast.R`), which turn the same formula used for `ergm` into the term list.
