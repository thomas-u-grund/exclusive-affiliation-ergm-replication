# ============================================================
# Simulation study, step 1: the synthetic world (fixed across replications)
#
# 1,000 players, 40 clubs of 25 players in 5 countries, agencies with portfolio sizes resampled from the real
# network (so the skew is realistic), 10 nationalities (1-5 = nationals of countries 1-5, 6-10 = foreign).
# Each agency has a home country (or none). Writes, in the same format as output/tp/ (ergm_tp_data.R):
#   sim/world/nodes.csv, edges.csv (a random starting arrangement), cov_pc_home.csv, cov_ac_home.csv
# Run from the project root: Rscript enrichment/analysis/sim/sim_world.R
# ============================================================
suppressMessages(library(tidyverse))
set.seed(20261109)
out <- "enrichment/analysis/sim/world"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

NP <- 1000; NC <- 40; NCTY <- 5; NNAT <- 10
# portfolio sizes: resample real agency sizes until they add up to NP
real <- read_csv("enrichment/analysis/output/tp/edges.csv", show_col_types = FALSE)
nodes_real <- read_csv("enrichment/analysis/output/tp/nodes.csv", show_col_types = FALSE)
real_sizes <- real |> filter(x %in% nodes_real$node[nodes_real$kind == 2]) |> count(x) |> pull(n)
sizes <- integer(0)
while (sum(sizes) < NP) sizes <- c(sizes, sample(real_sizes, 1))
sizes[length(sizes)] <- sizes[length(sizes)] - (sum(sizes) - NP)
sizes <- sizes[sizes > 0]
NA_ <- length(sizes)
N <- NP + NA_ + NC
cat(sprintf("players %d, agencies %d (sizes: mean %.2f, max %d, share of size 1: %.2f), clubs %d, nodes %d\n",
            NP, NA_, mean(sizes), max(sizes), mean(sizes == 1), NC, N))

club_cty <- rep(seq_len(NCTY), each = NC / NCTY)
# nationalities: 70% nationals of one of the five countries (equal shares), 30% foreign (5 groups)
nat <- sample(seq_len(NNAT), NP, replace = TRUE, prob = c(rep(0.14, 5), rep(0.06, 5)))
# agency home country: 75% have one of the five, 25% none (coded 0)
ag_home <- sample(0:NCTY, NA_, replace = TRUE, prob = c(0.25, rep(0.15, 5)))

nodes <- tibble(
  node = seq_len(N),
  kind = c(rep(1L, NP), rep(2L, NA_), rep(3L, NC)),
  nat = c(nat, rep(NA, NA_ + NC)),
  club_country = c(rep(NA, NP + NA_), club_cty),
  ag_home = c(rep(NA, NP), ag_home, rep(NA, NC)),
  type = c(rep(NA, NP), rep(1L, NA_), rep(2L, NC)),
  key_cty = c(rep(0L, NP), rep(-1L, NA_), club_cty),
  cls = c(rep(0L, NP), rep(1L, NA_ + NC)),
  isag = as.integer(kind == 2L), iscl = as.integer(kind == 3L),
  bmax = ifelse(kind == 1L, 1L, NP), bmin = ifelse(kind == 1L, 1L, 0L)
)
write_csv(nodes, file.path(out, "nodes.csv"), na = "")

# random starting arrangement with the given sizes (the sampler's burn-in removes it)
a_of <- sample(rep(NP + seq_len(NA_), sizes))
c_of <- sample(rep(NP + NA_ + seq_len(NC), each = NP / NC))
write_csv(tibble(p = c(seq_len(NP), seq_len(NP)), x = c(a_of, c_of)), file.path(out, "edges.csv"))

# P-C covariate: player is a national of the club's country
pc <- expand.grid(p = seq_len(NP), c = seq_len(NC)) |>
  mutate(v = as.numeric(nat[p] == club_cty[c])) |> filter(v != 0)
write_csv(tibble(i = pc$p, j = NP + NA_ + pc$c, v = pc$v), file.path(out, "cov_pc_home.csv"))
# A-C covariate W: agency's home country is the club's country
ac <- expand.grid(a = seq_len(NA_), c = seq_len(NC)) |>
  mutate(v = as.numeric(ag_home[a] == club_cty[c])) |> filter(v != 0)
write_csv(tibble(i = NP + ac$a, j = NP + NA_ + ac$c, v = ac$v), file.path(out, "cov_ac_home.csv"))
cat(sprintf("pc_home nonzero %d, ac_home nonzero %d\n", nrow(pc), nrow(ac)))
