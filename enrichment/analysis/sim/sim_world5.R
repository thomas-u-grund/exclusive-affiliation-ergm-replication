# ============================================================
# Simulation study, scenario 5 (2026-11-13): agency heterogeneity.
# Same world as sim_world.R; adds two agency-level covariates that mirror the real-data mechanisms that hid
# stacking in the homogeneous model:
#   cov_ap_spec.csv  agency-player: player's nationality is the agency's focus nationality (one per agency, random)
#   cov_ac_mkt.csv   agency-club:   club is in the agency's market country (one per agency with < 30 clients;
#                    agencies with 30+ clients have no single market, as in the data)
# Run from the project root: Rscript enrichment/analysis/sim/sim_world5.R
# ============================================================
suppressMessages(library(tidyverse))
set.seed(20261113)
w <- "enrichment/analysis/sim/world"
nodes <- read_csv(file.path(w, "nodes.csv"), show_col_types = FALSE)
edges <- read_csv(file.path(w, "edges.csv"), show_col_types = FALSE)
NP <- sum(nodes$kind == 1); ag <- nodes$node[nodes$kind == 2]; cl <- nodes$node[nodes$kind == 3]
size <- tabulate(edges$x, max(nodes$node))[ag]
focus <- sample(1:10, length(ag), replace = TRUE, prob = c(rep(0.14, 5), rep(0.06, 5)))
mkt <- ifelse(size < 30, sample(1:5, length(ag), replace = TRUE), 0L)
nat <- nodes$nat[seq_len(NP)]
ccty <- nodes$club_country[cl]

ap <- expand.grid(a = seq_along(ag), p = seq_len(NP)) |> filter(nat[p] == focus[a])
write_csv(tibble(i = ag[ap$a], j = ap$p, v = 1), file.path(w, "cov_ap_spec.csv"))
ac <- expand.grid(a = seq_along(ag), c = seq_along(cl)) |> filter(mkt[a] == ccty[c])
write_csv(tibble(i = ag[ac$a], j = cl[ac$c], v = 1), file.path(w, "cov_ac_mkt.csv"))
cat(sprintf("agencies %d, with a market %d (%.0f%% of players' agencies); ap_spec nonzero %d, ac_mkt nonzero %d\n",
            length(ag), sum(mkt > 0), 100 * sum(size[mkt > 0]) / NP, nrow(ap), nrow(ac)))
