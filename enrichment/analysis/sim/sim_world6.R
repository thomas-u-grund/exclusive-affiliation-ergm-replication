# Scenario 6 world: copy of world/ (scenario 5 covariates) plus node attribute key_clubs for a size-dependent stacking term:
#   tpcoloc(attr=key_clubs) counts pairs of clients of the same club, over agencies coded -1 (small: < 30 clients);
#   large agencies are coded -2 (not centres), clubs have unique keys 1..NC, players 0.
# Run from the project root: Rscript enrichment/analysis/sim/sim_world6.R
suppressMessages(library(tidyverse))
dir.create("enrichment/analysis/sim/world6", showWarnings = FALSE)
for (f in list.files("enrichment/analysis/sim/world", full.names = TRUE)) file.copy(f, "enrichment/analysis/sim/world6", overwrite = TRUE)
nodes <- read_csv("enrichment/analysis/sim/world6/nodes.csv", show_col_types = FALSE)
edges <- read_csv("enrichment/analysis/sim/world6/edges.csv", show_col_types = FALSE)
size <- tabulate(edges$x, nrow(nodes))
nodes$key_clubs <- case_when(nodes$kind == 1 ~ 0L, nodes$kind == 2 ~ ifelse(size[nodes$node] < 30, -1L, -2L),
                             TRUE ~ match(nodes$node, nodes$node[nodes$kind == 3]))
write_csv(nodes, "enrichment/analysis/sim/world6/nodes.csv", na = "")
cat("small agencies", sum(nodes$key_clubs == -1), "large", sum(nodes$key_clubs == -2), "\n")
