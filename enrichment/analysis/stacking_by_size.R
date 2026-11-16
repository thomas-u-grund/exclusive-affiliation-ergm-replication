# ============================================================
# Diagnostic for the mismatch between the joint model (stacking theta_S ~ 0) and the real-data stacking
# permutation test (Q = 1.22). Hypothesis: stacking is heterogeneous by agency size; S (pairs) is dominated by
# large agencies, the "stacked players" statistic by small ones.
# Null: agency labels permuted among players within club country x nationality x value tercile (as in
# ergm_followups.R 2b). Under the stage-1 joint model with theta_S = 0 this null is exact.
# Output: output/stacking_by_size.csv (observed, null mean, Q, p) per statistic and agency-size band.
# Run from the project root: Rscript enrichment/analysis/stacking_by_size.R
# ============================================================
suppressMessages(library(tidyverse))
tp <- "enrichment/analysis/output/tp"
nodes <- read_csv(file.path(tp, "nodes.csv"), show_col_types = FALSE)
edges <- read_csv(file.path(tp, "edges.csv"), show_col_types = FALSE)
P <- nodes |> filter(kind == 1)
NP <- nrow(P)
kind_of <- nodes$kind[edges$x]
ag <- integer(NP); cl <- integer(NP)
ag[edges$p[kind_of == 2]] <- edges$x[kind_of == 2]
cl[edges$p[kind_of == 3]] <- edges$x[kind_of == 3]
ag <- match(ag, sort(unique(ag))); cl <- match(cl, sort(unique(cl)))
NA_ <- max(ag); NC <- max(cl)
cty <- tapply(P$club_country, cl, first)[as.character(cl)]   # club country per player (constant within club)
v3 <- ntile(P$lmv, 3); v3 <- coalesce(as.character(v3), "na")
strata <- paste(P$club_country, P$nat, v3)
asize <- tabulate(ag, NA_)
band_of_ag <- cut(asize, c(0, 1, 3, 9, 29, Inf), labels = c("1", "2-3", "4-9", "10-29", "30+"))
cat("agencies per band:\n"); print(table(band_of_ag))
cat("players per band:\n"); print(tapply(asize, band_of_ag, sum))

stats <- function(a) {
  w <- matrix(tabulate((a - 1L) * NC + cl, nbins = NA_ * NC), NA_, NC, byrow = TRUE)
  c2 <- rowSums(w * (w - 1) / 2); st <- rowSums(w) - rowSums(w > 0)
  c(C2_all = sum(c2), stacked_all = sum(st),
    setNames(tapply(c2, band_of_ag, sum), paste0("C2_", levels(band_of_ag))),
    setNames(tapply(st, band_of_ag, sum), paste0("stacked_", levels(band_of_ag))))
}
obs <- stats(ag)
idx <- Filter(function(i) length(i) > 1, split(seq_len(NP), strata))
NPERM <- as.integer(Sys.getenv("NPERM", "2000"))
set.seed(20261113)
null <- replicate(NPERM, { a <- ag; for (i in idx) a[i] <- a[i][sample.int(length(i))]; stats(a) })
res <- tibble(stat = names(obs), observed = obs, null_mean = rowMeans(null),
              Q = obs / rowMeans(null), p_upper = (rowSums(null >= obs) + 1) / (NPERM + 1),
              p_lower = (rowSums(null <= obs) + 1) / (NPERM + 1))
print(res, n = Inf)
write_csv(res, "enrichment/analysis/output/stacking_by_size.csv")
