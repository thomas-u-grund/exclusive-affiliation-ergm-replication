# ============================================================
# Simulation study, step 4b: stacking permutation tests stratified by club country,
# as used for the real network (ergm_followups.R, section 2b).
# Players stay at their clubs; agency labels are reshuffled among players within strata:
#   cty      club country (keeps each agency's number of clients per country fixed)
#   cty_nat  club country x player nationality (the real-data test without the value tercile;
#            the simulated world has no market values)
# Two statistics: C2 = sum over agency-club cells of choose(w, 2) (as in sim_fit_r.R),
#                 stacked = players - distinct agency-club ties (as in the real-data test).
# Usage: Rscript enrichment/analysis/sim/sim_perm_country.R [scenarios, default 1 2 3 4]
# Writes sim/results/perm_country.csv (one row per scenario x replication; overwritten).
# ============================================================
suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
scens <- if (length(args)) as.integer(args) else 1:4
sim <- "enrichment/analysis/sim"
nodes <- read_csv(file.path(sim, "world/nodes.csv"), show_col_types = FALSE, guess_max = 1e6)
NP <- sum(nodes$kind == 1); NA_ <- sum(nodes$kind == 2); NC <- sum(nodes$kind == 3)
nat <- nodes$nat[nodes$kind == 1]
clubcty <- nodes$club_country[nodes$kind == 3]
NPERM <- 1000

stats <- function(ag, cl) {
  w <- tabulate((ag - 1L) * NC + cl, nbins = NA_ * NC)
  c(C2 = sum(w * (w - 1) / 2), stacked = NP - sum(w > 0))
}

rows <- list()
for (s in scens) {
  nets <- read_csv(file.path(sim, sprintf("nets/s%d_edges.csv", s)), show_col_types = FALSE) |>
    rename(sim = `_sim`, p = `_tail`, x = `_head`)
  for (r in sort(unique(nets$sim))) {
    e <- nets |> filter(sim == r)
    ag <- integer(NP); cl <- integer(NP)
    isag <- e$x <= NP + NA_
    ag[e$p[isag]] <- e$x[isag] - NP
    cl[e$p[!isag]] <- e$x[!isag] - NP - NA_
    cty <- clubcty[cl]
    obs <- stats(ag, cl)
    set.seed(1000 * s + r)
    out <- tibble(scen = s, rep = r)
    for (st in c("cty", "cty_nat")) {
      g <- if (st == "cty") cty else paste(cty, nat)
      idx <- Filter(function(i) length(i) > 1, split(seq_len(NP), g))
      null <- replicate(NPERM, {
        a <- ag
        for (i in idx) a[i] <- a[i][sample.int(length(i))]
        stats(a, cl)
      })
      for (k in names(obs)) {
        out[[paste0("Q_", k, "_", st)]] <- obs[[k]] / mean(null[k, ])
        out[[paste0("p_", k, "_", st)]] <- (sum(null[k, ] >= obs[[k]]) + 1) / (NPERM + 1)
      }
    }
    rows[[length(rows) + 1]] <- out
    cat(sprintf("scenario %d rep %d done\n", s, r))
  }
}
res <- bind_rows(rows)
write_csv(res, file.path(sim, "results/perm_country.csv"))
res |> group_by(scen) |>
  summarise(across(starts_with("Q_"), mean), across(starts_with("p_"), ~ mean(.x < .05), .names = "rej_{.col}")) |>
  print(width = 200)
