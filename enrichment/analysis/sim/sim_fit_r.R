# ============================================================
# Simulation study, step 4: the conditional model and the projection analyses (R ergm).
# For each simulated network of a scenario:
#   cond   agency-player ERGM given the clubs (players x agencies, constraints = ~degrees):
#          b1nodematch(club) = S, b1nodematch(nat) = H_A, b1nodematch(cty) = L, edgecov(w_ap) = T_W
#   perm   stacking permutation tests on the agency-club counts: Q = observed / null mean of sum C(w_ac, 2)
#          plain (clubs permuted across all players) and stratified (within nationality)
#   proj   binary agency-club projection ERGM (constraints = ~degrees):
#          b2nodematch(country) = clubs sharing an agency in the same country, edgecov(W) = agency home country
# Usage: Rscript enrichment/analysis/sim/sim_fit_r.R <scenario> <first rep> <last rep>
# Appends to sim/results/r_s<k>.csv (one row per replication).
# ============================================================
suppressMessages({library(tidyverse); library(ergm)})
args <- commandArgs(trailingOnly = TRUE)
s <- as.integer(args[1]); r1 <- as.integer(args[2]); r2 <- as.integer(args[3])
sim <- "enrichment/analysis/sim"
nodes <- read_csv(file.path(sim, "world/nodes.csv"), show_col_types = FALSE, guess_max = 1e6)
NP <- sum(nodes$kind == 1); NA_ <- sum(nodes$kind == 2); NC <- sum(nodes$kind == 3)
nat <- nodes$nat[nodes$kind == 1]
aghome <- nodes$ag_home[nodes$kind == 2]
clubcty <- nodes$club_country[nodes$kind == 3]
Wac <- outer(aghome, clubcty, "==") * 1                       # agencies x clubs
nets <- read_csv(file.path(sim, sprintf("nets/s%d_edges.csv", s)), show_col_types = FALSE) |>
  rename(sim = `_sim`, p = `_tail`, x = `_head`)
outf <- file.path(sim, sprintf("results/r_s%d.csv", s))
NPERM <- 1000

coefrow <- function(fit, prefix) {
  if (is.null(fit)) return(tibble())
  b <- coef(fit); se <- sqrt(diag(vcov(fit)))
  nm <- gsub("[^A-Za-z0-9]+", "_", names(b))
  bind_cols(as_tibble(setNames(as.list(b), paste0(prefix, "b_", nm))),
            as_tibble(setNames(as.list(se), paste0(prefix, "se_", nm))))
}

for (r in r1:r2) {
  t0 <- Sys.time()
  e <- nets |> filter(sim == r)
  ag <- integer(NP); cl <- integer(NP)
  isag <- e$x <= NP + NA_
  ag[e$p[isag]] <- e$x[isag] - NP
  cl[e$p[!isag]] <- e$x[!isag] - NP - NA_
  cty <- clubcty[cl]

  # ---- conditional model: players x agencies given clubs
  B <- matrix(0, NP, NA_); B[cbind(seq_len(NP), ag)] <- 1
  net_ap <- network(B, bipartite = NP, directed = FALSE, matrix.type = "adjacency")
  net_ap %v% "club" <- c(cl, rep(NA, NA_))
  net_ap %v% "nat" <- c(nat, rep(NA, NA_))
  net_ap %v% "cty" <- c(cty, rep(NA, NA_))
  w_ap <- Wac[, cl] |> t()                                    # players x agencies: W[agency, club(player)]
  fit_c <- tryCatch(
    ergm(net_ap ~ b1nodematch("club") + b1nodematch("nat") + b1nodematch("cty") + edgecov(w_ap),
         constraints = ~degrees,
         control = control.ergm(seed = r, MCMLE.maxit = 30, parallel = 0)),
    error = function(err) { message("cond failed: ", conditionMessage(err)); NULL })

  # ---- stacking permutation tests on agency-club counts
  stack_stat <- function(cl_) { w <- table(ag, cl_); sum(choose(w, 2)) }
  obs <- stack_stat(cl)
  set.seed(r)
  null_plain <- replicate(NPERM, stack_stat(sample(cl)))
  null_strat <- replicate(NPERM, { c2 <- cl; for (k in unique(nat)) { i <- which(nat == k); c2[i] <- c2[i][sample.int(length(i))] }; stack_stat(c2) })

  # ---- binary agency-club projection ERGM
  A <- (table(factor(ag, levels = seq_len(NA_)), factor(cl, levels = seq_len(NC))) > 0) * 1
  keep_ag <- rowSums(A) > 0
  net_ac <- network(unclass(A[keep_ag, ]), bipartite = sum(keep_ag), directed = FALSE, matrix.type = "adjacency")
  net_ac %v% "country" <- c(rep(NA, sum(keep_ag)), clubcty)
  W_k <- Wac[keep_ag, ]
  fit_p <- tryCatch(
    ergm(net_ac ~ b2nodematch("country") + edgecov(W_k), constraints = ~degrees,
         control = control.ergm(seed = r, MCMLE.maxit = 30, parallel = 0)),
    error = function(err) { message("proj failed: ", conditionMessage(err)); NULL })

  row <- bind_cols(
    tibble(scen = s, rep = r, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
           stack_obs = obs,
           Q_plain = obs / mean(null_plain), p_plain = mean(null_plain >= obs),
           Q_strat = obs / mean(null_strat), p_strat = mean(null_strat >= obs)),
    coefrow(fit_c, "cond_"), coefrow(fit_p, "proj_"))
  write_csv(row, outf, append = file.exists(outf))
  cat(sprintf("scenario %d rep %d: %.0f s, Q_plain %.2f (p %.3f), Q_strat %.2f (p %.3f)\n",
              s, r, row$seconds, row$Q_plain, row$p_plain, row$Q_strat, row$p_strat))
}
