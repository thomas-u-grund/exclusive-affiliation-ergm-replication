# What would a projection analysis conclude on the football network? Same analyses as the simulation study (sim/sim_fit_r.R), applied to the real data:
# stacking permutation tests (plain; clubs within nationality; agencies within club country; within country x nationality) and the binary agency-club projection ERGM.
suppressMessages({library(ergm); library(tidyverse)}); source("enrichment/analysis/ergm_excl_data.R")
B <- as.matrix(fb$net); n1 <- fb$n1; ag <- apply(B[, fb$isag == 1], 1, which.max); cl <- apply(B[, fb$iscl == 1], 1, which.max); nat <- fb$nat
cty <- fb$key_cty[fb$iscl == 1][cl]
stat <- function(a, c) sum(choose(table(a, c), 2)); obs <- stat(ag, cl); NP <- 1000; set.seed(42)
perm <- function(v, strata) { r <- v; for (k in unique(strata)) { i <- which(strata == k); r[i] <- r[i][sample.int(length(i))] }; r }
null <- list(plain = replicate(NP, stat(ag, sample(cl))),
             club_in_nat = replicate(NP, stat(ag, perm(cl, nat))),
             ag_in_cty = replicate(NP, stat(perm(ag, cty), cl)),
             ag_in_cty_nat = replicate(NP, stat(perm(ag, interaction(cty, nat)), cl)))
res <- tibble(test = names(null), observed = obs, null_mean = sapply(null, mean), Q = obs / sapply(null, mean), p = sapply(null, \(x) (1 + sum(x >= obs)) / (NP + 1)))
print(res)
# projection ERGM
A <- (table(factor(ag, 1:sum(fb$isag)), factor(cl, 1:sum(fb$iscl))) > 0) * 1; keep <- rowSums(A) > 0
net <- network(unclass(A[keep, ]), bipartite = sum(keep), directed = FALSE, matrix.type = "adjacency")
net %v% "country" <- c(rep(NA, sum(keep)), fb$key_cty[fb$iscl == 1])
fit <- ergm(net ~ b2nodematch("country"), constraints = ~degrees, control = control.ergm(seed = 1, MCMLE.maxit = 30, parallel = 0))
print(summary(fit)$coefficients)
dir.create("enrichment/analysis/output/tp/proj", showWarnings = FALSE)
write.csv(res, "enrichment/analysis/output/tp/proj/football_perm.csv", row.names = FALSE)
write.csv(summary(fit)$coefficients, "enrichment/analysis/output/tp/proj/football_projergm.csv")
