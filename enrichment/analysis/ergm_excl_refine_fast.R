# Moment-matching refinement with the compiled fast sampler (ergm.exclusive::ex_fast_simulate); same Newton scheme as ergm_excl_refine.R.
# Env MODEL, START (comma list, formula order), OUT, ITER, NS, IV, DAMP, CAP, SEED0, BURN; MKT_MIN as in ergm_excl_data.R.
suppressMessages({library(ergm); library(ergm.exclusive)}); source("enrichment/analysis/ergm_excl_data.R"); source("enrichment/analysis/ergm_excl_terms.R"); source("enrichment/analysis/ergm_excl_fast.R")
model <- Sys.getenv("MODEL"); th <- as.numeric(strsplit(Sys.getenv("START"), ",")[[1]]); out <- Sys.getenv("OUT", model)
niter <- as.integer(Sys.getenv("ITER", "30")); ns <- as.integer(Sys.getenv("NS", "400")); iv <- as.numeric(Sys.getenv("IV", "5e5"))
damp <- as.numeric(Sys.getenv("DAMP", "0.5")); seed0 <- as.integer(Sys.getenv("SEED0", "100")); cap <- as.numeric(Sys.getenv("CAP", "0.2")); burn <- as.numeric(Sys.getenv("BURN", "3e6"))
f <- excl_formula(model, fb); gobs <- summary(f); K <- length(gobs); stopifnot(length(th) == K); tm <- fast_terms(f)
dir.create("enrichment/analysis/output/tp/rrefine", showWarnings = FALSE)
tr <- matrix(NA, niter, 2 * K + 1); colnames(tr) <- c("it", paste0("th_", names(gobs)), paste0("z_", names(gobs)))
for (it in 1:niter) {
  s <- ex_fast_simulate(fb$net, fb$types, fb$nat, tm, th, nsim = ns, burnin = burn, interval = iv, seed = seed0 + it)$stats
  m <- colMeans(s); V <- cov(s); d <- gobs - m; z <- d / (sqrt(diag(V)) / sqrt(ns))
  step <- solve(V, d); if (max(abs(step)) > cap) step <- step * cap / max(abs(step))
  th <- th + damp * step; tr[it, ] <- c(it, th, z); cat(sprintf("it %d max|z| %.1f\n", it, max(abs(z))), file = stderr())
  write.csv(tr, sprintf("enrichment/analysis/output/tp/rrefine/trace_%s.csv", out), row.names = FALSE)
}
thA <- colMeans(tr[(niter %/% 2 + 1):niter, 2:(K + 1), drop = FALSE])
write.csv(data.frame(name = names(gobs), theta = thA), sprintf("enrichment/analysis/output/tp/rrefine/theta_%s.csv", out), row.names = FALSE); print(round(thA, 4))
