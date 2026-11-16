# Moment-matching refinement in R (same Newton scheme as refine_tp.do). Env MODEL, START (comma list, formula order), ITER, NS, IV, DAMP, OUT, SEED0.
suppressMessages({library(ergm); library(ergm.exclusive)}); source("enrichment/analysis/ergm_excl_data.R"); source("enrichment/analysis/ergm_excl_terms.R")
model <- Sys.getenv("MODEL"); th <- as.numeric(strsplit(Sys.getenv("START"), ",")[[1]]); out <- Sys.getenv("OUT", model)
niter <- as.integer(Sys.getenv("ITER", "20")); ns <- as.integer(Sys.getenv("NS", "40")); iv <- as.numeric(Sys.getenv("IV", "5e5"))
damp <- as.numeric(Sys.getenv("DAMP", "0.5")); seed0 <- as.integer(Sys.getenv("SEED0", "100")); cap <- as.numeric(Sys.getenv("CAP", "0.1"))
f <- excl_formula(model, fb); gobs <- summary(f); K <- length(gobs); stopifnot(length(th) == K)
dir.create("enrichment/analysis/output/tp/rrefine", showWarnings = FALSE)
tr <- matrix(NA, niter, 2 * K + 1); colnames(tr) <- c("it", paste0("th_", names(gobs)), paste0("z_", names(gobs)))
for (it in 1:niter) {
  s <- simulate(f, coef = th, constraints = ~exclusive(fb$types), nsim = ns, output = "stats",
                control = control.simulate.formula(MCMC.burnin = 3e6, MCMC.interval = iv), seed = seed0 + it)
  m <- colMeans(s); V <- cov(s); d <- gobs - m; z <- d / (sqrt(diag(V)) / sqrt(ns))
  step <- solve(V, d); if (max(abs(step)) > cap) step <- step * cap / max(abs(step))
  th <- th + damp * step
  tr[it, ] <- c(it, th, z); cat(sprintf("it %d max|z| %.1f\n", it, max(abs(z))), file = stderr())
  write.csv(tr, sprintf("enrichment/analysis/output/tp/rrefine/trace_%s.csv", out), row.names = FALSE)
}
thA <- colMeans(tr[(niter %/% 2 + 1):niter, 2:(K + 1), drop = FALSE])
write.csv(data.frame(name = names(gobs), theta = thA), sprintf("enrichment/analysis/output/tp/rrefine/theta_%s.csv", out), row.names = FALSE)
print(round(thA, 4))
