# Empirical-scale recovery: simulate a network from the fitted model (3s) on the real football margins and covariates, then re-estimate
# by moment matching (as ergm_excl_refine.R) from a perturbed start. Env REP (replicate id, also the seed), ITER (10), NS (40), IV (5e5).
suppressMessages({library(ergm); library(ergm.exclusive)}); source("enrichment/analysis/ergm_excl_data.R"); source("enrichment/analysis/ergm_excl_terms.R")
rep <- as.integer(Sys.getenv("REP")); niter <- as.integer(Sys.getenv("ITER", "14")); ns <- as.integer(Sys.getenv("NS", "40")); iv <- as.numeric(Sys.getenv("IV", "5e5"))
model <- Sys.getenv("MODEL", "3ms2")
if (model == "3s3") { tt <- read.csv("enrichment/analysis/output/tp/rrefine/theta_3s3b.csv")$theta } else {
  tt <- read.csv("enrichment/analysis/output/tp/refine/final_3ms22bR2.csv")$theta[-1]; tt <- tt[c(1:6, 9:12, 7, 8)] }
f0 <- excl_formula(model, fb); stopifnot(length(tt) == length(summary(f0)))
sim <- simulate(f0, coef = tt, constraints = ~exclusive(fb$types), nsim = 1, output = "network",
  control = control.simulate.formula(MCMC.burnin = 2e7, MCMC.interval = 1), seed = 70000 + rep)
fb$net <- sim; f <- excl_formula(model, fb); gobs <- summary(f); K <- length(gobs)
th <- 0.9 * tt; damp <- 0.8; cap <- 0.3; TH <- matrix(NA, niter, K)
for (it in 1:niter) {
  s <- simulate(f, coef = th, constraints = ~exclusive(fb$types), nsim = ns, output = "stats",
    control = control.simulate.formula(MCMC.burnin = 3e6, MCMC.interval = iv), seed = 71000 + 100 * rep + it)
  d <- gobs - colMeans(s); step <- solve(cov(s), d); if (max(abs(step)) > cap) step <- step * cap / max(abs(step))
  th <- th + damp * step; TH[it, ] <- th
}
thA <- colMeans(TH[(niter %/% 2 + 1):niter, , drop = FALSE])
s <- simulate(f, coef = thA, constraints = ~exclusive(fb$types), nsim = 2 * ns, output = "stats",
  control = control.simulate.formula(MCMC.burnin = 3e6, MCMC.interval = iv), seed = 72000 + rep)
se <- sqrt(diag(solve(cov(s)))); zm <- (gobs - colMeans(s)) / (sqrt(diag(cov(s))) / sqrt(nrow(s)))
dir.create(file.path("enrichment/analysis/output/tp", Sys.getenv("RECDIR", "recover")), showWarnings = FALSE)
write.csv(data.frame(rep = rep, name = names(gobs), truth = tt, est = thA, se = se, zmoment = zm),
  file.path("enrichment/analysis/output/tp", Sys.getenv("RECDIR", "recover"), sprintf("rec_%03d.csv", rep)), row.names = FALSE)
