# Recovery at application scale with the fast sampler: draw one network state from the fitted model (3c), take its statistics as "observed",
# re-estimate by the same moment-matching procedure as the application (start 0.9*truth), SE from the inverse covariance at the average.
# Env REP, MODEL (3s3), TRUTH (csv), RECDIR, ITER (24), NS (400), IV (5e5).
suppressMessages({library(ergm); library(ergm.exclusive)}); source("enrichment/analysis/ergm_excl_data.R"); source("enrichment/analysis/ergm_excl_terms.R"); source("enrichment/analysis/ergm_excl_fast.R")
rep <- as.integer(Sys.getenv("REP")); model <- Sys.getenv("MODEL", "3s3"); niter <- as.integer(Sys.getenv("ITER", "24")); ns <- as.integer(Sys.getenv("NS", "400")); iv <- as.numeric(Sys.getenv("IV", "5e5"))
tt <- read.csv(Sys.getenv("TRUTH", "enrichment/analysis/output/tp/rrefine/theta_3s3F.csv"))$theta
f <- excl_formula(model, fb); K <- length(summary(f)); stopifnot(length(tt) == K); tm <- fast_terms(f); nmx <- names(summary(f))
gobs <- ex_fast_simulate(fb$net, fb$types, fb$nat, tm, tt, nsim = 1, burnin = 2e7, interval = 1, seed = 70000 + rep)$stats[1, ]
th <- 0.9 * tt; damp <- 0.5; cap <- 0.2; TH <- matrix(NA, niter, K)
for (it in 1:niter) {
  s <- ex_fast_simulate(fb$net, fb$types, fb$nat, tm, th, nsim = ns, burnin = 3e6, interval = iv, seed = 71000 + 100 * rep + it)$stats
  d <- gobs - colMeans(s); step <- solve(cov(s), d); if (max(abs(step)) > cap) step <- step * cap / max(abs(step))
  th <- th + damp * step; TH[it, ] <- th
}
thA <- colMeans(TH[(niter %/% 2 + 1):niter, , drop = FALSE])
s <- ex_fast_simulate(fb$net, fb$types, fb$nat, tm, thA, nsim = 2 * ns, burnin = 3e6, interval = iv, seed = 72000 + rep)$stats
se <- sqrt(diag(solve(cov(s)))); zm <- (gobs - colMeans(s)) / (sqrt(diag(cov(s))) / sqrt(nrow(s)))
d <- file.path("enrichment/analysis/output/tp", Sys.getenv("RECDIR", "recover_fast")); dir.create(d, showWarnings = FALSE)
write.csv(data.frame(rep = rep, name = nmx, truth = tt, est = thA, se = se, zmoment = zm), file.path(d, sprintf("rec_%03d.csv", rep)), row.names = FALSE)
