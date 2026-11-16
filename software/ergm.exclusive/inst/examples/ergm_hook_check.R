# Check of the ergm() integration: simulate from  net ~ b1dsp(2)  with  constraints = ~exclusive(types)  and compare with
# exact enumeration on the toy network. Run: Rscript inst/examples/ergm_hook_check.R
suppressMessages({library(ergm); library(ergm.exclusive)})
toy <- ex_toy(); n <- toy$n; nA <- toy$nA; nC <- toy$nC
B <- matrix(0, n, nA + nC); B[cbind(1:n, toy$alpha)] <- 1; B[cbind(1:n, nA + toy$gamma)] <- 1
net <- network::network(B, bipartite = n, directed = FALSE, matrix.type = "bipartite")
types <- c(rep(1, nA), rep(2, nC))
one <- function(theta, nchain = 24, nsim = 5000, seed = 100) {
  ex <- ex_exact(toy, c(stacking = theta))
  M <- sapply(1:nchain, function(i) {
    s <- suppressMessages(simulate(net ~ b1dsp(2), coef = theta, nsim = nsim, output = "stats", seed = seed + i,
                  constraints = ~exclusive(types), control = control.simulate.formula(MCMC.burnin = 2000, MCMC.interval = 25)))
    c(mean(s), table(factor(s, 0:4)) / nsim) })
  m <- rowMeans(M); se <- apply(M, 1, sd) / sqrt(nchain)
  pmf <- ex$pmf[as.character(0:4)]; pmf[is.na(pmf)] <- 0
  data.frame(stat = c("mean S", paste0("P(S=", 0:4, ")")), exact = c(ex$mean_value, pmf), ergm_exclusive = m, mc_se = se,
             z = (m - c(ex$mean_value, pmf)) / se)
}
cat("theta_S = 0\n");   print(one(0), digits = 4)
cat("\ntheta_S = 0.8\n"); print(one(0.8), digits = 4)
