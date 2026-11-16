# Benchmark of the typed-swap sampler against exact enumeration on the 8-player toy network (Table "toy" of the paper).
# Run: Rscript inst/examples/toy_benchmark.R
suppressMessages(library(ergm.exclusive))
net <- ex_toy(); print(net)
bench <- function(theta, nchain = 10, nsim = 5000, seed = 1) {
  set.seed(seed); ex <- ex_exact(net, theta)
  runs <- lapply(seq_len(nchain), function(i) ex_simulate(net, theta, nsim = nsim, burnin = 2000, interval = 10))
  M <- t(sapply(runs, colMeans))
  data.frame(stat = names(theta), exact = ex$mean[names(theta)], sampler = colMeans(M),
             mc_se = apply(M, 2, sd) / sqrt(nchain), z = (colMeans(M) - ex$mean[names(theta)]) / (apply(M, 2, sd) / sqrt(nchain)))
}
cat("\nTheta = 0:\n");   print(bench(c(stacking = 0, hom_agency = 0, hom_club = 0, coloc = 0)), digits = 4)
cat("\nTheta != 0:\n");  print(bench(c(stacking = 0.8, hom_agency = 0.5, hom_club = -0.3, coloc = 0.4)), digits = 4)
