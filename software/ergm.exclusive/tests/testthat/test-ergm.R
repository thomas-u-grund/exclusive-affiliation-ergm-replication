test_that("simulate() with constraints = ~exclusive(types) reproduces the exact toy distribution", {
  skip_on_cran(); skip_if_not_installed("ergm")
  toy <- ex_toy(); n <- toy$n
  B <- matrix(0, n, toy$nA + toy$nC); B[cbind(1:n, toy$alpha)] <- 1; B[cbind(1:n, toy$nA + toy$gamma)] <- 1
  net <- network::network(B, bipartite = n, directed = FALSE, matrix.type = "bipartite")
  types <- c(rep(1, toy$nA), rep(2, toy$nC))
  for (theta in c(0, 0.8)) {
    ex <- ex_exact(toy, c(stacking = theta))
    M <- sapply(1:20, function(i) {
      s <- suppressMessages(simulate(net ~ b1dsp(2), coef = theta, nsim = 4000, output = "stats", seed = 40 + i,
                                     constraints = ~exclusive(types),
                                     control = ergm::control.simulate.formula(MCMC.burnin = 2000, MCMC.interval = 25)))
      mean(s) })
    expect_true(abs(mean(M) - ex$mean_value) < 4 * sd(M) / sqrt(length(M)))
  }
})

test_that("the exclusive constraint rejects a network that violates it", {
  skip_on_cran(); skip_if_not_installed("ergm")
  B <- matrix(0, 3, 4); B[1, 1] <- 1; B[2, 1] <- 1; B[3, 2] <- 1; B[1, 3] <- 1; B[2, 4] <- 1; B[3, 4] <- 1
  B[1, 2] <- 1                       # actor 1 now has two agencies
  net <- network::network(B, bipartite = 3, directed = FALSE, matrix.type = "bipartite")
  expect_error(suppressMessages(simulate(net ~ b1dsp(2), coef = 0, nsim = 2, output = "stats",
                                         constraints = ~exclusive(c(1, 1, 2, 2)))))
})
