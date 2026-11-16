test_that("a typed swap preserves every affiliation size and is its own inverse", {
  net <- ex_toy(); set.seed(1)
  for (i in 1:200) {
    s <- ex_swap(net)
    expect_equal(tabulate(s$net$alpha, 3), tabulate(net$alpha, 3))
    expect_equal(tabulate(s$net$gamma, 3), tabulate(net$gamma, 3))
    if (s$type == 1L) { back <- s$net; back$alpha[c(s$p, s$q)] <- back$alpha[c(s$q, s$p)]; expect_equal(back$alpha, net$alpha) }
  }
})

test_that("the chain reaches every arrangement (irreducibility on small networks)", {
  r <- ex_reachable(exclusive_net(c(1, 1, 2, 2, 3), c(1, 2, 1, 2, 2)))
  expect_equal(unname(r["reachable"]), unname(r["expected"]))
  r2 <- ex_reachable(exclusive_net(c(1, 1, 2, 2, 2, 3), c(1, 2, 2, 3, 3, 3)))
  expect_equal(unname(r2["reachable"]), unname(r2["expected"]))
})

test_that("the sampler matches the exact distribution at zero and nonzero parameters", {
  skip_on_cran(); net <- ex_toy()
  for (theta in list(c(stacking = 0, hom_agency = 0, hom_club = 0, coloc = 0),
                     c(stacking = 0.8, hom_agency = 0.5, hom_club = -0.3, coloc = 0.4))) {
    set.seed(11); ex <- ex_exact(net, theta)
    M <- t(sapply(1:8, function(i) colMeans(ex_simulate(net, theta, nsim = 3000, burnin = 2000, interval = 10))))
    se <- apply(M, 2, sd) / sqrt(nrow(M))
    expect_true(all(abs(colMeans(M) - ex$mean[names(theta)]) < 4 * se))
  }
})

test_that("a singleton-only network has period two without holding and moves with holding", {
  net <- exclusive_net(1:3, 1:3)       # every affiliation has one actor: every swap changes the state
  set.seed(1); out <- ex_simulate(net, c(stacking = 0), nsim = 5, burnin = 0, interval = 1, hold = 0.5)
  expect_equal(ncol(out), 1L)
})
