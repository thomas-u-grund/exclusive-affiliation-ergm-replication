toy_net2 <- function(x) {
  n <- x$n; B <- matrix(0, n, x$nA + x$nC); B[cbind(1:n, x$alpha)] <- 1; B[cbind(1:n, x$nA + x$gamma)] <- 1
  network::network(B, bipartite = n, directed = FALSE, matrix.type = "bipartite")
}
fast_terms <- function(x) {
  n2 <- x$nA + x$nC; W <- matrix(c(0.3, -0.2, 0.5), n2, n2, byrow = TRUE) * outer(1:n2, 1:n2, "+") / n2; W <- (W + t(W)) / 2
  M <- matrix(sin(seq_len(x$n * n2)), x$n, n2)
  list(list(type = "stack"), list(type = "homcentre", centre = c(rep(1, x$nA), rep(0, x$nC)), label = "hom_a"),
       list(type = "homcentre", centre = c(rep(0, x$nA), rep(1, x$nC)), label = "hom_c"),
       list(type = "coloc", centre = c(0.7, -1.2, 2, rep(0, x$nC)), key = c(rep(0, x$nA), x$club_cty), label = "coloc"),
       list(type = "acov", W = W), list(type = "edgecov", M = M))
}

test_that("fast sampler: initial statistics equal the ergm terms and tracked statistics equal a recomputation", {
  skip_on_cran(); toy <- ex_toy(); net <- toy_net2(toy); tm <- fast_terms(toy); types <- c(rep(1, toy$nA), rep(2, toy$nC))
  f <- net ~ stack() + homcentre(c(rep(1, 3), rep(0, 3)), toy$nat) + homcentre(c(rep(0, 3), rep(1, 3)), toy$nat) +
    coloc(c(0.7, -1.2, 2, 0, 0, 0), c(0, 0, 0, toy$club_cty)) + acov(tm[[5]]$W) + edgecov(tm[[6]]$M)
  r0 <- ex_fast_simulate(net, types, toy$nat, tm, rep(0, 6), nsim = 1, burnin = 0, interval = 1)
  expect_equal(unname(r0$stats[1, ]), unname(summary(f)))
  r <- ex_fast_simulate(net, types, toy$nat, tm, c(0.3, 0.2, -0.1, 0.4, 0.5, 0.3), nsim = 3, burnin = 2000, interval = 500, seed = 4)
  B <- matrix(0, toy$n, toy$nA + toy$nC); B[cbind(1:toy$n, r$alpha)] <- 1; B[cbind(1:toy$n, r$gamma)] <- 1
  net2 <- network::network(B, bipartite = toy$n, directed = FALSE, matrix.type = "bipartite")
  expect_equal(unname(r$stats[3, ]), unname(summary(update(f, net2 ~ .))), tolerance = 1e-8)
})

test_that("fast sampler matches the exact distribution (toy network, nonzero parameters)", {
  skip_on_cran(); toy <- ex_toy(); net <- toy_net2(toy); types <- c(rep(1, toy$nA), rep(2, toy$nC))
  tm <- list(list(type = "stack", label = "stacking"), list(type = "homcentre", centre = c(rep(1, 3), rep(0, 3)), label = "hom_agency"),
             list(type = "homcentre", centre = c(rep(0, 3), rep(1, 3)), label = "hom_club"),
             list(type = "coloc", centre = c(1, 1, 1, 0, 0, 0), key = c(0, 0, 0, toy$club_cty), label = "coloc"))
  theta <- c(stacking = 0.8, hom_agency = 0.5, hom_club = -0.3, coloc = 0.4); ex <- ex_exact(ex_toy(), theta)
  M <- t(sapply(1:8, function(i) colMeans(ex_fast_simulate(net, types, toy$nat, tm, theta, nsim = 3000, burnin = 2000, interval = 10, seed = i)$stats)))
  se <- apply(M, 2, sd) / sqrt(nrow(M))
  expect_true(all(abs(colMeans(M) - ex$mean[names(theta)]) < 4 * se))
})
