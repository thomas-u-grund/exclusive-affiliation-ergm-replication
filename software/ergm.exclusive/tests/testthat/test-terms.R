toy_net <- function(x) {
  n <- x$n; B <- matrix(0, n, x$nA + x$nC); B[cbind(1:n, x$alpha)] <- 1; B[cbind(1:n, x$nA + x$gamma)] <- 1
  network::network(B, bipartite = n, directed = FALSE, matrix.type = "bipartite")
}
toy_masks <- function(x) list(isag = c(rep(1, x$nA), rep(0, x$nC)), iscl = c(rep(0, x$nA), rep(1, x$nC)),
                              key = c(rep(0, x$nA), x$club_cty), types = c(rep(1, x$nA), rep(2, x$nC)))

test_that("full statistics from the C change statistics equal the direct computation on random arrangements", {
  skip_on_cran(); toy <- ex_toy(); m <- toy_masks(toy); set.seed(5)
  W <- matrix(runif((toy$nA + toy$nC)^2), toy$nA + toy$nC); W <- (W + t(W)) / 2; diag(W) <- 0
  for (i in 1:25) {
    x <- toy
    x$alpha <- sample(toy$alpha); x$gamma <- sample(toy$gamma)
    net <- toy_net(x); isag <- m$isag; iscl <- m$iscl; key <- m$key; nat <- x$nat
    s <- summary(net ~ b1dsp(2) + homcentre(isag, nat) + homcentre(iscl, nat) + coloc(isag, key) + acov(W))
    d <- ex_stats(x)
    expect_equal(unname(s[1:4]), unname(d[c("stacking", "hom_agency", "hom_club", "coloc")]))
    expect_equal(unname(s[5]), sum(W[cbind(x$alpha, x$nA + x$gamma)]))
  }
})

test_that("coloc with club identities as keys and a subset of centres counts stacking among those agencies", {
  skip_on_cran(); toy <- ex_toy(); m <- toy_masks(toy); set.seed(8)
  small <- c(1, 0, 1, 0, 0, 0); clubkey <- c(0, 0, 0, 1, 2, 3)
  for (i in 1:10) {
    x <- toy; x$alpha <- sample(toy$alpha); x$gamma <- sample(toy$gamma); net <- toy_net(x)
    s <- summary(net ~ coloc(small, clubkey))
    pr <- which(upper.tri(matrix(0, x$n, x$n)), arr.ind = TRUE)
    direct <- sum(small[x$alpha[pr[, 1]]] == 1 & x$alpha[pr[, 1]] == x$alpha[pr[, 2]] & x$gamma[pr[, 1]] == x$gamma[pr[, 2]])
    expect_equal(unname(s), direct)
  }
})

test_that("simulate() with all terms reproduces the exact means at nonzero parameters", {
  skip_on_cran(); toy <- ex_toy(); m <- toy_masks(toy); net <- toy_net(toy)
  isag <- m$isag; iscl <- m$iscl; key <- m$key; types <- m$types; nat <- toy$nat
  theta <- c(stacking = 0.6, hom_agency = 0.4, hom_club = -0.3, coloc = 0.35)
  ex <- ex_exact(toy, theta)
  M <- t(sapply(1:20, function(i) {
    s <- suppressMessages(simulate(net ~ b1dsp(2) + homcentre(isag, nat) + homcentre(iscl, nat) + coloc(isag, key),
                                   coef = unname(theta), nsim = 4000, output = "stats", seed = 90 + i, constraints = ~exclusive(types),
                                   control = ergm::control.simulate.formula(MCMC.burnin = 2000, MCMC.interval = 25)))
    colMeans(s) }))
  se <- apply(M, 2, sd) / sqrt(nrow(M))
  expect_true(all(abs(colMeans(M) - ex$mean[names(theta)]) < 4 * se))
})

test_that("coloc with centre weights weights every co-located pair by its centre's weight", {
  skip_on_cran(); toy <- ex_toy(); set.seed(9)
  w <- c(0.7, -1.2, 2, 0, 0, 0); clubkey <- c(0, 0, 0, 1, 2, 3)
  for (i in 1:10) {
    x <- toy; x$alpha <- sample(toy$alpha); x$gamma <- sample(toy$gamma); net <- toy_net(x)
    s <- summary(net ~ coloc(w, clubkey))
    pr <- which(upper.tri(matrix(0, x$n, x$n)), arr.ind = TRUE)
    same <- x$alpha[pr[, 1]] == x$alpha[pr[, 2]] & x$gamma[pr[, 1]] == x$gamma[pr[, 2]]
    expect_equal(unname(s), sum(w[x$alpha[pr[, 1]]][same]))
  }
})

test_that("stack equals b1dsp(2) and matches a brute-force count, in the stationary chain and from scratch", {
  skip_on_cran(); toy <- ex_toy(); set.seed(11)
  for (i in 1:10) {
    x <- toy; x$alpha <- sample(toy$alpha); x$gamma <- sample(toy$gamma); net <- toy_net(x)
    pr <- which(upper.tri(matrix(0, x$n, x$n)), arr.ind = TRUE)
    same <- x$alpha[pr[, 1]] == x$alpha[pr[, 2]] & x$gamma[pr[, 1]] == x$gamma[pr[, 2]]
    expect_equal(unname(summary(net ~ stack())), sum(same))
    expect_equal(unname(summary(net ~ stack())), unname(summary(net ~ b1dsp(2))))
  }
})

test_that("dstack counts the agency-club cells with at least two actors", {
  skip_on_cran(); toy <- ex_toy(); set.seed(13)
  for (i in 1:10) {
    x <- toy; x$alpha <- sample(toy$alpha); x$gamma <- sample(toy$gamma); net <- toy_net(x)
    n2 <- network::network.size(net) - network::get.network.attribute(net, "bipartite")
    cen <- c(rep(1, max(x$alpha)), rep(0, n2 - max(x$alpha)))
    expect_equal(unname(summary(net ~ dstack(cen))), sum(table(x$alpha, x$gamma) >= 2))
  }
})
