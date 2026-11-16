#' One typed swap
#'
#' Draw an actor p uniformly, a type uniformly (agency or club), and a second actor q uniformly among the others; exchange
#' the type-t partners of p and q. The proposal is symmetric: the same exchange is proposed from the new state with the
#' same probability. Margins (all agency and club sizes) and the one-of-each-type rule are preserved.
#' @param net an \code{exclusive_net}.
#' @return list with the proposed network \code{net}, and \code{type}, \code{p}, \code{q}.
#' @export
ex_swap <- function(net) {
  n <- net$n
  p <- sample.int(n, 1L); q <- sample.int(n - 1L, 1L); if (q >= p) q <- q + 1L
  type <- sample.int(2L, 1L)
  if (type == 1L) { net$alpha[c(p, q)] <- net$alpha[c(q, p)] } else { net$gamma[c(p, q)] <- net$gamma[c(q, p)] }
  list(net = net, type = type, p = p, q = q)
}

#' Simulate from the exclusive-affiliation ERGM with the typed-swap sampler
#'
#' Metropolis-Hastings with the symmetric typed-swap proposal; acceptance probability \code{min(1, exp(theta' (g' - g)))}.
#' This is a reference implementation (statistics recomputed from scratch at each step); use it for small networks and for
#' checking other implementations.
#' @param net starting \code{exclusive_net}.
#' @param theta named numeric vector of coefficients (names from \code{ex_stats}); omitted terms have coefficient 0.
#' @param nsim number of draws.
#' @param burnin,interval number of proposals before the first draw and between draws.
#' @param hold probability that a step does nothing (makes the chain aperiodic when all affiliations are singletons).
#' @return matrix with one row per draw and one column per statistic in \code{names(theta)}.
#' @export
ex_simulate <- function(net, theta, nsim = 1000, burnin = 1000, interval = 10, hold = 0) {
  terms <- names(theta)
  g <- ex_stats(net, terms)
  out <- matrix(NA_real_, nsim, length(terms), dimnames = list(NULL, terms))
  total <- burnin + (nsim - 1) * interval
  k <- 0L
  for (s in 0:total) {
    if (s >= burnin && (s - burnin) %% interval == 0) { k <- k + 1L; out[k, ] <- g }
    if (s == total) break
    if (hold > 0 && runif(1) < hold) next
    prop <- ex_swap(net)
    g2 <- ex_stats(prop$net, terms)
    if (log(runif(1)) < sum(theta * (g2 - g))) { net <- prop$net; g <- g2 }
  }
  out
}

#' States reachable by typed swaps (small networks)
#'
#' Breadth-first search over arrangements; returns the number of distinct states reachable from \code{net} and the number
#' of states expected from the multinomial formula, so that irreducibility can be checked directly on small examples.
#' @param net an \code{exclusive_net} with a small number of actors.
#' @export
ex_reachable <- function(net) {
  key <- function(a, g) paste(c(a, g), collapse = ",")
  seen <- new.env(); queue <- list(list(a = net$alpha, g = net$gamma)); seen[[key(net$alpha, net$gamma)]] <- TRUE
  n <- net$n; head <- 1L
  while (head <= length(queue)) {
    s <- queue[[head]]; head <- head + 1L
    for (p in 1:(n - 1)) for (q in (p + 1):n) for (type in 1:2) {
      a <- s$a; g <- s$g
      if (type == 1L) a[c(p, q)] <- a[c(q, p)] else g[c(p, q)] <- g[c(q, p)]
      k <- key(a, g)
      if (is.null(seen[[k]])) { seen[[k]] <- TRUE; queue[[length(queue) + 1L]] <- list(a = a, g = g) }
    }
  }
  expected <- exp(2 * lfactorial(n) - sum(lfactorial(tabulate(net$alpha))) - sum(lfactorial(tabulate(net$gamma))))
  c(reachable = length(queue), expected = expected)
}
