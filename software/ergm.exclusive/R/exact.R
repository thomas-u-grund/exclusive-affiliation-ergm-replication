# all distinct permutations of a multiset (rows)
.multiset_perms <- function(v) {
  if (length(v) <= 1L) return(matrix(v, 1L))
  out <- NULL
  for (x in unique(v)) {
    rest <- v[-match(x, v)]
    out <- rbind(out, cbind(x, .multiset_perms(rest)))
  }
  unname(out)
}

#' Exact enumeration of the sample space of a small network
#'
#' Enumerates every pair of assignments with the fibre sizes of \code{net} (\eqn{n!/\prod a_k! \times n!/\prod c_l!} states)
#' and returns the statistics of each. The statistics are computed by matrix products of pair-indicator matrices.
#' @param net an \code{exclusive_net}.
#' @param terms statistic names.
#' @return list with \code{stats} (matrix, one row per state, alpha varying fastest) and \code{n_alpha}, \code{n_gamma}.
#' @export
ex_enumerate <- function(net, terms = .stat_names) {
  PA <- .multiset_perms(sort(net$alpha)); PG <- .multiset_perms(sort(net$gamma))
  pr <- .pairs(net$n); p <- pr$p; q <- pr$q
  sn <- as.numeric(net$nat[p] == net$nat[q])
  A <- (PA[, p, drop = FALSE] == PA[, q, drop = FALSE]) * 1           # same-agency indicator per pair
  G <- (PG[, p, drop = FALSE] == PG[, q, drop = FALSE]) * 1           # same-club indicator
  K <- matrix(net$club_cty[PG[, p]] == net$club_cty[PG[, q]], nrow(PG)) * 1   # clubs in the same country
  out <- list()
  if ("stacking" %in% terms)   out$stacking   <- as.vector(A %*% t(G))
  if ("hom_agency" %in% terms) out$hom_agency <- rep(as.vector(A %*% sn), times = nrow(PG))
  if ("hom_club" %in% terms)   out$hom_club   <- rep(as.vector(G %*% sn), each = nrow(PA))
  if ("coloc" %in% terms)      out$coloc      <- as.vector(A %*% t(K))
  list(stats = do.call(cbind, out[terms]), n_alpha = nrow(PA), n_gamma = nrow(PG))
}

#' Exact distribution under the ERGM
#'
#' @param net an \code{exclusive_net}.
#' @param theta named coefficients (names from \code{ex_stats}).
#' @param value statistic whose distribution is returned (default \code{"stacking"}).
#' @return list: \code{mean} of every statistic in \code{names(theta)} and the probability mass function of \code{value}.
#' @export
ex_exact <- function(net, theta, value = "stacking") {
  terms <- unique(c(names(theta), value))
  e <- ex_enumerate(net, terms)
  lw <- as.vector(e$stats[, names(theta), drop = FALSE] %*% theta); w <- exp(lw - max(lw)); w <- w / sum(w)
  v <- e$stats[, value]
  list(mean = colSums(e$stats[, names(theta), drop = FALSE] * w),
       pmf = tapply(w, v, sum), mean_value = sum(w * v), kappa_ratio = log(sum(exp(lw - max(lw)))) + max(lw) - log(length(lw)))
}

#' The eight-actor toy network of the paper
#'
#' Eight players, three agencies of sizes 3, 2, 3 and three clubs of sizes 2, 3, 3 (560 x 560 arrangements). Players have a
#' binary nationality and clubs are in two countries, so that all four statistics are non-trivial.
#' @export
ex_toy <- function() {
  exclusive_net(alpha = c(1, 1, 1, 2, 2, 3, 3, 3), gamma = c(1, 2, 2, 1, 3, 3, 3, 2),
                nat = c(1, 1, 2, 2, 1, 2, 1, 2), club_cty = c(1, 1, 2))
}
