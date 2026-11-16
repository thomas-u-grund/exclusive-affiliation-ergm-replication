#' An exclusive multi-affiliation network
#'
#' Every actor (player) holds exactly one affiliation of each of two types (agency and club). The state is the pair of
#' assignments \code{alpha} (actor to agency) and \code{gamma} (actor to club); the sizes of all affiliations are fixed.
#'
#' @param alpha integer vector, agency of each actor (labels \code{1..nA}).
#' @param gamma integer vector, club of each actor (labels \code{1..nC}).
#' @param nat actor attribute used by the homophily statistics (any vector of length n).
#' @param club_cty club attribute (length nC) used by the co-location statistic.
#' @return an object of class \code{exclusive_net}.
#' @export
exclusive_net <- function(alpha, gamma, nat = NULL, club_cty = NULL) {
  n <- length(alpha)
  stopifnot(length(gamma) == n)
  if (is.null(nat)) nat <- rep(1L, n)
  if (is.null(club_cty)) club_cty <- rep(1L, max(gamma))
  stopifnot(length(nat) == n, length(club_cty) >= max(gamma))
  structure(list(alpha = as.integer(alpha), gamma = as.integer(gamma), n = n, nA = max(alpha), nC = max(gamma),
                 nat = nat, club_cty = club_cty), class = "exclusive_net")
}

#' @export
print.exclusive_net <- function(x, ...) {
  cat("Exclusive multi-affiliation network: ", x$n, " actors, ", x$nA, " agencies (sizes ",
      paste(tabulate(x$alpha, x$nA), collapse = ","), "), ", x$nC, " clubs (sizes ",
      paste(tabulate(x$gamma, x$nC), collapse = ","), ")\n", sep = "")
  invisible(x)
}

# pairs p < q as two index vectors
.pairs <- function(n) { m <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE); list(p = m[, 1], q = m[, 2]) }

.stat_names <- c("stacking", "hom_agency", "hom_club", "coloc")

#' Statistics
#'
#' \describe{
#'  \item{stacking}{number of actor pairs with the same agency and the same club (co-clients at the same club); in the
#'   bipartite representation, the number of actor pairs with exactly two shared partners.}
#'  \item{hom_agency}{pairs with the same agency and the same value of \code{nat}.}
#'  \item{hom_club}{pairs with the same club and the same value of \code{nat}.}
#'  \item{coloc}{pairs with the same agency whose clubs are in the same country (\code{club_cty}).}
#' }
#' @param net an \code{exclusive_net}.
#' @param terms character vector of statistic names.
#' @return named numeric vector.
#' @export
ex_stats <- function(net, terms = .stat_names) {
  pr <- .pairs(net$n); p <- pr$p; q <- pr$q
  sa <- net$alpha[p] == net$alpha[q]; sc <- net$gamma[p] == net$gamma[q]; sn <- net$nat[p] == net$nat[q]
  sk <- net$club_cty[net$gamma[p]] == net$club_cty[net$gamma[q]]
  all <- c(stacking = sum(sa & sc), hom_agency = sum(sa & sn), hom_club = sum(sc & sn), coloc = sum(sa & sk))
  all[terms]
}
