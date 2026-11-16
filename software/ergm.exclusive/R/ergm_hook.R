#' Exclusive-affiliation constraint for \code{ergm}
#'
#' Constrains the sample space of a bipartite network (actors = mode 1, affiliations = mode 2) to the arrangements in which
#' every actor has exactly one tie to an affiliation of each type and every affiliation keeps its size. Use as
#' \code{constraints = ~exclusive(types)}, where \code{types} is the type (1, 2, ...) of each mode-2 node. The sampler is the
#' typed swap (proposal \code{ExclusiveSwap}).
#' @param nw,arglist,... passed by \code{ergm}.
#' @export
InitErgmConstraint.exclusive <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = "types", vartypes = "numeric", defaultvalues = list(NULL), required = TRUE)
  if (!network::is.bipartite(nw)) ergm::ergm_Init_stop("The exclusive constraint needs a bipartite network (actors vs. affiliations).")
  nb2 <- network::network.size(nw) - network::get.network.attribute(nw, "bipartite")
  if (length(a$types) != nb2) ergm::ergm_Init_stop("types must have one entry per mode-2 node (", nb2, ").")
  B <- network::as.sociomatrix(nw)   # actors x affiliations for a bipartite network
  for (t in sort(unique(a$types))) {
    if (!all(rowSums(B[, a$types == t, drop = FALSE]) == 1))
      ergm::ergm_Init_stop("The network violates the exclusive constraint: every actor must have exactly one tie to an affiliation of type ", t, ".")
  }
  list(types = a$types, constrain = "exclusive", dependence = TRUE)
}

#' @rdname InitErgmConstraint.exclusive
#' @param arguments passed by \code{ergm}.
#' @export
InitErgmProposal.ExclusiveSwap <- function(arguments, nw) {
  types <- arguments$constraints$exclusive$types
  list(name = "ExclusiveSwap", inputs = as.double(types))
}

.onLoad <- function(libname, pkgname) {
  ergm:::ergm_proposal_table(Class = "c", Reference = "Bernoulli", Constraints = "&exclusive", Priority = 10,
                             Weights = "random", Proposal = "ExclusiveSwap")
}
