.bip_dims <- function(nw) {
  if (!network::is.bipartite(nw)) ergm::ergm_Init_stop("This term needs a bipartite network (actors vs. affiliations).")
  bip <- network::get.network.attribute(nw, "bipartite")
  list(n1 = bip, n2 = network::network.size(nw) - bip)
}

#' Cross-layer terms for exclusive multi-affiliation networks
#'
#' Terms for \code{ergm} formulas on a bipartite network whose mode 1 is the set of actors and whose mode 2 is the set of
#' affiliations of all types (agencies, clubs, ...).
#' \describe{
#'  \item{\code{coloc(centre, key)}}{For every \dQuote{centre} affiliation (\code{centre} = 1) and every pair of its actors,
#'    the number of pairs of the two actors' affiliations that have the same positive \code{key}. With agencies as centres
#'    and the country of a club as its key this is the number of co-client pairs whose clubs are in the same country; with
#'    every club its own key and a subset of the agencies as centres it is stacking among those agencies, so that stacking
#'    can differ by agency size. \\code{centre} may also hold weights: every pair of a centre's actors then counts with
#'    the centre's weight, which lets stacking vary with the size of the centre, e.g. weight \\code{log(size / 10)}.}
#'  \item{\code{stack(ntypes = 2)}}{Number of pairs of actors that share their affiliation of every type (for two types: the same
#'    agency and the same club). Equals \code{b1dsp(2)} on the constrained space and is much faster, since it scans the
#'    smaller of the partner sets.}
#'  \item{\code{dstack(centre)}}{Number of (centre, club) cells with at least two actors: stacked cells counted once, however large.}
#'  \item{\code{homcentre(centre, attr)}}{Number of pairs of actors with the same value of \code{attr} that share a centre
#'    affiliation (nationality homophily within agencies when the centres are the agencies, within clubs when they are the
#'    clubs).}
#'  \item{\code{acov(W)}}{\code{sum_p W[a(p), c(p)]}: an affiliation-affiliation covariate (matrix with one row and column per
#'    mode-2 node, symmetric) summed over the actors that link the two affiliations (for example the share of an agency's
#'    outside clients in the club's country).}
#' }
#' Actor-affiliation covariates (an agency's nationality profile, value gaps) are available through \code{ergm}'s
#' \code{edgecov()}.
#' @param nw,arglist,... passed by \code{ergm}.
#' @name ergm.exclusive-terms
NULL

#' @rdname ergm.exclusive-terms
#' @export
InitErgmTerm.coloc <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = c("centre", "key", "label"), vartypes = c("numeric", "numeric", "character"),
                            defaultvalues = list(NULL, NULL, NULL), required = c(TRUE, TRUE, FALSE))
  d <- .bip_dims(nw)
  if (length(a$centre) != d$n2 || length(a$key) != d$n2) ergm::ergm_Init_stop("centre and key need one entry per mode-2 node (", d$n2, ").")
  if (!is.numeric(a$centre) || anyNA(a$centre)) ergm::ergm_Init_stop("centre must be numeric: 0 = not a centre, 1 = centre, or a centre weight.")
  if (any(a$key < 0 | a$key != round(a$key))) ergm::ergm_Init_stop("key must be a non-negative integer (0 = no key).")
  if (any(a$centre != 0 & a$key > 0)) ergm::ergm_Init_stop("A node cannot be both a centre and have a key.")
  list(name = "coloc", coef.names = paste0("coloc", if (is.null(a$label)) "" else paste0(".", a$label)),
       inputs = c(d$n1, a$centre, a$key, .key_groups(a$key, d$n1)), dependence = TRUE, minval = 0)
}

# Key groups for the coloc term: number of keys G, offsets (G + 1) and the vertex ids (mode-2 nodes numbered after the
# n1 actors) of the nodes of each key, in key order. Lets the C code scan the nodes of a key instead of a large centre.
.key_groups <- function(key, n1) {
  G <- max(c(0L, key))
  ord <- order(key, seq_along(key)); ord <- ord[key[ord] > 0]
  cnt <- tabulate(key[key > 0], nbins = max(G, 1))[seq_len(G)]
  c(G, c(0, cumsum(cnt)), n1 + ord)
}

#' @rdname ergm.exclusive-terms
#' @export
InitErgmTerm.stack <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = c("ntypes"), vartypes = c("numeric"), defaultvalues = list(2), required = FALSE)
  d <- .bip_dims(nw)
  list(name = "stack", coef.names = "stack", inputs = c(d$n1, a$ntypes), dependence = TRUE, minval = 0)
}

#' @rdname ergm.exclusive-terms
#' @export
InitErgmTerm.dstack <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = c("centre", "label"), vartypes = c("numeric", "character"),
                            defaultvalues = list(NULL, NULL), required = c(TRUE, FALSE))
  d <- .bip_dims(nw)
  if (length(a$centre) != d$n2 || !all(a$centre %in% 0:1)) ergm::ergm_Init_stop("centre must be a 0/1 vector with one entry per mode-2 node (", d$n2, ").")
  list(name = "dstack", coef.names = paste0("dstack", if (is.null(a$label)) "" else paste0(".", a$label)),
       inputs = c(d$n1, a$centre), dependence = TRUE, minval = 0)
}

#' @rdname ergm.exclusive-terms
#' @export
InitErgmTerm.homcentre <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = c("centre", "attr", "label"), vartypes = c("numeric", "numeric,character,factor", "character"),
                            defaultvalues = list(NULL, NULL, NULL), required = c(TRUE, TRUE, FALSE))
  d <- .bip_dims(nw)
  if (length(a$centre) != d$n2 || length(a$attr) != d$n1) ergm::ergm_Init_stop("centre needs one entry per mode-2 node (", d$n2, "), attr one per actor (", d$n1, ").")
  if (!all(a$centre %in% 0:1)) ergm::ergm_Init_stop("centre must be 0/1.")
  at <- as.integer(factor(a$attr))
  list(name = "homcentre", coef.names = paste0("homcentre", if (is.null(a$label)) "" else paste0(".", a$label)),
       inputs = c(d$n1, a$centre, at), dependence = TRUE, minval = 0)
}

#' @rdname ergm.exclusive-terms
#' @export
InitErgmTerm.acov <- function(nw, arglist, ...) {
  a <- ergm::check.ErgmTerm(nw, arglist, varnames = c("W", "label"), vartypes = c("matrix", "character"),
                            defaultvalues = list(NULL, NULL), required = c(TRUE, FALSE))
  d <- .bip_dims(nw)
  if (!all(dim(a$W) == d$n2)) ergm::ergm_Init_stop("W must be ", d$n2, " x ", d$n2, " (one row and column per mode-2 node).")
  if (max(abs(a$W - t(a$W))) > 1e-12) ergm::ergm_Init_stop("W must be symmetric.")
  list(name = "acov", coef.names = paste0("acov", if (is.null(a$label)) "" else paste0(".", a$label)),
       inputs = c(d$n1, as.vector(t(a$W))), dependence = TRUE)
}
