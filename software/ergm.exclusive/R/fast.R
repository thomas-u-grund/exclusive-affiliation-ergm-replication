#' Dedicated fast sampler for exclusive two-affiliation networks
#'
#' A compiled sampler for the typed-swap chain that keeps count tables and updates every supported statistic in O(1) per term,
#' instead of going through the generic \code{ergm} machinery. Supported statistics: \code{stack}, \code{coloc}, \code{homcentre},
#' \code{acov} and \code{edgecov}-type actor-by-affiliation covariates. The proposal mixes the plain typed swap with a swap
#' between players that share a partner of the other type; both are symmetric, so the stationary distribution is the ERGM.
#'
#' @param net a bipartite \code{network} (mode 1 = players, mode 2 = affiliations) with one agency and one club per player.
#' @param types integer vector, type (1 = agency, 2 = club) of every mode-2 node.
#' @param nat integer or factor, attribute of the players used by all \code{homcentre} terms.
#' @param terms list of term specifications, in the order of the statistics: \code{list(type = "stack")},
#'   \code{list(type = "coloc", centre = <numeric, one per mode-2 node>, key = <integer, one per mode-2 node, 0 = none>)},
#'   \code{list(type = "homcentre", centre = <0/1 per mode-2 node>)}, \code{list(type = "acov", W = <n2 x n2 matrix>)},
#'   \code{list(type = "edgecov", M = <n1 x n2 matrix>)}.
#' @param theta coefficients, one per term.
#' @param nsim,burnin,interval number of draws, burn-in swaps and swaps between draws.
#' @param seed integer seed.
#' @return list with \code{stats} (nsim x number of terms), \code{alpha}, \code{gamma} (final agency and club of every player,
#'   mode-2 node indices) and \code{accept} (acceptance rate).
#' @export
ex_fast_simulate <- function(net, types, nat, terms, theta, nsim = 1, burnin = 1e6, interval = 1e5, seed = 1) {
  bip <- network::get.network.attribute(net, "bipartite"); n2 <- network::network.size(net) - bip
  B <- network::as.sociomatrix(net); stopifnot(nrow(B) == bip, ncol(B) == n2, length(types) == n2, length(theta) == length(terms))
  ag <- which(types == 1); cl <- which(types == 2); A <- length(ag); C <- length(cl)
  alpha <- as.integer(apply(B[, ag, drop = FALSE], 1, which.max) - 1L); gamma <- as.integer(apply(B[, cl, drop = FALSE], 1, which.max) - 1L)
  stopifnot(all(rowSums(B[, ag, drop = FALSE]) == 1), all(rowSums(B[, cl, drop = FALSE]) == 1))
  nat0 <- as.integer(factor(nat)) - 1L; V <- max(nat0) + 1L
  stack_col <- -1L; cols <- list(); homs <- list(); acv <- list(); edg <- list()
  for (j in seq_along(terms)) {
    tm <- terms[[j]]
    switch(tm$type,
      stack = { stack_col <- j - 1L },
      coloc = { key <- as.integer(factor(tm$key[cl], levels = sort(unique(tm$key[cl])))); kk <- tm$key[cl]; u <- sort(unique(kk[kk > 0]))
                kl <- match(kk, u); kl[is.na(kl)] <- 0L
                if (any(tm$key[ag] != 0)) stop("coloc: keys must be zero on agencies in the fast sampler")
                cols[[length(cols) + 1]] <- list(j - 1L, as.numeric(tm$centre[ag]), as.integer(kl), length(u)) },
      homcentre = { side <- if (all(tm$centre[cl] == 0)) 1L else if (all(tm$centre[ag] == 0)) 2L else stop("homcentre: centres must be agencies or clubs")
                    homs[[length(homs) + 1]] <- list(j - 1L, side, as.numeric(if (side == 1L) tm$centre[ag] else tm$centre[cl]), V) },
      acov = { acv[[length(acv) + 1]] <- list(j - 1L, tm$W) },
      edgecov = { edg[[length(edg) + 1]] <- list(j - 1L, tm$M) },
      stop("unknown term type ", tm$type))
  }
  inter <- function(lst, dims) {   # interleave matrices: element (i, k) of matrix j at ((i-1)*ncol + k-1)*m + j
    m <- length(lst); if (m == 0) return(list(integer(0), numeric(0)))
    arr <- array(0, c(m, dims[2], dims[1])); for (j in seq_len(m)) arr[j, , ] <- t(lst[[j]][[2]])
    list(vapply(lst, `[[`, 1L, 1L), as.vector(arr))
  }
  ac <- inter(acv, c(n2, n2)); ed <- inter(edg, c(bip, n2))
  spec <- list(as.integer(bip), as.integer(A), as.integer(C), length(terms), alpha, gamma, nat0, as.integer(ag - 1L), as.integer(cl - 1L),
               as.integer(stack_col), cols, homs, if (length(acv)) ac else list(), if (length(edg)) ed else list())
  # indices of club/agency in the agency/club-local numbering are positions within ag / cl
  r <- .Call("ex_fast_c", spec, as.numeric(theta), as.integer(nsim), as.numeric(burnin), as.numeric(interval), as.numeric(seed), PACKAGE = "ergm.exclusive")
  stats <- r[[1]]; colnames(stats) <- vapply(terms, function(t) if (is.null(t$label)) t$type else t$label, "")
  list(stats = stats, alpha = ag[r[[2]] + 1L], gamma = cl[r[[3]] + 1L], accept = r[[4]])
}
