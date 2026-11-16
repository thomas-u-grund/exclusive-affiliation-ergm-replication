# Fast sampler interface for the football models: converts an ergm formula from excl_formula() into the term list of
# ergm.exclusive::ex_fast_simulate (stack, coloc, homcentre, acov, edgecov), so both samplers use the same model specification.
fast_terms <- function(f) {
  e <- environment(f); rhs <- f[[3]]; calls <- list()
  while (is.call(rhs) && identical(rhs[[1]], as.name("+"))) { calls <- c(list(rhs[[3]]), calls); rhs <- rhs[[2]] }; calls <- c(list(rhs), calls)
  lapply(calls, function(cl) { nm <- as.character(cl[[1]]); a <- as.list(cl)[-1]; ev <- function(i) eval(a[[i]], e)
    switch(nm, stack = list(type = "stack", label = "stack"),
      coloc = list(type = "coloc", centre = ev(1), key = ev(2), label = paste0("coloc.", if (length(a) > 2) ev(3) else "")),
      homcentre = list(type = "homcentre", centre = ev(1), label = paste0("homcentre.", if (length(a) > 2) ev(3) else "")),
      acov = list(type = "acov", W = ev(1), label = paste0("acov.", if (length(a) > 1) ev(2) else "")),
      edgecov = list(type = "edgecov", M = ev(1), label = paste0("edgecov.", as.character(a[[1]]))), stop("term ", nm)) })
}
fast_sim <- function(f, fb, theta, ...) {
  tm <- fast_terms(f); stopifnot(length(tm) == length(theta))
  ex_fast_simulate(fb$net, fb$types, fb$nat, tm, theta, ...)
}
