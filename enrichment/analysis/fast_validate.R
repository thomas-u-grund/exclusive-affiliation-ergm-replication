# Validation of a fitted model with the fast sampler: 4 independent chains x NS draws at the estimate; moment match (d_sd), chain spread, ESS,
# SEs from the pooled covariance, and size-band co-client pair counts (same club S, same country L). Env MODEL, THFILE, TAG, NS (200), MKT_MIN.
suppressMessages({library(ergm); library(ergm.exclusive); library(coda)}); source("enrichment/analysis/ergm_excl_data.R"); source("enrichment/analysis/ergm_excl_terms.R"); source("enrichment/analysis/ergm_excl_fast.R")
model <- Sys.getenv("MODEL"); tag <- Sys.getenv("TAG", model); ns <- as.integer(Sys.getenv("NS", "200")); th <- read.csv(Sys.getenv("THFILE"))$theta
f <- excl_formula(model, fb); K <- length(th); stopifnot(length(summary(f)) == K); tm <- fast_terms(f); obs <- summary(f)
sz <- colSums(network::as.sociomatrix(fb$net)); band <- function(lo, hi) as.numeric(fb$isag == 1 & sz >= lo & sz <= hi); bn <- c("2-3", "4-9", "10-29", "30+")
bands <- list(band(2, 3), band(4, 9), band(10, 29), band(30, Inf))
mon <- c(lapply(bands, \(b) list(type = "coloc", centre = b, key = fb$club_id, label = "S")), lapply(bands, \(b) list(type = "coloc", centre = b, key = fb$key_cty, label = "L")))
# monitor terms with theta 0 are appended; they do not change the distribution
ch <- lapply(1:4, function(s) ex_fast_simulate(fb$net, fb$types, fb$nat, c(tm, mon), c(th, rep(0, 8)), nsim = ns, burnin = 1e7, interval = 1e6, seed = 500 + s)$stats)
M <- lapply(ch, \(x) x[, 1:K, drop = FALSE]); S <- do.call(rbind, M); nm <- names(obs); cm <- sapply(M, colMeans)
ess <- sapply(seq_len(K), \(j) sum(sapply(M, \(x) effectiveSize(x[, j]))))
r <- data.frame(obs = as.numeric(obs), sim = colMeans(S), dsd = (as.numeric(obs) - colMeans(S)) / apply(S, 2, sd), cmin = apply(cm, 1, min), cmax = apply(cm, 1, max), ess = ess,
                theta = th, se = sqrt(diag(solve(cov(S))))); rownames(r) <- nm; print(round(r, 3))
G <- do.call(rbind, lapply(ch, \(x) x[, K + 1:8])); colnames(G) <- c(paste0("S_", bn), paste0("L_", bn))
gb <- rbind(observed = summary(as.formula(paste("fb$net ~", paste(sprintf("coloc(bands[[%d]], fb$club_id)", 1:4), collapse = "+"), "+", paste(sprintf("coloc(bands[[%d]], fb$key_cty)", 1:4), collapse = "+")))), mean = colMeans(G), q025 = apply(G, 2, quantile, .025), q975 = apply(G, 2, quantile, .975))
colnames(gb) <- colnames(G); print(round(gb, 1))
cat("max |d_sd|:", round(max(abs(r$dsd)), 2), " correlation of stacking statistics:\n"); st <- grep("^coloc\\.(under10|10to29|large|small|zsize|hinge|size)|^stack", nm); if (length(st) > 1) print(round(cor(S[, st]), 2))
dir.create("enrichment/analysis/output/tp/fastval", showWarnings = FALSE)
write.csv(r, sprintf("enrichment/analysis/output/tp/fastval/val_%s.csv", tag)); write.csv(gb, sprintf("enrichment/analysis/output/tp/fastval/band_%s.csv", tag))
write.csv(data.frame(sim = 0:nrow(G), rbind(gb["observed", ], G)), sprintf("enrichment/analysis/output/tp/fastval/gof_%s.csv", tag), row.names = FALSE)
