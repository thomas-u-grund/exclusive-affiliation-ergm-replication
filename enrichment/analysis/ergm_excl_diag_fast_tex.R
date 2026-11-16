# Appendix convergence table from fast-sampler validations (fastval/val_<tag>.csv).
V <- "enrichment/analysis/output/tp/fastval/"
sp <- list("Core model" = "1F2", "Model 1" = "3mF2", "Model 2" = "3msF2", "Model 3" = "3s3F")
lab <- c(stack = "Stacking, all agencies", coloc.small = "Stacking, agencies $<30$ clients", coloc.under10 = "Stacking, agencies $<10$ clients", coloc.10to29 = "Stacking, agencies 10--29 clients",
  coloc.large = "Stacking, agencies $\\ge30$ clients", homcentre.agency = "Homophily in agencies", homcentre.club = "Homophily in clubs", coloc.country = "Co-location by country",
  acov.home = "Corridor", acov.val = "Agency--club value gap", acov.mkt = "Market presence", edgecov.ap_corr = "Agency--player profile", edgecov.ap_val = "Agency--player value gap",
  edgecov.pc_home = "Player--club home", edgecov.pc_val = "Player--club value gap")
f0 <- function(v) sapply(v, function(x) formatC(x, format = "f", digits = if (abs(x) >= 100) 0 else if (abs(x) >= 10) 1 else 2, big.mark = ","))
hdr <- c("\\begin{tabular}{l rrrr r}\\toprule", "Statistic & Observed & Simulated mean & (Obs$-$sim)/SD & Range of chain means & ESS\\\\\\midrule")
blk <- function(m) { d <- read.csv(paste0(V, "val_", sp[[m]], ".csv"), row.names = 1)
  c(sprintf("\\multicolumn{6}{l}{\\emph{%s}}\\\\", m), sprintf("%s & %s & %s & %.2f & %s--%s & %s\\\\", lab[rownames(d)], f0(d$obs), f0(d$sim), d$dsd, f0(d$cmin), f0(d$cmax), formatC(d$ess, format = "f", digits = 0))) }
wr <- function(ms, f) writeLines(c(hdr, unlist(lapply(seq_along(ms), function(k) c(if (k > 1) "\\addlinespace", blk(ms[k])))), "\\bottomrule", "\\end{tabular}"), f)
wr(c("Core model", "Model 1"), "paper/tables/diag_a.tex"); wr(c("Model 2", "Model 3"), "paper/tables/diag_b.tex")
