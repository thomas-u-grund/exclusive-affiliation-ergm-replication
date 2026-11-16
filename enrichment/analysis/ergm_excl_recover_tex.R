suppressMessages(library(dplyr))
d <- do.call(rbind, lapply(list.files(Sys.getenv("RECDIR_IN", "enrichment/analysis/output/tp/recover"), "rec_.*csv", full.names = TRUE), read.csv))
lab <- c(coloc.under10 = "Stacking, agencies $<10$ clients", coloc.10to29 = "Stacking, agencies 10--29 clients", coloc.small = "Stacking, agencies $<30$ clients", coloc.large = "Stacking, agencies $\\ge30$ clients", homcentre.agency = "Nationality homophily, agencies",
  homcentre.club = "Nationality homophily, clubs", coloc.country = "Co-location by country", acov.home = "Corridor", acov.val = "Agency--club value gap",
  acov.mkt = "Market presence", edgecov.ap_corr = "Agency--player profile", edgecov.ap_val = "Agency--player value gap",
  edgecov.pc_home = "Player--club home", edgecov.pc_val = "Player--club value gap")
t <- d %>% group_by(name) %>% summarise(truth = first(truth), bias = mean(est - truth), sd = sd(est), se = mean(se),
  cover = mean(abs(est - truth) < 1.96 * se), ok = mean(abs(zmoment) < 3), n = n()) %>% arrange(match(name, names(lab)))
print(as.data.frame(t), digits = 3)
f <- function(x, k = 3) formatC(x, format = "f", digits = k)
L <- sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", lab[t$name], f(t$truth), f(t$bias), f(t$sd), f(t$se), f(t$cover, 2), f(t$ok, 2))
writeLines(c("\\begin{tabular}{lrrrrrr}\\hline", "Parameter & Truth & Bias & SD & Mean SE & Coverage & Share $|z|<3$ \\\\ \\hline", L, "\\hline\\end{tabular}"), Sys.getenv("TEXOUT", "paper/tables/recover.tex"))
