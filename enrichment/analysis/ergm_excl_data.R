# Football network (agency - player - club) as a bipartite network for ergm.exclusive.
# Source: enrichment/analysis/output/tp/ (nodes.csv, edges.csv, cov_*.csv written by ergm_tp_data.R / ergm_tp_mkt.R).
# Usage: source("enrichment/analysis/ergm_excl_data.R") from the project root; creates `fb` (a list).
suppressMessages({library(tidyverse); library(network)})
.tp <- "enrichment/analysis/output/tp"
.nodes <- read_csv(file.path(.tp, "nodes.csv"), show_col_types = FALSE)
.edges <- read_csv(file.path(.tp, "edges.csv"), show_col_types = FALSE)
.n1 <- sum(.nodes$kind == 1); .n2 <- nrow(.nodes) - .n1
.B <- matrix(0, .n1, .n2); .B[cbind(.edges$p, .edges$x - .n1)] <- 1
.net <- network(.B, bipartite = .n1, directed = FALSE, matrix.type = "bipartite")
.aff <- .nodes[.nodes$kind != 1, ]
.cov <- function(name) read_csv(file.path(.tp, paste0("cov_", name, ".csv")), show_col_types = FALSE)
# actor x affiliation covariates (agency-player, player-club)
.apm <- function(name) { d <- .cov(name); r <- pmin(d$i, d$j); c <- pmax(d$i, d$j) - .n1; m <- matrix(0, .n1, .n2); m[cbind(r, c)] <- d$v; m }
# affiliation x affiliation covariates (agency-club), symmetric
.acm <- function(name) { d <- .cov(name); m <- matrix(0, .n2, .n2); m[cbind(d$i - .n1, d$j - .n1)] <- d$v; m[cbind(d$j - .n1, d$i - .n1)] <- d$v; m }
fb <- list(
  net = .net, n1 = .n1, n2 = .n2,
  types = ifelse(.aff$kind == 2, 1, 2), isag = as.numeric(.aff$kind == 2), iscl = as.numeric(.aff$kind == 3),
  key_cty = ifelse(.aff$kind == 3, .aff$key_cty, 0),               # club country (0 for agencies)
  small = as.numeric(.aff$key_clubs == -1), large = as.numeric(.aff$key_clubl == -1),   # agencies < 30 clients / >= 30
  club_id = ifelse(.aff$kind == 3, .aff$key_clubs, 0),             # every club its own key
  wsize = ifelse(.aff$kind == 2, log(pmax(colSums(.B), 1) / 10), 0),  # centred log portfolio size (agencies), centre = 10 clients
  nat = .nodes$nat[.nodes$kind == 1],
  ap_corr = .apm("ap_corr"), ap_val = .apm("ap_val"), pc_home = .apm("pc_home"), pc_val = .apm("pc_val"),
  ac_home = .acm("ac_home"), ac_val = .acm("ac_val"), ac_mkt = .acm("ac_mkt"))
# robustness check: keep market presence only for agencies with at least MKT_MIN outside clients (env MKT_MIN)
.mm <- as.numeric(Sys.getenv("MKT_MIN", "0"))
if (.mm > 0) { .on <- read_csv(file.path(.tp, "agency_outside_n.csv"), show_col_types = FALSE); .drop <- .on$node[.on$n_out < .mm] - fb$n1
  fb$ac_mkt[.drop, ] <- 0; fb$ac_mkt[, .drop] <- 0; cat("market presence dropped for", length(.drop), "agencies with <", .mm, "outside clients\n") }
rm(.B, .nodes, .edges, .aff)
cat("football network:", fb$n1, "players,", sum(fb$isag), "agencies,", sum(fb$iscl), "clubs,", network.edgecount(fb$net), "ties\n")
