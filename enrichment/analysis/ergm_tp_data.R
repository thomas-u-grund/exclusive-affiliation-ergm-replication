# ============================================================
# Data for the three-layer (agency - player - club) ERGM in nwergm (tripartite branch)
#
# Two-mode network: mode 1 = players, mode 2 = agencies and clubs ("affiliations").
# Every player has exactly one agency and one club. Node order:
#   1 .. NP            players
#   NP+1 .. NP+NA      agencies
#   NP+NA+1 .. N       clubs
#
# Writes enrichment/analysis/output/tp/:
#   nodes.csv  one row per node: node, kind (1 player, 2 agency, 3 club), attributes
#   edges.csv  player node, affiliation node
#   cov_*.csv  sparse dyadic covariates (i, j, value; i < j, only nonzero), one file per covariate
#
# Covariates (agency profiles from clients OUTSIDE the network, as in ergm_data.R):
#   A-P  ap_corr   share of the agency's outside clients with the player's nationality
#   A-P  ap_val    |log10 player value - mean log10 value of agency's outside clients|
#   P-C  pc_home   player's nationality = club's country
#   P-C  pc_val    |log10 player value - club squad mean log10 value (leave the player out)|
#   A-C  ac_home   share of the agency's outside clients who are nationals of the club's country (tpstarcov)
#   A-C  ac_val    |agency outside value - club squad mean value| (tpstarcov)
# ============================================================
source("enrichment/analysis/ergm_data.R")
out <- "enrichment/analysis/output/tp"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

pp <- pl |>
  left_join(attr |> transmute(player_id, nat = nationality1, pos = position_group,
                              age = as.numeric(age_at_snapshot_years), mv = as.numeric(mv_at_snapshot_eur)),
            by = "player_id") |>
  arrange(country, club_id, agency, player_id)
NP <- nrow(pp); NA_ <- length(agencies); NC <- nrow(clubs); N <- NP + NA_ + NC
p_node <- seq_len(NP)
a_node <- NP + match(pp$agency, agencies)
c_node <- NP + NA_ + match(pp$club_id, clubs$club_id)
cat(sprintf("players %d, agencies %d, clubs %d, nodes %d, ties %d\n", NP, NA_, NC, N, 2 * NP))

lmv <- ifelse(pp$mv > 0, log10(pp$mv), NA_real_)
# club squad mean value leaving the player out (all squad players in the snapshot, with or without agent)
sq <- attr |>
  filter(market_id %in% first_tier, snapshot_club_id %in% clubs$club_id) |>
  mutate(l = ifelse(as.numeric(mv_at_snapshot_eur) > 0, log10(as.numeric(mv_at_snapshot_eur)), NA_real_)) |>
  group_by(club_id = snapshot_club_id) |>
  summarise(s = sum(l, na.rm = TRUE), k = sum(!is.na(l)), .groups = "drop")
clubs2 <- clubs |> left_join(sq, by = "club_id")

# ---------------------------------------------------------------- nodes
cty_code <- match(clubs$country, countries)
nat_levels <- sort(unique(na.omit(pp$nat)))
pos_levels <- sort(unique(na.omit(pp$pos)))
nodes <- tibble(
  node = seq_len(N),
  kind = c(rep(1L, NP), rep(2L, NA_), rep(3L, NC)),
  id = c(pp$player_id, agencies, clubs$club_id),
  label = c(pp$player_id, agencies, clubs$club_name),
  nat = c(match(pp$nat, nat_levels), rep(NA, NA_ + NC)),
  pos = c(match(pp$pos, pos_levels), rep(NA, NA_ + NC)),
  agegrp = c(as.integer(cut(pp$age, c(-Inf, 21, 25, 29, Inf), right = FALSE)), rep(NA, NA_ + NC)),
  lmv = c(lmv, ag$ag_mv, clubs$club_mv),
  club_country = c(match(pp$country, countries), rep(NA, NA_), cty_code),
  # affiliation type for bnodematch1(.., byattr(type)): 1 agency, 2 club (players: missing)
  type = c(rep(NA, NP), rep(1L, NA_), rep(2L, NC)),
  # tpcoloc key: 0 player, -1 agency (centre), club country code for clubs
  key_cty = c(rep(0L, NP), rep(-1L, NA_), cty_code),
  # tpstarcov class: 0 player, 1 affiliation
  cls = c(rep(0L, NP), rep(1L, NA_ + NC)),
  # bd(): alter groups and bounds
  isag = as.integer(kind == 2L), iscl = as.integer(kind == 3L),
  bmax = ifelse(kind == 1L, 1L, NP), bmin = ifelse(kind == 1L, 1L, 0L)
)
write_csv(nodes, file.path(out, "nodes.csv"), na = "")
write_csv(tibble(level = seq_along(nat_levels), nat = nat_levels), file.path(out, "nat_levels.csv"))
write_csv(tibble(level = seq_along(countries), country = countries), file.path(out, "country_levels.csv"))
write_csv(tibble(p = c(p_node, p_node), x = c(a_node, c_node)), file.path(out, "edges.csv"))

# ---------------------------------------------------------------- covariates
wcov <- function(name, i, j, v) {
  keep <- !is.na(v) & v != 0
  lo <- pmin(i, j); hi <- pmax(i, j)
  stopifnot(length(lo) == length(v), all(lo < hi))
  tibble(i = lo[keep], j = hi[keep], v = v[keep]) |>
    write_csv(file.path(out, paste0("cov_", name, ".csv")))
  cat(sprintf("cov_%s: %d nonzero entries\n", name, sum(keep)))
}
mean_or <- function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }

# A-P: all agency x player dyads
gp <- expand.grid(a = seq_len(NA_), p = seq_len(NP))
pn <- match(pp$nat, nats)
ap_corr <- P[cbind(gp$a, ifelse(is.na(pn[gp$p]), NA, pn[gp$p]))]
ap_corr[is.na(ap_corr)] <- 0
wcov("ap_corr", NP + gp$a, gp$p, ap_corr)
ap_val <- abs(ag$ag_mv[gp$a] - lmv[gp$p])
ap_val[ag$ag_mv_missing[gp$a]] <- NA
wcov("ap_val", NP + gp$a, gp$p, mean_or(ap_val))

# P-C: all player x club dyads
gc <- expand.grid(c = seq_len(NC), p = seq_len(NP))
wcov("pc_home", NP + NA_ + gc$c, gc$p, as.numeric(!is.na(pp$nat[gc$p]) & pp$nat[gc$p] == clubs$country[gc$c]))
own <- pp$club_id[gc$p] == clubs2$club_id[gc$c] & !is.na(lmv[gc$p])
loo <- (clubs2$s[gc$c] - ifelse(own, lmv[gc$p], 0)) / (clubs2$k[gc$c] - own)
wcov("pc_val", NP + NA_ + gc$c, gc$p, mean_or(abs(lmv[gc$p] - loo)))

# A-C: agency x club (cross-layer, enters through tpstarcov)
gac <- expand.grid(a = seq_len(NA_), c = seq_len(NC))
wcov("ac_home", NP + gac$a, NP + NA_ + gac$c, nat_aff[cbind(gac$a, gac$c)])
ac_val <- abs(ag$ag_mv[gac$a] - clubs$club_mv[gac$c])
ac_val[ag$ag_mv_missing[gac$a]] <- NA
wcov("ac_val", NP + gac$a, NP + NA_ + gac$c, mean_or(ac_val))

# ---------------------------------------------------------------- observed statistics (for checking the Stata side)
co_ag <- pp |> count(agency, name = "n") |> summarise(s = sum(choose(n, 2))) |> pull(s)
stack <- pp |> count(agency, club_id, name = "n") |> summarise(s = sum(choose(n, 2))) |> pull(s)
coloc <- pp |> count(agency, country, name = "n") |> summarise(s = sum(choose(n, 2))) |> pull(s)
conat_ag <- pp |> filter(!is.na(nat)) |> count(agency, nat, name = "n") |> summarise(s = sum(choose(n, 2))) |> pull(s)
conat_cl <- pp |> filter(!is.na(nat)) |> count(club_id, nat, name = "n") |> summarise(s = sum(choose(n, 2))) |> pull(s)
obs <- tibble(stat = c("co_client_pairs", "b1dsp2_stacking", "tpcoloc_country", "conat_agency", "conat_club",
                       "ap_corr", "ap_val", "pc_home", "pc_val", "ac_home", "ac_val"),
              value = c(co_ag, stack, coloc, conat_ag, conat_cl,
                        sum(P[cbind(match(pp$agency, agencies), ifelse(is.na(pn), NA, pn))], na.rm = TRUE),
                        NA, sum(!is.na(pp$nat) & pp$nat == pp$country), NA,
                        sum(nat_aff[cbind(match(pp$agency, agencies), match(pp$club_id, clubs$club_id))]), NA))
write_csv(obs, file.path(out, "observed_stats_R.csv"))
print(obs)
