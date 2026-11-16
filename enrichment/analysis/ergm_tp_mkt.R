# ============================================================
# Agency market presence covariate for the three-layer ERGM (added 2026-11-13).
#   A-C  ac_mkt  share of the agency's outside clients (not in the network, not at a network club) whose club is in
#                the club's country (any tier). Club countries from enrichment/output/club_country.csv
#                (club_country.py: propagation over domestic league games); clubs without a country count as abroad.
# Motivation: a single co-location parameter L treats all agencies alike; the data show small agencies concentrated in
# one market and large ones spread out. ac_mkt measures each agency's market from clients outside the network, as
# ac_home does for nationality, so it cannot absorb stacking at network clubs.
# Writes enrichment/analysis/output/tp/cov_ac_mkt.csv (same format and node order as ergm_tp_data.R).
# Run from the project root: Rscript enrichment/analysis/ergm_tp_mkt.R
# ============================================================
source("enrichment/analysis/ergm_data.R")
out <- "enrichment/analysis/output/tp"
NP <- nrow(pl); NA_ <- length(agencies); NC <- nrow(clubs)
stopifnot(NP + NA_ + NC == nrow(read_csv(file.path(out, "nodes.csv"), show_col_types = FALSE)))

cc <- read_csv("enrichment/output/club_country.csv", col_types = "ccic") |>
  filter(!is.na(country), country != "") |>
  mutate(country = stringr::str_to_title(country))
cat("club countries:", nrow(cc), "; network countries:", paste(sort(unique(clubs$country)), collapse = ", "), "\n")
stopifnot(all(unique(clubs$country) %in% unique(cc$country)))

oc <- out_clients |> left_join(cc |> select(club_id, ccountry = country), by = "club_id")
cat(sprintf("outside client rows %d; with a known club country among the network countries %d (%.1f%%)\n",
            nrow(oc), sum(!is.na(oc$ccountry)), 100 * mean(!is.na(oc$ccountry))))
mkt <- oc |>
  crossing(country = sort(unique(clubs$country))) |>
  group_by(source_agent_id, country) |>
  summarise(share = mean(coalesce(ccountry == country, FALSE)), .groups = "drop") |>
  pivot_wider(names_from = country, values_from = share) |>
  right_join(tibble(source_agent_id = agencies), by = "source_agent_id") |>
  arrange(match(source_agent_id, agencies))
M <- as.matrix(mkt[, sort(unique(clubs$country))]); M[is.na(M)] <- 0
mkt_ac <- M[, match(clubs$country, colnames(M))]
gac <- expand.grid(a = seq_len(NA_), c = seq_len(NC))
v <- mkt_ac[cbind(gac$a, gac$c)]
i <- NP + gac$a; j <- NP + NA_ + gac$c; keep <- v != 0
tibble(i = i[keep], j = j[keep], v = v[keep]) |> write_csv(file.path(out, "cov_ac_mkt.csv"))
cat(sprintf("cov_ac_mkt: %d nonzero entries; mean share at the agency's own clubs' countries %.3f\n",
            sum(keep), mean(mkt_ac[cbind(match(pl$agency, agencies), match(pl$club_id, clubs$club_id))])))
# by agency size: how concentrated is each agency's outside market (max share over the 7 countries)?
asize <- tabulate(match(pl$agency, agencies), NA_)
print(tibble(band = cut(asize, c(0, 1, 3, 9, 29, Inf), labels = c("1", "2-3", "4-9", "10-29", "30+")),
             maxshare = apply(M, 1, max), n_out = ag$n_out) |>
        group_by(band) |> summarise(agencies = n(), mean_max_share = mean(maxshare), median_out = median(n_out)))
