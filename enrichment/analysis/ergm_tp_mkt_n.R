# Number of outside clients per agency (node order of the tripartite network) and the share with a known club country, for the
# market-presence robustness check. Writes output/tp/agency_outside_n.csv (col node, n_out, n_known).
source("enrichment/analysis/ergm_data.R")
out <- "enrichment/analysis/output/tp"; NP <- nrow(pl)
cc <- read_csv("enrichment/output/club_country.csv", col_types = "ccic") |> filter(!is.na(country), country != "")
oc <- out_clients |> left_join(cc |> select(club_id, ccountry = country), by = "club_id")
n <- oc |> group_by(source_agent_id) |> summarise(n_out = n(), n_known = sum(!is.na(ccountry)), .groups = "drop")
d <- tibble(node = NP + seq_along(agencies), source_agent_id = agencies) |> left_join(n, by = "source_agent_id") |> mutate(across(c(n_out, n_known), \(x) coalesce(x, 0L)))
write_csv(d |> select(node, n_out, n_known), file.path(out, "agency_outside_n.csv"))
print(summary(d$n_out)); cat("agencies with >=5 outside clients:", sum(d$n_out >= 5), " >=10:", sum(d$n_out >= 10), " small (<30 clients) among them: see model\n")
