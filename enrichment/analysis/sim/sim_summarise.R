# ============================================================
# Simulation study, step 5: summarise.
#   recovery.csv   joint and conditional model: mean estimate, bias, RMSE, mean SE / SD, 95% coverage, rejection
#                  rate at 5% (power where the truth is nonzero, false positives where it is zero)
#   projection.csv stacking permutation tests and the binary projection ERGM: rejection rates per scenario
# Run from the project root: Rscript enrichment/analysis/sim/sim_summarise.R
# ============================================================
suppressMessages(library(tidyverse))
sim <- "enrichment/analysis/sim"
truth <- tribble(
  ~scen, ~S,   ~HA,  ~HC, ~L,   ~TW, ~HOME,
  1,     0,    0.06, 0.1, 0,    0,   1.0,
  2,     0.08, 0.06, 0.1, 0,    0,   1.0,
  3,     0,    0,    0.1, 0.05, 0,   1.0,
  4,     0.08, 0.05, 0.1, 0.02, 0.5, 1.0) |>
  pivot_longer(-scen, names_to = "param", values_to = "true")

read_all <- function(pattern) {
  f <- list.files(file.path(sim, "results"), pattern = pattern, full.names = TRUE)
  map_dfr(f, ~ read_csv(.x, show_col_types = FALSE))
}
joint_map <- c(b1dsp2 = "S", b1nodematch_nat_1 = "HA", b1nodematch_nat_2 = "HC", tpcoloc_key_cty = "L",
               tpstarcov_ac_home = "TW", edgecov_pc_home = "HOME")
cond_map <- c(b1nodematch_club = "S", b1nodematch_nat = "HA", b1nodematch_cty = "L", edgecov_w_ap = "TW")

long <- function(d, prefix, map, model) {
  b <- d |> select(scen, rep, starts_with(paste0(prefix, "b_"))) |>
    pivot_longer(-c(scen, rep), names_to = "term", values_to = "est") |>
    mutate(term = sub(paste0(prefix, "b_"), "", term))
  se <- d |> select(scen, rep, starts_with(paste0(prefix, "se_"))) |>
    pivot_longer(-c(scen, rep), names_to = "term", values_to = "se") |>
    mutate(term = sub(paste0(prefix, "se_"), "", term))
  inner_join(b, se, by = c("scen", "rep", "term")) |>
    filter(term %in% names(map)) |> mutate(param = map[term], model = model)
}

joint <- read_all("^joint_s") |> filter(fitted == 1)
rres <- read_all("^r_s")
est <- bind_rows(long(joint, "", joint_map, "joint"), long(rres, "cond_", cond_map, "conditional")) |>
  left_join(truth, by = c("scen", "param"))

recovery <- est |>
  group_by(model, scen, param, true) |>
  summarise(n = n(), mean = mean(est), bias = mean(est - true), rmse = sqrt(mean((est - true)^2)),
            sd = sd(est), mean_se = mean(se), se_ratio = mean(se) / sd(est),
            coverage = mean(abs(est - true) <= 1.96 * se),
            reject5 = mean(abs(est / se) > 1.96), .groups = "drop") |>
  arrange(scen, param, model)
write_csv(recovery, file.path(sim, "results/recovery.csv"))
print(recovery, n = Inf, width = 200)

projection <- rres |>
  group_by(scen) |>
  summarise(n = n(),
            Q_plain = mean(Q_plain), reject_plain = mean(p_plain < 0.05),
            Q_strat = mean(Q_strat), reject_strat = mean(p_strat < 0.05),
            proj_country = mean(proj_b_b2nodematch_country, na.rm = TRUE),
            proj_country_reject = mean(abs(proj_b_b2nodematch_country / proj_se_b2nodematch_country) > 1.96, na.rm = TRUE),
            proj_W = mean(proj_b_edgecov_W_k, na.rm = TRUE),
            proj_W_reject = mean(abs(proj_b_edgecov_W_k / proj_se_edgecov_W_k) > 1.96, na.rm = TRUE))
write_csv(projection, file.path(sim, "results/projection.csv"))
print(projection, width = 200)
cat(sprintf("joint fits: %d of %d converged; median %.0f s\n",
            sum(joint$converged == 1, na.rm = TRUE), nrow(joint), median(joint$seconds)))
