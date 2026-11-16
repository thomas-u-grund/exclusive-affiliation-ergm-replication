# Figure: three patterns in the same football network by agency size -- nationality sorting, national-market concentration and
# same-club stacking -- against what random assignment would give. Run from the project root.
suppressMessages({library(tidyverse); library(patchwork); library(network)})
source("enrichment/analysis/ergm_excl_data.R")
pal <- c(blue = "#2a78d6", ink = "#1a1a19")
theme_set(theme_minimal(base_size = 10, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3), axis.text = element_text(colour = "grey25"), axis.title = element_text(colour = "grey25"),
  strip.text = element_text(colour = "grey15", hjust = 0), legend.position = "top", legend.title = element_blank(), plot.tag = element_text(face = "bold")))
B <- network::as.sociomatrix(fb$net); n1 <- fb$n1
ag <- apply(B[, fb$isag == 1], 1, which.max); cl <- apply(B[, fb$iscl == 1], 1, which.max)
cty <- fb$key_cty[fb$iscl == 1]; nat <- fb$nat; size <- tabulate(ag, sum(fb$isag))
band <- cut(size[ag], c(1, 3, 9, 29, Inf), c("2–3", "4–9", "10–29", "30+")); ok <- !is.na(band)   # players in agencies with >= 2 clients
stat <- function(ag, cl, nat) {   # per band: share of co-client pairs with the same nationality / same club country-of-modal-share / same club
  d <- tibble(ag = ag, cl = cl, nat = nat, cty = cty[cl], band = band)[ok, ]
  bind_rows(
    d |> group_by(band, ag) |> summarise(p = sum(choose(table(nat), 2)), n = choose(n(), 2), .groups = "drop") |> group_by(band) |> summarise(v = sum(p) / sum(n)) |> mutate(stat = "A"),
    d |> group_by(band, ag) |> summarise(m = max(table(cty)) / n(), .groups = "drop") |> group_by(band) |> summarise(v = mean(m)) |> mutate(stat = "B"),
    d |> group_by(band, ag) |> summarise(p = sum(choose(table(cl), 2)), n = choose(n(), 2), .groups = "drop") |> group_by(band) |> summarise(v = sum(p) / sum(n)) |> mutate(stat = "C"))
}
obs <- stat(ag, cl, nat); set.seed(1)
nul <- map_dfr(1:200, \(i) { pp <- sample(n1); stat(ag, cl, nat[pp]) |> filter(stat == "A") |> mutate(rep = i) }) |>   # nationality shuffled over players
  bind_rows(map_dfr(1:200, \(i) { pp <- sample(n1); stat(ag, cl[pp], nat) |> filter(stat != "A") |> mutate(rep = i) }))   # clubs shuffled over players
nl <- nul |> group_by(stat, band) |> summarise(m = mean(v), lo = quantile(v, .025), hi = quantile(v, .975), .groups = "drop")
lab <- c(A = "Nationality sorting:\nco-client pairs of the\nsame nationality", B = "Market concentration:\nclients in the agency's\nmain club country", C = "Stacking:\nco-client pairs at the\nsame club")
pl <- \(s) ggplot() + geom_linerange(data = filter(nl, stat == s), aes(band, ymin = lo, ymax = hi), colour = "grey60", linewidth = 3.5, alpha = .6) +
  geom_point(data = filter(obs, stat == s), aes(band, v), colour = pal["ink"], size = 2) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) + labs(x = "Clients of the agency in the network", y = NULL, subtitle = lab[[s]]) +
  theme(plot.subtitle = element_text(colour = "grey15", size = 8))
g <- (pl("A") | pl("B") | pl("C")) + plot_annotation(tag_levels = "A", caption = "Points = observed; grey bars = 95% range under random reassignment of nationalities (A) or of clubs (B, C) among players, 200 draws.", theme = theme(plot.caption = element_text(colour = "grey35", hjust = 0)))
ggsave("paper/figures/fig_patterns.pdf", g, width = 6.5, height = 3.2, device = cairo_pdf); print(obs |> pivot_wider(names_from = band, values_from = v), width = 120)
