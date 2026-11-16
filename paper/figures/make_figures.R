# Figures for the paper. Run from the project root: Rscript paper/figures/make_figures.R
# fig_recovery.pdf  simulation recovery (scenarios 1-4) and the stacking trap (scenarios 5-6)
# fig_gof_size.pdf  goodness of fit by agency size (real network vs networks simulated from fitted models)
# Palette (validated with the dataviz validate_palette script, light surface): blue, orange, aqua; observed = near-black.
# Shapes differ by series so the figures survive grayscale printing; text uses text colours, never series colours.
suppressMessages({library(tidyverse); library(patchwork)})
pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", ink = "#1a1a19", violet = "#7a4fd1")
theme_set(theme_minimal(base_size = 8.5, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        axis.text = element_text(colour = "grey25"), axis.title = element_text(colour = "grey25"),
        strip.text = element_text(colour = "grey15", face = "plain", hjust = 0),
        legend.position = "top", legend.title = element_blank(), plot.tag = element_text(face = "bold")))
res <- "enrichment/analysis/sim/results"; tp <- "enrichment/analysis/output/tp"; out <- "paper/figures"

# ---- Figure: recovery, scenarios 1-4 ------------------------------------------------------------------------
rec <- read_csv(file.path(res, "recovery.csv"), show_col_types = FALSE) |> filter(model == "joint") |>
  mutate(param = factor(param, c("S", "L", "TW", "HA", "HC", "HOME"),
                        c("Stacking", "Co-location by country", "Corridor", "Nationality homophily, agencies",
                          "Nationality homophily, clubs", "Home country"))) |>
  mutate(lo = mean - 1.96 * sd, hi = mean + 1.96 * sd, scen = factor(scen, 1:4, c("1", "2", "3", "4")))
pa <- ggplot(rec, aes(scen)) +
  geom_errorbar(aes(ymin = true, ymax = true), width = 0.55, colour = pal["ink"], linewidth = 0.5) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = pal["blue"], linewidth = 0.6) +
  geom_point(aes(y = mean), colour = pal["blue"], size = 1.8) +
  facet_wrap(~param, scales = "free_y", nrow = 2) +
  labs(x = "Scenario", y = "Estimate: mean \u00b1 1.96 SD over 100 fits",
       subtitle = "Black bar = generating value") + theme(plot.subtitle = element_text(colour = "grey35"))

# ---- Figure: the stacking trap, scenarios 5 and 6 ---------------------------------------------------------------
rd <- function(pat, sc, mod) list.files(res, pat, full.names = TRUE) |> map_dfr(\(f) read_csv(f, show_col_types = FALSE)) |>
  distinct(rep, .keep_all = TRUE) |> mutate(scenario = sc, model = mod)
d <- bind_rows(rd("^joint5_simple_", "5: all agencies stack (0.08)", "Homogeneous agencies"),
               rd("^joint5_full_",   "5: all agencies stack (0.08)", "Agency terms, one stacking parameter"),
               rd("^joint6_simple_", "6: only small agencies stack (0.12)", "Homogeneous agencies"),
               rd("^joint6_full_",   "6: only small agencies stack (0.12)", "Agency terms, one stacking parameter"),
               rd("^joint6_size_",   "6: only small agencies stack (0.12)", "Agency terms + size-dependent stacking")) |>
  mutate(small = b_tpcoloc_key_clubs,
         model = factor(model, c("Homogeneous agencies", "Agency terms, one stacking parameter", "Agency terms + size-dependent stacking")))
long <- bind_rows(d |> transmute(scenario, model, rep, what = "Common stacking parameter", est = b_b1dsp2),
                  d |> filter(!is.na(small)) |> transmute(scenario, model, rep, what = "Extra stacking, small agencies", est = small))
truth <- tribble(~scenario, ~what, ~t,
  "5: all agencies stack (0.08)", "Common stacking parameter", 0.08,
  "6: only small agencies stack (0.12)", "Common stacking parameter", 0,
  "6: only small agencies stack (0.12)", "Extra stacking, small agencies", 0.12)
mods <- levels(d$model); cols <- setNames(pal[c("blue", "orange", "aqua")], mods); shp <- setNames(c(16, 17, 15), mods)
sm <- long |> group_by(scenario, what, model) |> summarise(m = mean(est), lo = quantile(est, .025), hi = quantile(est, .975), .groups = "drop")
sm <- sm |> mutate(panel = factor(paste0(scenario, "\n", what), unique(paste0(scenario, "\n", what))))
truth <- truth |> mutate(panel = factor(paste0(scenario, "\n", what), levels(sm$panel)))
pb <- ggplot(sm, aes(model, m, colour = model, shape = model)) +
  geom_hline(data = truth, aes(yintercept = t), colour = pal["ink"], linewidth = 0.4, linetype = "22") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = cols) + scale_shape_manual(values = shp) +
  scale_x_discrete(labels = NULL) + labs(x = NULL, y = "Estimate: mean and central 95% of fits",
       subtitle = "Dashed line = generating value") + theme(plot.subtitle = element_text(colour = "grey35"),
       legend.direction = "vertical", axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(out, "fig_recovery.pdf"), (pa / pb) + plot_layout(heights = c(1.4, 1)) + plot_annotation(tag_levels = "A"),
       width = 6.5, height = 7.6, device = cairo_pdf)

# ---- Figure: goodness of fit by agency size -----------------------------------------------------------------------
fits <- c("3mb" = "Model 1", "3ms22b" = "Model 2", "3s3" = "Model 3")
g <- imap_dfr(fits, \(lab, f) { p <- file.path(tp, "gof", paste0("gof_", f, ".csv")); if (file.exists(p)) read_csv(p, show_col_types = FALSE) |> mutate(model = lab) })
bands <- c("2-3", "4-9", "10-29", "30+")
long <- g |> filter(sim > 0) |> pivot_longer(c(starts_with("S_"), starts_with("L_")), names_to = "k", values_to = "v") |>
  separate(k, c("stat", "band"), sep = "_", extra = "merge") |> filter(band %in% bands)
obs <- g |> filter(sim == 0) |> slice(1) |> pivot_longer(c(starts_with("S_"), starts_with("L_")), names_to = "k", values_to = "obs") |>
  separate(k, c("stat", "band"), sep = "_", extra = "merge") |> filter(band %in% bands)
lab_stat <- c(S = "Co-client pairs at the same club", L = "Co-client pairs at clubs in the same country")
lab_band <- c("2-3" = "2–3 clients", "4-9" = "4–9 clients", "10-29" = "10–29 clients", "30+" = "30+ clients")
sg <- long |> group_by(stat, band, model) |> summarise(m = mean(v), lo = quantile(v, .025), hi = quantile(v, .975), .groups = "drop") |>
  left_join(obs |> select(stat, band, obs), by = c("stat", "band")) |>
  mutate(stat = factor(stat, names(lab_stat), lab_stat), band = factor(band, bands, lab_band[bands]),
         model = factor(model, fits[fits %in% unique(g$model)]))
mc <- setNames(pal[c("blue", "orange", "aqua", "violet")][seq_along(levels(sg$model))], levels(sg$model))
ms <- setNames(c(16, 17, 15, 18)[seq_along(levels(sg$model))], levels(sg$model))
one <- function(st, lab) ggplot(filter(sg, stat == st), aes(model, m, colour = model, shape = model)) +
  geom_hline(aes(yintercept = obs), colour = pal["ink"], linewidth = 0.45) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~band, nrow = 1, scales = "free_y") + scale_colour_manual(values = mc) + scale_shape_manual(values = ms) +
  scale_x_discrete(labels = NULL) + scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA), expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = NULL, y = NULL, subtitle = lab) + theme(plot.subtitle = element_text(colour = "grey15", size = 9),
       axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
pg <- (one(lab_stat[["S"]], "Co-client pairs at the same club") / one(lab_stat[["L"]], "Co-client pairs at clubs in the same country")) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
pg <- pg + plot_annotation(caption = "Black line = observed network; points and ranges = mean and 95% range of 800 networks simulated from each fitted model.",
                           theme = theme(plot.caption = element_text(colour = "grey35", hjust = 0)))
ggsave(file.path(out, "fig_gof_size.pdf"), pg, width = 6.5, height = 4.2, device = cairo_pdf)
cat("models in GOF figure:", levels(sg$model), "\n")
