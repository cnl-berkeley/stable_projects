###############################################################################
# centroid/centroid_distance_stats.R
#
# Distance from each participant's maximally face-selective PMC patch to the
# precuneal sulci (Figure 4), the positional control (prculs-d), and the
# geometry of the face-selective patch.
#
# Reproduces
#   Methods, "Centroid analysis...": number of connected components of the
#     top-5% face-selective vertices (median, range), mass fraction carried by
#     the highest-mass component, and the percentage of PMC it covers
#   Results, "The prcus-p is the nearest precuneal sulcus...":
#     - mean distances to prcus-p, prcus-i, prcus-a per hemisphere (5%)
#     - linear contrast prcus-p -> prcus-a (estimate, 95% CI, t, p), from a
#       linear mixed-effects model (nlme::lme, distance ~ sulcus * hemisphere
#       * group, random intercept per participant)
#     - Holm-corrected pairwise steps within hemisphere
#     - number of hemispheres in which the prcus-p is closer than the prcus-a
#     - group main effect and sulcus-by-group interaction at 5%, 2.5% and 10%
#     - prculs-d posterior to prcus-p (paired t on A-P position)
#     - prculs-d farther than prcus-p from the patch (paired t, 95% CI, count)
#       at 5%, 2.5% and 10%, within each group, and the group difference
#       (Welch's t) at each threshold
#   Figure 4A-D
#
# Inputs (data/)
#   centroid_distance.csv    geodesic distance (mm) on the native surface from
#                            the area-weighted centroid of the face patch to the
#                            area-weighted centroid of each sulcus
#   face_patch_geometry.csv  connected-component description of the top-k%
#                            face-selective PMC vertices
#   sulcal_ap_position.csv   mean anterior-posterior position of the prcus-p
#                            and prculs-d labels (0 = pos, 1 = mcgs)
#
# Outputs (outputs/centroid/)
#   centroid_distance_stats.txt   every statistic listed above
#   centroid_descriptives.csv     mean/SD distance per threshold, hemisphere, sulcus
#   prculsd_control.csv           prculs-d minus prcus-p tests, all thresholds
#   figures/fig4_centroid_distance.{png,pdf}
#
# The face patch is the highest-mass connected component of the top 5% (2.5%,
# 10%) of face-selective (faces - objects) vertices within the PMC territory.
# The distances were computed from per-participant FreeSurfer surfaces with the
# scripts in surface_producers/ (surfaces not shared).
#
# Run from the repository root:  Rscript centroid/centroid_distance_stats.R
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse); library(nlme); library(emmeans); library(patchwork)
})
emm_options(msg.interaction = FALSE)
options(width = 110)

out_dir <- file.path("outputs", "centroid")
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

THRESHOLDS <- c(2.5, 5, 10)
PRIMARY <- 5
LADDER <- c("prcus-p", "prcus-i", "prcus-a")          # posterior -> anterior

dist <- read.csv("data/centroid_distance.csv", stringsAsFactors = FALSE) %>%
  mutate(hemi = factor(hemi, levels = c("lh", "rh")),
         group = factor(group, levels = c("Controls", "DPs")),
         subject = factor(sub))
geom <- read.csv("data/face_patch_geometry.csv", stringsAsFactors = FALSE)
appos <- read.csv("data/sulcal_ap_position.csv", stringsAsFactors = FALSE)
stopifnot(n_distinct(dist$sub) == 47,
          all(count(dist, threshold_pct, hemi, sulcus)$n == 47))

sink(file.path(out_dir, "centroid_distance_stats.txt"), split = TRUE)

#### 1. Face-patch geometry ###################################################
cat("1. FACE-PATCH GEOMETRY (94 hemispheres per threshold)\n\n")
print(as.data.frame(geom %>% group_by(threshold_pct) %>% summarise(
  n_hemispheres = n(),
  components_median = median(n_components),
  components_range = sprintf("%d-%d", min(n_components), max(n_components)),
  top_component_mass_fraction_median = round(median(top_component_mass_fraction), 3),
  top_component_pct_of_PMC_median = round(100 * median(top_component_vertices / n_pmc_vertices), 2),
  .groups = "drop")), row.names = FALSE)

#### 2. Descriptives ##########################################################
desc <- dist %>% group_by(threshold_pct, hemi, sulcus) %>%
  summarise(n = n(), mean_mm = mean(geodesic_mm), sd_mm = sd(geodesic_mm),
            median_mm = median(geodesic_mm), .groups = "drop")
write.csv(desc, file.path(out_dir, "centroid_descriptives.csv"), row.names = FALSE)
cat("\n2. MEAN DISTANCE (mm), top 5% threshold\n\n")
print(as.data.frame(desc %>% filter(threshold_pct == PRIMARY) %>%
                      mutate(across(ends_with("_mm"), ~ round(.x, 1)))), row.names = FALSE)

#### 3. The precuneal gradient (lme), at each threshold #######################
fit_gradient <- function(pct) {
  d <- dist %>% filter(threshold_pct == pct, sulcus %in% LADDER) %>%
    mutate(sulcus = factor(sulcus, levels = LADDER))
  stopifnot(nrow(d) == 47 * 2 * 3)
  lme(geodesic_mm ~ sulcus * hemi * group, random = ~ 1 | subject, data = d, method = "REML")
}
models <- setNames(lapply(THRESHOLDS, fit_gradient), THRESHOLDS)

m5 <- models[[as.character(PRIMARY)]]
cat("\n3. PRECUNEAL GRADIENT, top 5%: lme(distance ~ sulcus * hemi * group, random = ~1 | subject)\n\n")
print(anova(m5))
cat("\nLinear contrast across prcus-p -> prcus-i -> prcus-a, pooled over hemispheres and groups\n")
cat("(coefficients -1, 0, 1: the prcus-a minus prcus-p difference)\n")
print(summary(contrast(emmeans(m5, ~ sulcus), method = "poly"), infer = c(TRUE, TRUE)))
cat("\nPairwise steps within hemisphere, Holm-adjusted\n")
steps <- as.data.frame(summary(contrast(emmeans(m5, ~ sulcus | hemi), method = "pairwise",
                                        adjust = "holm"), infer = c(TRUE, TRUE)))
print(steps)
cat("\nHemispheres in which the prcus-p is closer to the patch than the prcus-a:\n")
print(as.data.frame(dist %>% filter(threshold_pct == PRIMARY, sulcus %in% c("prcus-p", "prcus-a")) %>%
  select(sub, hemi, sulcus, geodesic_mm) %>%
  pivot_wider(names_from = sulcus, values_from = geodesic_mm) %>%
  group_by(hemi) %>% summarise(closer = sum(`prcus-p` < `prcus-a`), n = n(), .groups = "drop")),
  row.names = FALSE)

cat("\nGroup effects at each threshold (same model; sequential ANOVA terms)\n")
group_tab <- map_dfr(THRESHOLDS, function(pct) {
  a <- anova(models[[as.character(pct)]])
  tibble(threshold_pct = pct,
         group_F = a["group", "F-value"], group_df = sprintf("%d,%d", a["group", "numDF"], a["group", "denDF"]),
         group_p = a["group", "p-value"],
         sulcus_group_F = a["sulcus:group", "F-value"],
         sulcus_group_df = sprintf("%d,%d", a["sulcus:group", "numDF"], a["sulcus:group", "denDF"]),
         sulcus_group_p = a["sulcus:group", "p-value"],
         linear_trend_mm = summary(contrast(emmeans(models[[as.character(pct)]], ~ sulcus),
                                            method = "poly"))$estimate[1])
})
print(as.data.frame(group_tab %>% mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)

#### 4. Positional control: the prculs-d ######################################
cat("\n4a. A-P POSITION: prculs-d minus prcus-p (paired t; negative = prculs-d posterior)\n\n")
ap_tab <- appos %>% pivot_wider(names_from = sulcus, values_from = ap_position) %>%
  mutate(d = `prculs-d` - `prcus-p`) %>% group_by(hemi) %>%
  summarise(n = n(), prcus_p = mean(`prcus-p`), prculs_d = mean(`prculs-d`),
            t = t.test(d)$statistic, df = t.test(d)$parameter, p = t.test(d)$p.value,
            prculsd_posterior = sprintf("%d/%d", sum(d < 0), n()), .groups = "drop")
print(as.data.frame(ap_tab), row.names = FALSE)

paired_diff <- function(x) {
  tt <- t.test(x$d)
  tibble(n = nrow(x), prcus_p_mm = mean(x$`prcus-p`), prculs_d_mm = mean(x$`prculs-d`),
         diff_mm = mean(x$d), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
         t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
         farther_in = sprintf("%d/%d", sum(x$d > 0), nrow(x)))
}
ctrl <- dist %>% filter(sulcus %in% c("prcus-p", "prculs-d")) %>%
  select(sub, group, hemi, threshold_pct, sulcus, geodesic_mm) %>%
  pivot_wider(names_from = sulcus, values_from = geodesic_mm) %>%
  mutate(d = `prculs-d` - `prcus-p`)
ctrl_all <- ctrl %>% group_by(threshold_pct, hemi) %>% group_modify(~ paired_diff(.x)) %>%
  ungroup() %>% mutate(group = "both")
ctrl_grp <- ctrl %>% group_by(threshold_pct, hemi, group) %>% group_modify(~ paired_diff(.x)) %>%
  ungroup() %>% mutate(group = as.character(group))
ctrl_between <- ctrl %>% group_by(threshold_pct, hemi) %>%
  summarise(welch_t = t.test(d ~ group)$statistic, welch_df = t.test(d ~ group)$parameter,
            welch_p = t.test(d ~ group)$p.value, .groups = "drop")
write.csv(bind_rows(ctrl_all, ctrl_grp) %>% left_join(ctrl_between, by = c("threshold_pct", "hemi")),
          file.path(out_dir, "prculsd_control.csv"), row.names = FALSE)
cat("\n4b. DISTANCE: prculs-d minus prcus-p (paired t; positive = prculs-d farther), both groups\n\n")
print(as.data.frame(ctrl_all %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
cat("\n4c. Within each group\n\n")
print(as.data.frame(ctrl_grp %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
cat("\n4d. Group difference in (prculs-d - prcus-p), Welch's t\n\n")
print(as.data.frame(ctrl_between %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
sink()

#### 5. Figure 4 ##############################################################
PANEL_D <- c("mcgs", "prcus-a", "prcus-i", "prcus-p", "prculs-d", "pos")   # anterior -> posterior
LADDER_AP <- rev(LADDER)                                                   # anterior -> posterior
PAL <- c("mcgs" = "#0C314D", "prcus-a" = "#66B9B8", "prcus-i" = "#518DF0",
         "prcus-p" = "#B5CFFF", "prculs-d" = "#2F5952", "pos" = "#6850A0")
GRP <- c(Controls = "#0072B2", DPs = "#D55E00")
HEMI <- c(lh = "Left", rh = "Right")
stars <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "n.s.")
theme_pub <- theme_classic(base_size = 7.5) +
  theme(axis.title = element_text(face = "bold.italic", colour = "grey35"),
        axis.text = element_text(colour = "grey20"),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", colour = "grey35", size = 8),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7.5),
        legend.key.size = unit(3.2, "mm"), panel.spacing.x = unit(4, "mm"))
ylab_mm <- "Distance to face-selective patch (mm)"
prim <- dist %>% filter(threshold_pct == PRIMARY)

# A: the precuneal gradient, with Holm-corrected pairwise steps
bA <- steps %>%
  mutate(a = sub("^\\(?([^)]*)\\)? - .*$", "\\1", contrast),
         b = sub("^.* - \\(?([^)]*)\\)?$", "\\1", contrast),
         x1 = match(a, LADDER_AP), x2 = match(b, LADDER_AP),
         lo = pmin(x1, x2), hi = pmax(x1, x2), span = hi - lo) %>%
  group_by(hemi) %>% arrange(span, lo, .by_group = TRUE) %>%
  mutate(y = 72 + (row_number() - 1) * 6.5, lab = stars(p.value)) %>% ungroup()
pa_dat <- prim %>% filter(sulcus %in% LADDER) %>% mutate(sulcus = factor(sulcus, levels = LADDER_AP))
pA <- ggplot(pa_dat, aes(sulcus, geodesic_mm)) +
  geom_line(aes(group = sub), colour = "grey75", alpha = .5, linewidth = .2) +
  geom_boxplot(aes(fill = sulcus), width = .55, outlier.shape = NA, colour = "black", linewidth = .3, alpha = .9) +
  geom_point(position = position_jitter(width = .07, height = 0, seed = 1), size = .45,
             alpha = .55, colour = "grey15") +
  geom_segment(data = bA, aes(x = lo, xend = hi, y = y, yend = y), inherit.aes = FALSE, linewidth = .3) +
  geom_segment(data = bA, aes(x = lo, xend = lo, y = y, yend = y - 1.8), inherit.aes = FALSE, linewidth = .3) +
  geom_segment(data = bA, aes(x = hi, xend = hi, y = y, yend = y - 1.8), inherit.aes = FALSE, linewidth = .3) +
  geom_text(data = bA, aes(x = (lo + hi) / 2, y = y, label = lab), inherit.aes = FALSE,
            vjust = -.15, size = 2.6) +
  facet_wrap(~ hemi, labeller = labeller(hemi = HEMI)) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 80, 20), expand = expansion(0)) +
  labs(x = NULL, y = ylab_mm) + theme_pub

# B: robustness across face-patch thresholds (means +/- 1 SE)
pB <- dist %>% filter(sulcus %in% LADDER) %>%
  group_by(hemi, sulcus, threshold_pct) %>%
  summarise(m = mean(geodesic_mm), se = sd(geodesic_mm) / sqrt(n()), .groups = "drop") %>%
  mutate(sulcus = factor(sulcus, levels = LADDER_AP)) %>%
  ggplot(aes(factor(threshold_pct), m, group = sulcus)) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .15, linewidth = .3) +
  geom_line(aes(colour = sulcus), linewidth = .6) +
  geom_point(aes(fill = sulcus), shape = 21, size = 1.8, stroke = .3, colour = "black") +
  facet_wrap(~ hemi, labeller = labeller(hemi = HEMI)) +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
  scale_y_continuous(limits = c(10, 40), breaks = seq(10, 40, 10), expand = expansion(0)) +
  labs(x = "Face-patch threshold (top % of PMC)", y = "Mean distance (mm)") +
  theme_pub + theme(axis.text.x = element_text(angle = 0, hjust = .5))

# C: NT vs DP
pC <- ggplot(pa_dat, aes(sulcus, geodesic_mm, fill = group)) +
  geom_boxplot(width = .7, outlier.size = .35, colour = "black", linewidth = .3, alpha = .6,
               position = position_dodge(.8)) +
  facet_wrap(~ hemi, labeller = labeller(hemi = HEMI)) +
  scale_fill_manual(values = GRP, labels = c(Controls = "NTs", DPs = "DPs")) +
  scale_y_continuous(limits = c(0, 75), breaks = seq(0, 75, 25), expand = expansion(0)) +
  labs(x = NULL, y = ylab_mm) + theme_pub

# D: the precuneal sulci, the prculs-d control, and the bordering mcgs and pos
# (mcgs and pos are shown for illustration only; descriptive panel)
pD <- prim %>% filter(sulcus %in% PANEL_D) %>%
  mutate(sulcus = factor(sulcus, levels = PANEL_D), border = sulcus %in% c("mcgs", "pos")) %>%
  ggplot(aes(sulcus, geodesic_mm)) +
  geom_boxplot(aes(fill = sulcus, alpha = border), width = .6, outlier.size = .35, colour = "black", linewidth = .3) +
  facet_wrap(~ hemi, labeller = labeller(hemi = HEMI)) +
  scale_fill_manual(values = PAL, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = .9, `TRUE` = .35), guide = "none") +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 100, 25), expand = expansion(0)) +
  labs(x = NULL, y = ylab_mm) + theme_pub

fig <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))
for (ext in c("png", "pdf"))
  ggsave(file.path(out_dir, "figures", paste0("fig4_centroid_distance.", ext)), fig,
         width = 6.5, height = 5.6, dpi = 300, bg = "white")
cat("Outputs written to", out_dir, "\n")
