###############################################################################
# dice/dice_group_tests.R
#
# Group comparisons (NT vs DP) of the overlap (Dice similarity coefficient,
# DSC) between the eight consistently present PMC sulci and the maximally
# face-selective region, and Figure 5B.
#
# Reproduces
#   Results, "The prcus-p and its immediate vicinity overlap...", group
#     paragraph: group main effect and sulcus-by-group interaction
#     (permutation tests) and the prcus-p NT vs DP comparison (Wilcoxon
#     rank-sum), uncorrected and Holm-corrected across the six
#     window-by-hemisphere tests, at the primary 5% threshold
#   Same paragraph and Supplemental Material: the same tests at the 2.5% and 10%
#     thresholds, including the count of tests reaching uncorrected p < .05
#   Figure 5B: mean DSC of each sulcus as a function of dilation distance
#
# Method
#   Each participant's DSC per window (0 mm, 0-5 mm, 0-10 mm) is the
#   trapezoidal area under the DSC-by-distance curve divided by the window
#   width (the 0-mm window is the DSC of the sulcus itself).
#   Permutation test (group labels shuffled across participants): statistics
#   are |mean over sulci of the DP - NT difference| (group main effect) and
#   the variance across sulci of that difference (sulcus-by-group).
#   p = (number of permuted statistics >= observed + 1) / (N_PERM + 1).
#   Holm correction is applied across the six window-by-hemisphere tests,
#   separately for each test type and threshold.
#
# Inputs (data/)   dice_overlap.csv
# Outputs (outputs/dice/)
#   dice_group_tests.csv     every group test, uncorrected and Holm p
#   dice_group_tests.txt     the same as text, with the counts quoted in the paper
#   figures/fig5B_dice_curves.{png,pdf}
#
# Usage (from the repository root):
#   Rscript dice/dice_group_tests.R [N_PERM]
# N_PERM defaults to 10000, the number used in the paper. The random seed is
# reset to 20261001 at the start of each threshold, so the default run
# reproduces the reported permutation p-values exactly.
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
N_PERM <- if (length(args) >= 1) as.integer(args[1]) else 10000L
SEED <- 20261001

out_dir <- file.path("outputs", "dice")
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# Sulcus order as in Figure 5; "sbps" in the data is the manuscript's "spls".
SULCI <- c("prcus-p", "prcus-i", "prcus-a", "prculs-d", "pos", "spls", "mcgs", "ifrms")
LADDER <- c(0, 2.5, 5, 7.5, 10)
WINDOWS <- list("0 mm" = 0, "0-5 mm" = LADDER[LADDER <= 5], "0-10 mm" = LADDER)
HEMI <- c(lh = "Left", rh = "Right")
THRESHOLDS <- c(5, 2.5, 10)

dice_all <- read.csv("data/dice_overlap.csv", stringsAsFactors = FALSE) %>%
  mutate(sulcus = recode(sulcus, sbps = "spls")) %>%
  filter(sulcus %in% SULCI, radius_mm %in% LADDER) %>%
  mutate(sulcus = factor(sulcus, levels = SULCI), group = factor(group, c("Controls", "DPs")),
         sub = factor(sub, levels = unique(sub)))   # participants kept in file order (fixes the
                                                    # permutation sequence for a given seed)

trap <- function(x, y) if (length(x) == 1) y else
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
perm_stats <- function(w, g) {                  # w: participants x sulci, g: TRUE = DP
  diff <- colMeans(w[g, , drop = FALSE]) - colMeans(w[!g, , drop = FALSE])
  c(main = abs(mean(diff)), inter = sum((diff - mean(diff))^2))
}

results <- map_dfr(THRESHOLDS, function(pct) {
  d <- dice_all %>% filter(threshold_pct == pct)
  stopifnot(all(count(d, hemi, sulcus, radius_mm)$n == 47))
  agg <- map_dfr(names(WINDOWS), function(w) d %>% filter(radius_mm %in% WINDOWS[[w]]) %>%
    arrange(radius_mm) %>% group_by(sub, group, hemi, sulcus) %>%
    summarise(stat = trap(radius_mm, dice), .groups = "drop") %>% mutate(window = w))
  set.seed(SEED)
  map_dfr(names(WINDOWS), function(w) map_dfr(names(HEMI), function(h) {
    x <- filter(agg, window == w, hemi == h)
    wide <- x %>% select(sub, group, sulcus, stat) %>%
      pivot_wider(names_from = sulcus, values_from = stat)
    M <- as.matrix(wide[, SULCI]); g <- wide$group == "DPs"
    obs <- perm_stats(M, g)
    null <- replicate(N_PERM, perm_stats(M, sample(g)))
    pp <- (rowSums(null >= obs) + 1) / (N_PERM + 1)
    pc <- filter(x, sulcus == "prcus-p")
    tibble(threshold_pct = pct, window = w, hemi = h,
           p_group_perm = pp[["main"]], p_sulcus_by_group_perm = pp[["inter"]],
           prcusp_NT_mean = mean(pc$stat[pc$group == "Controls"]),
           prcusp_DP_mean = mean(pc$stat[pc$group == "DPs"]),
           p_prcusp_ranksum = suppressWarnings(wilcox.test(stat ~ group, data = pc, exact = FALSE)$p.value))
  }))
}) %>%
  group_by(threshold_pct) %>%
  mutate(holm_group = p.adjust(p_group_perm, "holm"),
         holm_sulcus_by_group = p.adjust(p_sulcus_by_group_perm, "holm"),
         holm_prcusp_ranksum = p.adjust(p_prcusp_ranksum, "holm")) %>% ungroup()
write.csv(results, file.path(out_dir, "dice_group_tests.csv"), row.names = FALSE)

sink(file.path(out_dir, "dice_group_tests.txt"), split = TRUE)
cat("DSC GROUP TESTS (NT n = 25, DP n = 22); permutations per test:", N_PERM, "\n\n")
print(as.data.frame(results %>% mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)
cat("\nRanges per threshold (uncorrected / Holm across the six window-by-hemisphere tests):\n")
print(as.data.frame(results %>% group_by(threshold_pct) %>% summarise(
  min_p_group = min(p_group_perm), min_holm_group = min(holm_group),
  min_p_sulcus_by_group = min(p_sulcus_by_group_perm), min_holm_sulcus_by_group = min(holm_sulcus_by_group),
  min_p_prcusp = min(p_prcusp_ranksum), min_holm_prcusp = min(holm_prcusp_ranksum),
  n_uncorrected_below_05 = sum(c(p_group_perm, p_sulcus_by_group_perm, p_prcusp_ranksum) < .05),
  n_tests = 3 * n(), .groups = "drop") %>% mutate(across(where(is.double), ~ round(.x, 4)))),
  row.names = FALSE)
sink()

#### Figure 5B: mean DSC by dilation distance (NTs and DPs combined), top 5% ####
PAL <- c("prcus-p" = "#B5CFFF", "prcus-i" = "#518DF0", "prcus-a" = "#66B9B8",
         "prculs-d" = "#2F5952", "spls" = "#858CC2", "mcgs" = "#0C314D", "pos" = "#6850A0",
         "ifrms" = "#992F27")
curve <- dice_all %>% filter(threshold_pct == 5) %>%
  group_by(hemi, sulcus, radius_mm) %>%
  summarise(m = mean(dice), se = sd(dice) / sqrt(n()), .groups = "drop") %>%
  mutate(hemi = factor(HEMI[hemi], HEMI))
p <- ggplot(curve, aes(radius_mm, m)) + facet_wrap(~ hemi) +
  geom_ribbon(aes(ymin = pmax(m - se, 0), ymax = m + se, fill = sulcus), alpha = .18, colour = NA) +
  geom_line(aes(colour = sulcus), linewidth = .6) +
  geom_point(aes(fill = sulcus), shape = 21, size = 1.1, stroke = .25, colour = "black") +
  scale_colour_manual(values = PAL) + scale_fill_manual(values = PAL) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = LADDER, expand = expansion(mult = c(.02, .03))) +
  scale_y_continuous(limits = c(0, max(curve$m + curve$se) * 1.03), expand = expansion(0)) +
  labs(x = "Dilation from the sulcal border (mm)", y = "Dice overlap with face-selective area") +
  theme_classic(base_size = 7.5) +
  theme(axis.title = element_text(face = "bold.italic", colour = "grey35"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", colour = "grey35", size = 8),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 7.5))
for (ext in c("png", "pdf"))
  ggsave(file.path(out_dir, "figures", paste0("fig5B_dice_curves.", ext)), p,
         width = 6.5, height = 2.9, dpi = 300, bg = "white")
cat("Outputs written to", out_dir, "\n")
