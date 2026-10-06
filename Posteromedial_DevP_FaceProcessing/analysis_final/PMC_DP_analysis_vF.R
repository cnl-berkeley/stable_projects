###############################################################################
# PMC_DP_analysis_vF.R
#
# Face selectivity and sulcal incidence of posteromedial cortex (PMC) sulci in
# neurotypical (NT) participants and participants with developmental
# prosopagnosia (DP), with a replication of the face-selectivity analyses in
# 71 Human Connectome Project (HCP) participants.
#
# Reproduces
#   Results, "The ventral sub-splenial sulcus (sspls-v) shows reduced incidence
#     and face selectivity in developmental prosopagnosia" (all statistics)
#   Results, "A sulcal-functional gradient of face selectivity along the
#     anterior-posterior axis of PMC" (all statistics)
#   Supplemental Material, "Face selectivity of all sulci in PMC"
#   Figure 2B (sspls-v incidence) and Figure 2C (sspls-v face selectivity)
#   Figure 3B (precuneal face selectivity, Dataset 1 and HCP)
#   Table 1 (incidence of the 12 PMC sulci)
#
# Statistical procedures (as stated in the Methods)
#   Omnibus: mixed-model ANOVA (rstatix::anova_test) on mean face selectivity
#     with sulcus, hemisphere and group (Dataset 1) or sulcus and hemisphere
#     (HCP), with hemisphere as the repeated measure.
#   Post hoc tests: linear mixed-effects model fit to every available
#     observation (lmerTest::lmer, random intercepts for participant and
#     hemisphere within participant), contrasts via emmeans.
#     - one-sample tests against zero, Bonferroni-adjusted in two families:
#         x2  : across the two hemispheres within group (used for the sspls-v)
#         x12 : across the 12 sulci within each hemisphere and group (used for
#               the tests of all PMC sulci in the Supplemental Material and for the
#               prcus-i / prcus-p statement in the Results)
#     - pairwise between-sulcus contrasts within hemisphere and group, Tukey
#     - NT vs DP and right vs left contrasts within each sulcus
#   Incidence: chi-squared tests with Yates' continuity correction.
#
# Inputs (data/)
#   face_selectivity.csv   mean face selectivity per participant, hemisphere
#                          and sulcus (Dataset 1: % signal change, faces -
#                          objects; HCP: z, faces - all other categories)
#   sulcus_presence.csv    presence (1) / absence (0) of each of the 12 PMC
#                          sulci in each hemisphere of the 82 participants of
#                          Datasets 1 and 2
#
# Outputs (outputs/selectivity_incidence/)
#   dataset1_omnibus_anova.txt, hcp_omnibus_anova.txt
#   dataset1_selectivity_vs_zero_bonferroni_x2_hemispheres.csv
#   dataset1_selectivity_vs_zero_bonferroni_x12_sulci.csv
#   hcp_selectivity_vs_zero_bonferroni_x2_hemispheres.csv
#   hcp_selectivity_vs_zero_bonferroni_x12_sulci.csv
#   dataset1_between_sulci_tukey.csv, hcp_between_sulci_tukey.csv
#   dataset1_NT_vs_DP.csv
#   dataset1_right_vs_left.csv, hcp_right_vs_left.csv
#   incidence_counts_table1.csv
#   incidence_chisq.csv       (sspls-v: reported in the Results; icgs-p,
#                              sspls-d, prculs-v: exploratory, not reported)
#   figures/fig2B_ssplsv_incidence.pdf, fig2C_ssplsv_selectivity.pdf,
#   figures/fig3B_precuneal_selectivity.pdf
#
# Run from the repository root:  Rscript PMC_DP_analysis_vF.R
#
# Note: the sulcus called "spls" in the manuscript is labelled "sbps" in the
# data files (its alternative name, the subparietal sulcus).
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, ggplot2, readr
  library(rstatix)     # anova_test
  library(lmerTest)    # lmer with Satterthwaite degrees of freedom
  library(emmeans)
  library(patchwork)
})
emm_options(msg.interaction = FALSE)

out_dir <- file.path("outputs", "selectivity_incidence")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# The 12 PMC sulci, in the order used for model factors and tables.
pmc_sulci <- c("mcgs", "sbps", "icgs-p", "ifrms", "sspls-d", "sspls-v",
               "prcus-a", "prcus-i", "prcus-p", "prculs-v", "prculs-d", "pos")
variable_sulci <- c("sspls-v", "icgs-p", "sspls-d", "prculs-v")

face <- read.csv("data/face_selectivity.csv") %>%
  filter(sulcus %in% pmc_sulci) %>%
  mutate(label = factor(sulcus, levels = pmc_sulci))
presence <- read.csv("data/sulcus_presence.csv", stringsAsFactors = FALSE)

emm_df <- function(x, ...) as.data.frame(summary(x, ...))

# One-sample tests against zero in both Bonferroni families.
#   by_hemi_group: emmeans grouping for the x2 family (hemispheres within group)
vs_zero_x2 <- function(mod, by) {
  bind_rows(lapply(pmc_sulci, function(lab)
    emm_df(emmeans(mod, by, at = list(label = lab)),
           infer = c(FALSE, TRUE), null = 0, adjust = "bonferroni") %>%
      mutate(label = lab))) %>%
    mutate(adjustment = "Bonferroni x2 (hemispheres within group)")
}
vs_zero_x12 <- function(mod, spec) {
  emm_df(emmeans(mod, spec), infer = c(FALSE, TRUE), null = 0,
         adjust = "bonferroni") %>%
    mutate(adjustment = "Bonferroni x12 (sulci within hemisphere and group)")
}


#### 1. Face selectivity, Dataset 1 (NT n = 25, DP n = 22) ####################
# Participant factor levels are set explicitly (Dataset 1: file order; HCP: ascending subject ID): the
# mixed-model optimizer is sensitive to factor-level order in the last reported digits.
d1 <- face %>% filter(group %in% c("Controls", "DPs")) %>%
  mutate(sub = factor(sub, levels = unique(sub)))

d1_anova <- d1 %>%
  anova_test(mean_face_activation ~ group * hemi * label + Error(sub / hemi))
capture.output(print(d1_anova), file = file.path(out_dir, "dataset1_omnibus_anova.txt"))
print(d1_anova)

d1_mod <- lmer(mean_face_activation ~ label * group * hemi + (1 | sub / hemi), data = d1)

write.csv(vs_zero_x2(d1_mod, ~ hemi | group),
          file.path(out_dir, "dataset1_selectivity_vs_zero_bonferroni_x2_hemispheres.csv"),
          row.names = FALSE)
write.csv(vs_zero_x12(d1_mod, ~ label | hemi | group),
          file.path(out_dir, "dataset1_selectivity_vs_zero_bonferroni_x12_sulci.csv"),
          row.names = FALSE)
write.csv(emm_df(contrast(emmeans(d1_mod, ~ label | hemi | group),
                          method = "pairwise", adjust = "tukey")),
          file.path(out_dir, "dataset1_between_sulci_tukey.csv"), row.names = FALSE)
d1_group <- emm_df(contrast(emmeans(d1_mod, ~ group | hemi | label), method = "pairwise"))
write.csv(d1_group, file.path(out_dir, "dataset1_NT_vs_DP.csv"), row.names = FALSE)
write.csv(emm_df(contrast(emmeans(d1_mod, ~ hemi | group | label), method = "pairwise")),
          file.path(out_dir, "dataset1_right_vs_left.csv"), row.names = FALSE)
print(d1_group %>% filter(label == "sspls-v"))


#### 2. Face selectivity, HCP replication (n = 71) ###########################
hcp <- face %>% filter(group == "HCP_Controls") %>%
  mutate(sub = factor(sub, levels = sort(unique(as.numeric(sub)))))

hcp_anova <- hcp %>%
  anova_test(mean_face_activation ~ hemi * label + Error(sub / hemi))
capture.output(print(hcp_anova), file = file.path(out_dir, "hcp_omnibus_anova.txt"))
print(hcp_anova)

hcp_mod <- lmer(mean_face_activation ~ label * hemi + (1 | sub / hemi), data = hcp)

write.csv(vs_zero_x2(hcp_mod, ~ hemi) %>%
            mutate(adjustment = "Bonferroni x2 (hemispheres)"),
          file.path(out_dir, "hcp_selectivity_vs_zero_bonferroni_x2_hemispheres.csv"),
          row.names = FALSE)
write.csv(vs_zero_x12(hcp_mod, ~ label | hemi) %>%
            mutate(adjustment = "Bonferroni x12 (sulci within hemisphere)"),
          file.path(out_dir, "hcp_selectivity_vs_zero_bonferroni_x12_sulci.csv"),
          row.names = FALSE)
write.csv(emm_df(contrast(emmeans(hcp_mod, ~ label | hemi), method = "pairwise",
                          adjust = "tukey")),
          file.path(out_dir, "hcp_between_sulci_tukey.csv"), row.names = FALSE)
write.csv(emm_df(contrast(emmeans(hcp_mod, ~ hemi | label), method = "pairwise")),
          file.path(out_dir, "hcp_right_vs_left.csv"), row.names = FALSE)


#### 3. Sulcal incidence, Datasets 1 and 2 (NT n = 43, DP n = 39) ############
# Table 1: counts and percentages of hemispheres in which each sulcus is present.
incidence <- presence %>%
  group_by(sulcus, group, hemi) %>%
  summarise(n_present = sum(present), n_hemispheres = n(),
            percent = round(100 * mean(present), 1), .groups = "drop") %>%
  mutate(sulcus = factor(sulcus, levels = c(variable_sulci,
                                            setdiff(pmc_sulci, variable_sulci)))) %>%
  arrange(sulcus, group, hemi)
write.csv(incidence, file.path(out_dir, "incidence_counts_table1.csv"), row.names = FALSE)
cat("\nTotal sulci labelled:", sum(presence$present), "in",
    nrow(distinct(presence, sub, hemi)), "hemispheres of",
    n_distinct(presence$sub), "participants\n")

# Chi-squared tests (Yates' continuity correction, the chisq.test default for 2 x 2).
chisq_row <- function(dat, factor_name, sulcus, comparison) {
  ct <- chisq.test(dat[[factor_name]], dat$present)
  tibble(sulcus = sulcus, comparison = comparison,
         chisq = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value)
}
incidence_tests <- bind_rows(lapply(variable_sulci, function(s) {
  x <- presence %>% filter(sulcus == s)
  bind_rows(
    chisq_row(filter(x, hemi == "rh"), "group", s, "NT vs DP, right hemisphere"),
    chisq_row(filter(x, hemi == "lh"), "group", s, "NT vs DP, left hemisphere"),
    chisq_row(filter(x, group == "DPs"), "hemi", s, "right vs left, DPs"),
    chisq_row(filter(x, group == "Controls"), "hemi", s, "right vs left, NTs"))
})) %>%
  mutate(reported = ifelse(sulcus == "sspls-v", "Results", "exploratory (not reported)"))
write.csv(incidence_tests, file.path(out_dir, "incidence_chisq.csv"), row.names = FALSE)
print(as.data.frame(incidence_tests %>% filter(sulcus == "sspls-v")))


#### 4. Figures ###############################################################
hemi_labs <- c(lh = "Left", rh = "Right")
grp_fill <- scale_fill_manual(labels = c("NTs", "DPs"), values = c("#0072B2", "#D55E00"))
hcp_fill <- scale_fill_manual(labels = c("HCP NTs"), values = c("#332288"))

mean_activation_plot <- function(data) {
  ggplot(data, aes(x = label, y = mean_face_activation, fill = group)) +
    stat_summary(fun = mean, position = position_dodge(width = 0.95), geom = "bar",
                 color = "black", alpha = .6) +
    stat_summary(fun.data = mean_se, position = position_dodge(0.95),
                 geom = "errorbar", width = 0.2) +
    facet_wrap(vars(hemi), labeller = labeller(hemi = hemi_labs)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          axis.title.x = element_blank(), legend.position = "bottom",
          legend.title = element_blank())
}

# Figure 3B: face selectivity of the precuneal sulci (left: Dataset 1, right: HCP)
precuneal <- c("prcus-a", "prcus-i", "prcus-p")
p3_d1 <- mean_activation_plot(d1 %>% filter(label %in% precuneal) %>% droplevels()) +
  ylab("Faces - Objects (%)") + grp_fill + coord_cartesian(ylim = c(-.053, .6))
p3_hcp <- mean_activation_plot(hcp %>% filter(label %in% precuneal) %>% droplevels()) +
  ylab("Faces - All Other (Z-score)") + hcp_fill + coord_cartesian(ylim = c(NA, 1.6))
ggsave(file.path(fig_dir, "fig3B_precuneal_selectivity.pdf"),
       p3_d1 + p3_hcp + plot_layout(widths = c(3.6, 2)), width = 8, height = 3.5)

# Figure 2B: sspls-v incidence; Figure 2C: sspls-v face selectivity
ssplsv_inc <- incidence %>% filter(sulcus == "sspls-v")
p2b <- ggplot(ssplsv_inc, aes(x = group, y = percent, fill = group)) +
  geom_col(aes(alpha = hemi), color = "black", position = "dodge") +
  grp_fill +
  scale_alpha_manual(values = c(.8, 0.4), labels = c("Left", "Right")) +
  scale_x_discrete(labels = c(Controls = "NTs", DPs = "DPs")) +
  scale_y_continuous(name = "% of Hemispheres", breaks = seq(0, 100, 20), limits = c(0, 100)) +
  theme_classic() +
  theme(legend.position = "right", legend.title = element_blank(),
        axis.title.x = element_blank())
p2c_d1 <- mean_activation_plot(d1 %>% filter(label == "sspls-v") %>% droplevels()) +
  ylab("Faces - Objects (%)") + grp_fill + coord_cartesian(ylim = c(-0.17, 0.65))
p2c_hcp <- mean_activation_plot(hcp %>% filter(label == "sspls-v") %>% droplevels()) +
  ylab("Faces - All Other (Z-score)") + hcp_fill + coord_cartesian(ylim = c(-0.44, 1.65))
ggsave(file.path(fig_dir, "fig2B_ssplsv_incidence.pdf"), p2b, width = 3.5, height = 3.5)
ggsave(file.path(fig_dir, "fig2C_ssplsv_selectivity.pdf"),
       p2c_d1 + p2c_hcp + plot_layout(widths = c(2, 1)), width = 6, height = 3.5)

cat("\nOutputs written to", out_dir, "\n")
