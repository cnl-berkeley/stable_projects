###############################################################################
# demographics_covariates.R
#
# Participant demographics, group matching, demographic covariates of sulcal
# presence, and the distribution of Cambridge Face Memory Test (CFMT) scores
# for Datasets 1 and 2 (NT n = 43, DP n = 39).
#
# Reproduces
#   Abstract: percentage of female participants (63.4%)
#   Methods, "Participants (Dataset 1)" and "Participants (Dataset 2)":
#     numbers of participants and females, age mean (SD) and range
#   Methods, "Group matching (Datasets 1 and 2)": age by group (Welch's t test,
#     Wilcoxon rank-sum test) and sex by group (chi-squared with Yates'
#     continuity correction, Fisher's exact test)
#   Methods, "Face processing behavioral tasks": CFMT mean (SD), range and N
#   Methods, "Behavioral analysis": demographics of the 42 NTs with a CFMT score
#   Methods, "Extracting incidence rates...": logistic regressions of the
#     presence of each variably present sulcus on age and sex, per hemisphere
#     (8 models, 16 covariate tests)
#   Figure S1: CFMT score distributions by group
#   Methods, "Participants (Dataset 3)": number of HCP participants, number of
#     females, age mean (SD) and range (only if data/hcp_demographics.csv is supplied;
#     not redistributed, see README)
#
# Inputs (data/)
#   demographics.csv       one row per participant of Datasets 1 and 2 (82)
#   sulcus_presence.csv    presence/absence of each PMC sulcus per hemisphere
#   hcp_demographics.csv   optional, not included: one row per HCP participant (71): sub, age, sex
#
# Outputs (outputs/demographics/)
#   demographics_summary.txt          all statistics listed above
#   logistic_age_sex_presence.csv     coefficients of the 8 logistic models
#   figures/figS1_cfmt_distributions.{png,pdf}
#
# Run from the repository root:  Rscript demographics_covariates.R
###############################################################################

suppressPackageStartupMessages(library(tidyverse))

out_dir <- file.path("outputs", "demographics")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

demo <- read.csv("data/demographics.csv", stringsAsFactors = FALSE)
presence <- read.csv("data/sulcus_presence.csv", stringsAsFactors = FALSE)
stopifnot(nrow(demo) == 82, setequal(demo$sex, c("F", "M")))

sink(file.path(out_dir, "demographics_summary.txt"), split = TRUE)
options(width = 110)

describe <- function(d) d %>%
  summarise(n = n(), n_female = sum(sex == "F"), pct_female = round(100 * mean(sex == "F"), 1),
            age_mean = round(mean(age), 2), age_sd = round(sd(age), 2),
            age_min = min(age), age_max = max(age), .groups = "drop")

cat("1. DEMOGRAPHICS\n\n")
cat("All participants, Datasets 1 and 2:\n")
print(as.data.frame(describe(demo)))
cat("\nBy dataset and group:\n")
print(as.data.frame(describe(group_by(demo, dataset, group))))
cat("\nCombined sample by group:\n")
print(as.data.frame(describe(group_by(demo, group))))

cat("\n2. GROUP MATCHING (combined sample, 43 NTs vs 39 DPs)\n\n")
print(t.test(age ~ group, data = demo))                       # Welch's t test
print(wilcox.test(age ~ group, data = demo))
sex_tab <- table(demo$sex, demo$group)
print(sex_tab)
print(chisq.test(sex_tab))                                    # Yates' correction
print(fisher.test(sex_tab))

cat("\n3. CFMT (participants with a score)\n\n")
cfmt <- demo %>% filter(!is.na(CFMT))
print(as.data.frame(cfmt %>% group_by(group) %>%
  summarise(n = n(), mean = round(mean(CFMT), 1), sd = round(sd(CFMT), 1),
            min = min(CFMT), max = max(CFMT), .groups = "drop")))
cat("\nLASSO sample (NTs with a CFMT score):\n")
print(as.data.frame(describe(cfmt %>% filter(group == "Controls"))))

cat("\n4. LOGISTIC REGRESSION: presence ~ age + sex, per variably present sulcus and hemisphere\n\n")
variable_sulci <- c("sspls-v", "prculs-v", "sspls-d", "icgs-p")
logit <- presence %>%
  filter(sulcus %in% variable_sulci) %>%
  left_join(demo %>% select(sub, age, sex), by = "sub") %>%
  group_by(sulcus, hemi) %>%
  group_modify(function(d, k) {
    co <- summary(glm(present ~ age + sex, data = d, family = binomial))$coefficients
    tibble(term = rownames(co), estimate = co[, 1], se = co[, 2], z = co[, 3], p = co[, 4],
           n = nrow(d), n_present = sum(d$present))
  }) %>% ungroup()
write.csv(logit, file.path(out_dir, "logistic_age_sex_presence.csv"), row.names = FALSE)
print(as.data.frame(logit %>% filter(term != "(Intercept)") %>%
                      mutate(across(c(estimate, se, z, p), ~ round(.x, 4)))))
cat("\nMinimum p, age: ", round(min(logit$p[logit$term == "age"]), 4),
    "   Minimum p, sex: ", round(min(logit$p[logit$term == "sexM"]), 4), "\n", sep = "")

cat("\n5. HCP PARTICIPANTS (Dataset 3)\n\n")
hcp_file <- "data/hcp_demographics.csv"
if (file.exists(hcp_file)) {
  hcp <- read.csv(hcp_file, stringsAsFactors = FALSE)
  hcp_ids <- read.csv("data/face_selectivity.csv", stringsAsFactors = FALSE) %>%
    filter(group == "HCP_Controls") %>% distinct(sub)
  stopifnot(nrow(hcp) == 71, setequal(as.character(hcp$sub), as.character(hcp_ids$sub)),
            setequal(hcp$sex, c("F", "M")), !anyNA(hcp$age))
  print(as.data.frame(describe(hcp)))
} else {
  cat("data/hcp_demographics.csv not found; HCP demographics not computed.\n")
}
sink()

#### Figure S1: CFMT distributions ############################################
brks <- seq(20, 75, by = 5)
bins <- cfmt %>%
  mutate(group = factor(group, levels = c("Controls", "DPs")),
         bin = cut(CFMT, breaks = brks, right = TRUE, include.lowest = TRUE)) %>%
  count(group, bin, .drop = TRUE) %>%
  mutate(lo = brks[as.integer(bin)], hi = lo + 5)

p <- ggplot(bins, aes(xmin = lo, xmax = hi, ymin = 0, ymax = n, fill = group)) +
  geom_rect(alpha = 0.6, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(Controls = "#0072B2", DPs = "#D55E00"), name = NULL,
                    labels = c(Controls = "NTs", DPs = "DPs")) +
  scale_x_continuous(breaks = seq(20, 75, by = 10), limits = c(20, 75), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 15, 5), expand = expansion(0)) +
  labs(x = "CFMT score", y = "Number of participants") +
  theme_classic(base_size = 10) +
  theme(axis.title = element_text(face = "bold.italic", colour = "grey35", size = 10),
        axis.text = element_text(colour = "grey20", size = 9),
        legend.position = "bottom", legend.text = element_text(face = "italic", size = 10),
        plot.margin = margin(6, 10, 4, 4))
for (ext in c("png", "pdf"))
  ggsave(file.path(fig_dir, paste0("figS1_cfmt_distributions.", ext)), p,
         width = 5.42, height = 3.16, dpi = 300, bg = "white")
cat("Outputs written to", out_dir, "\n")
