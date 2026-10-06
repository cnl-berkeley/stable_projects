library(tidyverse)
library(car)
library(cowplot)
library(PupillometryR)
library(rstatix)
library(DescTools)
library(corrplot)
library(RVAideMemoire)
library(viridis)
library(lmPerm)
library(ggpubr)
library(ggsignif)

options(scipen=999, # 0
        show.signif.stars = FALSE) 

theme_set(theme_classic())
project_theme2 <- theme(
  # set axis title size
  axis.title.x = element_text(color="#333333", size=12),
  axis.title.y = element_text(color="#333333"),
  
  # set axis text size
  axis.text = element_text(color='#333333', size = 12),
  
  # rotate x-axis 45 degrees
  #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  
  #change background pannel
  panel.grid.major  = element_blank(), 
  panel.grid.minor=element_blank(),
  # change lengend size
  legend.title = element_text(color="#333333", size = 4),
  legend.text = element_text(color = '#333333', size = 10),
  legend.position="top",
  
  strip.text.x = element_text(size = 12),
  strip.placement = "outside"
)

project_theme3 <- theme(
  # set axis title size
  axis.title.x = element_text(color="#333333", size=12),
  axis.title.y = element_text(color="#333333", angle=90, size=12),
  
  # set axis text size
  axis.text.x = element_text(color='#333333', size = 12),
  axis.text.y = element_text(color='#333333', size = 12, angle = 90, hjust = 0.5),
  
  # rotate x-axis 45 degrees
  #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  
  #change background pannel
  panel.grid.major  = element_blank(), 
  panel.grid.minor=element_blank(),
  # change lengend size
  legend.title = element_text(color="#333333", size = 12),
  legend.text = element_text(color = '#333333', size = 12)
)

# set working directory
setwd()

# get data
pcgs_labels <- read.csv("./data/pcgs_labels.csv")
urgency_metrics <-  read.csv("./data/3Factor_impulsivity.csv")
urgency_demo <- read.csv("./data/demographics_urgency.csv")
urgency_scid <- read.csv("./data/scid_approach.csv")
pcgs_length_lh <- read.csv("./data/PCGS_lh_sulcal_length.csv")
pcgs_length_rh <- read.csv("./data/rh_PCGS_sulcal_length.csv")
pcgs_length <- rbind(pcgs_length_lh, pcgs_length_rh)
acc_lgi <-  read.csv("./data/acc_lgi.csv")
pcgs_depth <- read.csv("./data/pcgs_depth.csv")

# merge sulci + metrics
pcgs3factor <- merge(pcgs_labels, urgency_metrics, by = "ID")

# take out only relevant metrics from urgency metrics 
pcgs_3factor.clean <- pcgs_3factor %>% subset(PCGS_LH %in% c("present", "absent")) %>% dplyr::select(ID, PCGS_LH, PCGS_RH, PCGS_Asymmetry, Factor_FeelingsTrigger, Factor_PervInf, Factor_LackFollowThrough)
pcgs_3factor.clean <- pcgs_3factor.clean %>% mutate(asymmetry = case_when(PCGS_LH == "present" & PCGS_RH == "present" ~ "symmetrical",
                                                                          PCGS_LH == "present" & PCGS_RH == "absent" ~ "leftward",
                                                                          PCGS_LH == "absent" & PCGS_RH == "present" ~ "rightward",
                                                                          PCGS_LH == "absent" & PCGS_RH == "absent" ~ "symmetrical"))

# also add demographic variables 
pcgs_3factor.clean1 <- merge(pcgs_3factor.clean, urgency_demo, by = "ID", all.x = T, all.y = F)
pcgs_3factor.clean2 <- merge(pcgs_3factor.clean, urgency_scid, by = "ID", all.x = T, all.y = F)
#### data for table 1 ####
pcgs_3factor.clean1
pcgs_3factor.clean2$Lifetime_psychosis
pcgs_3factor.clean1 %>% group_by(Gender) %>% summarise(n = n(),
                                                      "perc" = n/120*100)
pcgs_3factor.clean1 %>% group_by(Race) %>% summarise(n = n(),
                                                     "perc" = n/120*100)
pcgs_3factor.clean2 %>% summarise(mean = mean(Education, na.rm=T),
                                  sd = sd(Education, na.rm=T),
                                  range = range(Education, na.rm=T))
pcgs_3factor.clean2 %>% group_by(Hispanic) %>% summarise(n = n(),
                                                     "perc" = n/120*100)
pcgs_3factor.clean2 %>% group_by(lifetime_MDE_dx) %>% summarise(n = n(),
                                                         "perc" = n/120*100)
pcgs_3factor.clean2 %>% group_by(Lifetime_anxiety_dx) %>% summarise(n = n(),
                                                                "perc" = n/120*100)
pcgs_3factor.clean2 %>% group_by(lifetimeAUD) %>% summarise(n = n(),
                                                                    "perc" = n/120*100)
pcgs_3factor.clean2 %>% group_by(lifetimeSUD) %>% summarise(n = n(),
                                                            "perc" = n/120*100)
pcgs_3factor.clean2 %>% group_by(Lifetime_psychosis) %>% summarise(n = n(),
                                                            "perc" = n/120*100)

pcgs_3factor.clean2 <- pcgs_3factor.clean2 %>% mutate(more_than_1 = lifetime_MDE_dx + Lifetime_anxiety_dx + lifetimeAUD + lifetimeSUD + Lifetime_psychosis)
pcgs_3factor.clean2 %>% group_by(more_than_1) %>% summarise(n = n(),
                                                                   "perc" = n/120*100)

pcgs_3factor.clean1 %>% summarise(mean = mean(Factor_PervInf, na.rm=T),
                                  sd = sd(Factor_PervInf, na.rm=T),
                                  range = range(Factor_PervInf, na.rm=T))
pcgs_3factor.clean1 %>% summarise(mean = mean(Factor_FeelingsTrigger, na.rm=T),
                                  sd = sd(Factor_FeelingsTrigger, na.rm=T),
                                  range = range(Factor_FeelingsTrigger, na.rm=T))
pcgs_3factor.clean1 %>% summarise(mean = mean(Factor_LackFollowThrough, na.rm=T),
                                  sd = sd(Factor_LackFollowThrough, na.rm=T),
                                  range = range(Factor_LackFollowThrough, na.rm=T))

# Age related to ERI? 
cor.test(pcgs_3factor.clean$Age, pcgs_3factor.clean$Factor_PervInf)
cor.test(pcgs_3factor.clean$Age, pcgs_3factor.clean$Factor_FeelingsTrigger)
cor.test(pcgs_3factor.clean$Age, pcgs_3factor.clean$Factor_LackFollowThrough)

#### Gender and age related to ERI? ####
pcgs_3factor.clean.gender <- pcgs_3factor.clean1
pcgs_3factor.clean.gender$Gender <- as.factor(pcgs_3factor.clean.gender$Gender)

rstatix::anova_summary(aov(Factor_PervInf ~ Gender + Age, data = pcgs_3factor.clean.gender), effect.size = "pes")
p_i_ag1 <- c(0.484, 0.076)
p.adjust(p_i_ag1, method = "fdr")

rstatix::anova_summary(aov(Factor_FeelingsTrigger ~ Gender + Age, data = pcgs_3factor.clean.gender), effect.size = "pes")
p_i_ag2 <- c(0.903, 0.325)
p.adjust(p_i_ag2, method = "fdr")

rstatix::anova_summary(aov(Factor_LackFollowThrough ~ Gender + Age, data = pcgs_3factor.clean.gender), effect.size = "pes")
p_i_ag3 <- c(0.711, 0.500)
p.adjust(p_i_ag3, method = "fdr")

#### Check incidence of PCGS ####
table(pcgs_3factor.clean$PCGS_LH) 
#92/120 (76%)
prop.test(x = 92, n = 120, p = 0.5,correct = TRUE)

table(pcgs_3factor.clean$PCGS_RH)
#68/120 (56%)

prop.test(x = 68, n = 120, p = 0.5, correct = TRUE)

prop.test(x = c(92, 68), n = c(120, 120), correct = TRUE)

pcgs_hemi.tab <- table(pcgs_3factor.clean$PCGS_LH, pcgs_3factor.clean$PCGS_RH)
mcnemar.test(pcgs_hemi.tab)

pcgs_incidence_plot.df <- pcgs_3factor.clean %>% dplyr::select(ID, PCGS_LH, PCGS_RH, Gender, Age, asymmetry) %>% 
  pivot_longer(cols = c("PCGS_LH", "PCGS_RH", "asymmetry"), 
               names_to = "Hemisphere", values_to = "PCGS_presence") %>%
  group_by(Hemisphere, PCGS_presence) %>% 
  summarise(n = n(), 
            incidence = n/120*100) 
  # mutate(group = case_when(Hemisphere == "PCGS_LH" ~ "hemisphere", 
  #                          Hemisphere == "PCGS_RH" ~ "hemisphere", 
  #                          Hemisphere == "asymmetry" ~ ""))
pcgs_incidence_plot.df$Hemisphere <- gsub("PCGS_LH", "left", pcgs_incidence_plot.df$Hemisphere)
pcgs_incidence_plot.df$Hemisphere <- gsub("PCGS_RH", "right", pcgs_incidence_plot.df$Hemisphere)
pcgs_incidence_plot.df$PCGS_presence <- gsub("leftward", "lw", pcgs_incidence_plot.df$PCGS_presence)
pcgs_incidence_plot.df$PCGS_presence <- gsub("rightward", "rw", pcgs_incidence_plot.df$PCGS_presence)
pcgs_incidence_plot.df$PCGS_presence <- gsub("symmetrical", "sy", pcgs_incidence_plot.df$PCGS_presence)
pcgs_incidence_plot.df$PCGS_presence <- factor(pcgs_incidence_plot.df$PCGS_presence, levels = c("absent", "present","rw", "lw", "sy"))
pcgs_incidence_plot.df$Hemisphere <- factor(pcgs_incidence_plot.df$Hemisphere, levels = c("left", "right", "asymmetry"))

pcgs_incidence_plot <- pcgs_incidence_plot.df %>% 
ggplot(aes(x = Hemisphere, y = incidence, fill = PCGS_presence)) +
  
geom_col(color = "black", width = .75) +
  
# set color palette
scale_fill_manual(breaks = c("present", "absent","rw", "lw", "sy"),
                  values = c("#a50f15", "#fee5d9", "#FFFFFF", "#ca0020", "#fb6a4a")) +
  
  scale_y_continuous(n.breaks = 5, limits = c(0,100)) +
  
  # set theme
  project_theme2 +
  
  guides(fill = guide_legend(override.aes = list(size = .5))) +
  #guides(fill = "none") +
  
  # set axis labels
  labs(x = "Hemisphere",
       y = "Incidence") 

### view
pcgs_incidence_plot

ggplot2::ggsave(filename = "~/Desktop/pcgs_project/pcgs_incidence_plot2.png",
                plot = pcgs_incidence_plot,
                device = "png",
                width = 4,
                height = 3,
                units = "in",
                dpi = "retina")

#### Is age/gender related to PCGS incidence ####
pcgs_3factor.clean.gender2 <- pcgs_3factor.clean.gender %>% dplyr::select(ID, PCGS_LH, PCGS_RH, asymmetry, Gender, Age) %>% 
  mutate(PCGS_LH_presence_binary = case_when(PCGS_LH == "absent" ~ 0, 
                                             PCGS_LH == "present" ~ 1),
         PCGS_RH_presence_binary = case_when(PCGS_RH == "absent" ~ 0, 
                                             PCGS_RH == "present" ~ 1)
  )

pcgs_pres_age.aov <- aov(Age ~ PCGS_LH + PCGS_RH, pcgs_3factor.clean.gender2)
summary(pcgs_pres_age.aov)
rstatix::anova_summary(pcgs_pres_age.aov, effect.size = "pes")

pcgs_pattern_age.aov <- aov(Age ~ asymmetry, pcgs_3factor.clean.gender2)
summary(pcgs_pattern_age.aov)
rstatix::anova_summary(pcgs_pattern_age.aov, effect.size = "pes")

p_pcgs_a <- c(0.597, 0.535, 0.965)
p.adjust(p_pcgs_a, method = "fdr")


lh_pcgs.tab <- table(pcgs_3factor.clean.gender$PCGS_LH, pcgs_3factor.clean.gender$Gender)
lh_pcgs.tab
chisq.test(lh_pcgs.tab)

rh_pcgs.tab <- table(pcgs_3factor.clean.gender$PCGS_RH, pcgs_3factor.clean.gender$Gender)
rh_pcgs.tab
chisq.test(rh_pcgs.tab)

pcgs_asy.tab <- table(pcgs_3factor.clean.gender$asymmetry, pcgs_3factor.clean.gender$Gender)
pcgs_asy.tab
chisq.test(pcgs_asy.tab)


p_pcgs_g <- c(0.7753, 0.2973, 0.8735)
p.adjust(p_pcgs_g, method = "fdr")


#### Analysis of Pcgs presence ~ disorder status ####
pcgs_3factor.clean2 <- pcgs_3factor.clean2 %>% mutate(lifetime_MDE_dx2 = case_when(lifetime_MDE_dx == 0 ~ "no", 
                                                                                   lifetime_MDE_dx == 1 ~ "yes",
                                                                                   lifetime_MDE_dx == 2 ~ "no",
                                                                                   lifetime_MDE_dx == NA ~ "no"
))
pcgs_3factor.clean2$lifetime_MDE_dx2[is.na(pcgs_3factor.clean2$lifetime_MDE_dx2)] <- "no"
pcgs_3factor.clean2 %>% group_by(PCGS_LH, lifetime_MDE_dx2) %>% summarise(n = n(), "perc" = n/120*100)
chisq.test(pcgs_3factor.clean2$PCGS_LH, pcgs_3factor.clean2$lifetime_MDE_dx2)

pcgs_3factor.clean2 %>% group_by(PCGS_RH, lifetime_MDE_dx2) %>% summarise(n = n(), "perc" = n/120*100)
chisq.test(pcgs_3factor.clean2$PCGS_RH, pcgs_3factor.clean2$lifetime_MDE_dx2)

pcgs_3factor.clean2 %>% group_by(asymmetry, lifetime_MDE_dx2) %>% summarise(n = n(), "perc" = n/120*100)
chisq.test(pcgs_3factor.clean2$asymmetry, pcgs_3factor.clean2$lifetime_MDE_dx2)

p_pcgs_mdd_dx <- c(0.6346, 0.8029, 0.4171)
p.adjust(p_pcgs_mdd_dx, method = "fdr")

pcgs_3factor.clean2 <- pcgs_3factor.clean2 %>% mutate(Lifetime_anxiety_dx2 = case_when(Lifetime_anxiety_dx == 0 ~ "no", 
                                                                                       Lifetime_anxiety_dx == 1 ~ "yes"
))

pcgs_3factor.clean2$Lifetime_anxiety_dx2[is.na(pcgs_3factor.clean2$Lifetime_anxiety_dx2)] <- "no"
chisq.test(pcgs_3factor.clean2$PCGS_LH, pcgs_3factor.clean2$Lifetime_anxiety_dx2)
pcgs_3factor.clean2 %>% group_by(PCGS_LH, Lifetime_anxiety_dx2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$PCGS_RH, pcgs_3factor.clean2$Lifetime_anxiety_dx2)
pcgs_3factor.clean2 %>% group_by(PCGS_RH, Lifetime_anxiety_dx2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$asymmetry, pcgs_3factor.clean2$Lifetime_anxiety_dx2)
pcgs_3factor.clean2 %>% group_by(asymmetry, Lifetime_anxiety_dx2) %>% summarise(n = n())

pcgs_3factor.clean2 <- pcgs_3factor.clean2 %>% mutate(lifetimeAUD2 = case_when(lifetimeAUD == 0 ~ "no", 
                                                                               lifetimeAUD == 1 ~ "yes"
))

pcgs_3factor.clean2$lifetimeAUD2[is.na(pcgs_3factor.clean2$lifetimeAUD2)] <- "no"
chisq.test(pcgs_3factor.clean2$PCGS_LH, pcgs_3factor.clean2$lifetimeAUD2)
pcgs_3factor.clean2 %>% group_by(PCGS_LH, lifetimeAUD2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$PCGS_RH, pcgs_3factor.clean2$lifetimeAUD2)
pcgs_3factor.clean2 %>% group_by(PCGS_RH, lifetimeAUD2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$asymmetry, pcgs_3factor.clean2$lifetimeAUD2)
pcgs_3factor.clean2 %>% group_by(asymmetry, lifetimeAUD2) %>% summarise(n = n())


pcgs_3factor.clean2 <- pcgs_3factor.clean2 %>% mutate(lifetimeSUD2 = case_when(lifetimeSUD == 0 ~ "no", 
                                                                               lifetimeSUD == 1 ~ "yes"
))

pcgs_3factor.clean2$lifetimeSUD2[is.na(pcgs_3factor.clean2$lifetimeSUD2)] <- "no"

chisq.test(pcgs_3factor.clean2$PCGS_LH, pcgs_3factor.clean2$lifetimeSUD2)
pcgs_3factor.clean2 %>% group_by(PCGS_LH, lifetimeSUD2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$PCGS_RH, pcgs_3factor.clean2$lifetimeSUD2)
pcgs_3factor.clean2 %>% group_by(PCGS_RH, lifetimeSUD2) %>% summarise(n = n())

chisq.test(pcgs_3factor.clean2$asymmetry, pcgs_3factor.clean2$lifetimeSUD2)
pcgs_3factor.clean2 %>% group_by(asymmetry, lifetimeSUD2) %>% summarise(n = n())

#### PCGS pres/abs ~ ERI #####
pcgs_pa_PIF.lm <- lm(Factor_PervInf ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean)
pcgs_pa_PIF.lm.aov <- Anova(pcgs_pa_PIF.lm, type = "III")
pcgs_pa_PIF.lm.aov
rstatix::anova_summary(pcgs_pa_PIF.lm.aov, effect.size = "pes")

pcgs_pa_FTA.lm <- lm(Factor_FeelingsTrigger ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean)
pcgs_pa_FTA.lm.aov <- Anova(pcgs_pa_FTA.lm, type = "III")
pcgs_pa_FTA.lm.aov
rstatix::anova_summary(pcgs_pa_FTA.lm.aov, effect.size = "pes")

pcgs_pa_LFT.lm <- lm(Factor_LackFollowThrough ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean)
pcgs_pa_LFT.lm.aov <- Anova(pcgs_pa_LFT.lm, type = "III")
pcgs_pa_LFT.lm.aov
rstatix::anova_summary(pcgs_pa_LFT.lm.aov, effect.size = "pes")
observed_F_A <- pcgs_pa_LFT.lm.aov$`F value`[2] 

pcgs_pa_gender_LFT.lm <- lm(Factor_LackFollowThrough ~ PCGS_LH*Gender, data = pcgs_3factor.clean1)
pcgs_pa_gender_LFT.lm.aov <- Anova(pcgs_pa_gender_LFT.lm, type = "III")
pcgs_pa_gender_LFT.lm.aov
rstatix::anova_summary(pcgs_pa_gender_LFT.lm.aov, effect.size = "pes")

p_1 <- c(0.147, 0.9460000, 0.147, 0.991, 0.000003302, 0.9460000)
p.adjust(p_1, method = "fdr")

set.seed(1)
perm_model_1 <- lmp(Factor_LackFollowThrough ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean, maxIter = 1000, perm = "Prob")
anova(perm_model_1)

set.seed(1)
perm.anova(Factor_LackFollowThrough ~ PCGS_LH, data = pcgs_3factor.clean, nperm = 2000,
           progress = TRUE)

# Create a data frame for ggplot2
df <- data.frame(F_statistic = perm_F_A)

# Plotting with ggplot2
ggplot(df, aes(x = F_statistic)) +
  geom_histogram(binwidth = 0.2, fill = "lightblue", color = "white") +
  geom_vline(aes(xintercept = observed_F_A), color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = observed_F_A, y = max(hist(perm_F_A, plot = FALSE)$counts) + 1, 
           label = paste("Observed F =", round(observed_F_A, 2)), color = "red", hjust = -0.1) +
  labs(title = "Permutation Test for Factor A",
       x = "F-statistic",
       y = "Frequency") +
  theme_minimal()

pcgs_3factor.clean.stats <- pcgs_3factor.clean %>% 
  group_by(PCGS_LH) %>% 
  summarise(
    sd = sd(Factor_LackFollowThrough, na.rm = T), 
    n = n(),
    Factor_LackFollowThrough = mean(Factor_LackFollowThrough, na.rm = T)) %>% 
  mutate(se = sd/sqrt(n),
         )
pcgs_3factor.clean.stats
# Increase = New Number - Original Number
# % increase = Increase ÷ Original Number × 100.
# Read more at: https://www.skillsyouneed.com/num/percent-change.html

(3.29-2.52)/2.52

# plot 
pcgs_presence_plot <- pcgs_3factor.clean %>% 
  ggplot(aes(x = PCGS_LH, y = Factor_LackFollowThrough, fill = PCGS_LH)) +
  
  # geom_point(position = position_jitter(width = .2),
  #            shape = 21,
  #            size = 1) +
  # 
  # geom_violin(alpha = .5) +
  # 
  # geom_boxplot(alpha = .75, outlier.alpha = 0, width = .2) +
  
  geom_flat_violin(aes(fill = PCGS_LH),
                   position = position_nudge(x = .1, y = 0),
                   #adjust = 1.0,
                   alpha = 1,
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = PCGS_LH, y = Factor_LackFollowThrough, fill = PCGS_LH),
             position = position_jitter(width = .025),
             shape = 21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_boxplot(position = position_nudge(x = .1, y = 0),
               lwd = .5,
               width = .1, 
               alpha = 1,
               outlier.alpha = 0) +
  
  # geom_pointrange(data = pcgs_asy.clean.stats,
  #                 aes(asymmetry, Factor_LackFollowThrough,
  #                     ymin=Factor_LackFollowThrough-se,
  #                     ymax=Factor_LackFollowThrough+se,
  #                     fill = asymmetry),
  #                 position = position_nudge(x = .1, y = 0),
  #                 lwd= 1,
  #                 shape = 21,
  #                 size = 1) +
  
  # set color palette
scale_fill_manual(breaks = c("present", "absent"),
                  values = c("#ca0020", "#fee5d9")) +
  
  # geom_bar(data = pcgs_3factor.clean.stats,
  #          aes(x=PCGS_LH, y=Factor_LackFollowThrough),
  #          stat="identity", color = "black", alpha = .25, width = .75) +
  # 
  # geom_errorbar(data = pcgs_3factor.clean.stats,
  #               aes(x=PCGS_LH, ymin=Factor_LackFollowThrough-se, ymax=Factor_LackFollowThrough+se),
  #               width=0, color="black", stat="identity", alpha = 1, linewidth = 1.5) +
  
  scale_y_continuous(n.breaks = 5, limits = c(0.5,5)) +
  
  #scale_fill_brewer(palette = "Dark2") + 
  
  # set theme
  project_theme3 +
  
  coord_flip() +
  
  guides(fill = "none") +
  
  # set axis labels
  labs(x = "Left PCGS",
       y = "LFT Severity") 

### view
pcgs_presence_plot

ggplot2::ggsave(filename = "~/Desktop/pcgs_project/pcgs_presence_plot.png",
                plot = pcgs_presence_plot,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")


# asymmetry analyses 
# slight alterations to data based on how groups run these (a)symmetry analyses 
pcgs_asy_PIF.lm <- lm(Factor_PervInf ~ asymmetry, data = pcgs_3factor.clean)
pcgs_asy_PIF.lm.aov <- Anova(pcgs_asy_PIF.lm, type = "III")
pcgs_asy_PIF.lm.aov
rstatix::anova_summary(pcgs_asy_PIF.lm.aov, effect.size = "pes")

pcgs_asy_FTA.lm <- lm(Factor_FeelingsTrigger ~ asymmetry, data = pcgs_3factor.clean)
pcgs_asy_FTA.lm.aov <- Anova(pcgs_asy_FTA.lm, type = "III")
pcgs_asy_FTA.lm.aov
rstatix::anova_summary(pcgs_asy_FTA.lm.aov, effect.size = "pes")

pcgs_asy_LFT.lm <- lm(Factor_LackFollowThrough ~ asymmetry, data = pcgs_3factor.clean)
pcgs_asy_LFT.lm.aov <- Anova(pcgs_asy_LFT.lm, type = "III")
pcgs_asy_LFT.lm.aov
rstatix::anova_summary(pcgs_asy_LFT.lm.aov, effect.size = "pes")

pcgs_asy_gender_LFT.lm <- lm(Factor_LackFollowThrough ~ asymmetry*Gender, data = pcgs_3factor.clean1)
pcgs_asy_gender_LFT.lm.aov <- Anova(pcgs_asy_gender_LFT.lm, type = "III")
pcgs_asy_gender_LFT.lm.aov
rstatix::anova_summary(pcgs_asy_gender_LFT.lm.aov, effect.size = "pes")

p_2 <- c(0.457, 0.653, 0.00001032)
p.adjust(p_2, method = "fdr")

pcgs_pa_LFT.lm.AIC <- AIC(lm(Factor_LackFollowThrough ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean))
pcgs_pa_PIF.lm.AIC <- AIC(lm(Factor_PervInf ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean))
pcgs_pa_FTA.lm.AIC <- AIC(lm(Factor_FeelingsTrigger ~ PCGS_LH + PCGS_RH, data = pcgs_3factor.clean))

pcgs_pa_PIF.lm.AIC-pcgs_pa_LFT.lm.AIC
pcgs_pa_FTA.lm.AIC-pcgs_pa_LFT.lm.AIC

pcgs_asy_PIF.lm.AIC <- AIC(pcgs_asy_PIF.lm)
pcgs_asy_FTA.lm.AIC <- AIC(pcgs_asy_FTA.lm)
pcgs_asy_LFT.lm.AIC <- AIC(pcgs_asy_LFT.lm)

pcgs_asy_PIF.lm.AIC-pcgs_asy_LFT.lm.AIC
pcgs_asy_FTA.lm.AIC-pcgs_asy_LFT.lm.AIC

set.seed(1)
perm_model_2 <- lmp(Factor_LackFollowThrough ~ asymmetry, data = pcgs_3factor.clean, maxIter = 1000, perm = "Prob")
anova(perm_model_2)

set.seed(1)
perm.anova(Factor_LackFollowThrough ~ asymmetry, data = pcgs_3factor.clean, nperm = 2000,
           progress = TRUE)

pcgs_asy_LFT.lm.aov.m1 <- emmeans::emmeans(pcgs_asy_LFT.lm, ~ asymmetry)
emmeans::contrast(pcgs_asy_LFT.lm.aov.m1, method='pairwise', adjust = 'none')

pcgs_3factor.clean2 <- pcgs_3factor.clean
pcgs_3factor.clean2$asymmetry <- gsub("leftward", "lw", pcgs_3factor.clean2$asymmetry)
pcgs_3factor.clean2$asymmetry <- gsub("rightward", "rw", pcgs_3factor.clean2$asymmetry)
pcgs_3factor.clean2$asymmetry <- gsub("symmetrical", "sy", pcgs_3factor.clean2$asymmetry)

pcgs_asy.clean.stats <- pcgs_3factor.clean2 %>% 
  group_by(asymmetry) %>% 
  summarise(
    sd = sd(Factor_LackFollowThrough, na.rm = T), 
    n = n(),
    Factor_LackFollowThrough = mean(Factor_LackFollowThrough, na.rm = T)) %>% 
  mutate(se = sd/sqrt(n),
  )

pcgs_asy_plot <- pcgs_3factor.clean2 %>% 
  ggplot(aes(x = asymmetry, y = Factor_LackFollowThrough, fill = asymmetry)) +

  geom_flat_violin(aes(fill = asymmetry),
                   position = position_nudge(x = .1, y = 0),
                   #adjust = 1.0,
                   alpha = 1,
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = asymmetry, y = Factor_LackFollowThrough, fill = asymmetry),
             position = position_jitter(width = .025),
             shape = 21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_boxplot(position = position_nudge(x = .1, y = 0),
               lwd = .5,
               width = .1, 
               alpha = 1,
               outlier.alpha = 0) +
  # geom_pointrange(data = pcgs_asy.clean.stats,
  #                 aes(asymmetry, Factor_LackFollowThrough,
  #                     ymin=Factor_LackFollowThrough-se,
  #                     ymax=Factor_LackFollowThrough+se,
  #                     fill = asymmetry),
  #                 position = position_nudge(x = .1, y = 0),
  #                 lwd= 1,
  #                 shape = 21,
  #                 size = 1) +
  
  # set color palette
  scale_fill_manual(breaks = c("lw", "rw", "sy"),
                    values = c("#ca0020", "#FFFFFF", "#fb6a4a")) +
  
  # geom_point(position = position_jitter(width = .2),
  #            shape = 21,
  #            size = 1) +
  # 
  # geom_violin(alpha = .5) +
  # 
  # geom_boxplot(alpha = .75, outlier.alpha = 0, width = .2) +
  
  # geom_bar(data = pcgs_asy.clean.stats,
  #          aes(x = PCGS_Asymmetry, y = Factor_LackFollowThrough),
  #          stat="identity", color = "black", alpha = .25, width = .75) +
  # 
  # geom_errorbar(data = pcgs_asy.clean.stats,
  #               aes(x = PCGS_Asymmetry, ymin = Factor_LackFollowThrough-se, ymax = Factor_LackFollowThrough+se),
  #               width=0, color="black", stat="identity", alpha = 1, linewidth = 1.5) +
  
  scale_y_continuous(n.breaks = 5, limits = c(0.5,5)) +
  
  #scale_fill_brewer(palette = "Dark2") + 
  
  # set theme
  project_theme3 +
  
  coord_flip() +
  
  guides(fill = "none") +
   
  # set axis labels
  labs(x = "PCGS Asymmetry",
       y = "LFT Severity") 

### view
pcgs_asy_plot

ggplot2::ggsave(filename = "~/Desktop/pcgs_project/pcgs_asy_plot.png",
                plot = pcgs_asy_plot,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")

# see if sample size differences between groups alters result 
pcgs_3factor.clean.lw <- pcgs_3factor.clean %>% subset(asymmetry == 'leftward')
pcgs_3factor.clean.lw <- pcgs_3factor.clean.lw$Factor_LackFollowThrough

pcgs_3factor.clean.rw <- pcgs_3factor.clean %>% subset(asymmetry == 'rightward')
pcgs_3factor.clean.rw <- pcgs_3factor.clean.rw$Factor_LackFollowThrough

pcgs_3factor.clean.sy <- pcgs_3factor.clean %>% subset(asymmetry == 'symmetrical')
pcgs_3factor.clean.sy <- pcgs_3factor.clean.sy$Factor_LackFollowThrough

var.test(pcgs_3factor.clean.lw, pcgs_3factor.clean.rw)
t.test(pcgs_3factor.clean.lw, pcgs_3factor.clean.rw, var.equal = T)
cohens_pimfs1 <- cohen.d(pcgs_3factor.clean.lw, pcgs_3factor.clean.rw, method = "cohen's d")
cohens_pimfs1

var.test(pcgs_3factor.clean.sy, pcgs_3factor.clean.rw)
t.test(pcgs_3factor.clean.sy, pcgs_3factor.clean.rw, var.equal = T)
cohens_pimfs2 <- cohen.d(pcgs_3factor.clean.sy, pcgs_3factor.clean.rw, method = "cohen's d")
cohens_pimfs2

set.seed(1)
B <- 1000
t.vect <- vector(length = B)
p.vect <- vector(length = B)
c.vect <- vector(length = B)
for(i in 1:B){
  boot.c <- sample(pcgs_3factor.clean.lw, size = 15, replace = F)
  boot.p <- pcgs_3factor.clean.rw
  
  ttest <- t.test(boot.c, boot.p, var.equal = F)
  effect_size <- cohen.d(boot.c, boot.p, method = "cohen's d")
  t.vect[i] <- ttest$statistic
  p.vect[i] <- ttest$p.value
  c.vect[i] <- effect_size$estimate
}


## Plot iterated t and p values

### Prepare data
t.vect.df <- as.data.frame(t.vect)
p.vect.df <- as.data.frame(p.vect)
c.vect.df <- as.data.frame(c.vect)

### T-value plot
DescTools::MedianCI(c.vect.df$c.vect, conf.level = 0.95)

c_dist_plot <- c.vect.df %>% ggplot(aes(x = c.vect)) + 
  
  # view distribution of values
  geom_density(color = "black", fill = '#ca0020') + 
  
  # view upper bound as a vertical line
  geom_vline(xintercept = 1.31, linetype = "dashed") +
  
  # view mean as a vertical line
  geom_vline(xintercept = 1.33, linetype = "solid") +
  
  # view lower bound as a vertical line
  geom_vline(xintercept = 1.35, linetype = "dashed") +
  
  # view 0 as a vertical line
  #geom_vline(xintercept = 0, linetype = "solid", color = "red") +
  
  # set axis scales
  scale_x_continuous(n.breaks = 5, limits = c(0, 2.5)) +
  scale_y_continuous(n.breaks = 8) + 
  
  labs(x = "lw vs.rw PCGS Asymmetry Cohen's d",
       y = "Density") +

  
  # set theme
  project_theme2

### view
c_dist_plot

ggplot2::ggsave(filename = "~/Desktop/pcgs_project/lw_rw_c_dist_plot.png",
                plot = c_dist_plot,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")


set.seed(1)
B2 <- 1000
t.vect2 <- vector(length = B2)
p.vect2 <- vector(length = B2)
c.vect2 <- vector(length = B2)
for(i in 1:B2){
  boot.c2 <- sample(pcgs_3factor.clean.sy, size = 15, replace = F)
  boot.p2 <- pcgs_3factor.clean.rw
  
  ttest2 <- t.test(boot.c2, boot.p2, var.equal = F)
  effect_size2 <- cohen.d(boot.c2, boot.p2, method = "cohen's d")
  t.vect2[i] <- ttest2$statistic
  p.vect2[i] <- ttest2$p.value
  c.vect2[i] <- effect_size2$estimate
}


## Plot iterated t and p values

### Prepare data
t.vect.df2 <- as.data.frame(t.vect2)
p.vect.df2 <- as.data.frame(p.vect2)
c.vect.df2 <- as.data.frame(c.vect2)

### T-value plot
DescTools::MedianCI(c.vect.df2$c.vect2, conf.level = 0.95)

c_dist_plot2 <- c.vect.df2 %>% ggplot(aes(x = c.vect2)) + 
  
  # view distribution of values
  geom_density(color = "black", fill = '#fb6a4a') + 
  
  # view upper bound as a vertical line
  geom_vline(xintercept = 1.45, linetype = "dashed") +
  
  # view mean as a vertical line
  geom_vline(xintercept = 1.48, linetype = "solid") +
  
  # view lower bound as a vertical line
  geom_vline(xintercept = 1.50, linetype = "dashed") +
  
  # view 0 as a vertical line
  #geom_vline(xintercept = 0, linetype = "solid", color = "red") +
  
  # set axis scales
  scale_x_continuous(n.breaks = 5, limits = c(0, 2.5)) +
  scale_y_continuous(n.breaks = 8) + 
  
  labs(x = "sy vs.rw PCGS Cohen's d",
       y = "Density") +
  
  
  # set theme
  project_theme2

### view
c_dist_plot2

ggplot2::ggsave(filename = "~/Desktop/pcgs_project/sy_rw_c_dist_plot.png",
                plot = c_dist_plot2,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")


#### Length ####
pcgs_length.w <- pcgs_length %>% pivot_wider(names_from = hemi, values_from = max_path_length) %>% rename(pcgs_lh_length = lh, pcgs_rh_length = rh)
pcgs_length.w <- merge(pcgs_3factor.clean1, pcgs_length.w, by = "ID", all.x=T)
pcgs_length.w2 <- pcgs_length.w
pcgs_length.w2[is.na(pcgs_length.w2)] <- 0
pcgs_length.w2 <- pcgs_length.w2 %>% mutate(ratio = (pcgs_rh_length-pcgs_lh_length)/(pcgs_rh_length+pcgs_lh_length))
pcgs_length.w2[is.na(pcgs_length.w2)] <- 0

rstatix::anova_summary(aov(pcgs_lh_length ~ Age + as.factor(Gender), data = pcgs_length.w), effect.size = "pes")
p_l_ga1 <- c(0.698, 0.795)
p.adjust(p_l_ga1, method = "fdr")

rstatix::anova_summary(aov(pcgs_rh_length ~ Age + as.factor(Gender), data = pcgs_length.w), effect.size = "pes")
p_l_ga2 <- c(0.559, .793)
p.adjust(p_l_ga2, method = "fdr")

rstatix::anova_summary(aov(ratio ~ Age + as.factor(Gender), data = pcgs_length.w2), effect.size = "pes")
p_l_ga3 <- c(0.857, 0.671)
p.adjust(p_l_ga3, method = "fdr")

pcgs_length_demo2_gender <- merge(pcgs_length, urgency_demo, by = "ID")
pcgs_length_demo2_gender.mod <- nlme::lme(max_path_length ~ hemi*as.factor(Gender), random = (~1|ID/hemi), 
                                data = pcgs_length_demo2_gender)
anova(pcgs_length_demo2_gender.mod)

summary(lm(Factor_PervInf ~ pcgs_lh_length + pcgs_rh_length, pcgs_length.w))
summary(lm(Factor_FeelingsTrigger ~ pcgs_lh_length + pcgs_rh_length, pcgs_length.w))
summary(lm(Factor_LackFollowThrough ~ pcgs_lh_length + pcgs_rh_length, pcgs_length.w))

summary(lm(Factor_PervInf ~ ratio, pcgs_length.w2))
summary(lm(Factor_FeelingsTrigger ~ ratio, pcgs_length.w2))
summary(lm(Factor_LackFollowThrough ~ ratio, pcgs_length.w2))
p_length_hemi_imp <- c(.607, 0.839, 0.547, 0.732, 0.574, 0.300, 0.196, 0.419, 0.0124)
p.adjust(p_length_hemi_imp, method = "fdr")

#### lgi ####
acc_lgi <- acc_lgi %>% rename(ID = sub) 

# Assume your data is in a data frame called df, with continuous dependent variable 'y', 
# categorical variable 'cat1' (3 levels), and categorical variable 'cat2' (2 levels)
acc_lgi$label <- gsub("superiorfrontal", "SFC", acc_lgi$label)
acc_lgi$label <- gsub("rostralanteriorcingulate", "rACC", acc_lgi$label)
acc_lgi$label <- gsub("caudalanteriorcingulate", "cACC", acc_lgi$label)


t.test(cortical_thickness_mean ~ hemi, acc_lgi)
acc_lgi_lh <- acc_lgi %>% subset(hemi == "lh")
p_lh <- ggplot(acc_lgi_lh, aes(x = label, y = cortical_thickness_mean, fill = label)) +
  geom_jitter(width = .2, shape = 21) +
  geom_boxplot(width = .5, alpha = .75) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("SFC", "rACC"), 
                                        c("SFC", "cACC"),
                                        c("rACC", "cACC")),
                     label = "p.signif") +  
  scale_y_continuous(n.breaks = 5, limits = c(1.5,3)) +
  scale_fill_manual(breaks = c("SFC", "cACC", "rACC"),
                    values = c("yellow", "red", "blue")) +
  project_theme3 +
  
  guides(fill = "none") +
  labs(x = "Region",
       y = "lgi")
p_lh
ggplot2::ggsave(filename = "~/Desktop/DropBox/current_projects/pcgs_project/p_lh.png",
                plot = p_lh,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

acc_lgi_rh <- acc_lgi %>% subset(hemi == "rh")
p_rh <- ggplot(acc_lgi_rh, aes(x = label, y = cortical_thickness_mean, fill = label)) +
  geom_jitter(width = .2, shape = 21) +
  geom_boxplot(width = .5, alpha = .75) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("SFC", "rACC"), 
                                        c("SFC", "cACC"),
                                        c("rACC", "cACC")),
                     label = "p.signif") +  
  scale_y_continuous(n.breaks = 5, limits = c(1.5,3)) +
  scale_fill_manual(breaks = c("SFC", "cACC", "rACC"),
                    values = c("yellow", "red", "blue")) +
  project_theme3 +
  
  guides(fill = "none") +
  labs(x = "Region",
       y = "lgi")
p_rh

ggplot2::ggsave(filename = "~/Desktop/DropBox/current_projects/pcgs_project/p_rh.png",
                plot = p_rh,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

pcgs_length$hemi

length_plot <- ggplot(pcgs_length, aes(x = hemi, y = max_path_length, fill = hemi)) +
  geom_jitter(width = .2, shape = 21) +
  geom_boxplot(width = .5, alpha = .75, outlier.alpha = 0) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("lh", "rh")),
                     label = "p.signif") +  
  
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +

  #scale_y_continuous(n.breaks = 5, limits = c(1.5,3)) +
  scale_fill_manual(breaks = c("lh", "rh"),
                    values = c("black", "white")) +
  project_theme3 +
  
  coord_flip() +
  
  guides(fill = "none") +
  labs(x = "Length",
       y = "hemisphere")
length_plot

ggplot2::ggsave(filename = "~/Desktop/DropBox/current_projects/pcgs_project/length_plot.png",
                plot = length_plot,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")

depth_plot <- ggplot(pcgs_depth, aes(x = hemi, y = sulcal_depth_mm, fill = hemi)) +
  geom_jitter(width = .2, shape = 21) +
  geom_boxplot(width = .5, alpha = .75, outlier.alpha = 0) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("lh", "rh")),
                     label = "p.signif") +  
  
  geom_hline(yintercept = 4, linetype = "dashed", color = "red", linewidth = 1) +
  
  #scale_y_continuous(n.breaks = 5, limits = c(1.5,3)) +
  scale_fill_manual(breaks = c("lh", "rh"),
                    values = c("black", "white")) +
  project_theme3 +
  
  coord_flip() +
  
  guides(fill = "none") +
  labs(x = "hemisphere",
       y = "Depth")
depth_plot

ggplot2::ggsave(filename = "~/Desktop/DropBox/current_projects/pcgs_project/depth_plot.png",
                plot = depth_plot,
                device = "png",
                width = 5,
                height = 3,
                units = "in",
                dpi = "retina")

pcgs_depth %>% group_by(hemi) %>% summarise(mean = mean(sulcal_depth_mm), 
                                            sd = sd(sulcal_depth_mm),
                                            range = range(sulcal_depth_mm))


acc_lgi2 <- merge(acc_lgi, pcgs_labels, by = "ID")
acc_lgi2.lh <- acc_lgi2 %>% subset(hemi == "rh")

pplot_lgi_pcgs <- ggplot(acc_lgi2.lh, aes(x = PCGS_RH, y = cortical_thickness_mean, fill = PCGS_RH)) +
  geom_jitter(width = .2, shape = 21) +
  geom_boxplot(width = .5, alpha = .75) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("present", "absent")),
                     label = "p.signif") +  
  scale_y_continuous(n.breaks = 5, limits = c(1.5,3)) +
  scale_fill_manual(breaks = c("present", "absent"),
                    values = c("gray", "white")) +
  project_theme3 +
  
  facet_wrap(~ label) +
  
  guides(fill = "none") +
  labs(x = "Region",
       y = "lgi")
  pplot_lgi_pcgs

acc_lgi.wide <- acc_lgi %>% select(ID, hemi, label, cortical_thickness_mean) %>% pivot_wider(names_from = c('label', 'hemi'), values_from = cortical_thickness_mean) 
acc_lgi.wide <- merge(acc_lgi.wide, urgency_metrics,  by = "ID")
acc_lgi.wide <- merge(acc_lgi.wide, urgency_demo,  by = "ID")

acc_lgi.wide <- acc_lgi.wide %>% mutate(cACC_ratio = (cACC_rh-cACC_lh)/(cACC_rh+cACC_lh),
                                        sF_ratio = (SF_rh-SF_lh)/(SF_rh+SF_lh),
                                        rACC_ratio = (rACC_rh-rACC_lh)/(rACC_rh+rACC_lh)
) 

rstatix::anova_summary(aov(cACC_rh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_rh_cacc_adj <- c(0.000143, 0.059052)
p.adjust(p_rh_cacc_adj, method = "fdr")

rstatix::anova_summary(aov(cACC_lh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_lh_cacc_adj <- c(0.003, 0.019)
p.adjust(p_lh_cacc_adj, method = "fdr")

rstatix::anova_summary(aov(SF_rh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_rh_sf_adj <- c(0.000001, 0.053000)
p.adjust(p_rh_sf_adj, method = "fdr")

rstatix::anova_summary(aov(SF_lh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_lh_sf_adj <- c(0.000008, 0.034000)
p.adjust(p_lh_sf_adj, method = "fdr")

rstatix::anova_summary(aov(rACC_rh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_rh_dacc_adj <- c(0.000141, 0.123000)
p.adjust(p_rh_dacc_adj, method = "fdr")

rstatix::anova_summary(aov(rACC_lh ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_lh_dacc_adj <- c(0.000612, 0.314000)
p.adjust(p_lh_dacc_adj)

rstatix::anova_summary(aov(cACC_ratio ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_cACC_ratio <- c(0.228, 0.719)
p.adjust(p_cACC_ratio, method = "fdr")

rstatix::anova_summary(aov(sF_ratio ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_sF_ratio<- c(0.547, 0.495)
p.adjust(p_sF_ratio, method = "fdr")

rstatix::anova_summary(aov(rACC_ratio ~ Age + as.factor(Gender), acc_lgi.wide), effect.size = "pes")
p_rACC_ratio <- c(0.563, 0.664)
p.adjust(p_rACC_ratio, method = "fdr")

# behavioral analysis
summary_model_LFT <- summary(lm(
  Factor_LackFollowThrough ~ caudalanteriorcingulate_rh + caudalanteriorcingulate_lh + 
    superiorfrontal_rh + superiorfrontal_lh + 
    rostralanteriorcingulate_rh + rostralanteriorcingulate_lh,
  data = acc_lgi.wide
))
summary_model_LFT
summary_model_LFT <- tidy(summary_model_LFT) 
summary_model_LFT <-  summary_model_LFT %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant)
summary_model_LFT

# FTA ~ lGI
summary_model_FTA <- summary(lm(
  Factor_FeelingsTrigger ~ caudalanteriorcingulate_rh + caudalanteriorcingulate_lh + 
    superiorfrontal_rh + superiorfrontal_lh + 
    rostralanteriorcingulate_rh + rostralanteriorcingulate_lh,
  data = acc_lgi.wide
))
summary_model_FTA
summary_model_FTA <- tidy(summary_model_FTA) 
summary_model_FTA <-  summary_model_FTA %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant)
summary_model_FTA

# PIF ~ lGI
summary_model_PIF <- summary(lm(
  Factor_PervInf ~ caudalanteriorcingulate_rh + caudalanteriorcingulate_lh + 
    superiorfrontal_rh + superiorfrontal_lh + 
    rostralanteriorcingulate_rh + rostralanteriorcingulate_lh,
  data = acc_lgi.wide
))
summary_model_PIF
summary_model_PIF <- tidy(summary_model_PIF) 
summary_model_PIF <-  summary_model_PIF %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant)
summary_model_PIF

# LFT ~ lGI
summary_model_ratio_LFT <- summary(lm(
  Factor_LackFollowThrough ~ cACC_ratio + sF_ratio + rACC_ratio,
  data = acc_lgi.wide
))
summary_model_ratio_LFT
summary_model_ratio_LFT <- tidy(summary_model_ratio_LFT) 
summary_model_ratio_LFT <-  summary_model_ratio_LFT %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant) %>%
  mutate(fdr_p.value = format(fdr_p.value, scientific = F))
summary_model_ratio_LFT

# FTA ~ lGI
summary_model_ratio_FTA <- summary(lm(
  Factor_FeelingsTrigger ~ cACC_ratio + sF_ratio + rACC_ratio,
  data = acc_lgi.wide
))
summary_model_ratio_FTA
summary_model_ratio_FTA <- tidy(summary_model_ratio_FTA) 
summary_model_ratio_FTA <-  summary_model_ratio_FTA %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant) %>%
  mutate(fdr_p.value = format(fdr_p.value, scientific = F))
summary_model_ratio_FTA

# PIF ~ lGI
summary_model_ratio_PIF <- summary(lm(
  Factor_PervInf ~ cACC_ratio + sF_ratio + rACC_ratio,
  data = acc_lgi.wide
))
summary_model_ratio_PIF
summary_model_ratio_PIF <- tidy(summary_model_ratio_PIF) 
summary_model_ratio_PIF <-  summary_model_ratio_PIF %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(fdr_p.value < 0.05, T, F)) %>% dplyr::select(-p.value, -significant) %>%
  mutate(fdr_p.value = format(fdr_p.value, scientific = F))
summary_model_ratio_PIF
