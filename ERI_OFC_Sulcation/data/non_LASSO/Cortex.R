### Load Packages ###
library(stats) #anovas, correlations, etc etc
library(tidyverse) #data manipulation
library(ggplot2) #visualizations
library(afex) #rm-anovas if you want
library(emmeans) #post hoc pairwise comparisons
library(rstatix) #provides eta^2 (anova effect size) values and nice aov tables
library(RColorBrewer) #pretty colors for plots go here for actual color names: https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=5
library(plotly) #interactive plots
library(ggnewscale) #extra data vis
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(GGally)
library(repr)
library(ragg)
library(gridExtra)
library(readr)
library(purrr)
library(cli)


### Clean Data ###

setwd("C:/Users/willl/Documents/Urgency")

## Load Behavioral Data
factor <- read.csv("3Factor_impulsivity.csv")
names(factor)[names(factor) == 'ID'] <- 'sub'
names(factor)[names(factor) == 'Factor_FeelingsTrigger'] <- 'fta'
names(factor)[names(factor) == 'Factor_LackFollowThrough'] <- 'lft'
names(factor)[names(factor) == 'Factor_PervInf'] <- 'pif'
factor <- factor %>% select(c(1, 108, 113, 116))

## Load Sulcal Type and Count
urgency_type <- read.csv("urgency_type.csv")
names(urgency_type)[names(urgency_type) == 'ï..sub'] <- 'sub'
urgency_type$type <- as.factor(urgency_type$type)

## Sulcal Type Merge with Behavioral Data
factor_type <- merge(factor, urgency_type, by="sub")

# Clean data to check for Asymmetry

urgency_type.lh <- factor_type %>% subset(hemi=="lh") %>% 
  rename(pos.lh = pos, 
         ios.lh = ios, 
         sf.lh = sf,
         nak.lh = nak,
         tot.lh = tot)
urgency_type.rh <- factor_type %>% subset(hemi=="rh") %>% 
  rename(pos.rh = pos, 
         ios.rh = ios, 
         sf.rh = sf,
         nak.rh = nak,
         tot.rh = tot)

# Asymmetry between hemis
urgency_type.asym <- merge(urgency_type.lh,urgency_type.rh, by="sub")
urgency_type.asym <- urgency_type.asym %>% 
  mutate(pos_asym = pos.lh - pos.rh, 
         ios_asym = ios.lh - ios.rh, 
         sf_asym = sf.lh - sf.rh,)

## Importa Sulcal Measurements
sulc_metrics <- read.csv("urgency_metrics_new.csv")
sulc_metrics <- sulc_metrics %>% select(c(2:4,6, 8, 14:16))

sulc_metrics_colfs <- read.csv("urgency_metrics_colfs.csv")
sulc_metrics_colfs <- sulc_metrics_colfs %>% select(c(2:4,6, 8, 14:16))

## Behavioral & Sulcal Data
sulc_metrics_final <- merge(sulc_metrics, factor, by="sub")
sulc_metrics_final_colfs <- merge(sulc_metrics_colfs, factor, by="sub")


# Set plot window size
options(repr.plot.width=20, repr.plot.height=8)

# get current working directory
getwd()

# Set working Directory
setwd("~/LASSO")

### BEHAVIORAL ANALYSIS ###

# Except when comparing asymmetry, each analysis was done separately on each hemisphere #

## Pervasive Influence of Feelings

# ANOVA comparing Asymmetry in number of each sulci between hemispheres
summary(aov(pif.x ~ as.factor(pos_asym) + as.factor(ios_asym) + as.factor(sf_asym), data=urgency_type.asym))

# ANOVA comparing number of each Variable sulcus 
summary(aov(pif ~ as.factor(pos.lh) + as.factor(ios.lh) + as.factor(sf.lh), data=urgency_type.lh)) #sf almost and sig when ios and pos taken out
pif.sf <- aov(pif ~ as.factor(sf.lh), data=urgency_type.lh) #sf
emmeans.sf <- emmeans(pif.sf, ~ sf.lh)
contrast(emmeans.sf, method='pairwise', adjust='tukey') #1-2
summary(aov(pif ~ as.factor(pos.rh) + as.factor(ios.rh) + as.factor(sf.rh), data=urgency_type.rh))

# ANOVA of total number of variable sulci (Total/Nak)
summary(aov(pif ~ as.factor(tot.lh), data=urgency_type.lh))
summary(aov(pif ~ as.factor(nak.lh), data=urgency_type.lh))
summary(aov(pif ~ as.factor(tot.rh), data=urgency_type.rh))
summary(aov(pif ~ as.factor(nak.rh), data=urgency_type.rh))

# ANOVA of sulcogyral type 
summary(aov(pif ~ as.factor(type), data=urgency_type.lh))
summary(aov(pif ~ as.factor(type), data=urgency_type.rh))
pif.type <- aov(pif.x ~ as.factor(type.x)*as.factor(type.y), data=urgency_type.asym) # SIGNIFICANT INTERACTION EFFECT
summary(pif.type) 
emmeans.pif.type <- emmeans(pif.type, ~ as.factor(type.x)*as.factor(type.y))
contrast(emmeans.pif.type, method='pairwise', adjust='tukey')

## Feelings Trigger Action

# ANOVA comparing Asymmetry in number of each sulci between hemispheres
summary(aov(fta.x ~ as.factor(pos_asym) + as.factor(ios_asym) + as.factor(sf_asym), data=urgency_type.asym))

# ANOVA comparing number of each Variable sulcus 
summary(aov(fta ~ as.factor(pos.lh) + as.factor(ios.lh) + as.factor(sf.lh), data=urgency_type.lh))
summary(aov(fta ~ as.factor(pos.rh) + as.factor(ios.rh) + as.factor(sf.rh), data=urgency_type.rh)) # POS

(count.pif.aov.rh <- aov(fta ~ as.factor(pos.rh) + as.factor(ios.rh) + as.factor(sf.rh), data=urgency_type.rh)) #POS
summary(count.pif.aov.rh)
emmeans.pos <- emmeans(count.pif.aov.rh, ~ pos.rh)
contrast(emmeans.pos, method='pairwise', adjust='tukey') # DOESN'T SURVIVE TUKEY

# ANOVA of total number of variable sulci (Total/Nak)
summary(aov(fta ~ as.factor(tot.lh), data=urgency_type.lh))
summary(aov(fta ~ as.factor(nak.lh), data=urgency_type.lh))
summary(aov(fta ~ as.factor(tot.rh), data=urgency_type.rh))
summary(aov(fta ~ as.factor(nak.rh), data=urgency_type.rh))

# ANOVA of sulcogyral type 
fta.type <- aov(fta ~ as.factor(type), data=urgency_type.lh)
summary(fta.type)
emmeans.fta.type <- emmeans(fta.type, ~ type)
contrast(emmeans.fta.type, method='pairwise', adjust='tukey') # DOESN'T SURVIVE TUKEY

summary(aov(fta ~ as.factor(type), data=urgency_type.rh))
fta.type <- aov(fta.x ~ as.factor(type.x)*as.factor(type.y), data=urgency_type.asym)

## Lack Follow Through

# ANOVA comparing Asymmetry in number of each sulci between hemispheres
summary(aov(lft.x ~ as.factor(pos_asym) + as.factor(ios_asym) + as.factor(sf_asym), data=urgency_type.asym))

# ANOVA comparing number of each Variable sulcus 
summary(aov(lft ~ as.factor(pos.lh) + as.factor(ios.lh) + as.factor(sf.lh), data=urgency_type.lh))
summary(aov(lft ~ as.factor(pos.rh) + as.factor(ios.rh) + as.factor(sf.rh), data=urgency_type.rh))

# ANOVA of total number of variable sulci (Total/Nak)
summary(aov(lft ~ as.factor(tot.lh), data=urgency_type.lh))
summary(aov(lft ~ as.factor(nak.lh), data=urgency_type.lh))
summary(aov(lft ~ as.factor(tot.rh), data=urgency_type.rh))
summary(aov(lft ~ as.factor(nak.rh), data=urgency_type.rh))

# ANOVA of sulcogyral type 
summary(aov(lft ~ as.factor(type), data=urgency_type.lh))
summary(aov(lft ~ as.factor(type), data=urgency_type.rh))
summary(aov(lft.x ~ as.factor(type.x)*as.factor(type.y), data=urgency_type.asym))


## Refining data for Lasso analysis

# Code for each makes data wide and splits by hemisphere

# Regular analysis

# Standardized Max Depth
tertiary.depth.lasso.lh <-  sulc_metrics_final %>% subset(hemi == "lh" & label %in% c("olfs", "tolfs", "mosa", "mosp", "losa", "losp", 'tos', "sf")) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )

tertiary.depth.lasso.rh <- sulc_metrics_final %>% subset(hemi == "rh" & label %in% c("olfs", "tolfs", "mosa", "mosp", "losa", "losp", 'tos', 'sf')) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )

tertiary_depth_lasso <- rbind(tertiary.depth.lasso.lh, tertiary.depth.lasso.rh)
tertiary_depth_lasso <- na.omit(tertiary_depth_lasso)
write.csv(tertiary_depth_lasso, "C:/Users/willl/Documents/tertiary_depth_lasso.csv")


# Cortical Thickness

tertiary.ct.lasso.lh <-  sulc_metrics_final %>% subset(hemi == "lh" & label %in% c("colfs", "mosa", "mosp", "losa", "losp", 'tos', "sf")) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )

tertiary.ct.lasso.rh <- sulc_metrics_final %>% subset(hemi == "rh" & label %in% c("sf", "colfs", "mosa", "mosp", "losa", "losp", 'tos')) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )

tertiary_ct_lasso <- rbind(tertiary.ct.lasso.lh, tertiary.ct.lasso.rh)
tertiary_ct_lasso <- na.omit(tertiary_ct_lasso)
write.csv(tertiary_ct_lasso, "C:/Users/willl/Documents/tertiary_ct_lasso.csv")

# Colfs comparison

# Standardized Max Depth
tertiary.depth.lasso.colfs.lh <-  sulc_metrics_final_colfs %>% subset(hemi == "lh" & label %in% c("colfs", "mosa", "mosp", "losa", "losp", 'tos', "sf")) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )

tertiary.depth.lasso.colfs.rh <- sulc_metrics_final_colfs %>% subset(hemi == "rh" & label %in% c("colfs", "mosa", "mosp", "losa", "losp", 'tos', 'sf')) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )

tertiary_depth_lasso_colfs <- rbind(tertiary.depth.lasso.colfs.lh, tertiary.depth.lasso.colfs.rh)
tertiary_depth_lasso_colfs <- na.omit(tertiary_depth_lasso_colfs)
write.csv(tertiary_depth_lasso_colfs, "C:/Users/willl/Documents/tertiary_depth_lasso_colfs.csv")

# Cortical Thickness

tertiary.ct.lasso.colfs.lh <-  sulc_metrics_final %>% subset(hemi == "lh" & label %in% c("colfs", "mosa", "mosp", "losa", "losp", 'tos', "sf")) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )

tertiary.ct.lasso.colfs.rh <- sulc_metrics_final %>% subset(hemi == "rh" & label %in% c("sf", "colfs", "mosa", "mosp", "losa", "losp", 'tos')) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )

tertiary_ct_lasso_colfs <- rbind(tertiary.ct.lasso.colfs.lh, tertiary.ct.lasso.colfs.rh)
tertiary_ct_lasso_colfs <- na.omit(tertiary_ct_lasso_colfs)
write.csv(tertiary_ct_lasso_colfs, "C:/Users/willl/Documents/tertiary_ct_lasso_colfs.csv")

### Variable Sulci Quantitative ###

# Looking for trends between variable sulci and ERI, as LASSO is incompatible with variable sulci
# (Note: LASSO included SF, as it was present in the vast majority of hemispheres)

## Refining data for Lasso analysis

# Code for each makes data wide and splits by hemisphere

# Making Data Longer for Corr tests

# Standardized Max Depth
variable_sulci_depth.lh <-  sulc_metrics_final %>% subset(hemi == "lh" & label %in% c("ios", "pos", "sf")) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )

variable_sulci_depth.rh <-  sulc_metrics_final %>% subset(hemi == "rh" & label %in% c("ios", "pos", "sf")) %>%
  select(
    c("sub", "hemi", "label", "sulcal_depth_mean_pct", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = sulcal_depth_mean_pct, # value columns(s) that hold data for each category column
  )


# Cortical Thickness

variable_sulci_ct.lh <-  sulc_metrics_final %>% subset(hemi == "lh" & label %in% c("ios", "pos", "sf")) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>%
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )

variable_sulci_depth.rh <-  sulc_metrics_final %>% subset(hemi == "rh" & label %in% c("ios", "pos", "sf")) %>%
  select(
    c("sub", "hemi", "label", "cortical_thickness_mean", "pif", "lft", "fta")
  ) %>% 
  pivot_wider(
    names_from = label, # category column(s) to pivot from long to wide
    values_from = cortical_thickness_mean, # value columns(s) that hold data for each category column
  )


## Standardized Max Depth Mean

# pif
summary(lm(variable_sulci_depth.lh$pif ~ variable_sulci_depth.lh$ios)) #TREND
summary(lm(variable_sulci_depth.lh$pif ~ variable_sulci_depth.lh$pos))
summary(lm(variable_sulci_depth.lh$pif ~ variable_sulci_depth.lh$sf))

summary(lm(variable_sulci_depth.rh$pif ~ variable_sulci_depth.rh$ios))
summary(lm(variable_sulci_depth.rh$pif ~ variable_sulci_depth.rh$pos))
summary(lm(variable_sulci_depth.rh$pif ~ variable_sulci_depth.rh$sf))

# fta
summary(lm(variable_sulci_depth.lh$fta ~ variable_sulci_depth.lh$ios))
summary(lm(variable_sulci_depth.lh$fta ~ variable_sulci_depth.lh$pos))
summary(lm(variable_sulci_depth.lh$fta ~ variable_sulci_depth.lh$sf)) #SIG

summary(lm(variable_sulci_depth.rh$fta ~ variable_sulci_depth.rh$ios))
summary(lm(variable_sulci_depth.rh$fta ~ variable_sulci_depth.rh$pos))
summary(lm(variable_sulci_depth.rh$fta ~ variable_sulci_depth.rh$sf))

# lft
summary(lm(variable_sulci_depth.lh$lft ~ variable_sulci_depth.lh$ios))
summary(lm(variable_sulci_depth.lh$lft ~ variable_sulci_depth.lh$pos))
summary(lm(variable_sulci_depth.lh$lft ~ variable_sulci_depth.lh$sf))

summary(lm(variable_sulci_depth.rh$lft ~ variable_sulci_depth.rh$ios))
summary(lm(variable_sulci_depth.rh$lft ~ variable_sulci_depth.rh$pos))
summary(lm(variable_sulci_depth.rh$lft ~ variable_sulci_depth.rh$sf))

## Cortical Thickness

# pif
summary(lm(variable_sulci_ct.lh$pif ~ variable_sulci_ct.lh$ios))
summary(lm(variable_sulci_ct.lh$pif ~ variable_sulci_ct.lh$pos))
summary(lm(variable_sulci_ct.lh$pif ~ variable_sulci_ct.lh$sf))

summary(lm(variable_sulci_ct.rh$pif ~ variable_sulci_ct.rh$ios))
summary(lm(variable_sulci_ct.rh$pif ~ variable_sulci_ct.rh$pos))
summary(lm(variable_sulci_ct.rh$pif ~ variable_sulci_ct.rh$sf))

# fta
summary(lm(variable_sulci_ct.lh$fta ~ variable_sulci_ct.lh$ios))
summary(lm(variable_sulci_ct.lh$fta ~ variable_sulci_ct.lh$pos))
summary(lm(variable_sulci_ct.lh$fta ~ variable_sulci_ct.lh$sf))

summary(lm(variable_sulci_ct.rh$fta ~ variable_sulci_ct.rh$ios))
summary(lm(variable_sulci_ct.rh$fta ~ variable_sulci_ct.rh$pos))
lm(variable_sulci_ct.rh$fta ~ variable_sulci_ct.rh$sf)

# lft
summary(lm(variable_sulci_ct.lh$lft ~ variable_sulci_ct.lh$ios))
summary(lm(variable_sulci_ct.lh$lft ~ variable_sulci_ct.lh$pos))
summary(lm(variable_sulci_ct.lh$lft ~ variable_sulci_ct.lh$sf))

summary(lm(variable_sulci_ct.rh$lft ~ variable_sulci_ct.rh$ios))
summary(lm(variable_sulci_ct.rh$lft ~ variable_sulci_ct.rh$pos))
summary(lm(variable_sulci_ct.rh$lft ~ variable_sulci_ct.rh$sf))

#### Plot

## Plots comparing sulcogyral count of this study to Nakamura Meta-Analysis Counts across different groups

type <- read.csv("C:/Users/willl/Documents/Urgency/type.csv")
type$Type.IV <- type$Type.IV * -1
type$N <- rowSums(type[,c(5,9)], na.rm = TRUE)
type[,c("I", "II", "III")] = type[,6:8]/type$N
type$Type <- as.factor(type$Type)
type <- type %>%
  select(
    c("Study", "Population", "Hemi", "N", "I", "II", "III")
  ) %>%
  pivot_longer(
    cols = c("I", "II", "III"),
    names_to =  "Type", # category column(s) to pivot from long to wide
    values_to = "Percent", # value columns(s) that hold data for each category column
  )

# Split all participants into healthy and transdiagnostic sub groups

type.healthy <- type %>% subset(Population == "Healthy Control")
type.diverse <- type %>% subset(Population != "Healthy Control")
urgency.types <- read.csv("C:/Users/willl/Documents/Urgency/urgency.types.csv")
names(urgency.types)[names(urgency.types) == 'ï..Hemi'] <- 'Hemi'

# All participants plot

plot.type <- ggplot(data = type, aes(x = Type, y = Percent)) +
  geom_boxplot(data = type, aes(fill = Type), # show summary statistics 
               width = .1,
               alpha=.5, 
               outlier.shape = NA) +
  theme(legend.position = "none") +
  ggtitle("All Participants") +
  ylab("Incidence Rate") +
  geom_point(data = urgency.types, aes(x = Type),
             color = 'red') +
  facet_wrap( ~ as.factor(Hemi)) +
  theme_classic()

ggplotly(plot.type)

# Healthy Participants plot

plot.type.healthy <- ggplot(data = type.healthy, aes(x = Type, y = Percent)) +
  geom_boxplot(aes(fill = Type), # show summary statistics 
               width = .1,
               alpha=.5, 
               outlier.shape = NA) +
  theme(legend.position = "none") +
  ggtitle("Healthy Participants") +
  ylab("Incidence Rate") +
  geom_point(data = urgency.types, aes(x = Type),
             color = 'red') +
  facet_wrap( ~ Hemi) +
  theme_classic()

ggplotly(plot.type.healthy)

ggsave(filename = "~/plot.type.healthy.png",
       plot = plot.type.healthy,
       device = "png",
       width = 6.5,
       height = 3.25,
       units = "in",
       dpi = "retina")

# Transdiagnostic participants plot

plot.type.diverse <- ggplot(data = type.diverse, aes(x = Type, y = Percent)) +
  geom_boxplot(aes(fill = Type), # show summary statistics 
               width = .1,
               alpha=.5, 
               outlier.shape = NA) +
  facet_wrap(~ Hemi) +
  theme(legend.position = "none") +
  ggtitle("Transdiagnostic Participants") +
  ylab("Incidence Rate") +
  geom_point(data = urgency.types, aes(x = Type),
             color = 'red') +
  facet_wrap( ~ as.factor(Hemi)) + 
  theme_classic()


ggplotly(plot.type.diverse)

ggsave(filename = "~/plot.type.diverse.png",
       plot = plot.type.diverse,
       device = "png",
       width = 6.5,
       height = 3.25,
       units = "in",
       dpi = "retina")


## Violin plots of sulcal metrics across sulci

# Cortical thickness

plot.ct <- sulc_metrics_final %>% ggplot(aes(x = label, y = cortical_thickness_mean)) +
  geom_violin(aes(fill = label), # show distribution
              alpha = .3
  ) +
  #geom_jitter(aes(color = label), # show ind. data points
  #position=position_jitterdodge(jitter.width = .2, 
  #                              dodge.width = .1)
  #) +
  geom_boxplot(aes(fill = label), # show summary statistics 
               width = .1,
               alpha=.5, 
               outlier.shape = NA
               
  ) +
  facet_wrap(~ hemi) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", # edit axis names
       y = "Cortical Thickness (mm)")

scale_y_continuous(breaks = c(), # set where ticks are
                   limits = c()) + # set limits
  scale_x_discrete(labels = c()) + # can change names of x-axis levels shown 
  guides() + # remove key (e.g., fill = 'none')
  scale_fill_manual(breaks = c(), # set equal to level names 
                    values = c() # colors for each
  ) +
  scale_color_manual(breaks = c(), # set equal to level names 
                     values = c() # colors for each
  ) 

plot.ct
ggsave(filename = "~/ct_plot.png",
       plot = plot.ct,
       device = "png",
       width = 6.25,
       height = 4,
       units = "in",
       dpi = "retina")

# Mean Sulcal Depth

plot.sd <- sulc_metrics_final %>% ggplot(aes(x = label, y = sulcal_depth_mean_pct)) +
  geom_violin(aes(fill = label), # show distribution
              alpha = .3
  ) +
  #geom_jitter(aes(color = label), # show ind. data points
  #position=position_jitterdodge(jitter.width = .2, 
  #                              dodge.width = .1)
  #) +
  geom_boxplot(aes(fill = label), # show summary statistics 
               width = .1,
               alpha=.5, 
               outlier.shape = NA
               
  ) +
  facet_wrap(~ hemi) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", # edit axis names
       y = "Sulcal Depth (FSU)")

scale_y_continuous(breaks = c(), # set where ticks are
                   limits = c()) + # set limits
  scale_x_discrete(labels = c()) + # can change names of x-axis levels shown 
  
  
  guides() + # remove key (e.g., fill = 'none')
  scale_fill_manual(breaks = c(), # set equal to level names 
                    values = c() # colors for each
  ) +
  scale_color_manual(breaks = c(), # set equal to level names 
                     values = c() # colors for each
  ) 

plot.sd
ggplotly(plot.sd)
ggsave(filename = "~/sd_plot.png",
       plot = plot.sd,
       device = "png",
       width = 6.5,
       height = 4,
       units = "in",
       dpi = "retina")

### Tukey Tests for Label and Hemisphere Interactions

## Used to assess structural differences between tolfs and olfs

# Cortical Thickness

cortical_thickness <- (aov(cortical_thickness_mean ~ label*hemi, data=sulc_metrics_final)) #sf
summary(cortical_thickness)

emmeans.cortical_thickness <- emmeans(cortical_thickness, ~ label | hemi)
contrast(emmeans.cortical_thickness, method='pairwise', adjust='tukey')

emmeans.cortical_thickness_hemi <- emmeans(cortical_thickness, ~ hemi | label)
contrast(emmeans.cortical_thickness_hemi, method='pairwise', adjust='tukey')

# Sulcal Depth

sulcal_depth <- (aov(sulcal_depth_mean_pct ~ label*hemi, data=sulc_metrics_final)) #sf
summary(sulcal_depth)

emmeans.sulcal_depth <- emmeans(sulcal_depth, ~label | hemi | hemi)
contrast(emmeans.sulcal_depth, method='pairwise', adjust='tukey')

emmeans.sulcal_depth_hemi <- emmeans(sulcal_depth, ~ hemi | label)
contrast(emmeans.sulcal_depth_hemi, method='pairwise', adjust='tukey')

# Cortical Thickness

cortical_thickness <- (aov(cortical_thickness_mean ~ label*hemi, data=sulc_metrics_final)) #sf
summary(cortical_thickness)

emmeans.cortical_thickness <- emmeans(cortical_thickness, ~ label | hemi)
contrast(emmeans.cortical_thickness, method='pairwise', adjust='tukey')

emmeans.cortical_thickness_hemi <- emmeans(cortical_thickness, ~ hemi | label)
contrast(emmeans.cortical_thickness_hemi, method='pairwise', adjust='tukey')

# Surface Area

surface_area <- (aov(total_surface_area_.mm.2. ~ label*hemi, data=sulc_metrics_final)) #sf
summary(surface_area)

emmeans.surface_area <- emmeans(sulcal_depth, ~ label)
contrast(emmeans.surface_area, method='pairwise', adjust='tukey')

emmeans.surface_area_hemi <- emmeans(sulcal_depth, ~ hemi | label)
contrast(emmeans.surface_area_hemi, method='pairwise', adjust='tukey')

## Plot comparison of sulcal depth and cortical thickness of olfs to that of the tolfs

# Filter the data to include only the measurements for the right hemisphere and the two sulci
olfs_metrics <- sulc_metrics %>% filter(label %in% c("olfs", "tolfs"))

# Create ct plot
olfs_ct.plot <- ggplot(olfs_metrics, aes(x = label, y = cortical_thickness_mean)) +
  facet_wrap(~ hemi, scales = "free_y") +
  geom_violin() +
  geom_boxplot(width = .15) +
  geom_signif(comparisons = list(c("olfs", "tolfs")), 
              map_signif_level=TRUE) +
  labs(x = '', y = "Cortical Thickness (mm)") +
  theme_classic() +
  theme(legend.position = "none") 


ggsave(filename = "~/olfs_ct.plot.png",
       plot = olfs_rh.plot,
       device = "png",
       width = 3.25,
       height = 4,
       units = "in",
       dpi = "retina")

# Create depth plot
olfs_sd.plot <- ggplot(olfs_metrics, aes(x = label, y = sulcal_depth_mean_pct)) +
  facet_wrap(~ hemi, scales = "free_y") +
  geom_violin() +
  geom_boxplot(width = .15) +
  geom_signif(comparisons = list(c("olfs", "tolfs")), 
              map_signif_level=TRUE) +
  labs(x = '', y = "Sulcal Depth (FSU)") +
  theme_classic() +
  theme(legend.position = "none") 


ggsave(filename = "~/olfs_sd.plot.png",
       plot = olfs_rh.plot,
       device = "png",
       width = 3.25,
       height = 4,
       units = "in",
       dpi = "retina")

# Combine Plots

olfs_sd_ct_plot <- ggarrange(olfs_ct.plot, olfs_sd.plot, ncol = 2, common.legend = TRUE, legend = "right")

ggsave(filename = "~/olfs__sd_ct_plot.png",
       plot = olfs_sd_ct_plot,
       device = "png",
       width = 6.5,
       height = 4,
       units = "in",
       dpi = "retina")

olfs_sd_ct_plot

## Column Plots for Components

# ios

ios_col <- urgency_type %>%
  group_by(ios, hemi) %>%
  summarize(Subjects = n()) %>%
  mutate(Proportion = Subjects / 123)

# Plot the summarized data
ios_plot <- ggplot(ios_col, aes(x = hemi, y = Proportion, fill = ios)) +
  geom_col() +
  labs(x = "", y = "Proportion of Participants", fill = "#ios") +
  theme_classic()

ggsave(filename = "~/ios_plot.png",
       plot = ios_plot,
       device = "png",
       width = 3.25,
       height = 3.25,
       units = "in",
       dpi = "retina")

# pos

pos_col <- urgency_type %>%
  group_by(pos, hemi) %>%
  summarize(Subjects = n()) %>%
  mutate(Proportion = Subjects / 123)

# Plot the summarized data
pos_plot <- ggplot(pos_col, aes(x = hemi, y = Proportion, fill = pos)) +
  geom_col() +
  labs(x = "", y = "Proportion of Participants", fill = "#pos") +
  theme_classic()

ggsave(filename = "~/pos_plot.png",
       plot = pos_plot,
       device = "png",
       width = 3.25,
       height = 3.25,
       units = "in",
       dpi = "retina")

## sf

sf_col <- urgency_type %>%
  group_by(sf, hemi) %>%
  summarize(Subjects = n()) %>%
  mutate(Proportion = Subjects / 123)

# Plot the summarized data
sf_plot <- ggplot(sf_col, aes(x = hemi, y = Proportion, fill = sf)) +
  geom_col() +
  labs(x = "", y = "Proportion of Participants", fill = "#sf") +
  theme_classic()

ggsave(filename = "~/sf_plot.png",
       plot = sf_plot,
       device = "png",
       width = 3.25,
       height = 3.25,
       units = "in",
       dpi = "retina")

## Column Plots for Type

type_col <- urgency_type %>%
  group_by(type, hemi) %>%
  summarize(Subjects = n()) %>%
  mutate(Proportion = Subjects / 123)

# Plot the summarized data
type_plot <- ggplot(type_col, aes(x = hemi, y = Proportion, fill = type)) +
  geom_col() +
  labs(x = "", y = "Proportion of Participants", fill = "Type") +
  theme_classic()

ggsave(filename = "~/type_plot.png",
       plot = type_plot,
       device = "png",
       width = 3.25,
       height = 3.25,
       units = "in",
       dpi = "retina")

### Demographic Tables

## Gather Data 

# Import Diagnosis and Medication Data

tertiary_depth_lasso <- subset(tertiary_depth_lasso, select = -c(X))
scid <- read.csv("scid_approach.csv")
scid <- subset(scid, select = -c(X))
med <- read.csv("MedData_approach.csv")
med <- subset(med, select = -c(X))

# Merge depth and med/diagnosis data

final <- merge(tertiary_depth_lasso, scid, by = "sub") %>%
  merge(med, by = "sub")

# Remove irrelevant data

final <- subset(final, select = -c(Gender_2_TEXT, Race_5_TEXT,
                                   SSRI_start, SSRI_dose_equiv_mg, 
                                   SNRI_start, SNRI_Dose_equiv_mg, 
                                   Tricyclic_start, Tricyclic_dose_equiv_mg,
                                   Other_antidep_start, other_antidep_dose_equiv_mg,
                                   benzodiaz_start, benzodiaz_dose_equiv_mg,
                                   atypical_antipsych_start, atypical_antipsych_dose_equiv_mg,
                                   antipsych_start, antipsych_dose_equiv_mg,
                                   antiseizureantiepileptic_start, antiseizureantiepileptic_dose_equiv_mg
)
)

# Assign correct data types to variables

names <- c(1:2,15:22,24:26, 28)
final[,names] <- lapply(final[,names] , factor)
#final[,3:13] <- lapply(final[,3:13] , numeric)

# Divide by Hemisphere

final_lh <- subset(final, final$hemi == "lh")
final_rh <- subset(final, final$hemi == "rh")

# Create demographic data

demo <- final %>%
  group_by(sub) %>%
  filter(!(hemi == "rh" & "lh" %in% hemi)) %>%
  ungroup()

write.csv(final, "final.csv")

### Make Demographic Tables

table(demo['sub']) #use to find missing hemis
table(demo['Gender'])
table(demo['Race'])
table(demo['Hispanic'])

mean(demo$Age, na.rm = TRUE)
sd(demo$Age, na.rm = TRUE)
range(demo$Age, na.rm = TRUE)

mean(demo$Education, na.rm = TRUE)
sd(demo$Education, na.rm = TRUE)
range(demo$Education, na.rm = TRUE)

table(demo$lifetime_MDE_dx)
table(demo$Lifetime_anxiety_dx)
table(demo$lifetimeAUD)
table(demo$lifetimeSUD)
table(demo$Lifetime_psychosis)

demomultiple <- demo %>%
  mutate(multiple = ifelse(rowSums(select(., 16:22) != 0) >= 2, 1, 0))

table(demomultiple$multiple)

demo <- demo %>%
  mutate_at(vars(3:13), ~as.numeric(as.character(.)))

final <- final %>%
  mutate_at(vars(3:13), ~as.numeric(as.character(.)))

# PIF demographic data

mean(demo$pif, na.rm = TRUE)
sd(demo$pif, na.rm = TRUE)
range(demo$pif, na.rm = TRUE)

# FTA demographic data

mean(demo$fta, na.rm = TRUE)
sd(demo$fta, na.rm = TRUE)
range(demo$fta, na.rm = TRUE)

# LFT demographic data

mean(demo$lft, na.rm = TRUE)
sd(demo$lft, na.rm = TRUE)
range(demo$lft, na.rm = TRUE)

# Medication Data

responses <- demo[, 28:39]

# Function to calculate count and percentage
count_and_percent <- function(x) {
  tbl <- table(x)
  percent <- prop.table(tbl) * 100
  data.frame(count = tbl, percent = percent)
}

# Apply the function to each column
response_counts <- lapply(responses, count_and_percent)

# Print the tables
response_counts
table(final$Cur_Meds)


## Run corrrelation tests on Demographic data

# Age x ERI (RH)

cor.test(final_rh$Age, final_rh$pif)
cor.test(final_rh$Age, final_rh$fta)
cor.test(final_rh$Age, final_rh$lft)

# Age x ERI (LH)

cor.test(final_lh$Age, final_lh$pif) # borderline
cor.test(final_lh$Age, final_lh$fta)
cor.test(final_lh$Age, final_lh$lft)


# Age x Sulcus (RH)

cor.test(final_rh$Age, final_rh$tolfs) # 
cor.test(final_rh$Age, final_rh$olfs)
cor.test(final_rh$Age, final_rh$mosa)
cor.test(final_rh$Age, final_rh$mosp)
cor.test(final_rh$Age, final_rh$losa)
cor.test(final_rh$Age, final_rh$losp)
cor.test(final_rh$Age, final_rh$tos)
cor.test(final_rh$Age, final_rh$sf)

# Age x Sulcus (LH)

cor.test(final_lh$Age, final_lh$tolfs)
cor.test(final_lh$Age, final_lh$olfs)
cor.test(final_lh$Age, final_lh$mosa)
cor.test(final_lh$Age, final_lh$mosp)
cor.test(final_lh$Age, final_lh$losa)
cor.test(final_lh$Age, final_lh$losp)
cor.test(final_lh$Age, final_lh$tos)
cor.test(final_lh$Age, final_lh$sf)

# T-Test comparing meds vs no meds on ERI

# Initialize a list to store t-test results
t_test_results <- list()

# Perform t-tests for each behavioral measure

behavioral_measures <- c("pif", "fta", "lft")
for (measure in behavioral_measures) {
  # Subset data for individuals taking medication and those not taking medication
  meds_group <- subset(final, Cur_Meds == 1)[[measure]]
  no_meds_group <- subset(final, Cur_Meds == 0)[[measure]]
  
  # Perform t-test
  t_test_results[[measure]] <- t.test(meds_group, no_meds_group)
  
  # Print the t-test result
  cat("T-test result for", measure, ":\n")
  print(t_test_results[[measure]])
  cat("\n")
}

t_test_results <- list()

# ANOVA comparing gender differences on ERI

summary(aov(demo$pif ~ demo$Gender))
summary(aov(demo$fta ~ demo$Gender))
summary(aov(demo$lft ~ demo$Gender))

# Correlation Tests comparing Age with ERI

cor.test(demo$pif, demo$Age)
cor.test(demo$fta, demo$Age)
cor.test(demo$lft, demo$Age)