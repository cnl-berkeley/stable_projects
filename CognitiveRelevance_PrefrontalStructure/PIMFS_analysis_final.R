# Libraries
## Requirements: packages below must be installed to load. 

## Run code below if packages are not installed
project_packages <- c("tidyverse", "cowplot", "PupillometryR", "effsize", "DescTools")
install.packages(project_packages)

## Load packages
library(tidyverse)
library(cowplot)
library(PupillometryR)
library(effsize)
library(DescTools)


# Set options and themes
## options
options(scipen=999, # 0
        show.signif.stars = FALSE) 

dodge <- position_dodge(width = .9)

theme_set(theme_classic())

## themes
project_theme <- theme(
  # set axis title size
  axis.title.x = element_blank(), 
  axis.title.y = element_text(color="#333333", angle=90, size=12),
  
  # set axis text size
  axis.text = element_text(color='#333333', size = 10),
  
  # rotate x-axis 45 degrees
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  
  #change background pannel
  panel.grid.major  = element_blank(), 
  panel.grid.minor=element_blank(),
  # change lengend size
  legend.title = element_text(color="#333333", size = 12),
  legend.text = element_text(color = '#333333', size = 10)
  )

project_theme2 <- theme(
  # set axis title size
  axis.title.x = element_blank(), 
  axis.title.y = element_text(color="#333333", angle=90, size=12),
  
  # set axis text size
  axis.text = element_text(color='#333333', size = 10),
  
  # rotate x-axis 45 degrees
  #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  
  #change background pannel
  panel.grid.major  = element_blank(), 
  panel.grid.minor=element_blank(),
  # change lengend size
  legend.title = element_text(color="#333333", size = 12),
  legend.text = element_text(color = '#333333', size = 10)
)


# Data
## Set working directory to where files are located
setwd("./CognitiveRelevance_PrefrontalStructure/data/")

HCP_pimfs <- read.csv("./hcp_pimfs.csv")
HCP_behavior <- read.csv("./hcp_behavior.csv")
chimp_pimfs <- read.csv("./chimp_incidence.csv")
pimfs_IR <- read.csv("./pimfs_incidence_pediatric.csv")

## Minor data manipulation
HCP_pimfs <- merge(HCP_pimfs, 
                   HCP_behavior, 
                   by = "sub")

HCP_pimfs <- HCP_pimfs %>% mutate(number_components_rh = case_when(num_comp_rh == 2 ~ "two",
                                                                   num_comp_rh == 1 ~ "<two",
                                                                   num_comp_rh == 0 ~ "<two"
),
                                  number_components_lh = case_when(num_comp_lh == 2 ~ "two",
                                                                   num_comp_lh == 1 ~ "<two",
                                                                   num_comp_lh == 0 ~ "<two"
))

HCP_pimfs <- HCP_pimfs %>% mutate(component_type_rh = case_when(num_comp_rh == 1 & pimfs_d_rh == "present" ~ "d",
                                                                num_comp_rh == 1 & pimfs_v_rh == "present" ~ "v"
),
                                  component_type_lh = case_when(num_comp_lh == 1 & pimfs_d_lh == "present" ~ "d",
                                                                num_comp_lh == 1 & pimfs_v_lh == "present" ~ "v"
))


##### Incidence rates of pimfs components ##### 
# 1. Does the number of components vary between hemispheres in the adult sample?
# 2. Does the number of components vary within hemisphere in the adult sample?
# 3. Does component type vary with hemispheres in the adult sample?
# 4. Does the number of components vary between age groups and/or species?

# Analyses 

## 1. 

### Left hemisphere
lh.HCP <- table(HCP_pimfs$num_comp_lh)
lh.HCP
lh_chi_HCP <- chisq.test(lh.HCP)
lh_chi_HCP
# X-squared = 54.333, df = 2, p-value = 0.000000000001591

### Right hemisphere
rh.HCP <- table(HCP_pimfs$num_comp_rh)
rh.HCP
rh_chi_HCP <- chisq.test(rh.HCP)
rh_chi_HCP
# X-squared = 68.083, df = 2, p-value = 0.000000000000001644


## 2. 
HCP_pimfs_lh <- HCP_pimfs %>% dplyr::select(sub, num_comp_lh) %>% rename(number_components = num_comp_lh)
HCP_pimfs_lh$hemi <- "lh"
HCP_pimfs_rh <- HCP_pimfs %>% dplyr::select(sub, num_comp_rh) %>% rename(number_components = num_comp_rh)
HCP_pimfs_rh$hemi <- "rh"
HCP_pimfs_hemi <- rbind(HCP_pimfs_lh, HCP_pimfs_rh)

lh_rh_fe <- fisher.test(HCP_pimfs_hemi$hemi, HCP_pimfs_hemi$number_components)
lh_rh_fe
# p-value = 0.6644


## 3. 

### Left hemisphere
lh_component <- HCP_pimfs %>% dplyr::select(sub, num_comp_lh, pimfs_d_lh, pimfs_v_lh, component_type_lh) %>% subset(component_type_lh %in% c("d","v"))
table(lh_component$component_type_lh)
lh_d_v <- chisq.test(table(lh_component$component_type_lh))
lh_d_v
# X-squared = 2, df = 1, p-value = 0.1573

### Right hemisphere
rh_component <- HCP_pimfs %>% dplyr::select(sub, num_comp_rh, pimfs_d_rh, pimfs_v_rh, component_type_rh) %>% subset(component_type_rh %in% c("d","v"))
table(rh_component$component_type_rh)
rh_d_v <- chisq.test(table(rh_component$component_type_rh))
rh_d_v
# X-squared = 0.066667, df = 1, p-value = 0.7963


## 4. 

## Left hemisphere
lh_inc_ya <- HCP_pimfs %>% dplyr::select(sub, num_comp_lh)
lh_inc_ya$species <- "human"
lh_inc_ya$age_group <- "human_young_adult"

lh_inc_ca <- pimfs_IR %>% subset(hemi =='lh') %>% dplyr::select(sub, num_components) %>% rename(num_comp_lh = num_components)
lh_inc_ca$species <- "human"
lh_inc_ca$age_group <- "human_pediatric"

lh_inc_ch <- chimp_pimfs %>% dplyr::select(sub, num_comp_lh)
lh_inc_ch$species <- "chimpanzee"
lh_inc_ch$age_group <- "chimpanzee"

lh_inc <- rbind(lh_inc_ya, lh_inc_ca, lh_inc_ch)
lh_inc.tab <- table(lh_inc$species, lh_inc$num_comp_lh)
chisq.test(lh_inc.tab)
# X-squared = 142.56, df = 2, p-value < 0.00000000000000022

lh_inc.ag <- lh_inc %>% subset(species == "human")
lh_inc.ag.tab <- table(lh_inc.ag$age_group, lh_inc.ag$num_comp_lh)
lh_inc.ag.res <- fisher.test(lh_inc.ag$age_group, lh_inc.ag$num_comp_lh)
# p-value = 0.9407

## Right hemisphere
rh_inc_ya <- HCP_pimfs %>% dplyr::select(sub, num_comp_rh)
rh_inc_ya$species <- "human"
rh_inc_ya$age_group <- "human_young_adult"

rh_inc_ca <- pimfs_IR %>% subset(hemi =='rh') %>% dplyr::select(sub, num_components) %>% rename(num_comp_rh = num_components)
rh_inc_ca$species <- "human"
rh_inc_ca$age_group <- "human_pediatric"

rh_inc_ch <- chimp_pimfs %>% dplyr::select(sub, num_comp_rh)
rh_inc_ch$species <- "chimpanzee"
rh_inc_ch$age_group <- "chimpanzee"

rh_inc <- rbind(rh_inc_ya, rh_inc_ca, rh_inc_ch)
rh_inc.tab <- table(rh_inc$species, rh_inc$num_comp_rh)
chisq.test(rh_inc.tab)
# X-squared = 127.68, df = 2, p-value < 0.00000000000000022

rh_inc.ag <- rh_inc %>% subset(species == "human")
rh_inc.ag.tab <- table(rh_inc.ag$age_group, rh_inc.ag$num_comp_rh)
rh_inc.ag.res <- fisher.test(rh_inc.ag$age_group, rh_inc.ag$num_comp_rh)
# p-value = 0.11


# Visualization 

## Pre-edited data for plot
HCP_pimfs.plot <- read.csv("./pimfs_for_paper.csv")
HCP_pimfs.plot$groups_for_plot <- factor(HCP_pimfs.plot$groups_for_plot, 
                                   levels = c("none", "ventral only", "dorsal only", "two"))
HCP_pimfs.plot$hemi <- gsub("lh", "left hemisphere", HCP_pimfs.plot$hemi)
HCP_pimfs.plot$hemi <- gsub("rh", "right hemisphere", HCP_pimfs.plot$hemi)

## Generate incidence rate plot
HCP_pimfs_plot2 <- HCP_pimfs.plot %>%
  ggplot(aes(x = species, y = percent, fill = groups_for_plot)) + 
  
  # set colors
  geom_col(color = 'black') + 
  
  scale_fill_manual(breaks = c("none", "ventral only", "dorsal only", "two"),
                    values = c("#cccccc", "#fecc5c", "#fd8d3c", "#f03b20")
  ) +
  
  # set labels
  labs(x = "",
       y = "number of participants",
       fill = "# components") +
  
  xlim("humans (young adult)", "humans (pediatric)", "chimpanzees") +
  
  # split by hemisphere
  facet_wrap(~ hemi) +
  
  # set theme 
  project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10) 

### view 
HCP_pimfs_plot2

### save as png
# ggplot2::ggsave(filename = "~/Downloads/HCP_pimfs_plot2.png",
#                 plot = HCP_pimfs_plot2,
#                 device = "png",
#                 width = 8,
#                 height = 6, 
#                 units = "in",
#                 dpi = "retina")



##### Reasoning analyses #####

# Prepare data
HCP_pimfs_analyses <- HCP_pimfs %>% dplyr::select(sub, num_comp_lh, num_comp_rh,
                                           number_components_lh, number_components_rh, 
                                           pimfs_v_lh, pimfs_v_rh, 
                                           pimfs_d_lh, pimfs_d_rh,
                                           PMAT24_A_CR, ProcSpeed_Unadj, 
                                           ListSort_Unadj, CardSort_Unadj, Gender
) %>% drop_na()


# Model 1: 
# 1. Does the number of the pimfs components relate to reasoning scores? 

## Analyses
var.test(HCP_pimfs_analyses$PMAT24_A_CR ~ HCP_pimfs_analyses$number_components_lh) #  p-value = 0.1132
var.test(HCP_pimfs_analyses$PMAT24_A_CR ~ HCP_pimfs_analyses$number_components_rh) #  p-value = 0.4537

### Left pimfs
HCP_pimfs_analyses$number_components_lh <- factor(HCP_pimfs_analyses$number_components_lh, levels = c("two", "<two"))
t.test(PMAT24_A_CR ~ number_components_lh, HCP_pimfs_analyses, var.equal = T)
# t = 2.5444, df = 69, p-value = 0.01319
effsize::cohen.d(PMAT24_A_CR ~ number_components_lh, HCP_pimfs_analyses)
# d estimate: 0.6713098 (medium)

### Right pimfs
HCP_pimfs_analyses$number_components_rh <- factor(HCP_pimfs_analyses$number_components_rh, levels = c("two", "<two"))
t.test(PMAT24_A_CR ~ number_components_rh, HCP_pimfs_analyses, var.equal = T)
# t = 1.2411, df = 69, p-value = 0.2188
effsize::cohen.d(PMAT24_A_CR ~ number_components_rh, HCP_pimfs_analyses)
# d estimate: 0.3525187 (small)


## Visualization 

### Compute mean and sd for Reasoning and hemi
pimfs_components.test.stats <- HCP_pimfs_analyses %>% 
  group_by(number_components_lh) %>% 
  summarise(
    mean = mean(PMAT24_A_CR, na.rm = T),
    sd = sd(PMAT24_A_CR, na.rm = T), 
    n = n()) %>% 
  mutate(se = sd/sqrt(n))

### Generate plot 
pimfs_presence_plot <- HCP_pimfs_analyses %>% 
  ggplot(aes(x = number_components_lh, y = PMAT24_A_CR)) + 
  
  
  # view distribution as violin
  geom_flat_violin(aes(fill = number_components_lh),
                   position = position_nudge(x = .1, y = 0),
                   #adjust = 1.0,
                   alpha = 1,
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = number_components_lh, y = PMAT24_A_CR, fill = number_components_lh),
             position = position_jitter(width = .05),
             shape = 21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_pointrange(data = pimfs_components.test.stats,
                  aes(number_components_lh, mean,
                      ymin=mean-sd,
                      ymax=mean+sd,
                      fill = number_components_lh),
                  position = position_nudge(x = .1, y = 0),
                  lwd= 1,
                  shape = 21,
                  size = 1) +
  
  # set color palette
  scale_fill_manual(breaks = c("two", "<two"),
                    values = c("#f03b20", "#fd8d3c")) + 
  
  # set axis scale
  scale_y_continuous(n.breaks = 10) +
  
  # set theme
  project_theme2 +
  
  # set axis labels
  labs(x = "number components",
       y = "reasoning performance") 

### view
pimfs_presence_plot

### save as png
# ggplot2::ggsave(filename = "~/Desktop/pimfs_presence_plot.png",
#                 plot = pimfs_presence_plot,
#                 device = "png",
#                 width = 3,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")



# Model 2: Reasoning performance ~ Presence/absence of pimfs-d and pimfs-v
# 1. Does the presence of the dorsal and/or ventral component of the pimfs components relate to reasoning scores? 

## Analyses 

### Left pimfs-v
HCP_pimfs_analyses$pimfs_v_lh <- factor(HCP_pimfs_analyses$pimfs_v_lh, levels = c("present", "absent"))
var.test(PMAT24_A_CR ~ pimfs_v_lh, HCP_pimfs_analyses)

t.test(PMAT24_A_CR ~ pimfs_v_lh, HCP_pimfs_analyses, var.equal = T)
# t = 3.4366, df = 69, p-value = 0.001002
effsize::cohen.d(PMAT24_A_CR ~ pimfs_v_lh, HCP_pimfs_analyses)
# d estimate: 1.02508 (large)

### Right pimfs-v
HCP_pimfs_analyses$pimfs_v_rh <- factor(HCP_pimfs_analyses$pimfs_v_rh, levels = c("present", "absent"))
var.test(PMAT24_A_CR ~ pimfs_v_rh, HCP_pimfs_analyses)

t.test(PMAT24_A_CR ~ pimfs_v_rh, HCP_pimfs_analyses, var.equal = T)
# t = 0.81485, df = 69, p-value = 0.418
effsize::cohen.d(PMAT24_A_CR ~ pimfs_v_rh, HCP_pimfs_analyses)
# d estimate: 0.305837 (small)

### Left pimfs-d
HCP_pimfs_analyses$pimfs_d_lh <- factor(HCP_pimfs_analyses$pimfs_d_lh, levels = c("present", "absent"))
var.test(PMAT24_A_CR ~ pimfs_d_lh, HCP_pimfs_analyses)

t.test(PMAT24_A_CR ~ pimfs_d_lh, HCP_pimfs_analyses, var.equal = T)
# t = -1.0574, df = 69, p-value = 0.294
effsize::cohen.d(PMAT24_A_CR ~ pimfs_d_lh, HCP_pimfs_analyses)
# d estimate: -0.3968581 (small)

### Right pimfs-d
HCP_pimfs_analyses$pimfs_d_rh <- factor(HCP_pimfs_analyses$pimfs_d_rh, levels = c("present", "absent"))
var.test(PMAT24_A_CR ~ pimfs_d_rh, HCP_pimfs_analyses)

t.test(PMAT24_A_CR ~ pimfs_d_rh, HCP_pimfs_analyses, var.equal = T)
# t = 1.316, df = 69, p-value = 0.1925
effsize::cohen.d(PMAT24_A_CR ~ pimfs_d_rh, HCP_pimfs_analyses)
# d estimate: 0.4694092 (small)


## Visualization 

### Compute mean and sd for Reasoning and hemi
pimfs_v_components.test.stats <- HCP_pimfs_analyses %>% 
  group_by(pimfs_v_lh) %>% 
  summarise(
    mean = mean(PMAT24_A_CR, na.rm = T),
    sd = sd(PMAT24_A_CR, na.rm = T), 
    n = n()) %>% 
  mutate(se = sd/sqrt(n))

### Generate plot 
pimfs_v_plot <- HCP_pimfs_analyses %>% 
  ggplot(aes(x = pimfs_v_lh, y = PMAT24_A_CR)) + 
  
  
  # view distribution as violin
  geom_flat_violin(aes(fill = pimfs_v_lh),
                   position = position_nudge(x = .1, y = 0),
                   #adjust = 1.0, 
                   alpha = 1, 
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = pimfs_v_lh, y = PMAT24_A_CR, fill = pimfs_v_lh),
             position = position_jitter(width = .05),
             shape=21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_pointrange(data = pimfs_v_components.test.stats, 
                  aes(pimfs_v_lh, mean, 
                      ymin=mean-sd, 
                      ymax=mean+sd, 
                      fill = pimfs_v_lh), 
                  position = position_nudge(x = .1, y = 0),
                  lwd= 1,
                  shape = 21, 
                  size = 1) +
  
  # set color palette
  scale_fill_manual(breaks = c("present", "absent"),
                    values = c("#fecc5c", "#cccccc")) +
  
  # set labels
  labs(x = "pimfs-v presence",
       y = "reasoning performance",
       fill = "presence") +
  
  # set theme
  project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10) 

### view
pimfs_v_plot

### save as png
# ggplot2::ggsave(filename = "~/Desktop/pimfs_v_plot.png",
#                 plot = pimfs_v_plot,
#                 device = "png",
#                 width = 3,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")


# Bootstrapped pimfs-v and reasoning analysis

## Sub-sampled t-tests matching sample size 
HCP_pimfs.pimfs_v_present <- HCP_pimfs_analyses %>% subset(pimfs_v_lh == 'present')
pimfs_v_lh_pres <- HCP_pimfs.pimfs_v_present$PMAT24_A_CR

HCP_pimfs.pimfs_v_absent <- HCP_pimfs_analyses %>% subset(pimfs_v_lh == 'absent')
pimfs_v_lh_abs <- HCP_pimfs.pimfs_v_absent$PMAT24_A_CR

t.test(pimfs_v_lh_pres, pimfs_v_lh_abs, var.equal = T)
cohens_pimfs <- cohen.d(pimfs_v_lh_pres, pimfs_v_lh_abs, method = "cohen's d")

set.seed(1)
B <- 1000
t.vect <- vector(length = B)
p.vect <- vector(length = B)
c.vect <- vector(length = B)
for(i in 1:B){
  boot.c <- sample(pimfs_v_lh_pres, size = 14, replace = F)
  boot.p <- pimfs_v_lh_abs
  
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
DescTools::MedianCI(t.vect.df$t.vect, conf.level = 0.95)

t_dist_plot <- t.vect.df %>% ggplot(aes(x = t.vect)) + 
  
  # view distribution of values
  geom_histogram(color = "black", fill = 'darkorange', bins = 50) +
  
  # view upper bound as a vertical line
  geom_vline(xintercept = 2.491942, linetype = "dashed") +
  
  # view mean as a vertical line
  geom_vline(xintercept = 2.441196, linetype = "solid") +
  
  # view lower bound as a vertical line
  geom_vline(xintercept = 2.378433, linetype = "dashed") +
  
  # view 0 as a vertical line
  geom_vline(xintercept = 0, linetype = "solid", color = "red") +
  
  # set axis scales
  scale_x_continuous(n.breaks = 10, limits = c(-5, 5)) +
  scale_y_continuous(n.breaks = 10) + 
  
  # set theme
  project_theme2

### view
t_dist_plot

### save as png
# ggplot2::ggsave(filename = "~/Desktop/t_dist_plot.png",
#                 plot = t_dist_plot,
#                 device = "png",
#                 width = 3,#4
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")

### P-value plot
DescTools::MedianCI(p.vect.df$p.vect, conf.level = 0.95)
p.vect.df$round <- round(p.vect.df$p.vect, 2)

p_dist_plot <- p.vect.df %>% ggplot(aes(x = p.vect)) + 
  
  # view distribution of values
  geom_histogram(color = "black", fill = 'orange', bins = 50) + 
  
  # view upper bound as a vertical line
  geom_vline(xintercept = 0.02578759, linetype = "dashed") +
  
  # view mean as a vertical line
  geom_vline(xintercept = 0.02241553, linetype = "solid") +
  
  # view lower bound as a vertical line
  geom_vline(xintercept = 0.02002359, linetype = "dashed") +
  
  # set axis scales
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) + 
  
  # set theme
  project_theme2

### view 
p_dist_plot

### save as png
# ggplot2::ggsave(filename = "~/Desktop/p_dist_plot.png",
#                 plot = p_dist_plot,
#                 device = "png",
#                 width = 3, #4
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")

### Cohen's D plot
DescTools::MedianCI(c.vect.df$c.vect, conf.level = 0.95)

c_dist_plot <- c.vect.df %>% ggplot(aes(x = c.vect)) + 
  
  # view distribution of values
  geom_histogram(color = "black", fill = 'orange', bins = 50) + 
  
  # view upper bound as a vertical line
  geom_vline(xintercept = 0.9418654, linetype = "dashed") +
  
  # view mean as a vertical line
  geom_vline(xintercept = 0.9226854, linetype = "solid") +
  
  # view lower bound as a vertical line
  geom_vline(xintercept = 0.8989633, linetype = "dashed") +
  
  # view 0 as a vertical line
  geom_vline(xintercept = 0, linetype = "solid", color = "red") +
  
  # set axis scales
  scale_x_continuous(n.breaks = 10, limits = c(-2.1, 2.1)) +
  scale_y_continuous(n.breaks = 10) + 

  # set theme
  project_theme2

### view
c_dist_plot

### save as png
# ggplot2::ggsave(filename = "~/Downloads/c_dist_plot.png",
#                 plot = c_dist_plot,
#                 device = "png",
#                 width = 4,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")


##### Control analyses #####

## General notes:
## We used the Akaike Information Criterion (AIC) to quantitatively compare models.
## If the ∆AIC is >2, it suggests an interpretable difference between models. 
## If the ∆AIC is >10, it suggests a strong difference between models.
## The lower AIC value indicating the preferred model

# Model 3: Processing speed performance ~ Presence/absence of left pimfs-v
# 1. Does the presence of the left ventral component of the pimfs components relate to processing speed scores? 

## Analysis
var.test(HCP_pimfs_analyses$ProcSpeed_Unadj ~ HCP_pimfs_analyses$pimfs_v_lh) # p-value = 0.5921
t.test(HCP_pimfs_analyses$ProcSpeed_Unadj ~ HCP_pimfs_analyses$pimfs_v_lh, var.equal = T)
cohen.d(ProcSpeed_Unadj ~ pimfs_v_lh, HCP_pimfs_analyses)


## Visualization

### Compute mean and sd for processing speed
pimfs_v_components.test.stats3 <- HCP_pimfs_analyses %>% 
  group_by(pimfs_v_lh) %>% 
  summarise(
    mean = mean(ProcSpeed_Unadj, na.rm = T),
    sd = sd(ProcSpeed_Unadj, na.rm = T), 
    n = n()) %>% 
  mutate(se = sd/sqrt(n))

### Generate plot 
pimfs_v_plot3 <- HCP_pimfs_analyses %>% 
  ggplot(aes(x = pimfs_v_lh, y = ProcSpeed_Unadj)) + 
  
  
  # view distribution as violin
  geom_flat_violin(aes(fill = pimfs_v_lh),
                   position = position_nudge(x = .1, y = 0),
                   alpha = 1, 
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = pimfs_v_lh, y = ProcSpeed_Unadj, fill = pimfs_v_lh),
             position = position_jitter(width = .05),
             shape=21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_pointrange(data = pimfs_v_components.test.stats3, 
                  aes(pimfs_v_lh, mean, 
                      ymin=mean-sd, 
                      ymax=mean+sd, 
                      fill = pimfs_v_lh), 
                  position = position_nudge(x = .1, y = 0),
                  lwd= 1,
                  shape = 21, 
                  size = 1) +
  
  # set colors 
  scale_fill_manual(breaks = c("present", "absent"),
                    values = c("#fecc5c", "#cccccc")) +
  
  # set labels
  labs(x = "pimfs-v presence",
       y = "ProcSpeed_Unadj",
       fill = "presence") +
  
  # set theme
  project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10) 

### view
pimfs_v_plot3

### save as png
# ggplot2::ggsave(filename = "~/Desktop/pimfs_v_plot3.png",
#                 plot = pimfs_v_plot3,
#                 device = "png",
#                 width = 3,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")


# Model 4: Working memory performance ~ Presence/absence of left pimfs-v
# 1. Does the presence of the left ventral component of the pimfs components relate to working memory scores? 

## Analysis
var.test(HCP_pimfs_analyses$ListSort_Unadj ~ HCP_pimfs_analyses$pimfs_v_lh) # p-value = 0.8218
t.test(HCP_pimfs_analyses$ListSort_Unadj ~ HCP_pimfs_analyses$pimfs_v_lh, var.equal = T)
cohen.d(ListSort_Unadj ~ pimfs_v_lh, HCP_pimfs_analyses)


## AIC comparison

### Models to calculate AIC
model2 <- lm(PMAT24_A_CR ~ pimfs_v_lh, 
             data = HCP_pimfs_analyses) 

model3 <- lm(ListSort_Unadj ~ pimfs_v_lh, 
             data = HCP_pimfs_analyses) 

### Calculate AIC
AIC2 <- AIC(model2)
AIC3 <- AIC(model3)

### Calculate difference
AIC3 - AIC2


## Visualization

### Compute mean and sd for working memory 
pimfs_v_components.test.stats2 <- HCP_pimfs_analyses %>% 
  group_by(pimfs_v_lh) %>% 
  summarise(
    mean = mean(ListSort_Unadj, na.rm = T),
    sd = sd(ListSort_Unadj, na.rm = T), 
    n = n()) %>% 
  mutate(se = sd/sqrt(n))

### Generate plot 
pimfs_v_plot2 <- HCP_pimfs_analyses %>% 
  ggplot(aes(x = pimfs_v_lh, y = ListSort_Unadj)) + 
  
  
  # view distribution as violin
  geom_flat_violin(aes(fill = pimfs_v_lh),
                   position = position_nudge(x = .1, y = 0),
                   #adjust = 1.0, 
                   alpha = 1, 
                   colour = 'black'
  ) +
  
  # view individual data points
  geom_point(aes(x = pimfs_v_lh, y = ListSort_Unadj, fill = pimfs_v_lh),
             position = position_jitter(width = .05),
             shape=21,
             size = 1.5,
             alpha = .75) +
  
  # include summary stats
  geom_pointrange(data = pimfs_v_components.test.stats2, 
                  aes(pimfs_v_lh, mean, 
                      ymin=mean-sd, 
                      ymax=mean+sd, 
                      fill = pimfs_v_lh), 
                  position = position_nudge(x = .1, y = 0),
                  lwd= 1,
                  shape = 21, 
                  size = 1) +
  
  # set colors
  scale_fill_manual(breaks = c("present", "absent"),
                    values = c("#fecc5c", "#cccccc")) +
  
  # set labels
  labs(x = "pimfs-v presence",
       y = "WM performance",
       fill = "presence") +
  
  # set theme
  project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10) 

### view
pimfs_v_plot2

### save as png
# ggplot2::ggsave(filename = "~/Desktop/pimfs_v_plot2.png",
#                 plot = pimfs_v_plot2,
#                 device = "png",
#                 width = 3,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina")