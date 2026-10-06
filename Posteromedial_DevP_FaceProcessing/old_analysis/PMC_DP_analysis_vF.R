###### PMC Analysis vF ######
### 09/15/2025 ###

library(car)
library(PupillometryR)
library(ggplot2)
library(gplots)
library(arm)
library(GGally)
library(performance)
library(ggpubr)
library(cocor)
library(reghelper)
library(sjstats)
library(lme4)
library(broom)
library(effectsize)
library(afex)
library(lsr)
library(patchwork)
library(emmeans)
library(ggpattern)
library(nlme)
library(multcomp)
library(tidyverse)
library(tidylog)



# use standardize() function from the arm package only
standardize <- arm::standardize
select <- dplyr::select

# lists
full_lasso_lst <- c("ifrms", "sspls-v", "prcus-i", "prcus-p", 'prcus-a', "sbps", "mcgs", "pos", "prculs-d")
full_pmc_lst <- c("icgs-p", "sspls-d", "ifrms", "sspls-v", "prcus-a", "prcus-i", 'prcus-p', 
                  "sbps", "mcgs", "pos", "prculs-d", "prculs-v")
full_pmc_lst_EXTRA <- c("icgs-p", "sspls-d", "ifrms", "sspls-v", "prcus-a", "prcus-i", 'prcus-p', "sbps", "mcgs", "pos", "prculs-d", "prculs-v",
                        "sbps_a", "sbps_p", "pmc_p", "pmc_pX", "pmc_pY", "pmc_pXY")
pcc_list <- c("ifrms", "icgs-p", "sspls-d", "sspls-v")
prc_list <- c("prcus-a", "prcus-i", "prcus-p", "prculs-v", "prculs-d")
border_list <- c("mcgs", "pos", "sbps")
consistent_list <- c("ifrms", "prcus-i", "prcus-p", 'prcus-a', "sbps", "mcgs", "pos", "prculs-d")


# Import data
dp <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/all_sulc_04_16_23_vf.csv", stringsAsFactors = T)
dp_demog <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/dp_demog_2023_10_19.csv")
act_df_all <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/act_df_all_2023_10_19.csv")
all_morph <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/all_morph_2023_10_19.csv")

dp_depth <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/dp_depth_2023_10_19.csv")
dp_thick <- read.csv("~/Documents/Documents_Primary/Labs/weiner_lab/dp_pcc/data/dp_thick_2023_10_19.csv")


#### Plot Parameters ####
## Setting up all Plotting parameters for all future plots ##
my_colors <- c(
  "#0072B2", "#D55E00"
)
my_colors_all <- c(
  "#332288", "#0072B2", "#D55E00"
)

## set boxplot theme
set_theme <- function(){
  
  # Customize axis and legend
  set_theme = theme(
    # set axis title size
    # axis.title.x = element_text(color="#333333", size=12), 
    # 
    # # set axis text size
    # axis.text = element_text( color='#333333', size = 14 ),
    # 
    # rotate x-axis 45 degrees
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    
    #change background pannel
    panel.grid.major  = element_blank(), 
    panel.grid.minor=element_blank(),
    # change lengend size
    legend.title = element_text(color="#333333",size = 12),
    # legend.title = element_blank(),
    legend.text = element_text( color = '#333333',size = 12),
    
    # set aspect ratio
    aspect.ratio = 1/1
  )
  # replace classic theme setting with our setting
  theme = theme_classic()%+replace% set_theme
  
  return(theme) 
}
## custom scales
set_scale_hemi <- function (){
  scale = scale_color_manual(name="Hemi", 
                             
                             # manually set color scale
                             values=c( 
                               "#999999",
                               "#000000"
                             )
  )
  return (scale)
}
set_scale_group <- function (){
  scale = scale_fill_manual(name="Group", 
                            # manually set color scale
                            values= my_colors
  )
  return (scale)
}
set_scale_all_groups <- function (){
  scale = scale_fill_manual(name="Group", 
                            # manually set color scale
                            values= my_colors_all
  )
  return (scale)
}
set_scale_group_color <- function (){
  scale = scale_color_manual(name="Group", 
                             
                             # manually set color scale
                             values=my_colors
  )
  return (scale)
}
## set boxplot
boxplot <- function (ggplot_object, xlab, ylab, ymin, ymax){
  
  plot=ggplot_object + 
    
    # boxplot
    geom_boxplot() +
    
    # set y-axis limits
    ylim(ymin, ymax) +
    
    # set axis labels
    xlab(xlab) + 
    ylab(ylab)
  
  return(plot)
}




#### Descriptive Stats of Sample ####
summary(dp_demog)
head(dp_demog, 10)

dp_demog_controls <- dp_demog %>% 
  filter(group == 'Controls')
dp_demog_dps <- dp_demog %>% 
  filter(group == 'DPs')

summary(dp_demog_controls)
summary(dp_demog_dps)

# CFMT
mean(dp_demog_controls$CFMT, na.rm = T)
sd(dp_demog_controls$CFMT, na.rm = T)
range(dp_demog_controls$CFMT, na.rm = T)
hist(dp_demog_controls$CFMT, breaks = c(15))

mean(dp_demog_dps$CFMT, na.rm = T)
sd(dp_demog_dps$CFMT, na.rm = T)
range(dp_demog_dps$CFMT, na.rm = T)
hist(dp_demog_dps$CFMT, breaks = c(15))

# Age
mean(dp_demog_controls$Age)
sd(dp_demog_controls$Age)
hist(dp_demog_controls$Age, breaks = c(15))

mean(dp_demog_dps$Age)
sd(dp_demog_dps$Age)
hist(dp_demog_dps$Age, breaks = c(15))




#### Analysis: Face Selectivity setup ####
mean_activation_plot_GRP <- function(data){
  p_controls <- data %>%
    ggplot(aes(x = label, y = mean_face_activation, fill = group)) +
    stat_summary(fun=mean,
                 position=position_dodge(width=0.95),
                 geom="bar",
                 color = "black",
                 alpha = .6) +
    stat_summary(fun.data=mean_se,
                 position=position_dodge(0.95),
                 geom="errorbar",
                 width = 0.2) +
    
    # geom_line(aes(group = sub), alpha=0.05) +
    # geom_point(alpha = 0.1) +
    
    # facet_wrap(~group, scales = "free") +
    # ylim(-0.31, 0.61) +
    # scale_fill_manual(values = c('#377eb8','#e41a1c')) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank())
  
  return(p_controls)
}

hemi.labs <- c("Left", "Right")
names(hemi.labs) <- c("lh", "rh")

hcp_plt_params <- list(ylab = "Faces - All Other (Z-score)", 
                       group_labels = c("HCP NTs"), 
                       group_colors = c("#332288")
)
nt_dp_plt_params <- list(ylab = "Faces - Objects (%)", 
                         group_labels = c("NTs", "DPs"), 
                         group_colors = c("#0072B2", "#D55E00"), 
                         ylim_coords = c(NA, NA)
)
nt_plt_params <- list(ylab = "Faces - Objects (%)", 
                      group_labels = c("NTs"), 
                      group_colors = c("#0072B2")
)
dp_plt_params <- list(ylab = "Faces - Objects (%)", 
                      group_labels = c("DPs"), 
                      group_colors = c("#D55E00")
)
nt_dp_hcp_plt_params <- list(ylab = "Faces - Objects (%)", 
                             group_labels = c("NTs", "DPs", "HCP_Controls"), 
                             group_colors = c("#0072B2", "#D55E00", "#332288"), 
                             ylim_coords = c(NA, NA)
)


act_df_all$label <- factor(act_df_all$label, levels = c("mcgs", "sbps", "icgs-p", "ifrms", "sspls-d", "sspls-v", "prcus-a", "prcus-i", "prcus-p", "prculs-v", "prculs-d", "pos", "sbps_a", "sbps_p", "pmc_p", "pmc_pX", "pmc_pY", "pmc_pXY"))


#### Analysis: Face Selectivity (data from Jiahui et al. 2018) ####

# testing for differences in hemisphere and group
sulc_data <- act_df_all %>% filter(group %in% c("Controls", "DPs")) %>%
  filter(label %in% full_pmc_lst)

### A 3-way (hemi x sulcus x group) RM-ANOVA to test for hemi and group effects
mod_anova_summary <- sulc_data %>% 
  rstatix::anova_test(mean_face_activation ~ group*hemi*label + Error(sub/hemi))
mod_anova_summary

# Now structure it so you can get post-hoc tests with emmeans
mod <- lmerTest::lmer(mean_face_activation ~ label*group*hemi + (1|sub/hemi), 
                         data=sulc_data)
# Test face selectivity of sulci
tests_by_label <- emmeans::emmeans(mod, ~ label | group | hemi)
test(tests_by_label, null = 0, adjust = "bonferroni", as.df = T)
# Now pairwise contrasts
pairwise_tests <- emmeans::emmeans(mod, ~ label | hemi | group)
post_hoc_tests <- emmeans::contrast(pairwise_tests, method='pairwise')
# More pairwise contrasts
pairwise_tests <- emmeans::emmeans(mod, ~ hemi | label | group)
post_hoc_tests <- emmeans::contrast(pairwise_tests, method='pairwise')



#### Analysis: Face Selectivity (replication using HCP data) ####

# testing for differences in hemisphere and group
sulc_data <- act_df_all %>% filter(group %in% c("HCP_Controls")) %>%
  filter(label %in% full_pmc_lst)

### A 3-way (hemi x sulcus x group) RM-ANOVA to test for hemi and group effects
mod_anova_summary <- sulc_data %>% 
  rstatix::anova_test(mean_face_activation ~ hemi*label + Error(sub/hemi))
mod_anova_summary

# Now structure it so you can get post-hoc tests with emmeans
mod <- lmerTest::lmer(mean_face_activation ~ label*hemi + (1|sub/hemi), 
                      data=sulc_data)
# Test face selectivity of sulci
tests_by_label <- emmeans::emmeans(mod, ~ label | hemi)
test(tests_by_label, null = 0, adjust = "bonferroni", as.df = T)
# Now pairwise contrasts
pairwise_tests <- emmeans::emmeans(mod, ~ label | hemi)
post_hoc_tests <- emmeans::contrast(pairwise_tests, method='pairwise')
# More pairwise contrasts
pairwise_tests <- emmeans::emmeans(mod, ~ hemi | label)
post_hoc_tests <- emmeans::contrast(pairwise_tests, method='pairwise')


# sulci ordered by mean face activation
allsulc_selectivity_sorted <- act_df_all %>% filter(group %in% c("Controls", "DPs", "HCP_Controls")) %>%
  # filter(label %in% full_lasso_lst) %>%
  filter(label %in% full_pmc_lst) %>% 
  dplyr::group_by(label, hemi) %>% 
  summarise(mean = mean(mean_face_activation)) %>% 
  arrange(desc(mean))
allsulc_selectivity_sorted




#### Analysis: Sulcal Incidence ####

## Test for incidence differences in DPs in the sspls-v, which is located in the familiar-face area identified in previous studies ##

# create incidence df
dp_pres_df <- all_morph %>%
  group_by(sub, hemi) %>% 
  summarise(ssplsv_present = if_else((sum(label == 'sspls-v') == 1), 1, 0),
            prculsv_present = if_else((sum(label == 'prculs-v') == 1), 1, 0),
            ssplsd_present = if_else((sum(label == 'sspls-d') == 1), 1, 0),
            icgsp_present = if_else((sum(label == 'icgs-p') == 1), 1, 0),
            ifrms_present = if_else((sum(label == 'ifrms') == 1), 1, 0))


dp_pres_df <- left_join(all_morph, dp_pres_df)
dp_pres_df <- dp_pres_df[c('sub', 'hemi', 'group', 'CFMT', 'Famous_Faces_UKidentified', 'Age', 'Sex', 
                           'ssplsv_present', 'prculsv_present', 'ssplsd_present', 'icgsp_present',
                           'ifrms_present')] %>% 
  distinct() %>% 
  filter(group %in% c("Controls", "DPs"))

dp_pres_df %>% 
  group_by(group, hemi) %>% 
  summarise(ssplsv_present = sum(ssplsv_present) / length(ssplsv_present))

dp_pres_df %>% 
  group_by(group, hemi) %>% 
  summarise(ssplsv_present = sum(ssplsv_present))


# Visualize incidences
dp_pres_df %>% 
  group_by(group, hemi) %>%
  # group_by(hemi) %>%
  summarize(ssplsv_present = sum(ssplsv_present) / length(ssplsv_present)) %>% 
  
  ggplot(aes( x = hemi, y = ssplsv_present)) + 
  geom_bar( position = 'dodge', stat = 'identity') + 
  facet_wrap(vars(group)) +
  ylim(c(0, 1)) +
  set_theme() + 
  set_scale_group()

lh_pres <- dp_pres_df %>% filter(hemi == 'lh')
rh_pres <- dp_pres_df %>% filter(hemi == 'rh')
dp_pres <- dp_pres_df %>% filter(group == 'DPs')
ctrl_pres <- dp_pres_df %>% filter(group == 'Controls')

### Test for incidence differences between NTs and DPs in the sspls-v ###
ssplsv_incidence_df <- dp_pres_df %>%
  group_by(group, hemi) %>%
  summarize(ct = sum(ssplsv_present),
            total_hemis = length(ssplsv_present),
            pct = round((sum(ssplsv_present) / length(ssplsv_present)) *100, 2)) %>% 
  as.data.frame()
ssplsv_incidence_df

grp_lh_chi <- chisq.test(lh_pres$group, lh_pres$ssplsv_present)
grp_rh_chi <- chisq.test(rh_pres$group, rh_pres$ssplsv_present) # sig. group effect
hemi_dp_chi <- chisq.test(dp_pres$hemi, dp_pres$ssplsv_present) # sig. hemi diff. in DPs
hemi_ctrl_chi <- chisq.test(ctrl_pres$hemi, ctrl_pres$ssplsv_present)

grp_rh_chi
hemi_dp_chi

# as an exploratory analysis, we can test the other sulci too (none approach significance)
ssplsd_incidence_df <- dp_pres_df %>%
  group_by(group, hemi) %>%
  summarize(ct = sum(ssplsd_present),
            total_hemis = length(ssplsd_present),
            pct = round((sum(ssplsd_present) / length(ssplsd_present)) *100, 2)) %>% 
  as.data.frame()
ssplsd_incidence_df
grp_lh_chi <- chisq.test(lh_pres$group, lh_pres$ssplsd_present)
grp_rh_chi <- chisq.test(rh_pres$group, rh_pres$ssplsd_present) 
hemi_dp_chi <- chisq.test(dp_pres$hemi, dp_pres$ssplsd_present) 
hemi_ctrl_chi <- chisq.test(ctrl_pres$hemi, ctrl_pres$ssplsd_present)

icgsp_incidence_df <- dp_pres_df %>%
  group_by(group, hemi) %>%
  summarize(ct = sum(icgsp_present),
            total_hemis = length(icgsp_present),
            pct = round((sum(icgsp_present) / length(icgsp_present)) *100, 2)) %>% 
  as.data.frame()
icgsp_incidence_df
grp_lh_chi <- chisq.test(lh_pres$group, lh_pres$icgsp_present)
grp_rh_chi <- chisq.test(rh_pres$group, rh_pres$icgsp_present) 
hemi_dp_chi <- chisq.test(dp_pres$hemi, dp_pres$icgsp_present) 
hemi_ctrl_chi <- chisq.test(ctrl_pres$hemi, ctrl_pres$icgsp_present)

prculsv_incidence_df <- dp_pres_df %>%
  group_by(group, hemi) %>%
  summarize(ct = sum(prculsv_present),
            total_hemis = length(prculsv_present),
            pct = round((sum(prculsv_present) / length(prculsv_present)) *100, 2)) %>% 
  as.data.frame()
prculsv_incidence_df
grp_lh_chi <- chisq.test(lh_pres$group, lh_pres$prculsv_present)
grp_rh_chi <- chisq.test(rh_pres$group, rh_pres$prculsv_present) 
hemi_dp_chi <- chisq.test(dp_pres$hemi, dp_pres$prculsv_present) 
hemi_ctrl_chi <- chisq.test(ctrl_pres$hemi, ctrl_pres$prculsv_present)


#### Plots ####

### Face Activation: Plots ###

p_controls_dps <- mean_activation_plot_GRP(act_df_all %>% filter((group %in% c("Controls", "DPs")) & 
                                                                   label %in% c("prcus-a", "prcus-i", "prcus-p"))) +
  ylab("Faces - Objects (%)") +
  theme(axis.title.y = element_text(size = 9)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "bottom") +
  scale_fill_manual(labels=c("NTs", "DPs"), values = c("#0072B2", "#D55E00")) +
  facet_wrap(vars(hemi), labeller = labeller(hemi = hemi.labs)) +
  coord_cartesian(ylim = c(-.053, .6))

p_hcp_controls <- mean_activation_plot_GRP(act_df_all %>% filter((group %in% c("HCP_Controls")) & 
                                                                   label %in% c("prcus-a", "prcus-i", "prcus-p"))) +
  labs(y = "Faces - All Other (Z-score)") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 9)) +
  scale_fill_manual(labels=c("HCP NTs"), values = c("#332288")) +
  facet_wrap(vars(hemi), labeller = labeller(hemi = hemi.labs)) +
  coord_cartesian(ylim = c(NA, 1.6))

# FIGURE 2 ---> Face activation for precuneal sulci
p_controls_dps + p_hcp_controls + plot_layout(widths = c(3.6,2))


### Now for the sspls-v activation plot
p_controls_dps_ssplsv <- mean_activation_plot_GRP(act_df_all %>% filter((group %in% c("Controls", "DPs")) & label %in% c("sspls-v"))) +
  ylab("Faces - Objects (%)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "bottom") +
  scale_fill_manual(labels=c("NTs", "DPs"), values = c("#0072B2", "#D55E00")) +
  facet_wrap(vars(hemi), labeller = labeller(hemi = hemi.labs)) +
  coord_cartesian(ylim = c(-0.17, 0.65)) +
  set_theme() + 
  theme_classic() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =30, hjust = 1, vjust = 1),
  )

p_hcp_controls_ssplsv <- mean_activation_plot_GRP(act_df_all %>% filter((group %in% c("HCP_Controls")) & label %in% c("sspls-v"))) +
  labs(y = "Faces - All Other (Z-score)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "bottom") +
  scale_fill_manual(labels=c("HCP NTs"), values = c("#332288")) +
  facet_wrap(vars(hemi), labeller = labeller(hemi = hemi.labs)) +
  coord_cartesian(ylim = c(-0.44, 1.65)) +
  set_theme() + 
  theme_classic() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =30, hjust = 1, vjust = 1))

# p_controls_dps_ssplsv + p_hcp_controls_ssplsv + plot_layout(widths = c(2,3))

## sspls-v incidence plot
pmc_incidence.figure.ssplsv <- ssplsv_incidence_df %>%
  filter(group != "HCP_Controls") %>%
  
  ggplot(aes(x = group, y = pct, fill = group)) +
  
  geom_col(aes(alpha = hemi),
           color = "black",
           position = "dodge",
           # alpha=.75
  ) +
  scale_fill_manual(labels=c("NTs", "DPs"), values = c("#0072B2", "#D55E00")) +
  scale_alpha_manual(values = c(.8, 0.4), labels = c("Left", "Right")) +
  scale_x_discrete(labels = c("Controls" = "NTs", "DPs" = "DPs")) +
  scale_y_continuous(name = "% of Hemispheres", seq(0,100,20), limits = c(0,100)) +
  set_theme() + 
  theme_classic() + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =30, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 15)
  )
# pmc_incidence.figure.ssplsv

# FIGURE 3: sspls-v incidence + face selectivity plots
layout <- 
  "AAABBC"
pmc_incidence.figure.ssplsv + p_controls_dps_ssplsv + p_hcp_controls_ssplsv + plot_layout(design = layout)



