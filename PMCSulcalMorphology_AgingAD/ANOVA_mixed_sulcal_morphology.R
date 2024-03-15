## ---------------------------
## Title: LME analyses
##
## Code for LME analyses comparing sulcal morphology across groups (Fig. 2, Table 3)
## for the following manuscript: Maboudian et al., JNeurosci 2024: 
## Defining overlooked structures reveals new associations between cortex and 
## cognition in aging and Alzheimer’s disease (accepted).
##
## Author: Samira Maboudian
##
## Date Modified: 02-01-2024
## ---------------------------

## Relevant Libraries ##
library('dplyr') #v1.1.0
library(nlme) #v3.1.160
library(effectsize) #v0.8.3
library(emmeans) #v1.8.4.1
library(tidyverse) #v2.0.0


## load data ##
morphology_all = read.csv('data/df_SulcMorph_AD72_CN72_YA72_final.csv')
length(unique(morphology_all$sub))

# categorize sulci as deep or shallow
shallow_sulci = c('ifrms','sspls-d', 'sspls-v', 'icgs-p', 'prculs-v')
morphology_all <- morphology_all %>%  mutate(
  depth_group = ifelse(label %in% shallow_sulci, "shallow", "deep"))
#summary(morphology_all)



## MODEL 1a: sulcal cortical thickness LME ANOVA, by depth grouping ##
model1a <- lme( cortical_thickness_mean ~ group * hemi * depth_group, 
               random = ~ 1|sub/hemi/depth_group,
               data = morphology_all)
summary(model1a)
model1a.aov <- anova(model1a)
model1a.aov

### partial effect size
eta_squared(model1a.aov, partial = TRUE)

### post-hoc comparisons
EMM.m1.group <- emmeans::emmeans(model1a, ~ hemi)
emmeans::contrast(EMM.m1.group, method='pairwise')  

EMM.m1.group_sulctype <- emmeans::emmeans(model1a, ~ group | depth_group)
emmeans::contrast(EMM.m1.group_sulctype, method='pairwise') 
eff.m1.group_sulctype<-eff_size(EMM.m1.group_sulctype, 
                                sigma = sigma(model1a), edf = 213) #eff size
eff.m1.group_sulctype



## MODEL 1b: sulcal cortical thickness LME ANOVA, by sulcus (label) ##
# these results are reported in Table 3.
model1b <- lme( cortical_thickness_mean ~group * hemi * label, 
               random = ~ 1|sub/hemi/label,
               data = morphology_all)
summary(model1b)
model1b.aov <- anova(model1b)
model1b.aov

### partial effect size
eta_squared(model1b.aov, partial = FALSE)

### post-hoc comparisons
EMM.m1b.group <- emmeans::emmeans(model1b, ~ group)
emmeans::contrast(EMM.m1b.group, method='pairwise')  

EMM.m1b.group_label <- emmeans::emmeans(model1b, ~ group | label)
emmeans::contrast(EMM.m1b.group_label, method='pairwise') 
eff.m1b.group_label<-eff_size(EMM.m1b.group_label, 
               sigma = sigma(model1b), edf = 213) #effect size
eff.m1b.group_label

EMM.m1b.group_hemi <- emmeans::emmeans(model1b, ~ group | hemi)
emmeans::contrast(EMM.m1b.group_hemi, method='pairwise') 




## MODEL 2: sulcal depth LME ANOVA, by sulcus (label) ##
model2 <- lme( sulcal_depth_mm ~group * hemi * label, 
                random = ~ 1|sub/hemi/label,
               data = morphology_all)
summary(model2)
model2.aov <- anova(model2) 
model2.aov

eta_squared(model2.aov, partial = FALSE)

EMM.m2.group <- emmeans::emmeans(model2, ~ group)
emmeans::contrast(EMM.m2.group, method='pairwise')  

EMM.m2.hemi <- emmeans::emmeans(model2, ~ hemi)
emmeans::contrast(EMM.m2.hemi, method='pairwise')  

EMM.m2.grp_label <- emmeans::emmeans(model2, ~  group | label)
emmeans::contrast(EMM.m2.grp_label, method='pairwise') 

EMM.m2.grp_label <- emmeans::emmeans(model2, ~  hemi | label)
emmeans::contrast(EMM.m2.grp_label, method='pairwise') 





