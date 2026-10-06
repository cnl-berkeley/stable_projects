#### libraries ####
library(extrafont)
library(tidyverse)
library(effectsize)
library(nlme)
library(ggbeeswarm)
library(nlme)
library(lme4)
library(lmerTest)
library(mediation)

#### plot information ####
FFA_project_theme <- 
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", face = "bold", size = 14, hjust = 0.5),

    # set axis title
    axis.title = element_text(family = "Arial", color = "black", face = "bold", size = 14), 

    # set axis text
    axis.text = element_text(family = "Arial", color = "black", size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", face = "bold", size = 14),
    legend.text = element_text(family = "Arial", color = "black", size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", face = "bold", size = 12)
  ) 

FFA_project_theme2 <- 
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", face = "bold", size = 14, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", face = "bold", size = 14), 
    
    # set axis text
    axis.text = element_text(family = "Arial", color = "black", size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", face = "bold", size = 14),
    legend.text = element_text(family = "Arial", color = "black", size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", face = "bold", size = 14),
    strip.background =  element_blank()
  ) 


#### FFA patterns in DPs vs NTs ####
FFA_demo <- read.csv("./data/dp_handedness.csv")
revision_metrics <- read.csv("./data/revision_metrics.csv")
FFA_DP <- read.csv("./data/FG_patterns.csv")

# 0. 
FFA_DP_inc.tab <- table(FFA_DP$group, FFA_DP$FUS_incidence)
FFA_DP_inc.tab
fisher.test(FFA_DP$group, FFA_DP$FUS_incidence) # p-value = 0.1515

FFA_DP.dp <- FFA_DP %>% 
  filter(group =='DPs')
fisher.test(FFA_DP.dp$hemi, FFA_DP.dp$FUS_incidence) # p-value = 0.747

FFA_DP.dp.specROI <- FFA_DP.dp %>% subset(FUS_specific_pattern %in% c("mFUS", "pFUS"))
FFA_DP.dp.specROI.tab <- table(FFA_DP.dp.specROI$FUS_specific_pattern)
chisq.test(FFA_DP.dp.specROI.tab)


FFA_DP.nt <- FFA_DP %>% 
  filter(group =='Controls')
fisher.test(FFA_DP.nt$hemi, FFA_DP.nt$FUS_incidence) # p-value = 0.747

FFA_DP.nt.specROI <- FFA_DP.nt %>% subset(FUS_specific_pattern %in% c("mFUS", "pFUS"))
FFA_DP.nt.specROI.tab <- table(FFA_DP.nt.specROI$FUS_specific_pattern)
chisq.test(FFA_DP.nt.specROI.tab)

# 0.1
lh_FFA_DP_inc <- FFA_DP %>% 
  filter(hemi =='lh')
lh_FFA_DP_inc.tab <- table(lh_FFA_DP_inc$group, lh_FFA_DP_inc$FUS_incidence)
fisher.test(lh_FFA_DP_inc$group, lh_FFA_DP_inc$FUS_incidence) # p-value = 0.3279

# 0.2
rh_FFA_DP_inc <- FFA_DP %>% 
  filter(hemi =='rh')
rh_FFA_DP_inc.tab <- table(rh_FFA_DP_inc$group, rh_FFA_DP_inc$FUS_incidence)
fisher.test(rh_FFA_DP_inc$group, rh_FFA_DP_inc$FUS_incidence) # p-value = 0.4796

# 1.
FFA_DP.tab <- table(FFA_DP$group, FFA_DP$FUS_general_pattern)
FFA_DP.tab
fisher.test(FFA_DP$group, FFA_DP$FUS_general_pattern) # p-value = 0.2066

# 1.1 
lh_FFA_DP <- FFA_DP %>% 
  filter(hemi =='lh')
lh_FFA_DP.tab <- table(lh_FFA_DP$group, lh_FFA_DP$FUS_general_pattern)
fisher.test(lh_FFA_DP$group, lh_FFA_DP$FUS_general_pattern) # p-value = 0.4983

# 1.2
rh_FFA_DP <- FFA_DP %>% 
  filter(hemi =='rh')
rh_FFA_DP.tab <- table(rh_FFA_DP$group, rh_FFA_DP$FUS_general_pattern)
fisher.test(rh_FFA_DP$group, rh_FFA_DP$FUS_general_pattern) # p-value = 0.2643

# plot 
FFA_DP.pct <- read.csv("./data/FUS_patterns_pct.csv")

FFA_DP.pct$hemi <- gsub("lh", "left", FFA_DP.pct$hemi)
FFA_DP.pct$hemi <- gsub("rh", "right", FFA_DP.pct$hemi)
FFA_DP.pct$group <- gsub("Controls", "NTs", FFA_DP.pct$group)
FFA_DP.pct$pattern <- gsub("mFUS", "mFus", FFA_DP.pct$pattern)
FFA_DP.pct$pattern <- gsub("pFUS", "pFus", FFA_DP.pct$pattern)

FFA_DP.pct$pattern <- factor(FFA_DP.pct$pattern, 
                             levels = c("none", "pFus", "mFus", "continuous", "separate")
)

FFA_DP_pattern.plot <- FFA_DP.pct %>%
  ggplot(aes(x = hemi, y = amount, fill = pattern)
         ) + 
  
  # set bar and colors
  geom_col(color = 'black', alpha = 0.90) + 
  
  scale_fill_manual(breaks = c("separate", "continuous", "mFus", "pFus", "none"),
                     values = c("#d7191c", "#d95f02", "#08519c", "#33a02c", "white")
                    ) +
  
  # set labels
  labs(x = "hemisphere",
       y = "% of participants",
       fill = "Pattern") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, 
                     limits = c(0,100)
                     ) +
  
  facet_wrap(~ group)

# view 
FFA_DP_pattern.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_DP_pattern.plot.png",
                plot = FFA_DP_pattern.plot,
                device = "png",
                width = 5,
                height = 5,
                units = "in",
                dpi = "retina")


#### morph of FFAs in DPs vs NTs ####
FFA_morphology <- read.csv("./data/FFA_DP_metrics.csv")
FFA_demo <- read.csv("./data/dp_handedness.csv")
FFA_demo <- FFA_demo %>% dplyr::select(sub, handedness)
FFA_revision_metrics <- read.csv("./data/revision_metrics.csv")
FFA_revision_metrics.cortex <- FFA_revision_metrics %>% 
  subset(label == "cortex") %>% 
  dplyr::select(sub, hemi, total_surface_area_.mm.2.) %>% 
  rename(cortex_sa = total_surface_area_.mm.2.)

t.test(cortex_sa ~ group, FFA_morphology)
anova(lme(cortex_sa ~ group * hemi + Age + handedness + Sex, random = (~1|sub/hemi), 
                      data = FFA_morphology))

# view summary
dp.mfus.sa.mod.aov <- anova(dp.mfus.sa.mod)
dp.mfus.sa.mod.aov

FFA_DP2 <- merge(FFA_DP, FFA_demo, by = c("sub"))
FFA_DP2 %>% subset(hemi == "lh") %>% group_by(group) %>% summarise(mean = mean(Age), sd = sd(Age))
FFA_DP2 %>% subset(hemi == "lh") %>% group_by(group) %>% summarise(range(Age))
FFA_DP.lh <- FFA_DP2 %>% subset(hemi == "lh")
t.test(Age ~ group, FFA_DP.lh)
FFA_DP.lh.gender <- table(FFA_DP.lh$Sex, FFA_DP.lh$group)
chisq.test(FFA_DP.lh.gender)
FFA_DP.lh.handed <- table(FFA_DP.lh$group, FFA_DP.lh$handedness)
chisq.test(FFA_DP.lh.handed)
table(FFA_demo$group, FFA_demo$handedness)

FFA_morphology <- merge(FFA_morphology, FFA_revision_metrics.cortex , by = c("sub", "hemi")
)
FFA_morphology <- merge(FFA_morphology, FFA_DP , by = c("sub", "hemi")
)
FFA_morphology <- merge(FFA_morphology, FFA_demo , by = c("sub")
)

FFA_morphology$sub <- as.factor(FFA_morphology$sub)
FFA_morphology$label <- gsub("pFUS", "pFus", FFA_morphology$label)
FFA_morphology$label <- gsub("mFUS", "mFus", FFA_morphology$label)
FFA_morphology$group <- gsub("Controls", "NTs", FFA_morphology$group)
FFA_morphology$hemi <- gsub("lh", "left", FFA_morphology$hemi)
FFA_morphology$hemi <- gsub("rh", "right", FFA_morphology$hemi)


# surface area
FFA_morphology.mfus.sa.data <-  FFA_morphology %>% subset(label == "mFus") %>% 
  group_by(group, hemi) %>% summarise(mean = mean(total_surface_area_.mm.2., na.rm = T), 
                                                 sd = sd(total_surface_area_.mm.2., na.rm = T), 
                                                 n=n(), 
                                                 se=sd/sqrt(n)) %>% 
                                        rename(total_surface_area_.mm.2. = mean)

FFA_morphology.mFus <- FFA_morphology %>% subset(label == "mFus")



FFA_morphology_group_SA_mFus.plot <- ggplot(FFA_morphology.mFus, aes(x = hemi, 
                                                                     y = total_surface_area_.mm.2., 
                                                                     fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_morphology.mfus.sa.data,
           aes(x=hemi, y=total_surface_area_.mm.2., fill = group), 
           stat="identity", color = "black", alpha = .5) +

   geom_errorbar(data = FFA_morphology.mfus.sa.data,
                 aes(x=hemi, ymin=total_surface_area_.mm.2.-se, ymax=total_surface_area_.mm.2.+se), 
                 width=0, color="black", alpha = 1, linewidth = 1.5) +

   geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +

  #geom_boxplot(outlier.alpha = 0, alpha = .5, width = .3, position = dodge) + 
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#1f78b4", "#a6cee3")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "surface area",
       fill = "group",
       title = "mFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, 
                     limits = c(0,800)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 
  

FFA_morphology_group_SA_mFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_SA_mFus.plot2.png",
                plot = FFA_morphology_group_SA_mFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")
  
FFA_morphology.pFus.sa.data <-  FFA_morphology %>% subset(label == "pFus") %>% 
  group_by(group) %>% summarise(mean = mean(total_surface_area_.mm.2., na.rm = T), 
                                      sd = sd(total_surface_area_.mm.2., na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(total_surface_area_.mm.2. = mean)

FFA_morphology.pFus <- FFA_morphology %>% subset(label == "pFus")


FFA_morphology_group_SA_pFus.plot <- ggplot(FFA_morphology.pFus, aes(x = hemi, 
                                                       y = total_surface_area_.mm.2., 
                                                       fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_morphology.pFus.sa.data,
           aes(x=hemi, y=total_surface_area_.mm.2., fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_morphology.pFus.sa.data,
                aes(x=hemi, ymin=total_surface_area_.mm.2.-se, ymax=total_surface_area_.mm.2.+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#33a02c", "#b2df8a")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "surface area",
       fill = "group",
       title = "pFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, 
                     limits = c(0,800)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_morphology_group_SA_pFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_SA_pFus.plot.png",
                plot = FFA_morphology_group_SA_pFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

##### FFA sa model #####
FFA_morphology.na <- FFA_morphology %>% drop_na()
FFA_morphology.sa <- FFA_morphology %>% 
  select(sub, hemi, label, total_surface_area_.mm.2., group) %>% 
  pivot_wider(names_from = label, values_from = total_surface_area_.mm.2.)
FFA_morphology.sa$mFus2 <- FFA_morphology.sa$mFus
FFA_morphology.sa$pFus2 <- FFA_morphology.sa$pFus
FFA_morphology.sa["mFus2"][is.na(FFA_morphology.sa["mFus2"])] <- 0
FFA_morphology.sa["pFus2"][is.na(FFA_morphology.sa["pFus2"])] <- 0

FFA_morphology.sa <- FFA_morphology.sa %>% mutate(faces = mFus2 + pFus2)
t.test(faces ~ group, FFA_morphology.sa, var.equal = T)
FFA_morphology.sa %>% group_by(group) %>% summarise(mean = mean(faces), sd = sd(faces))


FFA_morphology.na.mfus <- FFA_morphology.na %>% subset(label == "mFus")
dp.mfus.sa.mod <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                     data = FFA_morphology.na.mfus)

# view summary
dp.mfus.sa.mod.aov <- anova(dp.mfus.sa.mod)
dp.mfus.sa.mod.aov
# extract p-values from ANOVA table
pvals.mfus.sa <- dp.mfus.sa.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.sa <- p.adjust(pvals.mfus.sa, method = "fdr")

# combine into a new table
results.mfus.sa <- cbind(dp.mfus.sa.mod.aov, pvals_fdr.mfus.sa)
results.mfus.sa

effectsize::eta_squared(dp.mfus.sa.mod.aov)

FFA_morphology.na.pfus <- FFA_morphology.na %>% subset(label == "pFus")
dp.pfus.sa.mod <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_morphology.na.pfus)

# view summary
dp.pfus.sa.mod.aov <- anova(dp.pfus.sa.mod)
dp.pfus.sa.mod.aov
effectsize::eta_squared(dp.pfus.sa.mod.aov)

# extract p-values from ANOVA table
pvals.pfus.sa <- dp.pfus.sa.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.sa <- p.adjust(pvals.pfus.sa, method = "fdr")

# combine into a new table
results.pfus.sa <- cbind(dp.pfus.sa.mod.aov, pvals_fdr.pfus.sa)
results.pfus.sa


dp.ffa.sa.mod <- lme(total_surface_area_.mm.2. ~ group * hemi * label, random = (~1|sub/hemi/label), 
                  data = FFA_morphology.na)

# view summary
dp.ffa.sa.mod.aov <- anova(dp.ffa.sa.mod)
dp.ffa.sa.mod.aov
effectsize::eta_squared(dp.ffa.sa.mod.aov)

FFA_morphology %>% group_by(label, group) %>% summarise(mean = mean(total_surface_area_.mm.2., na.rm = T), 
                                                 sd = sd(total_surface_area_.mm.2., na.rm = T), 
                                                 n=n(), 
                                                 se=sd/sqrt(n)
                                                 )


dp.ffa.sa.mod.aov.i <- emmeans::emmeans(dp.ffa.sa.mod, ~ group | label)
emmeans::contrast(dp.ffa.sa.mod.aov.i, method='pairwise', adjust="tukey")
# label = mFUS:
#   contrast       estimate   SE df t.ratio p.value
# Controls - DPs     77.4 25.8 43   3.002  0.0045
# 
# label = pFUS:
#   contrast       estimate   SE df t.ratio p.value
# Controls - DPs     16.3 28.7 43   0.567  0.5734


# thickness
FFA_morphology.mfus.ct.data <-  FFA_morphology.na %>% subset(label == "mFus") %>% 
  group_by(group, hemi) %>% summarise(mean = mean(cortical_thickness_mean, na.rm = T), 
                                      sd = sd(cortical_thickness_mean, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(cortical_thickness_mean = mean)

FFA_morphology_group_CT_mFus.plot <- ggplot(FFA_morphology.mFus, aes(x = hemi, 
                                                                     y = cortical_thickness_mean, 
                                                                     fill = group)) + 
  
  geom_bar(data = FFA_morphology.mfus.ct.data,
           aes(x=hemi, y=cortical_thickness_mean, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_morphology.mfus.ct.data,
                aes(x=hemi, ymin=cortical_thickness_mean-se, ymax=cortical_thickness_mean+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#1f78b4", "#a6cee3")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "cortical thickness mean",
       fill = "group",
       title = "mFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5,
                     limits = c(0, 4.5)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_morphology_group_CT_mFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_CT_mFus.plot.png",
                plot = FFA_morphology_group_CT_mFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_morphology.pfus.ct.data <-  FFA_morphology.na %>% subset(label == "pFus") %>% 
  group_by(group, hemi) %>% summarise(mean = mean(cortical_thickness_mean, na.rm = T), 
                                      sd = sd(cortical_thickness_mean, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(cortical_thickness_mean = mean)


FFA_morphology_group_CT_pFus.plot <- ggplot(FFA_morphology.pFus, aes(x = hemi, 
                                                                     y = cortical_thickness_mean, 
                                                                     fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_morphology.pfus.ct.data,
           aes(x=hemi, y=cortical_thickness_mean, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_morphology.pfus.ct.data,
                aes(x=hemi, ymin=cortical_thickness_mean-se, ymax=cortical_thickness_mean+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#33a02c", "#b2df8a")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "cortical thickness mean",
       fill = "group",
       title = "pFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5,
                     limits = c(0, 4.5)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_morphology_group_CT_pFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_CT_pFus.plot.png",
                plot = FFA_morphology_group_CT_pFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


dp.mfus.ct.mod <- lme(cortical_thickness_mean ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_morphology.na.mfus)

# view summary
dp.mfus.ct.mod.aov <- anova(dp.mfus.ct.mod)
dp.mfus.ct.mod.aov

# extract p-values from ANOVA table
pvals.mfus.ct <- dp.mfus.ct.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.ct <- p.adjust(pvals.mfus.ct, method = "fdr")

# combine into a new table
results.mfus.ct <- cbind(dp.mfus.ct.mod.aov, pvals_fdr.mfus.ct)
results.mfus.ct

effectsize::eta_squared(dp.mfus.ct.mod.aov)

dp.pfus.ct.mod <- lme(cortical_thickness_mean ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_morphology.na.pfus)

# view summary
dp.pfus.ct.mod.aov <- anova(dp.pfus.ct.mod)
dp.pfus.ct.mod.aov

# extract p-values from ANOVA table
pvals.pfus.ct <- dp.pfus.ct.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.ct <- p.adjust(pvals.pfus.ct, method = "fdr")

# combine into a new table
results.pfus.ct <- cbind(dp.mfus.ct.mod.aov, pvals_fdr.pfus.ct)
results.pfus.ct

effectsize::eta_squared(dp.pfus.ct.mod.aov)

dp.ffa.ct.mod <- lme(cortical_thickness_mean ~ group * hemi * label, random = (~1|sub/hemi/label), 
                     data = FFA_morphology.na)

# view summary
dp.ffa.ct.mod.aov <- anova(dp.ffa.ct.mod)
dp.ffa.ct.mod.aov
effectsize::eta_squared(dp.ffa.ct.mod.aov)
FFA_morphology %>% group_by(label) %>% summarise(mean = mean(cortical_thickness_mean, na.rm = T), 
                                                 sd = sd(cortical_thickness_mean, na.rm = T), 
                                                 n=n(), 
                                                 se=sd/sqrt(n)
)

#### FFA ~ cfmt ####
mfs_DPs <- read.csv("./data/mfs.csv")
DP_sample <- mfs_DPs %>% subset(hemi == "lh") %>% dplyr::select(sub, CFMT, Old.New, Famous.Face)
DP_sample <- DP_sample[1:47,]
DP_sample$sub <- as.factor(DP_sample$sub)

FFA_morphology2 <- merge(FFA_morphology, DP_sample, by = c("sub"))

FFA_morphology2_SA_mFus_CFMT.plot <- FFA_morphology2 %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = total_surface_area_.mm.2., y = CFMT)) + 
                      
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +

  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "surface area",
       y = "CFMT",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
                     #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_SA_mFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_SA_mFus_CFMT.plot.png",
                plot = FFA_morphology2_SA_mFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


FFA_morphology2_SA_pFus_CFMT.plot <- FFA_morphology2 %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = total_surface_area_.mm.2., y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "surface area",
       y = "CFMT",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_SA_pFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_SA_pFus_CFMT.plot.png",
                plot = FFA_morphology2_SA_pFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

##### running correlations #####
FFA_morphology2_mfus <- FFA_morphology2 %>% subset(label == "mFus")
FFA_morphology2_mfus.lh <- FFA_morphology2_mfus %>% subset(hemi == "left")
FFA_morphology2_mfus.rh <- FFA_morphology2_mfus %>% subset(hemi == "right")

FFA_morphology2_pfus <- FFA_morphology2 %>% subset(label == "pFus")
FFA_morphology2_pfus.lh <- FFA_morphology2_pfus %>% subset(hemi == "left")
FFA_morphology2_pfus.rh <- FFA_morphology2_pfus %>% subset(hemi == "right")

# both
FFA_morphology2_mfus.na <- FFA_morphology2_mfus %>% select(ID, total_surface_area_.mm.2., mean_face_activation, 
                                                           cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()
cor.test(FFA_morphology2_mfus.na$total_surface_area_.mm.2., FFA_morphology2_mfus.na$CFMT, method = "spearman", exact = FALSE) #


FFA_morphology2_pfus.na <- FFA_morphology2_pfus %>% select(total_surface_area_.mm.2., mean_face_activation, 
                                                           cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()
cor.test(FFA_morphology2_pfus.na$total_surface_area_.mm.2., FFA_morphology2_pfus.na$CFMT, method = "spearman", exact = FALSE)

# left
FFA_morphology2_mfus.lh.na <- FFA_morphology2_mfus.lh %>% select(total_surface_area_.mm.2., mean_face_activation, 
                                                           cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()

cor.test(FFA_morphology2_mfus.lh.na$total_surface_area_.mm.2., FFA_morphology2_mfus.lh.na$CFMT, method = "spearman", exact = FALSE) #

FFA_morphology2_pfus.lh.na <- FFA_morphology2_pfus.lh %>% select(total_surface_area_.mm.2., mean_face_activation, 
                                                           cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()
cor.test(FFA_morphology2_pfus.lh.na$total_surface_area_.mm.2., FFA_morphology2_pfus.lh.na$CFMT, method = "spearman", exact = FALSE)


# right
FFA_morphology2_mfus.rh.na <- FFA_morphology2_mfus.rh %>% select(total_surface_area_.mm.2., mean_face_activation, 
                                                                 cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()

cor.test(FFA_morphology2_mfus.rh.na$total_surface_area_.mm.2., FFA_morphology2_mfus.rh.na$CFMT, method = "spearman", exact = FALSE) # 

FFA_morphology2_pfus.rh.na <- FFA_morphology2_pfus.rh %>% select(total_surface_area_.mm.2., mean_face_activation, 
                                                                 cortical_thickness_mean, total_gray_matter_volume_.mm.3., CFMT) %>% drop_na()
cor.test(FFA_morphology2_pfus.rh.na$total_surface_area_.mm.2., FFA_morphology2_pfus.rh.na$CFMT, method = "spearman", exact = FALSE)


#### crown analysis ####
FFA_sulc_morphology <- read.csv("./data/FFA_crown.csv")
FFA_sulc_morphology <- merge(FFA_sulc_morphology, FFA_revision_metrics.cortex , by = c("sub", "hemi")
)
FFA_sulc_morphology <- merge(FFA_sulc_morphology, FFA_DP , by = c("sub", "hemi")
)
FFA_sulc_morphology <- merge(FFA_sulc_morphology, FFA_demo , by = c("sub")
)
FFA_sulc_morphology <- merge(FFA_sulc_morphology, FFA_revision_metrics.cortex , by = c("sub", "hemi")
)
FFA_sulc_morphology$sub <- as.factor(FFA_sulc_morphology$sub)

FFA_sulc_morphology$label <- gsub("pFUS", "pFus", FFA_sulc_morphology$label)
FFA_sulc_morphology$label <- gsub("mFUS", "mFus", FFA_sulc_morphology$label)
FFA_sulc_morphology$group <- gsub("Controls", "NTs", FFA_sulc_morphology$group)
FFA_sulc_morphology$hemi <- gsub("lh", "left", FFA_sulc_morphology$hemi)
FFA_sulc_morphology$hemi <- gsub("rh", "right", FFA_sulc_morphology$hemi)

FFA_sulc_morphology.na <- FFA_sulc_morphology %>% drop_na()

FFA_sulc_morphology.na.mfus <- FFA_sulc_morphology.na %>% subset(label == "mFus")
dp.mfus.gc.mod <- lme(gyral_crown ~ group * hemi + Age + Sex + handedness + cortex_sa.x, random = (~1|sub/hemi), 
                      data = FFA_sulc_morphology.na.mfus)

# view summary
dp.mfus.gc.mod.aov <- anova(dp.mfus.gc.mod)
dp.mfus.gc.mod.aov

# extract p-values from ANOVA table
pvals.mfus.gc <- dp.mfus.gc.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.gc <- p.adjust(pvals.mfus.gc, method = "fdr")

# combine into a new table
results.mfus.gc <- cbind(dp.mfus.gc.mod.aov, pvals_fdr.mfus.gc)
results.mfus.gc

effectsize::eta_squared(dp.mfus.gc.mod.aov)

FFA_sulc_morphology.na.pfus <- FFA_sulc_morphology.na %>% subset(label == "pFus")
dp.pfus.gc.mod <- lme(gyral_crown ~ group * hemi + Age + Sex + handedness + cortex_sa.x, random = (~1|sub/hemi), 
                      data = FFA_sulc_morphology.na.pfus)

# view summary
dp.pfus.gc.mod.aov <- anova(dp.pfus.gc.mod)
dp.pfus.gc.mod.aov

# extract p-values from ANOVA table
pvals.pfus.gc <- dp.pfus.gc.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.gc <- p.adjust(pvals.pfus.gc, method = "fdr")

# combine into a new table
results.pfus.gc <- cbind(dp.pfus.gc.mod.aov, pvals_fdr.pfus.gc)
results.pfus.gc

effectsize::eta_squared(dp.pfus.gc.mod.aov)

FFA_sulc_morphology %>% group_by(label, group) %>% summarise(mean = mean(gyral_crown, na.rm = T), 
                                                      sd = sd(gyral_crown, na.rm = T), 
                                                      n=n(), 
                                                      se=sd/sqrt(n)
)

dp.ffa.gc.mod.i <- emmeans::emmeans(dp.ffa.gc.mod, ~ group | label)
emmeans::contrast(dp.ffa.gc.mod.i, method='pairwise', adjust="tukey")
# label = mFus:
#   contrast  estimate    SE df t.ratio p.value
# DPs - NTs     1.49 0.540 43   2.767  0.0083
# 
# label = pFus:
#   contrast  estimate    SE df t.ratio p.value
# DPs - NTs     2.29 0.592 43   3.871  0.0004


FFA_sulc_morphology.GC.data <-  FFA_sulc_morphology.na %>% subset(label == "mFus") %>% 
  group_by(group, hemi) %>% summarise(mean = mean(gyral_crown, na.rm = T), 
                                      sd = sd(gyral_crown, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(gyral_crown = mean)

FFA_sulc_morphology.mFus <- FFA_sulc_morphology %>% subset(label == "mFus")



FFA_morphology_group_GC_mFus.plot <- ggplot(FFA_sulc_morphology.mFus, aes(x = hemi, 
                                                                     y = gyral_crown, 
                                                                     fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_sulc_morphology.GC.data,
           aes(x=hemi, y=gyral_crown, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_sulc_morphology.GC.data,
                aes(x=hemi, ymin=gyral_crown-se, ymax=gyral_crown+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  #geom_boxplot(outlier.alpha = 0, alpha = .5, width = .3, position = dodge) + 
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#1f78b4", "#a6cee3")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "gyral height",
       fill = "group",
       title = "mFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10, limits = c(-12, 12) ) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_morphology_group_GC_mFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_GC_mFus.plot.png",
                plot = FFA_morphology_group_GC_mFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_sulc_morphology.GC.data <-  FFA_sulc_morphology.na %>% subset(label == "pFus") %>% 
  group_by(group, hemi) %>% summarise(mean = mean(gyral_crown, na.rm = T), 
                                      sd = sd(gyral_crown, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(gyral_crown = mean)

FFA_sulc_morphology.pFus <- FFA_sulc_morphology %>% subset(label == "pFus")


FFA_morphology_group_GC_pFus.plot <- ggplot(FFA_sulc_morphology.pFus, aes(x = hemi, 
                                                                     y = gyral_crown, 
                                                                     fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_sulc_morphology.GC.data,
           aes(x=hemi, y=gyral_crown, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_sulc_morphology.GC.data,
                aes(x=hemi, ymin=gyral_crown-se, ymax=gyral_crown+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#33a02c", "#b2df8a")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "gyral crown",
       fill = "group",
       title = "pFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 10, limits = c(-12, 12)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_morphology_group_GC_pFus.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology_group_GC_pFus.plot.png",
                plot = FFA_morphology_group_GC_pFus.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


#### gyral crown ~ CFMT ####
FFA_sulc_morphology2 <-  FFA_sulc_morphology %>% dplyr::select(sub,hemi,label, group, gyral_crown, total_surface_area_.mm.2., cortical_thickness_mean)
DP_sample2 <- mfs_DPs %>% subset(hemi == "lh") %>% dplyr::select(sub, CFMT)
DP_sample2 <- DP_sample2[1:47,]
DP_sample2$sub <- as.factor(DP_sample2$sub)

FFA_sulc_morphology2 <- merge(FFA_sulc_morphology2, DP_sample2, by = c("sub"))
FFA_sulc_morphology2$group <- gsub("Controls", "NTs", FFA_sulc_morphology2$group)

FFA_sulc_morphology2_mFus <- FFA_sulc_morphology2 %>% subset(label == "mFus")
FFA_sulc_morphology2_pFus <- FFA_sulc_morphology2 %>% subset(label == "pFus")

cor.test(FFA_sulc_morphology2_mFus$gyral_crown, FFA_sulc_morphology2_mFus$CFMT, method = "spearman", exact = F) #

FFA_sulc_morphology2_mFus.lh <- FFA_sulc_morphology2_mFus %>% subset(hemi == "left")
cor.test(FFA_sulc_morphology2_mFus.lh$gyral_crown, FFA_sulc_morphology2_mFus.lh$CFMT, method = "spearman", exact = F) #
FFA_sulc_morphology2_mFus.rh <- FFA_sulc_morphology2_mFus %>% subset(hemi == "right")
cor.test(FFA_sulc_morphology2_mFus.rh$gyral_crown, FFA_sulc_morphology2_mFus.rh$CFMT, method = "spearman", exact = F) #

cor.test(FFA_sulc_morphology2_pFus$gyral_crown, FFA_sulc_morphology2_pFus$CFMT, method = "spearman", exact = F) #

FFA_sulc_morphology2_pFus.lh <- FFA_sulc_morphology2_pFus %>% subset(hemi == "left")
cor.test(FFA_sulc_morphology2_pFus.lh$gyral_crown, FFA_sulc_morphology2_pFus.lh$CFMT, method = "spearman", exact = F) #
FFA_sulc_morphology2_pFus.rh <- FFA_sulc_morphology2_pFus %>% subset(hemi == "right")
cor.test(FFA_sulc_morphology2_pFus.rh$gyral_crown, FFA_sulc_morphology2_pFus.rh$CFMT, method = "spearman", exact = F) #

FFA_sulc_morphology2_pFus <- FFA_sulc_morphology2 %>% subset(label == "pFus")

FFA_morphology2_GC_mFus_CFMT.plot <- FFA_sulc_morphology2 %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = gyral_crown, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "gyral height",
       y = "CFMT",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(-10, 5)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_GC_mFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_GC_mFus_CFMT.plot.png",
                plot = FFA_morphology2_GC_mFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_morphology2_GC_pFus_CFMT.plot <- FFA_sulc_morphology2 %>% 
  subset(label == "pFus") %>%
  
  ggplot(aes(x = gyral_crown, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "gyral height",
       y = "CFMT",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(-10, 5)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_GC_pFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_GC_pFus_CFMT.plot.png",
                plot = FFA_morphology2_GC_pFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

cor.test(FFA_sulc_morphology2_mFus$cortical_thickness_mean, FFA_sulc_morphology2_mFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_sulc_morphology2_pFus$cortical_thickness_mean, FFA_sulc_morphology2_pFus$CFMT, method = "spearman", exact = F) #

FFA_morphology2_CT_mFus_CFMT.plot <- FFA_sulc_morphology2 %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = cortical_thickness_mean, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "CFMT",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(1.5,4.5)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_CT_mFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_CT_mFus_CFMT.plot.png",
                plot = FFA_morphology2_CT_mFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_morphology2_CT_pFus_CFMT.plot <- FFA_sulc_morphology2 %>% 
  subset(label == "pFus") %>%
  
  ggplot(aes(x = cortical_thickness_mean, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical_thickness_mean",
       y = "CFMT",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(1.5, 4.5)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_morphology2_CT_pFus_CFMT.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_morphology2_CT_pFus_CFMT.plot.png",
                plot = FFA_morphology2_CT_pFus_CFMT.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

#### metric correlations ####
FFA_all_morph <- merge(FFA_morphology, FFA_sulc_morphology2, by = c("sub", "hemi", "label"))
FFA_all_morph_mFus <- FFA_all_morph %>% subset(label == "mFus")
FFA_all_morph_pFus <- FFA_all_morph %>% subset(label == "pFus")

# mfus
# across groups + hemis -- do separate too?
FFA_all_morph_mFus.lh <- FFA_all_morph_mFus %>% subset(hemi == "left")
FFA_all_morph_mFus.rh <- FFA_all_morph_mFus %>% subset(hemi == "right")

FFA_all_morph_mFus.lh.dp <- FFA_all_morph_mFus %>% subset(hemi == "left" & group == "DPs")
FFA_all_morph_mFus.lh.nt <- FFA_all_morph_mFus %>% subset(hemi == "left" & group == "NTs")

cor.test(FFA_all_morph_mFus$gyral_crown, FFA_all_morph_mFus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph_mFus$cortical_thickness_mean, FFA_all_morph_mFus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph_mFus.lh$gyral_crown, FFA_all_morph_mFus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #
cor.test(FFA_all_morph_mFus.rh$gyral_crown, FFA_all_morph_mFus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

FFA_all_morph_mFus2 <- FFA_all_morph_mFus %>% subset(!is.na(CFMT))

mod.D <- lm(CFMT ~ total_surface_area_.mm.2., FFA_all_morph_mFus2)
summary(mod.D)
mod.M <- lm(gyral_crown ~ total_surface_area_.mm.2., FFA_all_morph_mFus2)
summary(mod.M)
mod.Y <- lm(CFMT ~ total_surface_area_.mm.2. + gyral_crown, FFA_all_morph_mFus2)
summary(mod.Y)

set.seed(1)
results <- mediate(mod.M, mod.Y, treat = 'total_surface_area_.mm.2.', mediator = 'gyral_crown',
                   boot = TRUE, sims = 1000)
summary(results)

mod.D <- lm(CFMT ~ total_surface_area_.mm.2., FFA_all_morph_mFus2)
summary(mod.D)
mod.M <- lm(sulcal_pit ~ total_surface_area_.mm.2., FFA_all_morph_mFus2)
summary(mod.M)
mod.Y <- lm(CFMT ~ total_surface_area_.mm.2. + sulcal_pit, FFA_all_morph_mFus2)
summary(mod.Y)


# pfus
# across groups + hemis
cor.test(FFA_all_morph_pFus$cortical_thickness_mean, FFA_all_morph_pFus$gyral_crown) # 

cor.test(FFA_all_morph_pFus$cortical_thickness_mean, FFA_all_morph_pFus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph_pFus$gyral_crown, FFA_all_morph_pFus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

FFA_all_morph_pFus.lh <- FFA_all_morph_pFus %>% subset(hemi == "left")
FFA_all_morph_pFus.rh <- FFA_all_morph_pFus %>% subset(hemi == "right")

cor.test(FFA_all_morph_pFus.lh$gyral_crown, FFA_all_morph_pFus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) # 
cor.test(FFA_all_morph_pFus.rh$gyral_crown, FFA_all_morph_pFus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) # 

cor.test(FFA_all_morph_pFus.lh$cortical_thickness_mean, FFA_all_morph_pFus.lh$total_surface_area_.mm.2.) # 
cor.test(FFA_all_morph_pFus.rh$cortical_thickness_mean, FFA_all_morph_pFus.rh$total_surface_area_.mm.2.) # 


FFA_all_morph_mFus_gc_sa.plot <- FFA_all_morph %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = gyral_crown, y = total_surface_area_.mm.2.)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "gyral height",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 750)) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) + #  limits = c(-8, 6)
  guides(fill = "none", shape = "none", color = "none")

FFA_all_morph_mFus_gc_sa.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_all_morph_mFus_gc_sa.plot.png",
                plot = FFA_all_morph_mFus_gc_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_all_morph_mFus_ct_sa.plot <- FFA_all_morph %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = cortical_thickness_mean, y = total_surface_area_.mm.2.)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 750)) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 3, limits = c(1.5, 4)) + #  limits = c(-8, 6)
  guides(fill = "none", shape = "none", color = "none")

FFA_all_morph_mFus_ct_sa.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_all_morph_mFus_ct_sa.plot.png",
                plot = FFA_all_morph_mFus_ct_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_all_morph_pFus_gc_sa.plot <- FFA_all_morph %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = gyral_crown, y = total_surface_area_.mm.2.)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "gyral height",
       y = "surface area",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 750)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) + # , limits = c(-11, 5)
  guides(fill = "none", shape = "none", color = "none")

FFA_all_morph_pFus_gc_sa.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_all_morph_pFus_gc_sa.plot.png",
                plot = FFA_all_morph_pFus_gc_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_all_morph.na <- FFA_all_morph %>% drop_na()

FFA_all_morph_pFus_ct_sa.plot <- FFA_all_morph %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = cortical_thickness_mean, y = total_surface_area_.mm.2.)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "surface area",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 750)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 3, limits = c(1.5, 4)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_all_morph_pFus_ct_sa.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_all_morph_pFus_ct_sa.plot2.png",
                plot = FFA_all_morph_pFus_ct_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_all_morph_mFus_ct_sa.plot <- FFA_all_morph.na %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = cortical_thickness_mean, y = total_surface_area_.mm.2.)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 750)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(1.5, 4)) +
  guides(fill = "none", shape = "none", color = "none")

FFA_all_morph_mFus_ct_sa.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_all_morph_mFus_ct_sa.plot.png",
                plot = FFA_all_morph_mFus_ct_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


##### checking outlier relationships #####
# Identify outliers in x (±3 sd)
FFA_sulc_morphology2.1 <- FFA_sulc_morphology2 %>% dplyr::select(sub, group, hemi, label, gyral_crown, total_surface_area_.mm.2., CFMT)

# Compute IQR thresholds for mfus gc
FFA_sulc_morphology2.1.mfus <- FFA_sulc_morphology2.1 %>% subset(label == "mFus") %>% drop_na() 
Q1 <- quantile(FFA_sulc_morphology2.1.mfus$gyral_crown, 0.25)
Q3 <- quantile(FFA_sulc_morphology2.1.mfus$gyral_crown, 0.75)
IQR_val <- IQR(FFA_sulc_morphology2.1.mfus$gyral_crown)

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Flag outliers using IQR rule
FFA_sulc_morphology2.1.mfus <- FFA_sulc_morphology2.1.mfus %>%
  mutate(outlier = gyral_crown < lower_bound | gyral_crown > upper_bound)

# Plot: highlight outliers in red
mfus_gc_outlier.plot <- ggplot(FFA_sulc_morphology2.1.mfus, aes(x = gyral_crown, y = CFMT)) +
  geom_jitter(aes(fill = outlier, shape = group, alpha = hemi), size = 3) +
  
  geom_vline(xintercept = lower_bound, linetype = "dotted", color = "red") +
  geom_vline(xintercept = upper_bound, linetype = "dotted", color = "red") +
  
  scale_fill_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  scale_alpha_manual(breaks = c("left", "right"), values = c(.75, .5)) + 
  labs(x = "gyral height",
       y = "CFMT",
       fill = "hemisphere") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none", alpha = "none")

mfus_gc_outlier.plot

ggplot2::ggsave(filename = "~/Downloads/mfus_gc_outlier.plot.png",
                plot = mfus_gc_outlier.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

# Compute IQR thresholds for pfus gc
FFA_sulc_morphology2.1.pfus <- FFA_sulc_morphology2.1 %>% subset(label == "pFus") %>% drop_na() 
Q1 <- quantile(FFA_sulc_morphology2.1.pfus$gyral_crown, 0.25)
Q3 <- quantile(FFA_sulc_morphology2.1.pfus$gyral_crown, 0.75)
IQR_val <- IQR(FFA_sulc_morphology2.1.pfus$gyral_crown)

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Flag outliers using IQR rule
FFA_sulc_morphology2.1.pfus <- FFA_sulc_morphology2.1.pfus %>%
  mutate(outlier = gyral_crown < lower_bound | gyral_crown > upper_bound)

# Plot: highlight outliers in red
pfus_gc_outlier.plot <- ggplot(FFA_sulc_morphology2.1.pfus, aes(x = gyral_crown, y = CFMT)) +
  geom_jitter(aes(fill = outlier, shape = group, alpha = hemi), size = 3) +
  
  geom_vline(xintercept = lower_bound, linetype = "dotted", color = "red") +
  geom_vline(xintercept = upper_bound, linetype = "dotted", color = "red") +
  
  scale_fill_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  scale_alpha_manual(breaks = c("left", "right"), values = c(.75, .5)) + 
  labs(x = "gyral height",
       y = "CFMT",
       fill = "hemisphere") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none", alpha = "none")

pfus_gc_outlier.plot

ggplot2::ggsave(filename = "~/Downloads/pfus_gc_outlier.plot.png",
                plot = pfus_gc_outlier.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


# Identify outliers for surfact area
# Compute IQR thresholds for mfus gc
# Compute IQR thresholds for mfus gc
Q1 <- quantile(FFA_sulc_morphology2.1.mfus$total_surface_area_.mm.2., 0.25)
Q3 <- quantile(FFA_sulc_morphology2.1.mfus$total_surface_area_.mm.2., 0.75)
IQR_val <- IQR(FFA_sulc_morphology2.1.mfus$total_surface_area_.mm.2.)

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Flag outliers using IQR rule
FFA_sulc_morphology2.1.mfus <- FFA_sulc_morphology2.1.mfus %>%
  mutate(sa_outlier = total_surface_area_.mm.2. < lower_bound | total_surface_area_.mm.2. > upper_bound)

# Plot: highlight outliers in red
mfus_sa_outlier.plot <- ggplot(FFA_sulc_morphology2.1.mfus, aes(x = total_surface_area_.mm.2., y = CFMT)) +
  geom_jitter(aes(fill = sa_outlier, shape = group, alpha = hemi), size = 3) +
  
  geom_vline(xintercept = lower_bound, linetype = "dotted", color = "red") +
  geom_vline(xintercept = upper_bound, linetype = "dotted", color = "red") +
  
  scale_fill_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  scale_alpha_manual(breaks = c("left", "right"), values = c(.75, .5)) + 
  labs(x = "surface area",
       y = "CFMT",
       fill = "hemisphere") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none", alpha = "none")

mfus_sa_outlier.plot

ggplot2::ggsave(filename = "~/Downloads/mfus_sa_outlier.plot.png",
                plot = mfus_sa_outlier.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

# Compute IQR thresholds for pfus gc
Q1 <- quantile(FFA_sulc_morphology2.1.pfus$total_surface_area_.mm.2., 0.25)
Q3 <- quantile(FFA_sulc_morphology2.1.pfus$total_surface_area_.mm.2., 0.75)
IQR_val <- IQR(FFA_sulc_morphology2.1.pfus$total_surface_area_.mm.2.)

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Flag outliers using IQR rule
FFA_sulc_morphology2.1.pfus <- FFA_sulc_morphology2.1.pfus %>%
  mutate(sa_outlier = total_surface_area_.mm.2. < lower_bound | total_surface_area_.mm.2. > upper_bound)

# Plot: highlight outliers in red
pfus_sa_outlier.plot <- ggplot(FFA_sulc_morphology2.1.pfus, aes(x = total_surface_area_.mm.2., y = CFMT)) +
  geom_jitter(aes(fill = sa_outlier, shape = group, alpha = hemi), size = 3) +
  
  geom_vline(xintercept = lower_bound, linetype = "dotted", color = "red") +
  geom_vline(xintercept = upper_bound, linetype = "dotted", color = "red") +
  
  scale_fill_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  scale_alpha_manual(breaks = c("left", "right"), values = c(.75, .5)) + 
  labs(x = "surface area",
       y = "CFMT",
       fill = "hemisphere") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none", alpha = "none")

pfus_sa_outlier.plot

ggplot2::ggsave(filename = "~/Downloads/pfus_sa_outlier.plot.png",
                plot = pfus_sa_outlier.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

#### HCP ####
HCP_FFA_SA <- read.csv("./FFA_metrics/sulc-corr-others/FFA_va.csv")
HCP_FFA_SA <- HCP_FFA_SA %>% rename(lh_pFus_sa = lh_pFus,
                                    lh_mFus_sa = lh_mFus,
                                    rh_pFus_sa = rh_pFus,
                                    rh_mFus_sa = rh_mFus
  
) 
HCP_FFA_SA <- tibble::rowid_to_column(HCP_FFA_SA, "ID")
HCP_FFA_SA.lh <- HCP_FFA_SA %>% select(ID, lh_pFus_sa, lh_mFus_sa) %>% 
  pivot_longer(cols = c(lh_pFus_sa, lh_mFus_sa), names_to = "label", values_to = "surface_area") %>%
  mutate(label = recode(label, lh_pFus_sa = "pFus", lh_mFus_sa = "mFus"), hemi = "left")

HCP_FFA_SA.rh <- HCP_FFA_SA %>% select(ID, rh_pFus_sa, rh_mFus_sa) %>% 
  pivot_longer(cols = c(rh_pFus_sa, rh_mFus_sa), names_to = "label", values_to = "surface_area") %>%
  mutate(label = recode(label, rh_pFus_sa = "pFus", rh_mFus_sa = "mFus"), hemi = "right")

HCP_FFA_SA_long <- rbind(HCP_FFA_SA.lh, HCP_FFA_SA.rh)

HCP_FFA_gc <- read.csv("./FFA_metrics/sulc-corr-others/FFA_gyralCrown.csv")
HCP_FFA_gc <- HCP_FFA_gc %>% rename(lh_pFus_gc = lh_pFus,
                                    lh_mFus_gc = lh_mFus,
                                    rh_pFus_gc = rh_pFus,
                                    rh_mFus_gc = rh_mFus
                                    
)
HCP_FFA_gc <- tibble::rowid_to_column(HCP_FFA_gc, "ID")
HCP_FFA_gc.lh <- HCP_FFA_gc %>% select(ID, lh_pFus_gc, lh_mFus_gc) %>% 
  pivot_longer(cols = c(lh_pFus_gc, lh_mFus_gc), names_to = "label", values_to = "gyral_crown") %>%
  mutate(label = recode(label, lh_pFus_gc = "pFus", lh_mFus_gc = "mFus"), hemi = "left")

HCP_FFA_gc.rh <- HCP_FFA_gc %>% select(ID, rh_pFus_gc, rh_mFus_gc) %>% 
  pivot_longer(cols = c(rh_pFus_gc, rh_mFus_gc), names_to = "label", values_to = "gyral_crown") %>%
  mutate(label = recode(label, rh_pFus_gc = "pFus", rh_mFus_gc = "mFus"), hemi = "right")

HCP_FFA_gc_long <- rbind(HCP_FFA_gc.lh, HCP_FFA_gc.rh)


HCP_FFA_sp <- read.csv("./FFA_metrics/sulc-corr-others/FFA_sulcBtm.csv")
HCP_FFA_sp <- HCP_FFA_sp %>% rename(lh_pFus_sp = lh_pFus,
                                    lh_mFus_sp = lh_mFus,
                                    rh_pFus_sp = rh_pFus,
                                    rh_mFus_sp = rh_mFus
                                    
)
HCP_FFA_sp <- tibble::rowid_to_column(HCP_FFA_sp, "ID")

HCP_FFA_sp.lh <- HCP_FFA_sp %>% select(ID, lh_pFus_sp, lh_mFus_sp) %>% 
  pivot_longer(cols = c(lh_pFus_sp, lh_mFus_sp), names_to = "label", values_to = "sulcal_pit") %>%
  mutate(label = recode(label, lh_pFus_sp = "pFus", lh_mFus_sp = "mFus"), hemi = "left")

HCP_FFA_sp.rh <- HCP_FFA_sp %>% select(ID, rh_pFus_sp, rh_mFus_sp) %>% 
  pivot_longer(cols = c(rh_pFus_sp, rh_mFus_sp), names_to = "label", values_to = "sulcal_pit") %>%
  mutate(label = recode(label, rh_pFus_sp = "pFus", rh_mFus_sp = "mFus"), hemi = "right")

HCP_FFA_sp_long <- rbind(HCP_FFA_sp.lh, HCP_FFA_sp.rh)


HCP_FFA_ct <- read.csv("./FFA_metrics/sulc-corr-others/FFA_thickness.csv")
HCP_FFA_ct <- HCP_FFA_ct %>% rename(lh_pFus_ct = lh_pFus,
                                    lh_mFus_ct = lh_mFus,
                                    rh_pFus_ct = rh_pFus,
                                    rh_mFus_ct = rh_mFus
                                    
)
HCP_FFA_ct <- tibble::rowid_to_column(HCP_FFA_ct, "ID")

HCP_FFA_ct.lh <- HCP_FFA_ct %>% dplyr::select(ID, lh_pFus_ct, lh_mFus_ct) %>% 
  pivot_longer(cols = c(lh_pFus_ct, lh_mFus_ct), names_to = "label", values_to = "cortical_thickness") %>%
  mutate(label = recode(label, lh_pFus_ct = "pFus", lh_mFus_ct = "mFus"), hemi = "left")

HCP_FFA_ct.rh <- HCP_FFA_ct %>% dplyr::select(ID, rh_pFus_ct, rh_mFus_ct) %>% 
  pivot_longer(cols = c(rh_pFus_ct, rh_mFus_ct), names_to = "label", values_to = "cortical_thickness") %>%
  mutate(label = recode(label, rh_pFus_ct = "pFus", rh_mFus_ct = "mFus"), hemi = "right")

HCP_FFA_ct_long <- rbind(HCP_FFA_ct.lh, HCP_FFA_ct.rh)

HCP_FFA_metrics_wide_full <- merge(HCP_FFA_SA, HCP_FFA_gc, by = "ID")
HCP_FFA_metrics_wide_full <- merge(HCP_FFA_metrics_wide_full, HCP_FFA_sp, by = "ID")
HCP_FFA_metrics_wide_full <- merge(HCP_FFA_metrics_wide_full, HCP_FFA_ct, by = "ID")


HCP_FFA_metrics_long_full <- merge(HCP_FFA_SA_long, HCP_FFA_gc_long, by = c("ID", "label", "hemi"))
HCP_FFA_metrics_long_full <- merge(HCP_FFA_metrics_long_full, HCP_FFA_sp_long, by = c("ID", "label", "hemi"))
HCP_FFA_metrics_long_full <- merge(HCP_FFA_metrics_long_full, HCP_FFA_ct_long, by = c("ID", "label", "hemi"))

HCP_FFA_metrics_long_full.mfus <- HCP_FFA_metrics_long_full %>% subset(label == "mFus")

HCP_FFA_metrics_long_full.pfus <- HCP_FFA_metrics_long_full %>% subset(label == "pFus")

# gc ~ sa
cor.test(HCP_FFA_metrics_long_full.mfus$gyral_crown, HCP_FFA_metrics_long_full.mfus$surface_area, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$lh_mFus_gc, HCP_FFA_metrics_wide_full$lh_mFus_sa, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$rh_mFus_gc, HCP_FFA_metrics_wide_full$rh_mFus_sa, method = "spearman", exact = FALSE)

cor.test(HCP_FFA_metrics_long_full.pfus$gyral_crown, HCP_FFA_metrics_long_full.pfus$surface_area, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$lh_pFus_gc, HCP_FFA_metrics_wide_full$lh_pFus_sa, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$rh_pFus_gc, HCP_FFA_metrics_wide_full$rh_pFus_sa, method = "spearman", exact = FALSE)

# ct ~ sa
cor.test(HCP_FFA_metrics_long_full.mfus$cortical_thickness, HCP_FFA_metrics_long_full.mfus$surface_area, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$lh_mFus_ct, HCP_FFA_metrics_wide_full$lh_mFus_sa, method = "spearman", exact = FALSE)
cor.test(HCP_FFA_metrics_wide_full$rh_mFus_ct, HCP_FFA_metrics_wide_full$rh_mFus_sa, method = "spearman", exact = FALSE)

cor.test(HCP_FFA_metrics_long_full.pfus$cortical_thickness, HCP_FFA_metrics_long_full.pfus$surface_area)
cor.test(HCP_FFA_metrics_wide_full$lh_pFus_ct, HCP_FFA_metrics_wide_full$lh_pFus_sa)
cor.test(HCP_FFA_metrics_wide_full$rh_pFus_ct, HCP_FFA_metrics_wide_full$rh_pFus_sa)

HCP_FFA_metrics_long_full.mfus$surface_area

HCP_FFA_all_morph_mFus_gc_sa.plot <- HCP_FFA_metrics_long_full.mfus %>% 
  
  ggplot(aes(x = gyral_crown, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), #method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", #method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  # set labels
  labs(x = "gyral height",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) + #  limits = c(-8, 6)
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_mFus_gc_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_mFus_gc_sa.plot.png",
                plot = HCP_FFA_all_morph_mFus_gc_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


HCP_FFA_all_morph_pFus_gc_sa.plot <- HCP_FFA_metrics_long_full.pfus %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = gyral_crown, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  # set labels
  labs(x = "gyral height",
       y = "surface area",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) + # , limits = c(-11, 5)
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_pFus_gc_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_pFus_gc_sa.plot.png",
                plot = HCP_FFA_all_morph_pFus_gc_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


HCP_FFA_all_morph_mFus_sp_sa.plot <- HCP_FFA_metrics_long_full.mfus %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = sulcal_pit, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  # set labels
  labs(x = "sulcal pit",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) + # , limits = c(0, 750)
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) + # , limits = c(-8, 6)
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_mFus_sp_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_mFus_sp_sa.plot.png",
                plot = HCP_FFA_all_morph_mFus_sp_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


HCP_FFA_all_morph_pFus_sp_sa.plot <- HCP_FFA_metrics_long_full.pfus %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = sulcal_pit, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "sulcal pit",
       y = "surface area",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_pFus_sp_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_pFus_sp_sa.plot.png",
                plot = HCP_FFA_all_morph_pFus_sp_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


HCP_FFA_all_morph_mFus_ct_sa.plot <- HCP_FFA_metrics_long_full.mfus %>% 
  
  subset(label == "mFus") %>%
  
  ggplot(aes(x = cortical_thickness, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#1f78b4", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "surface area",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) + # , limits = c(0, 750)
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 4.5)) + # , limits = c(-8, 6)
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_mFus_ct_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_mFus_ct_sa.plot2.png",
                plot = HCP_FFA_all_morph_mFus_ct_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


HCP_FFA_all_morph_pFus_ct_sa.plot <- HCP_FFA_metrics_long_full.pfus %>% 
  
  subset(label == "pFus") %>%
  
  ggplot(aes(x = cortical_thickness, y = surface_area)) + 
  
  geom_jitter(aes(fill = hemi), shape = 21, alpha = .3) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#33a02c", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  # set labels
  labs(x = "cortical thickness mean",
       y = "surface area",
       fill = "hemisphere",
       title = "pFus") +
  
  theme_classic() +
  
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(0, 1200)) + 
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5, limits = c(0, 4.5)) +
  guides(fill = "none", color = "none")

HCP_FFA_all_morph_pFus_ct_sa.plot

ggplot2::ggsave(filename = "~/Downloads/HCP_FFA_all_morph_pFus_ct_sa.plot2.png",
                plot = HCP_FFA_all_morph_pFus_ct_sa.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


#### revision analyses ####

##### removing DP that did did not reach −2 SD #####
FFA_morphology.rev <-FFA_morphology %>% filter(sub != "45")

FFA_morphology.rev.na <- FFA_morphology.rev %>% drop_na()
FFA_morphology.rev.na.mfus <- FFA_morphology.rev.na %>% subset(label == "mFus")
dp.mfus.sa.mod2 <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_morphology.rev.na.mfus)

# view summary
dp.mfus.sa.mod2.aov <- anova(dp.mfus.sa.mod2)
dp.mfus.sa.mod2.aov
# extract p-values from ANOVA table
pvals.mfus.sa2 <- dp.mfus.sa.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.sa2 <- p.adjust(pvals.mfus.sa2, method = "fdr")

# combine into a new table
results.mfus.sa2 <- cbind(dp.mfus.sa.mod2.aov, pvals_fdr.mfus.sa2)
results.mfus.sa2

effectsize::eta_squared(dp.mfus.sa.mod2.aov)

dp.mfus.ct.mod2 <- lme(cortical_thickness_mean ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                       data = FFA_morphology.rev.na.mfus)

# view summary
dp.mfus.ct.mod2.aov <- anova(dp.mfus.ct.mod2)
dp.mfus.ct.mod2.aov
# extract p-values from ANOVA table
pvals.mfus.ct2 <- dp.mfus.ct.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.ct2 <- p.adjust(pvals.mfus.ct2, method = "fdr")

# combine into a new table
results.mfus.ct2 <- cbind(dp.mfus.ct.mod2.aov, pvals_fdr.mfus.ct2)
results.mfus.ct2

effectsize::eta_squared(dp.mfus.ct.mod2.aov)

FFA_morphology.rev.na.pfus <- FFA_morphology.rev.na %>% subset(label == "pFus")
dp.pfus.sa.mod2 <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_morphology.rev.na.pfus)

# view summary
dp.pfus.sa.mod2.aov <- anova(dp.pfus.sa.mod2)
dp.pfus.sa.mod2.aov
effectsize::eta_squared(dp.pfus.sa.mod2.aov)

# extract p-values from ANOVA table
pvals.pfus.sa2 <- dp.pfus.sa.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.sa2 <- p.adjust(pvals.pfus.sa2, method = "fdr")

# combine into a new table
results.pfus.sa2 <- cbind(dp.pfus.sa.mod2.aov, pvals_fdr.pfus.sa2)
results.pfus.sa2

dp.pfus.ct.mod2 <- lme(cortical_thickness_mean ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                       data = FFA_morphology.rev.na.pfus)

# view summary
dp.pfus.ct.mod2.aov <- anova(dp.pfus.ct.mod2)
dp.pfus.ct.mod2.aov
effectsize::eta_squared(dp.pfus.ct.mod2.aov)

# extract p-values from ANOVA table
pvals.pfus.ct2 <- dp.pfus.ct.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.ct2 <- p.adjust(pvals.pfus.ct2, method = "fdr")

# combine into a new table
results.pfus.ct2 <- cbind(dp.pfus.ct.mod2.aov, pvals_fdr.pfus.ct2)
results.pfus.ct2

FFA_sulc_morphology.rev <- FFA_sulc_morphology %>% filter(sub != "45")
FFA_sulc_morphology.rev.na <- FFA_sulc_morphology.rev %>% drop_na()

FFA_sulc_morphology.rev.na.mfus <- FFA_sulc_morphology.rev.na %>% subset(label == "mFus")
dp.mfus.gc.mod2 <- lme(gyral_crown ~ group * hemi + Age + Sex + handedness + cortex_sa.x, random = (~1|sub/hemi), 
                      data = FFA_sulc_morphology.rev.na.mfus)

# view summary
dp.mfus.gc.mod2.aov <- anova(dp.mfus.gc.mod2)
dp.mfus.gc.mod2.aov

# extract p-values from ANOVA table
pvals.mfus.gc2 <- dp.mfus.gc.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.mfus.gc2 <- p.adjust(pvals.mfus.gc2, method = "fdr")

# combine into a new table
results.mfus.gc2 <- cbind(dp.mfus.gc.mod2.aov, pvals_fdr.mfus.gc2)
results.mfus.gc2

effectsize::eta_squared(dp.mfus.gc.mod2.aov)

FFA_sulc_morphology.rev.na.pfus <- FFA_sulc_morphology.rev.na %>% subset(label == "pFus")
dp.pfus.gc.mod2 <- lme(gyral_crown ~ group * hemi + Age + Sex + handedness + cortex_sa.x, random = (~1|sub/hemi), 
                      data = FFA_sulc_morphology.rev.na.pfus)

# view summary
dp.pfus.gc.mod2.aov <- anova(dp.pfus.gc.mod2)
dp.pfus.gc.mod2.aov

# extract p-values from ANOVA table
pvals.pfus.gc2 <- dp.pfus.gc.mod2.aov[["p-value"]]

# apply FDR correction
pvals_fdr.pfus.gc2 <- p.adjust(pvals.pfus.gc2, method = "fdr")

# combine into a new table
results.pfus.gc2 <- cbind(pvals_fdr.pfus.gc2, pvals_fdr.pfus.gc2)
results.pfus.gc2

effectsize::eta_squared(dp.pfus.gc.mod2.aov)


FFA_all_morph.rev_mfus <- FFA_sulc_morphology.rev.na %>% subset(label == "mFus")
FFA_all_morph.rev_pfus <- FFA_sulc_morphology.rev.na %>% subset(label == "pFus")

# mfus
FFA_all_morph.rev_mfus.lh <- FFA_all_morph.rev_mfus %>% subset(hemi == "left")
FFA_all_morph.rev_mfus.rh <- FFA_all_morph.rev_mfus %>% subset(hemi == "right")

cor.test(FFA_all_morph.rev_mfus$gyral_crown, FFA_all_morph.rev_mfus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_mfus$cortical_thickness_mean, FFA_all_morph.rev_mfus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_mfus.lh$gyral_crown, FFA_all_morph.rev_mfus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #
cor.test(FFA_all_morph.rev_mfus.rh$gyral_crown, FFA_all_morph.rev_mfus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_mfus.lh$cortical_thickness_mean, FFA_all_morph.rev_mfus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) 
cor.test(FFA_all_morph.rev_mfus.rh$cortical_thickness_mean, FFA_all_morph.rev_mfus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) 

# pfus
FFA_all_morph.rev_pfus.lh <- FFA_all_morph.rev_pfus %>% subset(hemi == "left")
FFA_all_morph.rev_pfus.rh <- FFA_all_morph.rev_pfus %>% subset(hemi == "right")

cor.test(FFA_all_morph.rev_pfus$gyral_crown, FFA_all_morph.rev_pfus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_pfus$cortical_thickness_mean, FFA_all_morph.rev_pfus$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_pfus.lh$gyral_crown, FFA_all_morph.rev_pfus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #
cor.test(FFA_all_morph.rev_pfus.rh$gyral_crown, FFA_all_morph.rev_pfus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) #

cor.test(FFA_all_morph.rev_pfus.lh$cortical_thickness_mean, FFA_all_morph.rev_pfus.lh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) 
cor.test(FFA_all_morph.rev_pfus.rh$cortical_thickness_mean, FFA_all_morph.rev_pfus.rh$total_surface_area_.mm.2., method = "spearman", exact = FALSE) 

# cfmt 
DP_sample2 <- mfs_DPs %>% subset(hemi == "lh") %>% dplyr::select(sub, CFMT, Old.New, group)
DP_sample2 <- DP_sample2[1:47,]
DP_sample2$sub <- as.factor(DP_sample2$sub)

FFA_all_morph.rev_cfmt <- merge(FFA_sulc_morphology.rev.na, DP_sample2, by = c("sub"))

FFA_all_morph.rev_cfmt_mFus <- FFA_all_morph.rev_cfmt %>% subset(label == "mFus")
FFA_all_morph.rev_cfmt_pFus <- FFA_all_morph.rev_cfmt %>% subset(label == "pFus")
FFA_all_morph.rev_cfmt_mFus.lh <- FFA_all_morph.rev_cfmt_mFus %>% subset(hemi == "left")
FFA_all_morph.rev_cfmt_mFus.rh <- FFA_all_morph.rev_cfmt_mFus %>% subset(hemi == "right")
FFA_all_morph.rev_cfmt_pFus.lh <- FFA_all_morph.rev_cfmt_pFus %>% subset(hemi == "left")
FFA_all_morph.rev_cfmt_pFus.rh <- FFA_all_morph.rev_cfmt_pFus %>% subset(hemi == "right")

cor.test(FFA_all_morph.rev_cfmt_mFus$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_mFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.lh$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_mFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.rh$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_mFus.rh$CFMT, method = "spearman", exact = F) #

cor.test(FFA_all_morph.rev_cfmt_pFus$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_pFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.lh$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_pFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.rh$total_surface_area_.mm.2., FFA_all_morph.rev_cfmt_pFus.rh$CFMT, method = "spearman", exact = F) #


cor.test(FFA_all_morph.rev_cfmt_mFus$gyral_crown, FFA_all_morph.rev_cfmt_mFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.lh$gyral_crown, FFA_all_morph.rev_cfmt_mFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.rh$gyral_crown, FFA_all_morph.rev_cfmt_mFus.rh$CFMT, method = "spearman", exact = F) #

cor.test(FFA_all_morph.rev_cfmt_pFus$gyral_crown, FFA_all_morph.rev_cfmt_pFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.lh$gyral_crown, FFA_all_morph.rev_cfmt_pFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.rh$gyral_crown, FFA_all_morph.rev_cfmt_pFus.rh$CFMT, method = "spearman", exact = F) #


cor.test(FFA_all_morph.rev_cfmt_mFus$cortical_thickness_mean, FFA_all_morph.rev_cfmt_mFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.lh$cortical_thickness_mean, FFA_all_morph.rev_cfmt_mFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_mFus.rh$cortical_thickness_mean, FFA_all_morph.rev_cfmt_mFus.rh$CFMT, method = "spearman", exact = F) #

cor.test(FFA_all_morph.rev_cfmt_pFus$cortical_thickness_mean, FFA_all_morph.rev_cfmt_pFus$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.lh$cortical_thickness_mean, FFA_all_morph.rev_cfmt_pFus.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_all_morph.rev_cfmt_pFus.rh$cortical_thickness_mean, FFA_all_morph.rev_cfmt_pFus.rh$CFMT, method = "spearman", exact = F) #


##### lat fus size analyses #####
FFA_revision_metrics.fg <- FFA_revision_metrics %>% subset(label == "G_oc-temp_lat-fusifor") 
FFA_revision_metrics.fg <- merge(FFA_revision_metrics.fg, FFA_DP , by = c("sub", "hemi")
)
FFA_revision_metrics.fg <- merge(FFA_revision_metrics.fg, FFA_demo , by = c("sub")
)
FFA_revision_metrics.fg <- merge(FFA_revision_metrics.fg, FFA_revision_metrics.cortex , by = c("sub", "hemi")
)
FFA_revision_metrics.fg$sub <- as.factor(FFA_revision_metrics.fg$sub)
FFA_revision_metrics.fg$group <- gsub("Controls", "NTs", FFA_revision_metrics.fg$group)
FFA_revision_metrics.fg$hemi <- gsub("lh", "left", FFA_revision_metrics.fg$hemi)
FFA_revision_metrics.fg$hemi <- gsub("rh", "right", FFA_revision_metrics.fg$hemi)


##### Lateral FG area ~ DP? #####
FFA_revision_metrics.fg.sa.mod <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                      data = FFA_revision_metrics.fg)

# view summary
FFA_revision_metrics.fg.sa.mod.aov <- anova(FFA_revision_metrics.fg.sa.mod)
FFA_revision_metrics.fg.sa.mod.aov

# extract p-values from ANOVA table
pvals.latFG <- FFA_revision_metrics.fg.sa.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.latFG <- p.adjust(pvals.latFG, method = "fdr")

# combine into a new table
results.latFG <- cbind(FFA_revision_metrics.fg.sa.mod.aov, pvals_fdr.latFG)
results.latFG

effectsize::eta_squared(FFA_revision_metrics.fg.sa.mod.aov)


FFA_revision_metrics.fg %>% group_by(Sex) %>% summarise(mean = mean(total_surface_area_.mm.2., na.rm = T), 
                                                          sd = sd(total_surface_area_.mm.2., na.rm = T), 
                                                          n=n(), 
                                                          se=sd/sqrt(n))


FFA_revision_metrics.fg.rev <- FFA_revision_metrics.fg %>% filter(sub != "45")
FFA_revision_metrics.fg.sa.mod <- lme(total_surface_area_.mm.2. ~ group * hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                                      data = FFA_revision_metrics.fg.rev)

# view summary
FFA_revision_metrics.fg.sa.mod.aov <- anova(FFA_revision_metrics.fg.sa.mod)
FFA_revision_metrics.fg.sa.mod.aov

# extract p-values from ANOVA table
pvals.latFG <- FFA_revision_metrics.fg.sa.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.latFG <- p.adjust(pvals.latFG, method = "fdr")

# combine into a new table
results.latFG <- cbind(FFA_revision_metrics.fg.sa.mod.aov, pvals_fdr.latFG)
results.latFG

effectsize::eta_squared(FFA_revision_metrics.fg.sa.mod.aov)

FFA_revision_metrics.fg.data <-  FFA_revision_metrics.fg %>%
  group_by(group, hemi) %>% summarise(mean = mean(total_surface_area_.mm.2., na.rm = T), 
                                      sd = sd(total_surface_area_.mm.2., na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(total_surface_area_.mm.2. = mean)


FFA_revision_metrics.fg.plot <- ggplot(FFA_revision_metrics.fg, aes(x = hemi, 
                                                                     y = total_surface_area_.mm.2., 
                                                                     fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_revision_metrics.fg.data,
           aes(x=hemi, y=total_surface_area_.mm.2., fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_revision_metrics.fg.data,
                aes(x=hemi, ymin=total_surface_area_.mm.2.-se, ymax=total_surface_area_.mm.2.+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  #geom_boxplot(outlier.alpha = 0, alpha = .5, width = .3, position = dodge) + 
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#bd0026", "#f03b20")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "surface area",
       fill = "group",
       title = "mFus") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_revision_metrics.fg.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.fg.plot.png",
                plot = FFA_revision_metrics.fg.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_revision_metrics.fg2 <- merge(FFA_revision_metrics.fg, DP_sample2, by = c("sub"))
FFA_revision_metrics.fg2.lh <- FFA_revision_metrics.fg2 %>% subset(hemi == "left")
FFA_revision_metrics.fg2.rh <- FFA_revision_metrics.fg2 %>% subset(hemi == "right")

cor.test(FFA_revision_metrics.fg2$total_surface_area_.mm.2., FFA_revision_metrics.fg2$CFMT, method = "spearman", exact = FALSE) #
cor.test(FFA_revision_metrics.fg2.lh$total_surface_area_.mm.2., FFA_revision_metrics.fg2.lh$CFMT, method = "spearman", exact = FALSE) #
cor.test(FFA_revision_metrics.fg2.rh$total_surface_area_.mm.2., FFA_revision_metrics.fg2.rh$CFMT, method = "spearman", exact = FALSE) #

FFA_revision_metrics.fg2 <- merge(FFA_revision_metrics.fg.rev, DP_sample2, by = c("sub"))
FFA_revision_metrics.fg2.lh <- FFA_revision_metrics.fg2 %>% subset(hemi == "left")
FFA_revision_metrics.fg2.rh <- FFA_revision_metrics.fg2 %>% subset(hemi == "right")

cor.test(FFA_revision_metrics.fg2$total_surface_area_.mm.2., FFA_revision_metrics.fg2$CFMT, method = "spearman", exact = FALSE) #
cor.test(FFA_revision_metrics.fg2.lh$total_surface_area_.mm.2., FFA_revision_metrics.fg2.lh$CFMT, method = "spearman", exact = FALSE) #
cor.test(FFA_revision_metrics.fg2.rh$total_surface_area_.mm.2., FFA_revision_metrics.fg2.rh$CFMT, method = "spearman", exact = FALSE) #

FFA_revision_metrics.fg2.cfmt.plot <- FFA_revision_metrics.fg2 %>% 

  ggplot(aes(x = total_surface_area_.mm.2., y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#bd0026", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "surface area",
       y = "CFMT",
       fill = "hemisphere",
       title = "mFus") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none")

FFA_revision_metrics.fg2.cfmt.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.fg2.cfmt.plot.png",
                plot = FFA_revision_metrics.fg2.cfmt.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

##### PPA #####
FFA_revision_metrics.PPA <- FFA_revision_metrics %>% subset(label == "PPA") 
FFA_revision_metrics.PPA <- merge(FFA_revision_metrics.PPA, FFA_DP , by = c("sub", "hemi")
)
FFA_revision_metrics.PPA <- merge(FFA_revision_metrics.PPA, FFA_demo , by = c("sub")
)
FFA_revision_metrics.PPA <- merge(FFA_revision_metrics.PPA, FFA_revision_metrics.cortex , by = c("sub", "hemi")
)
FFA_revision_metrics.PPA$sub <- as.factor(FFA_revision_metrics.PPA$sub)
FFA_revision_metrics.PPA$group <- gsub("Controls", "NTs", FFA_revision_metrics.PPA$group)
FFA_revision_metrics.PPA$hemi <- gsub("lh", "left", FFA_revision_metrics.PPA$hemi)
FFA_revision_metrics.PPA$hemi <- gsub("rh", "right", FFA_revision_metrics.PPA$hemi)

FFA_revision_metrics.PPA.gc.mod <- lme(gyral_crown ~ group * hemi + Age + handedness + Sex, random = (~1|sub/hemi), 
                                       data = FFA_revision_metrics.PPA)
# view summary
FFA_revision_metrics.PPA.gc.mod.aov <- anova(FFA_revision_metrics.PPA.gc.mod)
FFA_revision_metrics.PPA.gc.mod.aov

# extract p-values from ANOVA table
pvals.ppa.gc <- FFA_revision_metrics.PPA.gc.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.ppa.gc <- p.adjust(pvals.ppa.gc, method = "fdr")

# combine into a new table
results.ppa.gc <- cbind(FFA_revision_metrics.PPA.gc.mod.aov, pvals_fdr.ppa.gc)
results.ppa.gc

effectsize::eta_squared(FFA_revision_metrics.PPA.gc.mod.aov)


FFA_revision_metrics.PPA.data <-  FFA_revision_metrics.PPA %>%
  group_by(group, hemi) %>% summarise(mean = mean(gyral_crown, na.rm = T), 
                                      sd = sd(gyral_crown, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(gyral_crown = mean)

FFA_revision_metrics.PPA.gc.plot <- ggplot(FFA_revision_metrics.PPA, aes(x = hemi, 
                                                                    y = gyral_crown, 
                                                                    fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_revision_metrics.PPA.data,
           aes(x=hemi, y=gyral_crown, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_revision_metrics.PPA.data,
                aes(x=hemi, ymin=gyral_crown-se, ymax=gyral_crown+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  #geom_boxplot(outlier.alpha = 0, alpha = .5, width = .3, position = dodge) + 
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#fed98e", "#ffffd4")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "gyral height",
       fill = "group",
       title = "PPA") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_revision_metrics.PPA.gc.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.PPA.gc.plot.png",
                plot = FFA_revision_metrics.PPA.gc.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")



FFA_revision_metrics.PPA.ct.mod <- lme(cortical_thickness_mean ~ group + hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                                       data = FFA_revision_metrics.PPA)
# view summary
FFA_revision_metrics.PPA.ct.mod.aov <- anova(FFA_revision_metrics.PPA.ct.mod)
FFA_revision_metrics.PPA.ct.mod.aov

# extract p-values from ANOVA table
pvals.ppa.ct <- FFA_revision_metrics.PPA.ct.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.ppa.gc <- p.adjust(pvals.ppa.ct, method = "fdr")

# combine into a new table
results.pfus.gc2 <- cbind(FFA_revision_metrics.PPA.ct.mod.aov, pvals_fdr.ppa.gc)
results.pfus.gc2

effectsize::eta_squared(FFA_revision_metrics.PPA.ct.mod.aov)

FFA_revision_metrics.PPA.data.ct <-  FFA_revision_metrics.PPA %>%
  group_by(group, hemi) %>% summarise(mean = mean(cortical_thickness_mean, na.rm = T), 
                                      sd = sd(cortical_thickness_mean, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(cortical_thickness_mean = mean)

FFA_revision_metrics.PPA.ct.plot <- ggplot(FFA_revision_metrics.PPA, aes(x = hemi, 
                                                                         y = cortical_thickness_mean, 
                                                                         fill = group)) + 
  
  #geom_violin(alpha = .5) +
  geom_bar(data = FFA_revision_metrics.PPA.data.ct,
           aes(x=hemi, y=cortical_thickness_mean, fill = group), 
           stat="identity", color = "black", alpha = .5) +
  
  geom_errorbar(data = FFA_revision_metrics.PPA.data.ct,
                aes(x=hemi, ymin=cortical_thickness_mean-se, ymax=cortical_thickness_mean+se), 
                width=0, color="black", alpha = 1, linewidth = 1.5) +
  
  geom_beeswarm(aes(shape = group), cex = 4, alpha = .75) +
  
  #geom_boxplot(outlier.alpha = 0, alpha = .5, width = .3, position = dodge) + 
  
  scale_fill_manual(breaks = c("NTs", "DPs"),
                    values = c("#fed98e", "#ffffd4")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "hemisphere",
       y = "cortical thickness",
       fill = "group",
       title = "PPA") +
  
  # set theme customization
  theme_classic() +
  FFA_project_theme +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5, limits = c(2,4)) +
  
  facet_wrap(~ group) +
  
  guides(fill = "none", shape = "none") 


FFA_revision_metrics.PPA.ct.plot  

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.PPA.ct.plot.png",
                plot = FFA_revision_metrics.PPA.ct.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


FFA_revision_metrics.PPA.cfmt <- merge(FFA_revision_metrics.PPA, DP_sample2, by = c("sub"))
FFA_revision_metrics.PPA.cfmt.lh <- FFA_revision_metrics.PPA.cfmt %>% subset(hemi == "left")
FFA_revision_metrics.PPA.cfmt.rh <- FFA_revision_metrics.PPA.cfmt %>% subset(hemi == "right")

cor.test(FFA_revision_metrics.PPA.cfmt$gyral_crown, FFA_revision_metrics.PPA.cfmt$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.lh$gyral_crown, FFA_revision_metrics.PPA.cfmt.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.rh$gyral_crown, FFA_revision_metrics.PPA.cfmt.rh$CFMT, method = "spearman", exact = F) #

FFA_revision_metrics.PPA.cfmt.gc.plot <- FFA_revision_metrics.PPA.cfmt %>% 
  
  ggplot(aes(x = gyral_crown, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#fed98e", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "gyral height",
       y = "CFMT",
       fill = "hemisphere",
       title = "PPA") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none")

FFA_revision_metrics.PPA.cfmt.gc.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.PPA.cfmt.gc.plot.png",
                plot = FFA_revision_metrics.PPA.cfmt.gc.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")


cor.test(FFA_revision_metrics.PPA.cfmt$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.lh$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.rh$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt.rh$CFMT, method = "spearman", exact = F) #
ppa2 <- c(0.3724, 0.7968, 0.3699)
p.adjust(ppa2, method = "fdr")

FFA_revision_metrics.PPA.cfmt.ct.plot <- FFA_revision_metrics.PPA.cfmt %>% 
  
  ggplot(aes(x = cortical_thickness_mean, y = CFMT)) + 
  
  geom_jitter(aes(shape = group, fill = hemi)) + 
  
  geom_smooth(aes(color = hemi), method = "lm", 
              se = FALSE, linetype = "dashed", fullrange = T) +
  
  geom_smooth(color = "#fed98e", method = "lm", 
              se = FALSE, linetype = "solid", fullrange = T) +
  
  
  scale_fill_manual(breaks = c("left", "right"),
                    values = c("black", "gray")
  ) +
  
  scale_color_manual(breaks = c("left", "right"),
                     values = c("black", "gray")
  ) +
  
  scale_shape_manual(breaks = c("NTs", "DPs"),
                     values = c(21, 24)) +
  
  # set labels
  labs(x = "cortical thickness",
       y = "CFMT",
       fill = "hemisphere",
       title = "PPA") +
  
  theme_classic() +
  FFA_project_theme2 +
  
  # set axis scale
  scale_y_continuous(n.breaks = 5) +
  #limits = c(0,800)) +
  scale_x_continuous(n.breaks = 5) +
  guides(fill = "none", shape = "none", color = "none")

FFA_revision_metrics.PPA.cfmt.ct.plot

ggplot2::ggsave(filename = "~/Downloads/FFA_revision_metrics.PPA.cfmt.ct.plot.png",
                plot = FFA_revision_metrics.PPA.cfmt.ct.plot,
                device = "png",
                width = 3,
                height = 5,
                units = "in",
                dpi = "retina")

FFA_revision_metrics.PPA.rev <- FFA_revision_metrics.PPA %>% filter(sub != "45")

FFA_revision_metrics.PPA.gc.mod <- lme(gyral_crown ~ group * hemi + Age + handedness + Sex, random = (~1|sub/hemi), 
                                       data = FFA_revision_metrics.PPA.rev)
# view summary
FFA_revision_metrics.PPA.gc.mod.aov <- anova(FFA_revision_metrics.PPA.gc.mod)
FFA_revision_metrics.PPA.gc.mod.aov

# extract p-values from ANOVA table
pvals.ppa.gc <- FFA_revision_metrics.PPA.gc.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.ppa.gc <- p.adjust(pvals.ppa.gc, method = "fdr")

# combine into a new table
results.ppa.gc <- cbind(FFA_revision_metrics.PPA.gc.mod.aov, pvals_fdr.ppa.gc)
results.ppa.gc

effectsize::eta_squared(FFA_revision_metrics.PPA.gc.mod.aov)

FFA_revision_metrics.PPA.ct.mod <- lme(cortical_thickness_mean ~ group + hemi + Age + handedness + Sex + cortex_sa, random = (~1|sub/hemi), 
                                       data = FFA_revision_metrics.PPA.rev)
# view summary
FFA_revision_metrics.PPA.ct.mod.aov <- anova(FFA_revision_metrics.PPA.ct.mod)
FFA_revision_metrics.PPA.ct.mod.aov

# extract p-values from ANOVA table
pvals.ppa.ct <- FFA_revision_metrics.PPA.ct.mod.aov[["p-value"]]

# apply FDR correction
pvals_fdr.ppa.gc <- p.adjust(pvals.ppa.ct, method = "fdr")

# combine into a new table
results.pfus.gc2 <- cbind(FFA_revision_metrics.PPA.ct.mod.aov, pvals_fdr.ppa.gc)
results.pfus.gc2

effectsize::eta_squared(FFA_revision_metrics.PPA.ct.mod.aov)

FFA_revision_metrics.PPA.cfmt <- merge(FFA_revision_metrics.PPA.rev, DP_sample2, by = c("sub"))
FFA_revision_metrics.PPA.cfmt.lh <- FFA_revision_metrics.PPA.cfmt %>% subset(hemi == "left")
FFA_revision_metrics.PPA.cfmt.rh <- FFA_revision_metrics.PPA.cfmt %>% subset(hemi == "right")

cor.test(FFA_revision_metrics.PPA.cfmt$gyral_crown, FFA_revision_metrics.PPA.cfmt$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.lh$gyral_crown, FFA_revision_metrics.PPA.cfmt.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.rh$gyral_crown, FFA_revision_metrics.PPA.cfmt.rh$CFMT, method = "spearman", exact = F) #

cor.test(FFA_revision_metrics.PPA.cfmt$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.lh$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt.lh$CFMT, method = "spearman", exact = F) #
cor.test(FFA_revision_metrics.PPA.cfmt.rh$cortical_thickness_mean, FFA_revision_metrics.PPA.cfmt.rh$CFMT, method = "spearman", exact = F) #
