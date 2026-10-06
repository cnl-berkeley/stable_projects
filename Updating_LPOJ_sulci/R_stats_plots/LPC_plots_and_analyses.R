#### libraries ####
# Unmark and run code below if packages are not installed
# project_packages <- c("tidyverse", "afex", "car")
# install.packages(project_packages)

library(tidyverse)
library(afex)
library(car)

#### Plot code ####
set_theme <- function(){
  
  # Customize axis and legend
  set_theme = theme(
    # set axis title size
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    
    # set axis text size
    axis.text = element_text(color = "#333333", size = 14),
    
    # rotate x-axis 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    
    # set strip text
    strip.text = element_text(size = 14),
    
    # remove strip itself
    strip.background = element_blank(),
    
    #change background pannel
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # change legend size
    legend.title = element_text(color = "#333333", size = 12),
    legend.text = element_text(color = "#333333", size = 12)
    #panel.border = element_rect(colour = "black", 
    #                            fill=NA, size = 1),

  )
  # replace classic theme setting with our setting
  theme = theme_classic()%+replace% set_theme
  
  return(theme) 
}

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

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

#### Morphology dataframes ####
lpc_morph_final <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/HCP_lpc_anatomy_72_subs_final.csv")
lpc_morph_final$label <- gsub("SLOS4", "pAngs-d", lpc_morph_final$label)
lpc_morph_final$label <- gsub("SLOS3", "pAngs-v",lpc_morph_final$label)
lpc_morph_final$label <- gsub("SLOS2", "slocs-d", lpc_morph_final$label )
lpc_morph_final$label <- gsub("SLOS", "slocs-v", lpc_morph_final$label )
lpc_morph_final$label <- gsub("iTOS", "lTOS", lpc_morph_final$label)
lpc_morph_final$label <- gsub("STS", "1.STS", lpc_morph_final$label)
lpc_morph_final$label <- gsub("cSTS3", "2.cSTS3", lpc_morph_final$label)
lpc_morph_final$label <- gsub("cSTS2", "3.cSTS2", lpc_morph_final$label)
lpc_morph_final$label <- gsub("cSTS1", "4.cSTS1", lpc_morph_final$label)
lpc_morph_final$label <- gsub("SmgS", '5.SmgS', lpc_morph_final$label)
lpc_morph_final$label <- gsub("slocs-v", '6.slocs-v', lpc_morph_final$label)
lpc_morph_final$label <- gsub("slocs-d", '7.slocs-d', lpc_morph_final$label)
lpc_morph_final$label <- gsub("pAngs-v", '8.pAngs-v', lpc_morph_final$label)
lpc_morph_final$label <- gsub("pAngs-d", '9.pAngs-d', lpc_morph_final$label)
lpc_morph_final$label <- gsub("lTOS", '10.lTOS', lpc_morph_final$label)
lpc_morph_final$label <- gsub("mTOS", '11.mTOS', lpc_morph_final$label)
lpc_morph_final$label <- gsub("IPS-PO", '12.IPS-PO', lpc_morph_final$label)
lpc_morph_final$label <- gsub("pips", '13.pips', lpc_morph_final$label)
lpc_morph_final$label <- gsub("sB", '14.sB', lpc_morph_final$label)
lpc_morph_final$label <- gsub("aipsJ", '15.aipsJ', lpc_morph_final$label)
lpc_morph_final$label <- gsub("IPS", '16.IPS', lpc_morph_final$label)
lpc_morph_final$label <- gsub("SPS", '17.SPS', lpc_morph_final$label)
lpc_morph_final$label <- gsub("c1.STS1", '4.cSTS1', lpc_morph_final$label)
lpc_morph_final$label <- gsub("c1.STS2", '3.cSTS2', lpc_morph_final$label)
lpc_morph_final$label <- gsub("c1.STS3", '2.cSTS3', lpc_morph_final$label)
lpc_morph_final$label <- gsub("12.16.IPS-PO", '12.IPS-PO', lpc_morph_final$label)

lpc_morph_final$label <- factor(lpc_morph_final$label, 
                                                levels = c('1.STS',
                                                           '2.cSTS3',
                                                           '3.cSTS2', 
                                                           '4.cSTS1',
                                                           '5.SmgS',
                                                           '6.slocs-v',
                                                           '7.slocs-d',
                                                           '8.pAngs-v',
                                                           '9.pAngs-d',
                                                           '10.lTOS',
                                                           '11.mTOS',
                                                           '12.IPS-PO',
                                                           '13.pips',
                                                           '14.sB',
                                                           '15.aipsJ',
                                                           '16.IPS',
                                                           '17.SPS'
                                                           ))

unique(lpc_morph_final$label)


#### Incidence analysis ####
lpc_morph_final <- lpc_morph_final %>% mutate(pres_abs = case_when(
  vertices > 0 ~ "present",
  TRUE ~ "absent"))  

lpc_morph_final.chisq <- lpc_morph_final %>% mutate(pres_abs2 = case_when(pres_abs == "absent" ~ 0, pres_abs == "present" ~ 1))

lpc_morph_final.chisq.slocs.pangs <- lpc_morph_final.chisq %>% subset(label %in% c("8.pAngs-d", "7.pAngs-v", "6.slocs-d", "5.slocs-v"))

lpc_morph_final.chisq.slocs.pangs.glm <- glm(pres_abs2 ~ label*hemi, data=lpc_morph_final.chisq.slocs.pangs,
                                                 family=binomial(link = "logit"))
lpc_morph_final.chisq.slocs.pangs.glm.aov <- car::Anova(lpc_morph_final.chisq.slocs.pangs.glm, test="LR", type="III")
lpc_morph_final.chisq.slocs.pangs.glm.aov

lpc_morph_final.chisq.slocs.pangs.glm.m1 <- emmeans::emmeans(lpc_morph_final.chisq.slocs.pangs.glm, ~ label)
emmeans::contrast(lpc_morph_final.chisq.slocs.pangs.glm.m1, method='pairwise')

lpc_morph_final.chisq.slocs.pangs %>% 
  group_by(label, hemi, pres_abs) %>% 
  subset(pres_abs == "present") %>% 
  summarise(n = n(), pct = n/72*100)


#### Morphology and architecture analyses #### 

# Normalization
lobe_morphology <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lobe_morphology.csv")
lobe_morphology.human <- lobe_morphology %>% subset(species == "human")
lpc_morph_final_analysis <- merge(lpc_morph_final_analysis, lobe_morphology.human, by = c("sub", "hemi"))
lpc_morph_final_analysis <- lpc_morph_final_analysis %>% mutate(normalized_SA = total_surface_area_.mm.2./cortex_sa)

lpc_morph_final_analysis$sub <- as.factor(lpc_morph_final_analysis$sub)

lpc_morph_final_analysis.plot <- lpc_morph_final_analysis %>% subset(label %in% c(
                                           '2.cSTS3',
                                           '3.cSTS2', 
                                           '4.cSTS1',
                                           '5.SmgS',
                                           '6.slocs-v',
                                           '7.slocs-d',
                                           '8.pAngs-v',
                                           '9.pAngs-d',
                                           '10.lTOS',
                                           '11.mTOS',
                                           '12.IPS-PO',
                                           '13.pips',
                                           '14.sB',
                                           '15.aipsJ',
                                           '16.IPS',
                                           '17.SPS'))
                               

##### Supplemental Figure 7 plots #####

lpc_morph_final_analysis.sd.data <-  lpc_morph_final_analysis.plot %>% drop_na() %>%
  group_by(label, hemi) %>% summarise(mean = mean(sulcal_depth_mean_pct, na.rm = T), 
                                      sd = sd(sulcal_depth_mean_pct, na.rm = T), 
                                      n=n(), 
                                      se=sd/sqrt(n)) %>% 
  rename(sulcal_depth_mean_pct = mean)


# Depth
lpc_morph_final_sd <-  ggplot(lpc_morph_final_analysis.plot, 
                              aes(x = label, 
                                  y = sulcal_depth_mean_pct, 
                                  fill = hemi)) +
  
  geom_point(shape = 21, color = 'black', alpha = .25, 
             position = position_jitterdodge(jitter.width = .2, dodge.width = .5, 
                                             jitter.height = 0)) +
  
  geom_boxplot(alpha = 0.8, outlier.alpha = 0, width = 0.5) +
  
  set_theme() +
  scale_fill_manual(breaks = c("lh", "rh"),
                      values = c("darkgray", "white"))

lpc_morph_final_sd

# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_sulcdepth2.png",
#                   plot = lpc_morph_final_sd,
#                   device = "png",
#                   width = 10,
#                   height = 4, 
#                   units = "in",
#                   dpi = "retina")


# Surface area
lpc_morph_final_sa <- lpc_morph_final_analysis.plot %>% 

  ggplot(aes(x = label, y = normalized_SA, fill = hemi)) +
  
  geom_point(shape = 21, color = 'black', alpha = .25, 
             position = position_jitterdodge(jitter.width = .2, dodge.width = .5, 
                                             jitter.height = 0)) +
  
  geom_boxplot(alpha = 0.8, outlier.alpha = 0, width = 0.5) +
  
  set_theme() +
  scale_fill_manual(breaks = c("lh", "rh"),
                    values = c("darkgray", "white"))
  
lpc_morph_final_sa
  
# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_sa2.png",
#                   plot = lpc_morph_final_sa,
#                   device = "png",
#                   width = 10,
#                   height = 4, 
#                   units = "in",
#                   dpi = "retina")


# Cortical thickness
lpc_morph_final_ct <- lpc_morph_final_analysis.plot %>% 
  
  ggplot(aes(x = label, y = cortical_thickness_mean, fill = hemi)) +
  
  geom_point(shape = 21, color = 'black', alpha = .25, 
             position = position_jitterdodge(jitter.width = .2, dodge.width = .5, 
                                             jitter.height = 0)) +
  
  geom_boxplot(alpha = 0.8, outlier.alpha = 0, width = 0.5) +
  
  set_theme() +
  scale_fill_manual(breaks = c("lh", "rh"),
                    values = c("darkgray", "white"))

lpc_morph_final_ct

# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_ct.png",
#                 plot = lpc_morph_final_ct,
#                 device = "png",
#                 width = 10,
#                 height = 4, 
#                 units = "in",
#                 dpi = "retina")


# Myelination
lpc_morph_final_my <- lpc_morph_final_analysis.plot %>% 
  
  ggplot(aes(x = label, y = myelin, fill = hemi)) +
  
  geom_point(shape = 21, color = 'black', alpha = .25, 
             position = position_jitterdodge(jitter.width = .2, dodge.width = .5, 
                                             jitter.height = 0)) +
  
  geom_boxplot(alpha = 0.8, outlier.alpha = 0, width = 0.5) +
  
  set_theme() +
  scale_fill_manual(breaks = c("lh", "rh"),
                    values = c("darkgray", "white"))

lpc_morph_final_my

# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_my.png",
#                 plot = lpc_morph_final_my,
#                 device = "png",
#                 width = 10,
#                 height = 4, 
#                 units = "in",
#                 dpi = "retina")


#### Morphology and architecture statistics ####
lpc_morph_final_depth_ratio_analysis <- lpc_morph_final_analysis %>% 
  
  select(sub, hemi, label, sulcal_depth_mean_pct, normalized_SA, cortical_thickness_mean, myelin) %>%
  
  mutate(scaled_sulcal_depth_mean_pct = scale(sulcal_depth_mean_pct), 
         scaled_normalized_SA = scale(normalized_SA),
         scaled_cortical_thickness_mean = scale(cortical_thickness_mean),
         scaled_myelin = scale(myelin)
  ) %>% 

  select(sub, hemi, label, scaled_sulcal_depth_mean_pct, scaled_normalized_SA, scaled_cortical_thickness_mean, scaled_myelin) %>%
  
  pivot_longer(cols = c("scaled_sulcal_depth_mean_pct", "scaled_normalized_SA",
                        "scaled_cortical_thickness_mean", "scaled_myelin"),
               names_to = "metric",
               values_to = "value") %>%
  
  mutate(metric = recode(metric, "scaled_sulcal_depth_mean_pct" = "Sulcal Depth"),
         metric = recode(metric, "scaled_normalized_SA" = "Surface Area"),
         metric = recode(metric, "scaled_cortical_thickness_mean" = "Cortical Thickness"),
         metric = recode(metric, "scaled_myelin" = "Myelination")
  )


###### Main analysis: slocs-v vs cSTS3 vs lTOS ######
lpc_morph_final_depth_ratio.main_analysis <- lpc_morph_final_depth_ratio_analysis %>%
  
  subset(label %in% c("2.cSTS3", "6.slocs-v", "10.lTOS"))

model_depth_ratio.main2 <- aov_ez('sub', dv = 'value', data = lpc_morph_final_depth_ratio.main_analysis, within = c('hemi','label','metric'),
                                      anova_table = list(es = "pes"))
model_depth_ratio.main2

model_depth_ratio.main.ixn <- emmeans::emmeans(model_depth_ratio.main2, ~ label | metric)
emmeans::contrast(model_depth_ratio.main.ixn, method='pairwise')

model_depth_ratio.main.ixn3 <- emmeans::emmeans(model_depth_ratio.main2, ~ hemi | label | metric)
emmeans::contrast(model_depth_ratio.main.ixn3, method='pairwise')


####### Figure 2a ####### 
lpc_morph_final_depth_ratio.main_analysis.mean_sd <- lpc_morph_final_depth_ratio.main_analysis %>% 
  group_by(label, metric, hemi) %>%
  summarise(n = n(), 
            mean = mean(value), 
            sd = sd(value), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(value)/sqrt(length(value)),
            upper = mean + 1.96*sd(value)/sqrt(length(value)))

lpc_morph_final_depth_ratio.main_analysis.mean_sd$metric <- factor(lpc_morph_final_depth_ratio.main_analysis.mean_sd$metric, 
                                                                    levels = c("Surface Area", "Cortical Thickness", "Myelination", "Sulcal Depth"))

lpc_morph_final_depth_ratio.main_analysis.mean_sd.plot <- lpc_morph_final_depth_ratio.main_analysis.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = metric, y = mean+se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = metric, y = mean-se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = metric, y = mean, 
                   group = label, 
                   color = label, fill = label),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = metric, y = mean, 
                 fill = label), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(-1.50, 1.50), n.breaks = 5) +
  
  labs(x = "Metric",
       y = "Standardized Units",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("2.cSTS3", "6.slocs-v", "10.lTOS"),
                    values = c("#0000ff", "#969696", "#00ff00")) +
  
  scale_color_manual(breaks = c("2.cSTS3", "6.slocs-v", "10.lTOS"),
                     values = c("#0000ff", "#969696", "#00ff00")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

lpc_morph_final_depth_ratio.main_analysis.mean_sd.plot


# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_depth_ratio.main_analysis.mean_sd.plot.png",
#                 plot = lpc_morph_final_depth_ratio.main_analysis.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


##### Supplemental analysis: cSTS ###### 
lpc_morph_final_depth_ratio_analysis_cSTS <- lpc_morph_final_depth_ratio_analysis %>%
  
  subset(label %in% c("4.cSTS1", "3.cSTS2", "2.cSTS3"))

model_depth_ratio.supp.cSTS <- aov_ez('sub', dv = 'value', 
                                      data = lpc_morph_final_depth_ratio_analysis_cSTS, 
                                      within = c('hemi','label','metric'),
                                  anova_table = list(es = "pes"))
model_depth_ratio.supp.cSTS

model_depth_ratio_csts.ixn <- emmeans::emmeans(model_depth_ratio.supp.cSTS, ~ label | metric)
emmeans::contrast(model_depth_ratio_csts.ixn, method='pairwise')

model_depth_ratio_csts.ixn2 <- emmeans::emmeans(model_depth_ratio.supp.cSTS, ~ hemi | label | metric)
emmeans::contrast(model_depth_ratio_csts.ixn2, method='pairwise')


##### Supplemental Figure 8a ##### 
lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd <- lpc_morph_final_depth_ratio_analysis_cSTS %>% 
  group_by(label, metric, hemi) %>%
  summarise(n = n(), 
            mean = mean(value), 
            sd = sd(value), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(value)/sqrt(length(value)),
            upper = mean + 1.96*sd(value)/sqrt(length(value)))

lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd$metric <- factor(lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd$metric, 
                                                                   levels = c("Surface Area", "Cortical Thickness", "Myelination", "Sulcal Depth"))


lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd.plot <- lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = metric, y = mean+se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = metric, y = mean-se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = metric, y = mean, 
                   group = label, 
                   color = label, fill = label),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = metric, y = mean, 
                 fill = label), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(-1.50, 1.50), n.breaks = 5) +
  
  labs(x = "Metric",
       y = "Standardized Units",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("2.cSTS3", "3.cSTS2", "4.cSTS1"),
                    values = c("#0000ff", "#0066ff", "#3399ff")) +
  
  scale_color_manual(breaks = c("2.cSTS3", "3.cSTS2", "4.cSTS1"),
                     values = c("#0000ff", "#0066ff", "#3399ff")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd.plot

# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd.plot.png",
#                 plot = lpc_morph_final_depth_ratio_analysis_cSTS.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


###### Supplemental analysis: pips vs aipsJ ###### 
lpc_morph_final_depth_ratio_analysis_ips <- lpc_morph_final_depth_ratio_analysis %>%
  
  subset(label %in% c("13.pips","15.aipsJ"))

model_depth_ratio.supp.ips <- aov_ez('sub', dv = 'value', data = lpc_morph_final_depth_ratio_analysis_ips, within = c('hemi','label','metric'),
                                      anova_table = list(es = "pes"))
model_depth_ratio.supp.ips

model_depth_ratio.supp.ips.ixn <- emmeans::emmeans(model_depth_ratio.supp.ips, ~ label | metric)
emmeans::contrast(model_depth_ratio.supp.ips.ixn, method='pairwise')


###### Supplemental Figure 8c ######
lpc_morph_final_depth_ratio_analysis_ips.mean_sd <- lpc_morph_final_depth_ratio_analysis_ips %>% 
  group_by(label, metric, hemi) %>%
  summarise(n = n(), 
            mean = mean(value), 
            sd = sd(value), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(value)/sqrt(length(value)),
            upper = mean + 1.96*sd(value)/sqrt(length(value)))

lpc_morph_final_depth_ratio_analysis_ips.mean_sd$metric <- factor(lpc_morph_final_depth_ratio_analysis_ips.mean_sd$metric, 
                                                                   levels = c("Surface Area", "Cortical Thickness", "Myelination", "Sulcal Depth"))


lpc_morph_final_depth_ratio_analysis_ips.mean_sd.plot <- lpc_morph_final_depth_ratio_analysis_ips.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = metric, y = mean+se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = metric, y = mean-se, 
                   group = label, 
                   color = label), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = metric, y = mean, 
                   group = label, 
                   color = label, fill = label),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = metric, y = mean, 
                 fill = label), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(-1.50, 1.50), n.breaks = 5) +
  
  labs(x = "Metric",
       y = "Standardized Units",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("13.pips", "15.aipsJ"),
                    values = c("#ff99ff", "#ff0000")) +
  
  scale_color_manual(breaks = c("13.pips", "15.aipsJ"),
                     values = c("#ff99ff", "#ff0000")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

lpc_morph_final_depth_ratio_analysis_ips.mean_sd.plot

# ggplot2::ggsave(filename = "~/Downloads/lpc_morph_final_depth_ratio_analysis_ips.mean_sd.plot.png",
#                 plot = lpc_morph_final_depth_ratio_analysis_ips.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


#### Functional fingerprint analyses ####
slocs_v_dice_lh <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_v_dice_kong_69_subs.csv")
slocs_v_dice_rh <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_v_dice_kong_69_subs.csv")
slocs_v_dice <- rbind(slocs_v_dice_lh, slocs_v_dice_rh)
slocs_v_dice$label_name <- gsub("SLOS", 'slocs-v', slocs_v_dice$label_name)

slocs_d_dice_lh <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_d_dice_kong_48_subs.csv")
slocs_d_dice_rh <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_d_dice_kong_46_subs.csv")
slocs_d_dice <- rbind(slocs_d_dice_lh, slocs_d_dice_rh)
slocs_d_dice$label_name <- gsub("SLOS2", 'slocs-d', slocs_d_dice$label_name)

slocs_dice <- rbind(slocs_v_dice, slocs_d_dice)
slocs_dice$sub <- as.factor(slocs_dice$sub)


##### Main analysis: slocs-v, cSTS3, lTOS #####
Consistent_LPC_sulci_dice <-  read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/consistent_lpc_sulci_dice_kong_70_subs.csv")
Consistent_LPC_sulci_dice$sub <- as.factor(Consistent_LPC_sulci_dice$sub)

all_LPC_func <- rbind(slocs_dice, Consistent_LPC_sulci_dice)
all_LPC_func_analysis <- all_LPC_func %>% subset(label_name %in% c("slocs-v", "cSTS3", "iTOS"))

main_lpc_func_analysis.aov2 <- aov_ez('sub', dv = 'dice_coeff', 
                                      data = all_LPC_func_analysis, 
                                      within = c('hemi','label_name','network_name'),
                                      anova_table = list(es = "pes"))
main_lpc_func_analysis.aov2

main_lpc_func_analysis.aov2.ixn <- emmeans::emmeans(main_lpc_func_analysis.aov2, ~ label_name | network_name)
emmeans::contrast(main_lpc_func_analysis.aov2.ixn, method='pairwise')

main_lpc_func_analysis.aov2.ixn2 <- emmeans::emmeans(main_lpc_func_analysis.aov2, ~ hemi | label_name | network_name)
emmeans::contrast(main_lpc_func_analysis.aov2.ixn2, method='pairwise')


##### Figure 2b #####
all_LPC_func_analysis.mean_sd <- all_LPC_func_analysis %>% 
  group_by(label_name, network_name, hemi) %>%
  summarise(n = n(), 
            mean = mean(dice_coeff), 
            sd = sd(dice_coeff), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)),
            upper = mean + 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)))

all_LPC_func_analysis.mean_sd.plot <- all_LPC_func_analysis.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = network_name, y = mean+se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = network_name, y = mean-se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = network_name, y = mean, 
                   group = label_name, 
                   color = label_name, fill = label_name),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = network_name, y = mean, 
                 fill = label_name), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(0, 1), n.breaks = 5) +
  
  labs(x = "Network Name",
       y = "Dice Coeff",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("cSTS3", "slocs-v", "iTOS"),
                    values = c("#0000ff", "#969696", "#00ff00")) +
  
  scale_color_manual(breaks = c("cSTS3", "slocs-v", "iTOS"),
                     values = c("#0000ff", "#969696", "#00ff00")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

all_LPC_func_analysis.mean_sd.plot

# ggplot2::ggsave(filename = "~/Downloads/all_LPC_func_analysis.mean_sd.plot.png",
#                 plot = all_LPC_func_analysis.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


##### Supplemental analysis: compare cSTS components #####
Consistent_LPC_sulci_dice.cSTS <- all_LPC_func %>% 
  subset(label_name %in% c("cSTS1", "cSTS2", "cSTS3"))

cSTS_func_analysis.aov <- aov_ez('sub', dv = 'dice_coeff', 
                                 data = Consistent_LPC_sulci_dice.cSTS, 
                                 within = c('hemi','label_name','network_name'),
                                      anova_table = list(es = "pes"))
cSTS_func_analysis.aov

cSTS_func_analysis.aov.ixn <- emmeans::emmeans(cSTS_func_analysis.aov, ~ label_name | network_name)
emmeans::contrast(cSTS_func_analysis.aov.ixn, method='pairwise')

cSTS_func_analysis.aov.ixn2 <- emmeans::emmeans(cSTS_func_analysis.aov, ~ hemi | label_name | network_name)
emmeans::contrast(cSTS_func_analysis.aov.ixn2, method='pairwise')


##### Supplemental Figure 8b #####
Consistent_LPC_sulci_dice.cSTS.mean_sd <- Consistent_LPC_sulci_dice.cSTS %>% 
  group_by(label_name, network_name, hemi) %>%
  summarise(n = n(), 
            mean = mean(dice_coeff), 
            sd = sd(dice_coeff), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)),
            upper = mean + 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)))

Consistent_LPC_sulci_dice.cSTS.mean_sd.plot <- Consistent_LPC_sulci_dice.cSTS.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = network_name, y = mean+se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = network_name, y = mean-se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = network_name, y = mean, 
                   group = label_name, 
                   color = label_name, fill = label_name),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = network_name, y = mean, 
                 fill = label_name), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(0, 1), n.breaks = 5) +
  
  labs(x = "Network Name",
       y = "Dice Coeff",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("cSTS1", "cSTS2", "cSTS3"),
                    values = c("#3399ff", "#0066ff", "#0000ff")) +
  
  scale_color_manual(breaks = c("cSTS1", "cSTS2", "cSTS3"),
                     values = c("#3399ff", "#0066ff", "#0000ff")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

Consistent_LPC_sulci_dice.cSTS.mean_sd.plot

# ggplot2::ggsave(filename = "~/Downloads/Consistent_LPC_sulci_dice.cSTS.mean_sd.plot.png",
#                 plot = Consistent_LPC_sulci_dice.cSTS.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


##### Supplemental analysis: pips vs aipsJ #####
Consistent_LPC_sulci_dice.ips <- all_LPC_func %>% 
  subset(label_name %in% c("pips", "aipsJ"))

ips_func_analysis.aov <- aov_ez('sub', dv = 'dice_coeff', 
                                data = Consistent_LPC_sulci_dice.ips, 
                                within = c('hemi','label_name','network_name'),
                                anova_table = list(es = "pes"))
ips_func_analysis.aov

ips_func_analysis.aov.ixn <- emmeans::emmeans(ips_func_analysis.aov, ~ label_name | network_name)
emmeans::contrast(ips_func_analysis.aov.ixn, method='pairwise')

ips_func_analysis.aov.ixn2 <- emmeans::emmeans(ips_func_analysis.aov, ~ hemi | label_name | network_name)
emmeans::contrast(ips_func_analysis.aov.ixn2, method='pairwise')


##### Supplemental Figure 8d #####
Consistent_LPC_sulci_dice.ips.mean_sd <- Consistent_LPC_sulci_dice.ips %>% 
  group_by(label_name, network_name, hemi) %>%
  summarise(n = n(), 
            mean = mean(dice_coeff), 
            sd = sd(dice_coeff), 
            se=sd/sqrt(n),
            lower = mean - 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)),
            upper = mean + 1.96*sd(dice_coeff)/sqrt(length(dice_coeff)))

Consistent_LPC_sulci_dice.ips.mean_sd.plot <- Consistent_LPC_sulci_dice.ips.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = network_name, y = mean+se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") + 
  geom_polygon(aes(x = network_name, y = mean-se, 
                   group = label_name, 
                   color = label_name), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = network_name, y = mean, 
                   group = label_name, 
                   color = label_name, fill = label_name),
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = network_name, y = mean, 
                 fill = label_name), 
             size = 3, 
             shape = 21) + 
  
  coord_radar() +
  
  scale_y_continuous(limits = c(0, 1), n.breaks = 5) +
  
  labs(x = "Network Name",
       y = "Dice Coeff",
       
       color = "Sulcus",
       fill = "Sulcus") +  
  
  scale_fill_manual(breaks = c("pips", "aipsJ"),
                    values = c("#ff99ff", "#ff0000")) +
  
  scale_color_manual(breaks = c("pips", "aipsJ"),
                     values = c("#ff99ff", "#ff0000")) +
  
  guides(color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) +
  
  facet_wrap(~ hemi)

Consistent_LPC_sulci_dice.ips.mean_sd.plot

# ggplot2::ggsave(filename = "~/Downloads/Consistent_LPC_sulci_dice.ips.mean_sd.plot.png",
#                 plot = Consistent_LPC_sulci_dice.ips.mean_sd.plot,
#                 device = "png",
#                 width = 10,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


#### SLOCS-V and SLOCS-D HCPMMP and Amunts ####

##### Glasser et al., 2016 - HCPMMP ##### 
lh_slocs_v_HCPMMP <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_v_HCPMMP_ROIs_updated_71_subs.csv")
rh_slocs_v_HCPMMP <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_v_HCPMMP_ROIs_updated_71_subs.csv")

lh_slocs_d_HCPMMP <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_d_HCPMMP_ROIs_updated_50_subs.csv")
rh_slocs_d_HCPMMP <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_d_HCPMMP_ROIs_updated_48_subs.csv")

slocs_HCPMMP <- rbind(lh_slocs_v_HCPMMP, rh_slocs_v_HCPMMP, lh_slocs_d_HCPMMP, rh_slocs_d_HCPMMP)
unique(slocs_HCPMMP$network_name)

slocs_HCPMMP$label_name <- gsub("SLOS_fsavg", "slocs-v", slocs_HCPMMP$label_name)
slocs_HCPMMP$label_name <- gsub("SLOS2_fsavg", "slocs-d", slocs_HCPMMP$label_name)
slocs_HCPMMP$hemi <- gsub("lh", "Left", slocs_HCPMMP$hemi)
slocs_HCPMMP$hemi <- gsub("rh", "Right", slocs_HCPMMP$hemi)

slocs_HCPMMP$network_name <- gsub("L_PGp_ROI", "PGp-HCPMMP", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_PGp_ROI", "PGp-HCPMMP", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_LO3_ROI", "LO3", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_LO3_ROI", "LO3", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_TPOJ1_ROI", "TPOJ1", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_TPOJ1_ROI", "TPOJ1", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_TPOJ2_ROI", "TPOJ2", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_TPOJ2_ROI", "TPOJ2", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_TPOJ3_ROI", "TPOJ3", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_TPOJ3_ROI", "TPOJ3", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_PGi_ROI", "PGi", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_PGi_ROI", "PGi", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_PGs_ROI", "PGs-HCPMMP", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_PGs_ROI", "PGs-HCPMMP", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_STV_ROI", "STV", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_STV_ROI", "STV", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_PFm_ROI", "PFm", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_PFm_ROI", "PFm", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("L_IP0_ROI", "IP0", slocs_HCPMMP$network_name)
slocs_HCPMMP$network_name <- gsub("R_IP0_ROI", "IP0", slocs_HCPMMP$network_name)

##### Amunts et al. 2020 - JulichBrain #####
lh_slocs_v_amunts <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_v_Amunts_ROIs_updated_71_subs.csv")
rh_slocs_v_amunts <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_v_Amunts_ROIs_updated_71_subs.csv")

lh_slocs_d_amunts <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/lh_slocs_d_Amunts_ROIs_updated_50_subs.csv")
rh_slocs_d_amunts <- read.csv("~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/rh_slocs_d_Amunts_ROIs_updated_48_subs.csv")

slocs_amunts <- rbind(lh_slocs_v_amunts, rh_slocs_v_amunts, lh_slocs_d_amunts, rh_slocs_d_amunts)
unique(slocs_amunts$network_name)

slocs_amunts$label_name <- gsub("SLOS_fsavg", "slocs-v", slocs_amunts$label_name)
slocs_amunts$label_name <- gsub("SLOS2_fsavg", "slocs-d", slocs_amunts$label_name)
slocs_amunts$hemi <- gsub("lh", "Left", slocs_amunts$hemi)
slocs_amunts$hemi <- gsub("rh", "Right", slocs_amunts$hemi)

slocs_amunts$network_name <- gsub("Area_hOc4lp_LOC", "hOc4lp-Amunts", slocs_amunts$network_name)
slocs_amunts$network_name <- gsub("Area_hIP4_IPS", "hIP4-Amunts", slocs_amunts$network_name)
slocs_amunts$network_name <- gsub("Area_hIP5_IPS", "hIP5-Amunts", slocs_amunts$network_name)
slocs_amunts$network_name <- gsub("Area_PGp_IPL", "PGp-Amunts", slocs_amunts$network_name)

slocs_HCPMMP_amunts <- rbind(slocs_HCPMMP, slocs_amunts)
slocs_HCPMMP_amunts2 <- slocs_HCPMMP_amunts %>% subset(label_name == "slocs-v" & 
                                                       network_name %in% c("PGs-HCPMMP", "PGp-HCPMMP", "hIP4-Amunts", "PGp-Amunts"))


##### Figure 4 #####
slocs_HCPMMP_amunts2.mean_sd <- slocs_HCPMMP_amunts2 %>% 
  group_by(network_name, label_name, hemi) %>% 
  summarise(n = n(), mean = mean(dice_coeff), sd = sd(dice_coeff), se=sd/sqrt(n))

slocs_HCPMMP_amunts2.mean_sd.radar.plot <- slocs_HCPMMP_amunts2.mean_sd %>%
  
  ggplot() +  
  
  geom_polygon(aes(x = network_name, y = mean+se, 
                   group = hemi, 
                   color = hemi), 
               fill = NA,
               linetype = "dashed") + 
  
  geom_polygon(aes(x = network_name, y = mean-se, 
                   group = hemi, 
                   color = hemi), 
               fill = NA,
               linetype = "dashed") +
  
  geom_polygon(aes(x = network_name, y = mean, 
                   group = hemi, 
                   color = hemi, fill = hemi), 
               fill = NA,
               size = 1,
               alpha = .2) +
  
  geom_point(aes(x = network_name, y = mean, 
                 fill = hemi), 
             size = 3, 
             shape = 21) + 
  
  coord_radar(start = pi/2) +
  
  scale_y_continuous(limits = c(0, 1), n.breaks = 5) +
  
  labs(x = "Metric",
       y = "Dice Coefficient",
       
       color = "Hemisphere",
       fill = "Hemisphere") +  
  
  scale_fill_manual(breaks = c("Left", "Right"),
                    values = c("#404040", "white")) +
  
  scale_color_manual(breaks = c("Left", "Right"),
                     values = c("#404040", "#969696")) +
  
  guides( 
    color = "none") +
  
  theme_minimal() +
  
  theme( 
    # set plot title
    plot.title = element_text(family = "Arial", color = "black", 
                              size = 0, hjust = 0.5),
    
    # set axis title
    axis.title = element_text(family = "Arial", color = "black", 
                              face = "italic", size = 14), 
    
    # set axis text
    axis.text.y = element_text(family = "Arial", color = "black",  size = 12),
    axis.text.x = element_text(family = "Arial", color = "black",  size = 12),
    
    # set legend
    legend.title = element_text(family = "Arial", color = "black", 
                                size = 14),
    legend.text = element_text(family = "Arial", color = "black", 
                               size = 12),
    
    # facet
    strip.text = element_text(family = "Arial", color = "black", 
                              face = "italic", 
                              size = 12, hjust = 0.5),
    
    panel.grid.minor = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.major = element_line(linewidth = 0.5, color = "gray")) 

slocs_HCPMMP_amunts2.mean_sd.radar.plot

# ggplot2::ggsave(filename = "~/Downloads/slocs_HCPMMP_amunts2.mean_sd.radar.plot.png",
#                 plot = slocs_HCPMMP_amunts2.mean_sd.radar.plot,
#                 device = "png",
#                 width = 5,
#                 height = 5,
#                 units = "in",
#                 dpi = "retina") 


#### LASSO plots ####

##### Coefficient plot #####
lasso_coefs_lpc <- read.csv('~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/LPC_coeffs_lasso.csv')
lasso_coefs_lpc$Sulci <- factor(lasso_coefs_lpc$Sulci, 
                                levels = c('1.STS',
                                           '2.cSTS3',
                                           '3.cSTS2', 
                                           '4.cSTS1',
                                           '5.SmgS',
                                           '6.slocs-v',
                                           '10.lTOS',
                                           '11.mTOS',
                                           '12.IPS-PO',
                                           '13.pips',
                                           '14.sB',
                                           '15.aipsJ',
                                           '16.IPS',
                                           '17.SPS'))

lasso_coefs_lpc_plt <- ggplot(data = lasso_coefs_lpc, aes(x=as.factor(alpha_level), y=as.factor(Sulci))) +
  # tile plot
  geom_tile(aes(fill = beta), linewidth=0.3) +
  
  # set coordinates
  coord_equal() +
  scale_y_discrete(limits = rev(levels(lasso_coefs_lpc$sulcus))) +
  
  # set color gradient
  scale_fill_gradient2(low = scales::muted("#1b9e77"), mid = "white", high = scales::muted("#d95f02"),
                       midpoint = 0, space = "Lab", na.value = "grey50",
                       guide = "colourbar", aesthetics = "fill") +
  
  
  # set text object and axis labels
  geom_text(aes(label=round(beta, digits = 2), 
                size = 2), size = 4)+
  
  xlab("alpha level") +
  ylab('lateral parietal sulci') +
  labs(fill = "beta") +
  set_theme()
# set theme
# theme(
#   # set axis title
#   axis.title.x = element_text(color="#333333",
#                               size=12), 
#   axis.title.y = element_text(color="#333333", 
#                               angle=90,size=14),
#   # set axis text
#   axis.text = element_text( color='#333333',
#                             size = 12),
#   # set pannel
#   
#   #panel.border = element_rect(colour = "black", 
#   #                fill=NA, size=1),
#   # set legend
#   legend.title = element_text(color="#333333", size = 12),
#   legend.text = element_text( color = '#333333', size = 12),
#   axis.line = element_blank(),
#   axis.ticks = element_blank(),
#   # set aspect ratio
#   aspect.ratio = 1/1
# )
lasso_coefs_lpc_plt

# ggplot2::ggsave(filename = "~/Downloads/lasso_coefs_lpc_plt.png",
#                 plot = lasso_coefs_lpc_plt,
#                 device = "png",
#                 width = 7,
#                 height = 7, 
#                 units = "in",
#                 dpi = "retina")


##### Smoothed MSE plot #####
lasso_mse_lpc <- read.csv('~/Desktop/Updating_LPOJ_sulci/R_stats_plots/dataframes/LPC_MSE_lasso.csv')

spline_int_lpc <- as.data.frame(spline(log(lasso_mse_lpc$alpha), lasso_mse_lpc$MSE))
lasso_mse_lpc_plt <- ggplot(spline_int_lpc, aes(x=x, y=y)) +
  # add smooth line
  # geom_smooth(color = "gray" ,
  #            span = .99,
  #            se = FALSE,
  #            size = 1) +
  geom_line(color = "gray", size = 1) +
  
  # vertical line to indiciate lowest MSE at best alpha level
  geom_vline(xintercept = log(0.05), 
             linetype="dotted", 
             size = 2) +
  # axis labels
  xlab ('log(alpha)')+
  ylab(' Cross Validated MSE') +
  
  #set_theme() +
  # set theme
  # theme_classic() +
  # 
  # # customize theme
  theme(
    # set axis title size
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    
    # set axis text size
    axis.text = element_text(color = "#333333", size = 14),
    
    # rotate x-axis 45 degrees
    # axis.text.x = element_text(angle = 0),
    # axis.text.y = element_text(angle = 0),
    
    # set strip text
    strip.text = element_text(size = 14),
    
    # remove strip itself
    strip.background = element_blank(),
    
    #change background pannel
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # change legend size
    legend.title = element_text(color = "#333333", size = 12),
    legend.text = element_text(color = "#333333", size = 12),
    panel.border = element_rect(colour = "black", 
                                fill=NA, size = 1)) +
  scale_x_continuous(limits = c(log(.0133), 0))  
#scale_y_continuous(limits = c(25, 29))
lasso_mse_lpc_plt

# ggplot2::ggsave(filename = "~/Downloads/lasso_mse_lpc_plt3.png",
#                 plot = lasso_mse_lpc_plt,
#                 device = "png",
#                 width = 7,
#                 height = 3, 
#                 units = "in",
#                 dpi = "retina")


##### Measured ~ Predicted plot #####
model_lpc_1a <- read.csv("~/Desktop/Updating_LPOJ_sulci/LASSO/plots/model_1a.csv")
cor.test(model_lpc_1a$Measured, model_lpc_1a$Predicted, method = "spearman")

model_lpc_1a.plot <- model_lpc_1a %>% 
  ggplot(aes(x = Measured, y = Predicted)) +
  
  geom_jitter(shape = 21, fill = "white", color = 'black') +
  
  
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  scale_x_continuous(n.breaks = 6) + 
  scale_y_continuous(n.breaks = 5) +
  theme(
    # set axis title size
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    
    # set axis text size
    axis.text = element_text(color = "#333333", size = 14),
    
    # rotate x-axis 45 degrees
    axis.text.x = element_text(angle = 0),
    axis.text.y = element_text(angle = 0),
    
    # set strip text
    strip.text = element_text(size = 14),
    
    # remove strip itself
    strip.background = element_blank(),
    
    #change background pannel
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    # change legend size
    legend.title = element_text(color = "#333333", size = 12),
    legend.text = element_text(color = "#333333", size = 12)) 
# panel.border = element_rect(colour = "black", 
#                             fill=NA, size = 1)) +

model_lpc_1a.plot

# ggplot2::ggsave(filename = "~/Downloads/model_lpc_1a.plot3.png",
#                 plot = model_lpc_1a.plot,
#                 device = "png",
#                 width = 7,
#                 height = 5, 
#                 units = "in",
#                 dpi = "retina")