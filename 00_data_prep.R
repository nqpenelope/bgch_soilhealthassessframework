library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(cowplot)
library(emmeans)
library(modelr)
library(effectsize)
library(ggpubr)
library(ggpmisc)

custom_theme <- theme_cowplot()
custom_theme <- custom_theme + theme_bw() + theme(strip.background = element_blank(), panel.border = element_rect(colour = "black"), strip.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.25)), axis.text = element_text(size = rel(1.25), colour = 1), legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5)), legend.key = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  panel.background = element_rect(colour = "black", linewidth = 1, fill = NA), plot.margin=unit(rep(0.5, 4),"cm"))
theme_set(custom_theme)

options(scipen = 999)
set.seed(42)

# management palette / labels 
pal_mgt <- c("cropland" = "#c5a05f", "reference" = "#134b46")
labl_mgt <- c("Row Crop", "Reference")

# reference-composition palette (04_references.R)
pal_ref <- c("#008F85", "#134b46", "black")
labl_ref <- c("Managed Reference", "Unmanaged Reference", "All References")

# SOM is converted to SOC throughout using the SHAP conversion factor
som_factor <- 1.9

dat # REDACTED: unfortunately we are not able to publish the raw data at this time due to the terms of farmer contracts and our shared data management plan 

dat_r <- dat %>% filter(!is.na(shmgtcat) & shmgtcat != 'other')
dat_r$shmgtcat <- factor(dat_r$shmgtcat, levels = c("baseline", "shms", "reference"))

dat_x <- dat_r %>% filter(covercrop.y == 'MISSING')
dat_r <- dat_r %>% filter(!slocid %in% dat_x$slocid)

dat_r <- dat_r %>% dplyr::select(-covercrop.x, -rotation.x)
dat_r <- dat_r %>% dplyr::rename(covercrop1 = covercrop.y, rotation = rotation.y)

# management practice variables - collapse any number of cover crops into binary has cover crops variable
# create cover crop / winter wheat variable
dat_r[which(dat_r$covercrop1 == '2' | dat_r$covercrop1 == '3'), 'covercrop1'] <- '1'
dat_r <- dat_r %>% mutate(covercrop2 = covercrop1)
dat_r[which(dat_r$rotation == '1'), 'covercrop2'] <- '1'

dat_r$am <- factor(dat_r$am)
levels(dat_r$am) <- c('No Manure', 'Manure')

dat_r$covercrop1 <- factor(dat_r$covercrop1)
levels(dat_r$covercrop1) <- c('No WW+CC', 'WW+CC')

dat_r$covercrop2 <- factor(dat_r$covercrop2)
levels(dat_r$covercrop2) <- c('No WW', 'WW')

dat_r %>% group_by(am) %>% summarize(n())
dat_r %>% group_by(covercrop1) %>% summarize(n())
dat_r %>% group_by(covercrop2) %>% summarize(n())

dat_r <- dat_r %>% mutate(Mgmt_binary = case_match(shmgtcat,
                                                   c('shms', 'baseline') ~ 'cropland',
                                                   'reference' ~ 'reference'),
                          Mgmt_binary = factor(Mgmt_binary, levels = c('cropland', 'reference')))

# tillage classes from GTIR
dat_r <- dat_r %>% mutate(tillage = case_when(gtir_max < 6 ~ 'no-till',
                                              gtir_max < 12 ~ 'min till',
                                              gtir_max >= 12 ~ 'conventional',
                                              .default = NA))

# remove one PMN outlier
dat_r[which(dat_r$pmn > 90 & dat_r$shmgtcat != 'reference'), 'pmn'] <- NA

dat_r %>% group_by(shmgtcat) %>% summarize(n())

dat_pm <- dat_r %>% filter(!is.na(pmn))
dat_sg <- dat_r %>% filter(!is.na(stab10_gmean))
dat_crop <- dat_r %>% filter(shmgtcat != 'reference')
dat_ref <- dat_r %>% filter(shmgtcat == 'reference')

# Reference-based (SHI) models

oc_mod <- lm(log(soc_0_15) ~ Mgmt_binary + clay_0_15, data = dat_r)
summary(oc_mod)
plot(fitted(oc_mod), resid(oc_mod)); abline(0, 0)

poxc_mod <- lm(log(pox_c) ~ Mgmt_binary + clay_0_15 + mat_daymet, data = dat_r)
summary(poxc_mod)
plot(fitted(poxc_mod), resid(poxc_mod)); abline(0, 0)

# AS-Image (stab10_gmean)
as_mod <- lm(log(stab10_gmean) ~ Mgmt_binary + clay_0_15, data = dat_sg)
summary(as_mod)
plot(fitted(as_mod), resid(as_mod)); abline(0, 0)

# AS-Yoder (stab_macroagg)
sm_mod <- lm(log(stab_macroagg) ~ Mgmt_binary + clay_0_15, data = dat_r)
summary(sm_mod)
plot(fitted(sm_mod), resid(sm_mod)); abline(0, 0)

# PMC follows Bower et al. (2026): separate cropland and reference models.
# Cropland depends on pH (management term dropped, cropland pooled);
# reference depends on clay.
cm_mod_crop <- lm(log(cminp_96hr_raw) ~ soil_ph, data = dat_crop)
summary(cm_mod_crop)
plot(fitted(cm_mod_crop), resid(cm_mod_crop)); abline(0, 0)

cm_mod_ref <- lm(log(cminp_96hr_raw) ~ clay_0_15, data = dat_ref)
summary(cm_mod_ref)
plot(fitted(cm_mod_ref), resid(cm_mod_ref)); abline(0, 0)

pmn_mod <- lm(log1p(pmn) ~ Mgmt_binary, data = dat_pm)
summary(pmn_mod)
plot(fitted(pmn_mod), resid(pmn_mod)); abline(0, 0)

# Proportion of reference

dat_r$oc_prop_potential <- dat_r$soc_0_15 / exp(predict(oc_mod, newdata = dat_r %>%
                                                          mutate(Mgmt_binary = 'reference')))

dat_r$px_prop_potential <- dat_r$pox_c / exp(predict(poxc_mod, newdata = dat_r %>%
                                                       mutate(Mgmt_binary = 'reference')))

dat_r$sm_prop_potential <- dat_r$stab_macroagg / exp(predict(sm_mod, newdata = dat_r %>%
                                                               mutate(Mgmt_binary = 'reference')))

dat_r$cm_prop_potential <- dat_r$cminp_96hr_raw / exp(predict(cm_mod_ref, newdata = dat_r %>%
                                                                mutate(Mgmt_binary = 'reference')))

# AS-Image and PMN are fitted on reduced frames, so join back on slocid
dat_sg$as_prop_potential <- dat_sg$stab10_gmean / exp(predict(as_mod, newdata = dat_sg %>%
                                                               mutate(Mgmt_binary = 'reference')))
dat_sg <- dat_sg %>% dplyr::select(slocid, as_prop_potential)

dat_pm$pmn_prop_potential <- dat_pm$pmn / expm1(predict(pmn_mod, newdata = dat_pm %>%
                                                          mutate(Mgmt_binary = 'reference')))
dat_pm <- dat_pm %>% dplyr::select(slocid, pmn_prop_potential)

dat_r <- dat_r %>% left_join(dat_sg, by = 'slocid') %>% left_join(dat_pm, by = 'slocid')

# refresh the subsets now that the proportions are attached
dat_crop <- dat_r %>% filter(shmgtcat != 'reference')
dat_ref <- dat_r %>% filter(shmgtcat == 'reference')

# Join distribution-based (SHAP) scores

# keep only the scores and texture group
dat_shap_scores #redacted

dat_both <- dat_r %>% left_join(dat_shap_scores, by = 'slocid')

# row-cropped fields only; this frame is in Figures 3 and 4 and Table 3
dat_both_crop <- dat_both %>%
  filter(shmgtcat != 'reference') %>%
  filter(!is.na(som_score))

dat_both_crop %>% group_by(shmgtcat) %>% summarise(n())
