# 01 - indicator relationships and sensitivity to management

source('final/code/00_data_prep.R')

library(GGally)
library(factoextra)
library(car)

# Figure 1. Pairwise correlations among indicators

dat_corr <- dat_r %>% dplyr::select(soc_0_15, pox_c, stab10_gmean, stab_macroagg, cminp_96hr_raw, pmn)
colnames(dat_corr) <- c("SOC", "POXC", "AS-Image", "AS-Yoder", "PMC", "PMN")

p_corr <- ggpairs(
  dat_corr,
  columns = 1:6,
  upper = list(continuous = wrap("cor", size = 6, digits = 2)),
  lower = list(continuous = wrap("points", alpha = 0.5)),
  diag  = list(continuous = wrap("densityDiag"))
)

for (i in 1:6) {
  for (j in 1:6) {
    p_corr[i, j] <- p_corr[i, j] +
      theme_bw(base_size = 10) +
      theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))
    if (i >= j) p_corr[i, j] <- p_corr[i, j] + scale_x_continuous(limits = c(0, NA))
    if (i >  j) p_corr[i, j] <- p_corr[i, j] + scale_y_continuous(limits = c(0, NA))
  }
}

p_corr <- p_corr + theme(text = element_text(size = 12))
p_corr

ggsave('final/figs/fig1_indicators_corr.png', p_corr, width = 9, height = 9, dpi = 300)

# Figure 2. PCA of indicators by management

dat_pca <- dat_r %>%
  dplyr::select(soc_0_15, pox_c, stab10_gmean, stab_macroagg, cminp_96hr_raw, pmn,
                newsand_0_15, Mgmt_binary) %>%
  drop_na()

colnames(dat_pca)[1:6] <- c("SOC", "POXC", "AS-Image", "AS-Yoder", "PMC", "PMN")

dat_pca$Mgmt_binary <- factor(dat_pca$Mgmt_binary,
                              levels = c("cropland", "reference"),
                              labels = labl_mgt)

pca <- prcomp(dat_pca %>% dplyr::select(SOC, POXC, `AS-Image`, `AS-Yoder`, PMC, PMN), scale = TRUE)

sand_std <- scale(dat_pca$newsand_0_15)
sand_arrow <- data.frame(x = cor(sand_std, pca$x[, 1]), y = cor(sand_std, pca$x[, 2]))
sand_arrow

# PCA coordinates in the same scaling factoextra plots in
scores <- as.data.frame(factoextra::get_pca_ind(pca)$coord[, 1:2])
colnames(scores) <- c("Dim1", "Dim2")
scores$Group <- dat_pca$Mgmt_binary

# 50% ellipse
ell_df <- data.frame()
for (g in levels(scores$Group)) {
  df_g <- scores %>% filter(Group == g)
  ell <- car::dataEllipse(x = df_g$Dim1, y = df_g$Dim2, levels = 0.5,
                          draw = FALSE, robust = TRUE, segments = 200)
  ell <- as.data.frame(ell)
  colnames(ell) <- c("Dim1", "Dim2")
  ell$Group <- g
  ell_df <- bind_rows(ell_df, ell)
}

p_base <- fviz_pca_biplot(
  pca,
  geom.ind    = "point",
  habillage   = dat_pca$Mgmt_binary,
  addEllipses = FALSE,
  col.var     = "#CE6650",
  repel       = TRUE
)

# tidy up the axis labels
xl <- round(as.numeric(sub(".*\\((.*)%\\).*", "\\1", p_base$labels$x)))
yl <- round(as.numeric(sub(".*\\((.*)%\\).*", "\\1", p_base$labels$y)))

p_pca <- p_base +
  geom_path(data = ell_df, aes(x = Dim1, y = Dim2, color = Group),
            linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_manual(values = unname(pal_mgt), labels = labl_mgt, name = "Management") +
  scale_shape_discrete(name = "Management") +
  geom_segment(data = sand_arrow, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.9,
               color = "#CE6650", inherit.aes = FALSE) +
  annotate("text", x = sand_arrow$x * 1.15, y = sand_arrow$y * 1.15,
           label = "Sand", color = "#CE6650", size = 4) +
  labs(title = NULL,
       x = paste0("Dimension 1 (", xl, "%)"),
       y = paste0("Dimension 2 (", yl, "%)"))

p_pca

ggsave('final/figs/fig2_pca.png', p_pca, width = 7, height = 5, dpi = 300)

# Management practice models

#   SOC       clay + Manure                    (orig paper)
#   POXC      clay + MAT + Manure              (lowest BIC)
#   AS-Image  clay + Tillage Intensity + WW+CC (orig paper)
#   AS-Yoder  clay                             (lowest BIC, no practice term)
#   PMC       pH + Manure + WW                 (orig paper)
#   PMN       Manure                           (lowest BIC)

dat_agg <- dat_crop %>% filter(!is.na(stab10_gmean))

# Table 2. Partial eta-squared for the raw indicators

oc_ancova   <- lm(log(soc_0_15) ~ clay_0_15 + am, data = dat_crop)
poxc_ancova <- lm(log(pox_c) ~ clay_0_15 + mat_daymet + am, data = dat_crop)
sl_ancova   <- lm(log(stab10_gmean) ~ clay_0_15 + gtir_max + covercrop1, data = dat_agg)
as_ancova   <- lm(log(stab_macroagg) ~ clay_0_15, data = dat_crop)
pmc_ancova  <- lm(log(cminp_96hr_raw) ~ soil_ph + am + covercrop2, data = dat_crop)
pmn_ancova  <- lm(log1p(pmn) ~ am, data = dat_crop)

summary(oc_ancova)
summary(poxc_ancova)
summary(sl_ancova)
summary(as_ancova)
summary(pmc_ancova)
summary(pmn_ancova)

# Type II sums of squares, so each term is adjusted for every other
# term and the result does not depend on the order terms are entered
eta_raw <- bind_rows(
  eta_squared(Anova(pmc_ancova, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMC"),
  eta_squared(Anova(sl_ancova, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "AS-Image"),
  eta_squared(Anova(pmn_ancova, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMN"),
  eta_squared(Anova(oc_ancova, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "SOC"),
  eta_squared(Anova(poxc_ancova, type = 2)) %>% as.data.frame() %>% mutate(Indicator = "POXC"),
  eta_squared(Anova(as_ancova, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "AS-Yoder")
) %>%
  mutate(Parameter = dplyr::recode(Parameter,
                                   clay_0_15  = "Clay",
                                   mat_daymet = "Mean Annual Temperature",
                                   soil_ph    = "Soil pH",
                                   am         = "Manure",
                                   covercrop1 = "WW+CC",
                                   covercrop2 = "WW",
                                   gtir_max   = "Tillage Intensity")) %>%
  # one-way models return Eta2 rather than Eta2_partial; they are the same thing
  mutate(Eta2_partial = coalesce(Eta2_partial, Eta2)) %>%
  dplyr::select(Indicator, Parameter, Eta2_partial) %>%
  pivot_wider(names_from = Parameter, values_from = Eta2_partial) %>%
  dplyr::select(Indicator, Clay, `Mean Annual Temperature`, `Soil pH`,
                Manure, `WW+CC`, WW, `Tillage Intensity`) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

eta_raw

write.csv(eta_raw, 'final/outputs/table2_eta2_raw_indicators.csv', row.names = FALSE)

# Table S1. Partial eta-squared after distribution-based scoring

oc_score_mod   <- lm(som_score ~ clay_0_15 + am, data = dat_both_crop)
poxc_score_mod <- lm(active_c_score ~ clay_0_15 + mat_daymet + am, data = dat_both_crop)
as_score_mod   <- lm(agg_stab_score ~ clay_0_15, data = dat_both_crop)
pmc_score_mod  <- lm(respiration_score ~ soil_ph + am + covercrop2, data = dat_both_crop)
pmn_score_mod  <- lm(pmn_score ~ am, data = dat_both_crop)

summary(oc_score_mod)
summary(poxc_score_mod)
summary(as_score_mod)
summary(pmc_score_mod)
summary(pmn_score_mod)

eta_shap <- bind_rows(
  eta_squared(Anova(pmc_score_mod, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMC"),
  eta_squared(Anova(pmn_score_mod, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMN"),
  eta_squared(Anova(oc_score_mod, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "SOC"),
  eta_squared(Anova(poxc_score_mod, type = 2)) %>% as.data.frame() %>% mutate(Indicator = "POXC"),
  eta_squared(Anova(as_score_mod, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "AS-Yoder")
) %>%
  mutate(Parameter = dplyr::recode(Parameter,
                                   clay_0_15  = "Clay",
                                   mat_daymet = "Mean Annual Temperature",
                                   soil_ph    = "Soil pH",
                                   am         = "Manure",
                                   covercrop1 = "WW+CC",
                                   covercrop2 = "WW",
                                   gtir_max   = "Tillage Intensity")) %>%
  mutate(Eta2_partial = coalesce(Eta2_partial, Eta2)) %>%
  dplyr::select(Indicator, Parameter, Eta2_partial) %>%
  pivot_wider(names_from = Parameter, values_from = Eta2_partial) %>%
  dplyr::select(Indicator, Clay, `Mean Annual Temperature`, `Soil pH`, Manure, WW) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

eta_shap

write.csv(eta_shap, 'final/outputs/tableS1_eta2_after_scores.csv', row.names = FALSE)

# Table S2. Partial eta-squared after reference-based proportions

oc_prop_mod   <- lm(oc_prop_potential ~ clay_0_15 + am, data = dat_crop)
poxc_prop_mod <- lm(px_prop_potential ~ clay_0_15 + mat_daymet + am, data = dat_crop)
sl_prop_mod   <- lm(as_prop_potential ~ clay_0_15 + gtir_max + covercrop1, data = dat_crop)
as_prop_mod   <- lm(sm_prop_potential ~ clay_0_15, data = dat_crop)
pmc_prop_mod  <- lm(cm_prop_potential ~ soil_ph + am + covercrop2, data = dat_crop)
pmn_prop_mod  <- lm(pmn_prop_potential ~ am, data = dat_crop)

summary(oc_prop_mod)
summary(poxc_prop_mod)
summary(sl_prop_mod)
summary(as_prop_mod)
summary(pmc_prop_mod)
summary(pmn_prop_mod)

eta_shi <- bind_rows(
  eta_squared(Anova(pmc_prop_mod, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMC"),
  eta_squared(Anova(sl_prop_mod, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "AS-Image"),
  eta_squared(Anova(pmn_prop_mod, type = 2))  %>% as.data.frame() %>% mutate(Indicator = "PMN"),
  eta_squared(Anova(oc_prop_mod, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "SOC"),
  eta_squared(Anova(poxc_prop_mod, type = 2)) %>% as.data.frame() %>% mutate(Indicator = "POXC"),
  eta_squared(Anova(as_prop_mod, type = 2))   %>% as.data.frame() %>% mutate(Indicator = "AS-Yoder")
) %>%
  mutate(Parameter = dplyr::recode(Parameter,
                                   clay_0_15  = "Clay",
                                   mat_daymet = "Mean Annual Temperature",
                                   soil_ph    = "Soil pH",
                                   am         = "Manure",
                                   covercrop1 = "WW+CC",
                                   covercrop2 = "WW",
                                   gtir_max   = "Tillage Intensity")) %>%
  mutate(Eta2_partial = coalesce(Eta2_partial, Eta2)) %>%
  dplyr::select(Indicator, Parameter, Eta2_partial) %>%
  pivot_wider(names_from = Parameter, values_from = Eta2_partial) %>%
  dplyr::select(Indicator, Clay, `Mean Annual Temperature`, `Soil pH`,
                Manure, `WW+CC`, WW, `Tillage Intensity`) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

eta_shi

write.csv(eta_shi, 'final/outputs/tableS2_eta2_after_proportions.csv', row.names = FALSE)
