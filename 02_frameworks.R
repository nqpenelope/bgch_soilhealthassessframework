# 02 - outcomes of the two interpretive frameworks

source('final/code/00_data_prep.R')

library(grid)
library(gridExtra)
library(gtable)

# Figure 3A. Reference-based proportions
# Six indicators, in the order used across both panels

shi_levels <- c("SOC", "POXC", "PMC", "PMN", "AS-Yoder", "AS-Image")

shi_long <- dat_both_crop %>%
  dplyr::select(oc_prop_potential, px_prop_potential, sm_prop_potential,
                cm_prop_potential, pmn_prop_potential, as_prop_potential) %>%
  pivot_longer(cols = everything(), names_to = "indicator", values_to = "value") %>%
  mutate(indicator = dplyr::recode(indicator,
                                   oc_prop_potential  = "SOC",
                                   px_prop_potential  = "POXC",
                                   sm_prop_potential  = "AS-Yoder",
                                   cm_prop_potential  = "PMC",
                                   pmn_prop_potential = "PMN",
                                   as_prop_potential  = "AS-Image"),
         indicator = factor(indicator, levels = shi_levels))

p_shi <- ggplot(shi_long, aes(x = value)) +
  geom_density(linewidth = 1, na.rm = TRUE, colour = pal_mgt[["cropland"]]) +
  facet_wrap(~ indicator, scales = "fixed", ncol = 6, dir = "h", drop = FALSE) +
  coord_cartesian(xlim = c(0, 1.25)) +
  labs(x = "Reference-based score", y = "Density") +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5, size = 14),
    panel.grid = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

p_shi

# Figure 3B. Distribution-based scores

shap_levels <- c("SOC", "POXC", "PMC", "PMN", "AS-Yoder", "placeholder")

shap_long <- dat_both_crop %>%
  dplyr::select(som_score, active_c_score, agg_stab_score,
                respiration_score, pmn_score) %>%
  pivot_longer(cols = everything(), names_to = "indicator", values_to = "value") %>%
  mutate(indicator = dplyr::recode(indicator,
                                   som_score         = "SOC",
                                   active_c_score    = "POXC",
                                   agg_stab_score    = "AS-Yoder",
                                   respiration_score = "PMC",
                                   pmn_score         = "PMN"))

shap_long <- bind_rows(
  shap_long,
  tibble(indicator = "placeholder",
         value = NA_real_)
) %>%
  mutate(indicator = factor(indicator, levels = shap_levels))

p_shap <- ggplot(shap_long, aes(x = value)) +
  geom_density(linewidth = 1, na.rm = TRUE, colour = pal_mgt[["cropland"]]) +
  facet_wrap(~ indicator, scales = "free_x", ncol = 6, dir = "h", drop = FALSE,
             labeller = labeller(indicator = function(x) ifelse(grepl("placeholder", x), "", x))) +
  labs(x = "Distribution-based score", y = "Density") +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0.5, size = 14),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12)
  )

p_shap

# Figure 3. Stack the two panels

p_shi_grob  <- ggplotGrob(p_shi)
p_shap_grob <- ggplotGrob(p_shap)

# strip the axes off the empty sixth panel in the bottom row
for (nm in grep("axis-b-6|axis-l-6|axis-b-1-6|axis-l-6-1", p_shap_grob$layout$name, value = TRUE)) {
  p_shap_grob$grobs[[which(p_shap_grob$layout$name == nm)]] <- nullGrob()
}

fig3 <- arrangeGrob(
  p_shi_grob,
  p_shap_grob,
  ncol = 1,
  heights = unit(c(1, 1), c("null", "null")),
  padding = unit(0.5, "line")
)

grid.newpage()
grid.draw(fig3)
grid.text("A", x = 0.02, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("B", x = 0.02, y = 0.48, gp = gpar(fontsize = 18, fontface = "bold"))

png("final/figs/fig3_density.png", width = 12, height = 6, units = "in", res = 300)
grid.draw(fig3)
grid.text("A", x = 0.02, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("B", x = 0.02, y = 0.48, gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

# Figure 4. Rank comparison between the two approaches

a_cols <- c("oc_prop_potential", "px_prop_potential", "sm_prop_potential",
            "cm_prop_potential", "pmn_prop_potential")
b_cols <- c("som_score", "active_c_score", "agg_stab_score",
            "respiration_score", "pmn_score")

df_ranked <- dat_both_crop %>%
  mutate(across(all_of(c(a_cols, b_cols)),
                ~ rank(.x, ties.method = "average", na.last = "keep"),
                .names = "{.col}_rank"))

plot_data <- tibble()
cor_results <- tibble(pair = character(), spearman_rho = numeric())

for (i in seq_along(a_cols)) {
  pair_label <- paste(a_cols[i], "vs", b_cols[i])

  temp_df <- df_ranked %>%
    transmute(x_rank = .data[[paste0(a_cols[i], "_rank")]],
              y_rank = .data[[paste0(b_cols[i], "_rank")]],
              pair = pair_label)

  plot_data <- bind_rows(plot_data, temp_df)

  cor_results <- bind_rows(cor_results, tibble(
    pair = pair_label,
    spearman_rho = cor(temp_df$x_rank, temp_df$y_rank, method = "spearman", use = "complete.obs")
  ))
}

facet_order <- c(
  "oc_prop_potential vs som_score",
  "px_prop_potential vs active_c_score",
  "cm_prop_potential vs respiration_score",
  "pmn_prop_potential vs pmn_score",
  "sm_prop_potential vs agg_stab_score"
)

facet_labels <- c(
  "oc_prop_potential vs som_score" = "SOC",
  "px_prop_potential vs active_c_score" = "POXC",
  "cm_prop_potential vs respiration_score" = "PMC",
  "pmn_prop_potential vs pmn_score" = "PMN",
  "sm_prop_potential vs agg_stab_score" = "AS-Yoder"
)

plot_data$pair   <- factor(plot_data$pair, levels = facet_order, labels = facet_labels)
cor_results$pair <- factor(cor_results$pair, levels = facet_order, labels = facet_labels)

p_rank <- ggplot(plot_data, aes(x = x_rank, y = y_rank)) +
  geom_point(colour = pal_mgt[["cropland"]]) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 1) +
  facet_wrap(~ pair, scales = "free") +
  labs(x = "Rank of reference-based score", y = "Rank of distribution-based score") +
  geom_text(data = cor_results,
            aes(x = Inf, y = -Inf, label = paste0("ρ = ", round(spearman_rho, 2))),
            hjust = 1.1, vjust = -1.1, inherit.aes = FALSE, size = 4, color = "black") +
  theme_minimal()

p_rank

ggsave('final/figs/fig4_rank_comparison.png', p_rank, width = 5.5, height = 4.5, dpi = 300)

cor_results

write.csv(cor_results, 'final/outputs/fig4_spearman_rho.csv', row.names = FALSE)

# Model-based ratings; PMC uses its cropland-only pH model (cm_mod_crop), following Bower et al. (2026).

# every field predicted as though it were cropland
dat_asif_crop <- dat_both_crop %>% mutate(Mgmt_binary = factor('cropland', levels = levels(dat_r$Mgmt_binary)))

oc_pi60  <- exp(predict(oc_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.60))
px_pi60  <- exp(predict(poxc_mod, newdata = dat_asif_crop, interval = 'prediction', level = 0.60))
as_pi60  <- exp(predict(as_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.60))
sm_pi60  <- exp(predict(sm_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.60))
cm_pi60  <- exp(predict(cm_mod_crop,  newdata = dat_asif_crop, interval = 'prediction', level = 0.60))
pmn_pi60 <- expm1(predict(pmn_mod, newdata = dat_asif_crop, interval = 'prediction', level = 0.60))

oc_pi20  <- exp(predict(oc_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.20))
px_pi20  <- exp(predict(poxc_mod, newdata = dat_asif_crop, interval = 'prediction', level = 0.20))
as_pi20  <- exp(predict(as_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.20))
sm_pi20  <- exp(predict(sm_mod,   newdata = dat_asif_crop, interval = 'prediction', level = 0.20))
cm_pi20  <- exp(predict(cm_mod_crop,  newdata = dat_asif_crop, interval = 'prediction', level = 0.20))
pmn_pi20 <- expm1(predict(pmn_mod, newdata = dat_asif_crop, interval = 'prediction', level = 0.20))

dat_ratings <- dat_both_crop %>%
  mutate(
    oc_crop_rate = case_when(is.na(soc_0_15) ~ NA_character_,
                             soc_0_15 <= oc_pi60[, 'lwr'] ~ 'very low',
                             soc_0_15 <= oc_pi20[, 'lwr'] ~ 'low',
                             soc_0_15 <= oc_pi20[, 'upr'] ~ 'medium',
                             soc_0_15 <= oc_pi60[, 'upr'] ~ 'high',
                             .default = 'very high'),
    px_crop_rate = case_when(is.na(pox_c) ~ NA_character_,
                             pox_c <= px_pi60[, 'lwr'] ~ 'very low',
                             pox_c <= px_pi20[, 'lwr'] ~ 'low',
                             pox_c <= px_pi20[, 'upr'] ~ 'medium',
                             pox_c <= px_pi60[, 'upr'] ~ 'high',
                             .default = 'very high'),
    as_crop_rate = case_when(is.na(stab10_gmean) ~ NA_character_,
                             stab10_gmean <= as_pi60[, 'lwr'] ~ 'very low',
                             stab10_gmean <= as_pi20[, 'lwr'] ~ 'low',
                             stab10_gmean <= as_pi20[, 'upr'] ~ 'medium',
                             stab10_gmean <= as_pi60[, 'upr'] ~ 'high',
                             .default = 'very high'),
    sm_crop_rate = case_when(is.na(stab_macroagg) ~ NA_character_,
                             stab_macroagg <= sm_pi60[, 'lwr'] ~ 'very low',
                             stab_macroagg <= sm_pi20[, 'lwr'] ~ 'low',
                             stab_macroagg <= sm_pi20[, 'upr'] ~ 'medium',
                             stab_macroagg <= sm_pi60[, 'upr'] ~ 'high',
                             .default = 'very high'),
    cm_crop_rate = case_when(is.na(cminp_96hr_raw) ~ NA_character_,
                             cminp_96hr_raw <= cm_pi60[, 'lwr'] ~ 'very low',
                             cminp_96hr_raw <= cm_pi20[, 'lwr'] ~ 'low',
                             cminp_96hr_raw <= cm_pi20[, 'upr'] ~ 'medium',
                             cminp_96hr_raw <= cm_pi60[, 'upr'] ~ 'high',
                             .default = 'very high'),
    pmn_crop_rate = case_when(is.na(pmn) ~ NA_character_,
                              pmn <= pmn_pi60[, 'lwr'] ~ 'very low',
                              pmn <= pmn_pi20[, 'lwr'] ~ 'low',
                              pmn <= pmn_pi20[, 'upr'] ~ 'medium',
                              pmn <= pmn_pi60[, 'upr'] ~ 'high',
                              .default = 'very high'),
    across(ends_with('_crop_rate'),
           ~ factor(.x, levels = c('very low', 'low', 'medium', 'high', 'very high'))))

table(dat_ratings$oc_crop_rate, useNA = 'ifany')

ratings_nomanure <- dat_ratings %>%
  dplyr::select(slocid, ends_with('_crop_rate')) %>%
  pivot_longer(cols = -slocid, names_to = 'var', values_to = 'crop_rate') %>%
  mutate(var = sub('_crop_rate', '', var))

write.csv(ratings_nomanure, 'final/outputs/shi_ratings_nomanure.csv', row.names = FALSE)

rate_bands <- c("very low", "low", "medium", "high", "very high")

soc_rated <- dat_ratings %>%
  filter(!is.na(som_score), !is.na(oc_crop_rate)) %>%
  mutate(
    dist_class = factor(case_when(som_score < 20 ~ "very low",
                                  som_score < 40 ~ "low",
                                  som_score < 60 ~ "medium",
                                  som_score < 80 ~ "high",
                                  .default = "very high"),
                        levels = rate_bands),
    ref_class = factor(oc_crop_rate, levels = rate_bands))

soc_xtab <- soc_rated %>%
  count(dist_class, ref_class, .drop = FALSE) %>%
  pivot_wider(names_from = ref_class, values_from = n, values_fill = 0) %>%
  arrange(dist_class)

soc_xtab

write.csv(soc_xtab, 'final/outputs/table3_rating_crosstab.csv', row.names = FALSE)

# row percentages
soc_xtab_props <- soc_xtab %>%
  pivot_longer(cols = -dist_class, names_to = "ref_class", values_to = "n") %>%
  group_by(dist_class) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup()

soc_xtab_props

# how often the two frameworks agree
mean(soc_rated$dist_class == soc_rated$ref_class)
mean(abs(as.integer(soc_rated$dist_class) - as.integer(soc_rated$ref_class)) <= 1)

# Table 4. Four example fields, all loam

tab4_slocids <- c(2513, 26321, 26265, 26376)

table4 <- dat_ratings %>%
  filter(slocid %in% tab4_slocids) %>%
  mutate(
    Management = paste0(
      if_else(tillage == 'no-till', 'No-till', 'Annual tillage'),
      if_else(am == 'Manure', ', manure', ''),
      if_else(rotation == '1', ', corn-soy-wheat', ', corn-soy'),
      if_else(covercrop1 == 'WW+CC', ', cover crop', '')),
    across(c(som_score, active_c_score, respiration_score, pmn_score, agg_stab_score),
           ~ case_when(.x < 20 ~ 'very low',
                       .x < 40 ~ 'low',
                       .x < 60 ~ 'medium',
                       .x < 80 ~ 'high',
                       .default = 'very high'),
           .names = "{.col}_band")) %>%
  arrange(match(slocid, tab4_slocids)) %>%
  dplyr::select(slocid, sitecode, mucode, Management,
                ref_SOC = oc_crop_rate, ref_POXC = px_crop_rate, ref_PMC = cm_crop_rate,
                ref_PMN = pmn_crop_rate, ref_ASYoder = sm_crop_rate, ref_ASImage = as_crop_rate,
                dist_SOC = som_score_band, dist_POXC = active_c_score_band,
                dist_PMC = respiration_score_band, dist_PMN = pmn_score_band,
                dist_ASYoder = agg_stab_score_band)

table4

write.csv(table4, 'final/outputs/table4_example_fields.csv', row.names = FALSE)

# the underlying values, for checking the bands
dat_ratings %>%
  filter(slocid %in% tab4_slocids) %>%
  arrange(match(slocid, tab4_slocids)) %>%
  dplyr::select(slocid, newtexture_0_15, clay_0_15, soc_0_15, oc_prop_potential, som_score,
                cminp_96hr_raw, cm_prop_potential, respiration_score,
                stab10_gmean, as_prop_potential)
