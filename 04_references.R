# 04 - influence of reference site composition

source('final/code/00_data_prep.R')

# reference site descriptions
refs #redacted

dat_refcls <- dat_r %>%
  dplyr::select(-any_of('perennialdesc')) %>%
  left_join(refs %>% dplyr::select(slocid, perennialdesc), by = 'slocid') %>%
  mutate(shmgtcat_binref = case_when(
    shmgtcat == 'baseline' ~ 'cropland',
    shmgtcat == 'shms' ~ 'cropland',
    shmgtcat == 'reference' & perennialdesc %in% c('pasture', 'hayland', 'orchard', 'lawn') ~ 'managed',
    shmgtcat == 'reference' & perennialdesc %in% c('fencerow', 'forest', 'woodlot', 'meadow') ~ 'unmanaged'))

dat_refcls$shmgtcat_binref <- factor(dat_refcls$shmgtcat_binref,
                                     levels = c('cropland', 'managed', 'unmanaged'))

dat_refcls %>% group_by(shmgtcat_binref) %>% summarise(n())

dat_refcls_ref <- dat_refcls %>% filter(shmgtcat == 'reference')

# Models split by reference type

oc_mod_bin   <- lm(log(soc_0_15) ~ shmgtcat_binref + clay_0_15, data = dat_refcls)
poxc_mod_bin <- lm(log(pox_c) ~ shmgtcat_binref + clay_0_15 + mat_daymet, data = dat_refcls)
cm_mod_bin   <- lm(log(cminp_96hr_raw) ~ shmgtcat_binref + clay_0_15, data = dat_refcls_ref)
pmn_mod_bin  <- lm(log1p(pmn) ~ shmgtcat_binref, data = dat_refcls)
sm_mod_bin   <- lm(log(stab_macroagg) ~ shmgtcat_binref + clay_0_15, data = dat_refcls)
as_mod_bin   <- lm(log(stab10_gmean) ~ shmgtcat_binref + clay_0_15, data = dat_refcls)

summary(oc_mod_bin)
summary(poxc_mod_bin)
summary(cm_mod_bin)
summary(pmn_mod_bin)
summary(sm_mod_bin)
summary(as_mod_bin)

plot(fitted(oc_mod_bin), resid(oc_mod_bin)); abline(0, 0)
plot(fitted(poxc_mod_bin), resid(poxc_mod_bin)); abline(0, 0)
plot(fitted(cm_mod_bin), resid(cm_mod_bin)); abline(0, 0)
plot(fitted(pmn_mod_bin), resid(pmn_mod_bin)); abline(0, 0)
plot(fitted(sm_mod_bin), resid(sm_mod_bin)); abline(0, 0)
plot(fitted(as_mod_bin), resid(as_mod_bin)); abline(0, 0)

# Managed vs unmanaged contrasts

pairs_oc   <- pairs(emmeans(oc_mod_bin, 'shmgtcat_binref', data = dat_refcls), type = "response", reverse = TRUE)
pairs_poxc <- pairs(emmeans(poxc_mod_bin, 'shmgtcat_binref', data = dat_refcls), type = "response", reverse = TRUE)
pairs_cm   <- pairs(emmeans(cm_mod_bin, 'shmgtcat_binref', data = dat_refcls_ref), type = "response", reverse = TRUE)
pairs_pmn  <- pairs(emmeans(pmn_mod_bin, 'shmgtcat_binref', data = dat_refcls), reverse = TRUE)
pairs_sm   <- pairs(emmeans(sm_mod_bin, 'shmgtcat_binref', data = dat_refcls), type = "response", reverse = TRUE)
pairs_as   <- pairs(emmeans(as_mod_bin, 'shmgtcat_binref', data = dat_refcls), type = "response", reverse = TRUE)

pvals <- bind_rows(
  as.data.frame(pairs_oc)   %>% mutate(indicator = "SOC"),
  as.data.frame(pairs_poxc) %>% mutate(indicator = "POXC"),
  as.data.frame(pairs_cm)   %>% mutate(indicator = "PMC"),
  as.data.frame(pairs_pmn)  %>% mutate(indicator = "PMN"),
  as.data.frame(pairs_sm)   %>% mutate(indicator = "AS-Yoder"),
  as.data.frame(pairs_as)   %>% mutate(indicator = "AS-Image")
) %>%
  filter(grepl("unmanaged / managed|unmanaged - managed", contrast)) %>%
  mutate(label = paste0("p = ", signif(p.value, 1)))

pvals

# panel order used
ind_levels <- c("SOC", "POXC", "PMC", "PMN", "AS-Yoder", "AS-Image")

# How much greater are unmanaged references than managed ones

emm_oc   <- as.data.frame(emmeans(oc_mod_bin, 'shmgtcat_binref', type = 'response', data = dat_refcls))
emm_poxc <- as.data.frame(emmeans(poxc_mod_bin, 'shmgtcat_binref', type = 'response', data = dat_refcls))
emm_cm   <- as.data.frame(emmeans(cm_mod_bin, 'shmgtcat_binref', type = 'response', data = dat_refcls_ref))
emm_pmn  <- as.data.frame(regrid(emmeans(pmn_mod_bin, 'shmgtcat_binref', data = dat_refcls)))
emm_sm   <- as.data.frame(emmeans(sm_mod_bin, 'shmgtcat_binref', type = 'response', data = dat_refcls))
emm_as   <- as.data.frame(emmeans(as_mod_bin, 'shmgtcat_binref', type = 'response', data = dat_refcls))

ref_means <- bind_rows(
  emm_oc   %>% mutate(indicator = "SOC"),
  emm_poxc %>% mutate(indicator = "POXC"),
  emm_cm   %>% mutate(indicator = "PMC"),
  emm_pmn  %>% mutate(indicator = "PMN"),
  emm_sm   %>% mutate(indicator = "AS-Yoder"),
  emm_as   %>% mutate(indicator = "AS-Image")
) %>%
  filter(shmgtcat_binref != 'cropland') %>%
  dplyr::select(indicator, shmgtcat_binref, response) %>%
  pivot_wider(names_from = shmgtcat_binref, values_from = response)

ref_magnitude <- ref_means %>%
  left_join(pvals %>% dplyr::select(indicator, p.value), by = "indicator") %>%
  mutate(ratio = unmanaged / managed,
         pct_greater = round((ratio - 1) * 100),
         across(c(managed, unmanaged), ~ signif(.x, 3)),
         ratio = round(ratio, 3),
         p.value = signif(p.value, 2)) %>%
  mutate(indicator = factor(indicator, levels = ind_levels)) %>%
  arrange(indicator)

ref_magnitude

write.csv(ref_magnitude, 'final/outputs/fig7_reference_magnitudes.csv', row.names = FALSE)

# Clay-adjusted predictions

clay_seq <- seq(min(dat_refcls$clay_0_15), max(dat_refcls$clay_0_15), length.out = 100)
mat_med <- median(dat_refcls$mat_daymet, na.rm = TRUE)

new_bin <- expand.grid(clay_0_15 = clay_seq,
                       shmgtcat_binref = c('managed', 'unmanaged'),
                       mat_daymet = mat_med)

new_all <- data.frame(clay_0_15 = clay_seq, Mgmt_binary = 'reference', mat_daymet = mat_med)

pred_bin <- bind_rows(
  new_bin %>% mutate(indicator = "SOC",      fit = exp(predict(oc_mod_bin, newdata = .))),
  new_bin %>% mutate(indicator = "POXC",     fit = exp(predict(poxc_mod_bin, newdata = .))),
  new_bin %>% mutate(indicator = "PMC",      fit = exp(predict(cm_mod_bin, newdata = .))),
  new_bin %>% mutate(indicator = "PMN",      fit = expm1(predict(pmn_mod_bin, newdata = .))),
  new_bin %>% mutate(indicator = "AS-Yoder", fit = exp(predict(sm_mod_bin, newdata = .))),
  new_bin %>% mutate(indicator = "AS-Image", fit = exp(predict(as_mod_bin, newdata = .)))
)

pred_all <- bind_rows(
  new_all %>% mutate(indicator = "SOC",      fit = exp(predict(oc_mod, newdata = .))),
  new_all %>% mutate(indicator = "POXC",     fit = exp(predict(poxc_mod, newdata = .))),
  new_all %>% mutate(indicator = "PMC",      fit = exp(predict(cm_mod_ref, newdata = .))),
  new_all %>% mutate(indicator = "PMN",      fit = expm1(predict(pmn_mod, newdata = .))),
  new_all %>% mutate(indicator = "AS-Yoder", fit = exp(predict(sm_mod, newdata = .))),
  new_all %>% mutate(indicator = "AS-Image", fit = exp(predict(as_mod, newdata = .)))
)

ind_labels <- c("SOC" = "SOC (%)",
                "POXC" = "POXC (mg/kg)",
                "PMC" = "PMC (mg/kg CO2-C)",
                "PMN" = "PMN (mg/kg)",
                "AS-Yoder" = "AS-Yoder (%)",
                "AS-Image" = "AS-Image")

pred_bin$indicator <- factor(pred_bin$indicator, levels = ind_levels)
pred_all$indicator <- factor(pred_all$indicator, levels = ind_levels)
pvals$indicator    <- factor(pvals$indicator, levels = ind_levels)

# Figure 7

fig7 <- ggplot() +
  geom_line(data = pred_bin,
            aes(x = clay_0_15, y = fit, color = shmgtcat_binref), linewidth = 1.5) +
  geom_line(data = pred_all,
            aes(x = clay_0_15, y = fit, color = "all"), linewidth = 0.9, linetype = 'dashed') +
  geom_text(data = pvals, aes(x = Inf, y = -Inf, label = label),
            hjust = 1.12, vjust = -0.9, size = 4.2, color = pal_ref[2], inherit.aes = FALSE) +
  facet_wrap(~ indicator, scales = "free_y", ncol = 3,
             labeller = labeller(indicator = ind_labels)) +
  scale_color_manual(name = "",
                     breaks = c("managed", "unmanaged", "all"),
                     values = c(managed = pal_ref[1], unmanaged = pal_ref[2], all = pal_ref[3]),
                     labels = labl_ref) +
  scale_x_continuous(limits = c(0, 50)) +
  expand_limits(y = 0) +
  labs(x = "Clay (%)", y = "Indicator value") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 13, hjust = 0.5),
        panel.grid = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(1.1, "cm"))

fig7

ggsave('final/figs/fig7_reference_composition.png', fig7, width = 8.5, height = 6, dpi = 300)
