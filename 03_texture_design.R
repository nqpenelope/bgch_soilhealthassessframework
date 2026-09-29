# 03 - influence of texture category and of experimental design

source('final/code/00_data_prep.R')

library(purrr)
library(scales)

# SHAP SOM scoring curves by texture group

z_break <- qnorm(c(0.2, 0.4, 0.6, 0.8))

som_breaks <- tibble(
  texture = rep(c("Coarse", "Medium", "Fine"), each = 4),
  score   = rep(c(20, 40, 60, 80), 3),
  som     = c(2.45, 3.18, 3.77, 4.50,
              3.45, 4.15, 4.74, 5.42,
              4.00, 4.73, 5.35, 6.05)
)

som_params <- tibble()
for (tex in c("Coarse", "Medium", "Fine")) {
  b <- som_breaks %>% filter(texture == tex)
  fit <- lm(b$som ~ z_break)
  som_params <- bind_rows(som_params, tibble(texture = tex,
                                             mu = coef(fit)[1],
                                             sigma = coef(fit)[2]))
}

som_params

# Table 5. Influence of texture category assignment

dat_regA <- dat_crop %>%
  filter(newsand_0_15 >= 49, newsand_0_15 <= 55, clay_0_15 >= 12, clay_0_15 <= 18) %>%
  mutate(region = "A")

dat_regB <- dat_crop %>%
  filter(clay_0_15 >= 25, clay_0_15 <= 31, newsand_0_15 > 24) %>%
  mutate(region = "B")

dat_regions <- bind_rows(dat_regA, dat_regB)

region_summary <- dat_regions %>%
  group_by(region) %>%
  summarise(n = n(),
            clay_mean = round(mean(clay_0_15), 0),
            clay_sd = round(sd(clay_0_15), 0),
            sand_mean = round(mean(newsand_0_15), 0),
            sand_sd = round(sd(newsand_0_15), 0),
            median_soc = median(soc_0_15),
            reference_proportion = round(median(oc_prop_potential, na.rm = TRUE), 2),
            .groups = "drop") %>%
  mutate(median_som = median_soc * som_factor)

region_summary

# score each region's median SOM under the two texture groups it could be assigned to
table5 <- bind_rows(
  region_summary %>% filter(region == "A") %>% mutate(texture = "Coarse (sandy loam)"),
  region_summary %>% filter(region == "A") %>% mutate(texture = "Medium (loam)"),
  region_summary %>% filter(region == "B") %>% mutate(texture = "Medium (loam)"),
  region_summary %>% filter(region == "B") %>% mutate(texture = "Fine (clay loam)")
) %>%
  mutate(texture_group = sub(" .*", "", texture)) %>%
  left_join(som_params, by = c("texture_group" = "texture")) %>%
  mutate(score = round(pnorm((median_som - mu) / sigma) * 100),
         score_rating = case_when(score < 20 ~ "very low",
                                  score < 40 ~ "low",
                                  score < 60 ~ "medium",
                                  score < 80 ~ "high",
                                  .default = "very high")) %>%
  mutate(median_soc = round(median_soc, 1), median_som = round(median_som, 1)) %>%
  dplyr::select(region, n, clay_mean, clay_sd, sand_mean, sand_sd, median_soc,
                reference_proportion, median_som, texture, score, score_rating)

table5

write.csv(table5, 'final/outputs/table5_texture_boundary_scores.csv', row.names = FALSE)

# Figure 5 (left). SOM scoring curves by texture group, drawn on a SOC axis, over the five score bands.

soc_grid <- seq(0, 9 / som_factor, length.out = 400)

som_curves <- tibble()
for (tex in c("Coarse", "Medium", "Fine")) {
  p <- som_params %>% filter(texture == tex)
  som_curves <- bind_rows(som_curves, tibble(
    texture = tex,
    soc = soc_grid,
    score = pnorm((soc_grid * som_factor - p$mu) / p$sigma) * 100
  ))
}

som_curves$texture <- factor(som_curves$texture, levels = c("Coarse", "Medium", "Fine"))

p_curves <- ggplot(som_curves, aes(x = soc, y = score, linetype = texture)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0,  ymax = 20),  fill = "#d73027", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 40),  fill = "#fc8d59", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 40, ymax = 60),  fill = "#fee08b", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 60, ymax = 80),  fill = "#d9ef8b", inherit.aes = FALSE) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 80, ymax = 100), fill = "#1a9850", inherit.aes = FALSE) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(name = "Texture", values = c(Coarse = "solid", Medium = "dotted", Fine = "dashed")) +
  scale_x_continuous(limits = c(0, 9 / som_factor), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Soil Organic Carbon (%)", y = "Distribution-based score") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.78, 0.35),
        legend.background = element_rect(fill = "white", colour = "black"))

p_curves

ggsave('final/figs/fig5_left_som_curves_soc_axis.png', p_curves, width = 5, height = 4, dpi = 300)

# Figure 6 and Table 6. Influence of experimental design

dat_design <- dat_r %>%
  mutate(soiltex = ifelse(is.na(soiltexture_sieve_0_15), soiltexture_hyd_0_15, soiltexture_sieve_0_15),
         som = soc_0_15 * som_factor,
         shmgt_binary = case_when(shmgtcat %in% c("baseline", "shms") ~ "Cropland",
                                  shmgtcat == "reference" ~ "Reference"),
         shap_as = stab_macroagg^5,
         # the high PMN values blanked for the reference models in
         # 00_data_prep.R are kept here, because the scoring curves
         # are built from the distribution of what was measured
         pmn = dat$pmn[match(slocid, dat$slocid)],
         shap_pmn = pmn^0.5)

dat_design$shmgt_binary <- factor(dat_design$shmgt_binary, levels = c("Cropland", "Reference"))

medium_df <- dat_design %>% filter(soiltex %in% c('silt loam', 'loam', 'silt', 'sandy clay loam'))

medium_df %>% group_by(shmgt_binary) %>% summarise(n())

# OTSP scoring curves, supplied as breakpoints
shap_curves <- readxl::read_excel("data/shap_curves_loam.xlsx") %>%
  mutate(Score = as.numeric(Score_Level),
         Indicator_Value = case_when(Indicator == "agg_stab" ~ Value^5,
                                     Indicator == "pmn" ~ Value^0.5,
                                     .default = Value),
         Curve = "75% Row Crop") %>%
  dplyr::select(Indicator, Indicator_Value, Score, Curve) %>%
  arrange(Indicator, Indicator_Value)

indicator_list <- list(som = "som", respiration = "cminp_96hr_raw", agg_stab = "shap_as",
                       pmn = "shap_pmn", active_c = "pox_c")

indicator_labels <- c(som = "SOC (%)", pmn = "PMN", respiration = "PMC",
                      agg_stab = "AS-Yoder", active_c = "POXC")

xlims <- list(som = c(0, 15), respiration = c(0, 60), agg_stab = c(0, 10e+09),
              pmn = c(0, 10), active_c = c(0, 1250))

set.seed(123)
n_iter <- 100
target_prop <- 0.8

# resample the medium-textured soils
sim_params <- tibble()
for (cat_focus in c("Cropland", "Reference")) {
  for (sim_id in 1:n_iter) {
    n_total <- nrow(medium_df)
    n_target <- round(n_total * target_prop)
    n_others <- n_total - n_target

    df_target <- medium_df %>%
      filter(shmgt_binary == cat_focus) %>%
      sample_n(n_target, replace = n_target > sum(medium_df$shmgt_binary == cat_focus))

    df_others <- medium_df %>%
      filter(shmgt_binary != cat_focus) %>%
      sample_n(n_others, replace = n_others > sum(medium_df$shmgt_binary != cat_focus))

    df_sim <- bind_rows(df_target, df_others)

    for (ind_name in names(indicator_list)) {
      vals <- df_sim[[indicator_list[[ind_name]]]]
      sim_params <- bind_rows(sim_params, tibble(
        Category = cat_focus, Simulation = sim_id, Indicator = ind_name,
        mu = mean(vals, na.rm = TRUE), sigma = sd(vals, na.rm = TRUE)
      ))
    }
  }
}

# mean simulated curve per indicator, on a grid of indicator values
sims_summaries <- tibble()
for (ind_name in names(indicator_list)) {
  grid <- tibble(Indicator_Value = seq(xlims[[ind_name]][1], xlims[[ind_name]][2], length.out = 400))

  sims_summaries <- bind_rows(sims_summaries,
    sim_params %>%
      filter(Indicator == ind_name) %>%
      tidyr::crossing(grid) %>%
      mutate(Score = pnorm((Indicator_Value - mu) / sigma) * 100) %>%
      group_by(Category, Indicator_Value) %>%
      summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
      mutate(Indicator = ind_name,
             Curve = ifelse(Category == "Cropland", "80% Row Crop", "20% Row Crop")) %>%
      dplyr::select(Indicator, Curve, Indicator_Value, Score)
  )
}

# Greenbelt curve from the full medium-textured dataset
greenbelt_curves <- tibble()
for (ind_name in names(indicator_list)) {
  vals <- medium_df[[indicator_list[[ind_name]]]]
  grid <- seq(xlims[[ind_name]][1], xlims[[ind_name]][2], length.out = 400)

  greenbelt_curves <- bind_rows(greenbelt_curves, tibble(
    Indicator = ind_name,
    Curve = "57% Row Crop",
    Indicator_Value = grid,
    Score = pnorm((grid - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)) * 100
  ))
}

curve_tbl <- bind_rows(sims_summaries, greenbelt_curves, shap_curves)

curve_breaks <- c("20% Row Crop", "57% Row Crop", "80% Row Crop", "75% Row Crop")
curve_labels <- c("20% Row Crop (simulated)", "57% Row Crop (Greenbelt)",
                  "80% Row Crop (simulated)", "75% Row Crop (OTSP)")
curve_colors <- c("20% Row Crop" = "grey60", "57% Row Crop" = "grey40",
                  "80% Row Crop" = "grey20", "75% Row Crop" = "black")
curve_lines  <- c("20% Row Crop" = "solid", "57% Row Crop" = "solid",
                  "80% Row Crop" = "solid", "75% Row Crop" = "dotted")

plot_list <- list()
for (ind_name in names(indicator_list)) {
  curve_i <- curve_tbl %>% filter(Indicator == ind_name)

  # SOM is drawn on a SOC axis; the other indicators keep their own units
  if (ind_name == "som") {
    curve_i <- curve_i %>% mutate(Indicator_Value = Indicator_Value / som_factor)
    xlim_i <- xlims[[ind_name]] / som_factor
  } else {
    xlim_i <- xlims[[ind_name]]
  }

  plot_list[[ind_name]] <- ggplot(curve_i, aes(x = Indicator_Value, y = Score,
                                               linetype = Curve, color = Curve)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0,  ymax = 20),  fill = "#d73027", inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 40),  fill = "#fc8d59", inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 40, ymax = 60),  fill = "#fee08b", inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 60, ymax = 80),  fill = "#d9ef8b", inherit.aes = FALSE) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 80, ymax = 100), fill = "#1a9850", inherit.aes = FALSE) +
    geom_line(linewidth = 1.15) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(20, 100, by = 20), expand = c(0, 0)) +
    scale_x_continuous(limits = xlim_i, expand = c(0, 0),
                       labels = if (ind_name == "agg_stab") label_scientific(digits = 1) else label_number()) +
    scale_linetype_manual(name = "Row Crop Proportion", breaks = curve_breaks,
                          values = curve_lines, labels = curve_labels) +
    scale_color_manual(name = "Row Crop Proportion", breaks = curve_breaks,
                       values = curve_colors, labels = curve_labels) +
    labs(x = indicator_labels[[ind_name]], y = "Distribution-based score", title = "") +
    theme_bw() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.margin = margin(t = 6, r = 14, b = 8, l = 6))
}

fig6 <- ggarrange(plotlist = plot_list[c("som", "active_c", "respiration", "pmn", "agg_stab")],
                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "right", labels = "AUTO")

fig6

ggsave('final/figs/fig6_design_scoring_curves.png', fig6, height = 6, width = 11, dpi = 300)

# Table 6. Score given under each of the four curves

rowcrop_means <- medium_df %>%
  filter(shmgtcat == 'shms') %>%
  summarise(som = median(som, na.rm = TRUE),
            active_c = median(pox_c, na.rm = TRUE),
            agg_stab = median(shap_as, na.rm = TRUE),
            respiration = median(cminp_96hr_raw, na.rm = TRUE),
            pmn = median(shap_pmn, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Indicator", values_to = "Mean_Value")

table6 <- tibble()
for (ind_name in rowcrop_means$Indicator) {
  val <- rowcrop_means$Mean_Value[rowcrop_means$Indicator == ind_name]

  for (cv in curve_breaks) {
    curve_i <- curve_tbl %>% filter(Indicator == ind_name, Curve == cv) %>% arrange(Indicator_Value)

    table6 <- bind_rows(table6, tibble(
      Indicator = ind_name,
      Curve = cv,
      Mean_Value = val,
      Score = approx(curve_i$Indicator_Value, curve_i$Score, xout = val, rule = 2)$y
    ))
  }
}

# report the SOM mean as SOC, and back-transform AS-Yoder and PMN
table6 <- table6 %>%
  mutate(Mean_Value = case_when(Indicator == "som" ~ Mean_Value / som_factor,
                                Indicator == "agg_stab" ~ Mean_Value^(1 / 5),
                                Indicator == "pmn" ~ Mean_Value^2,
                                .default = Mean_Value),
         Indicator = dplyr::recode(Indicator, som = "SOC", active_c = "POXC",
                                   agg_stab = "AS-Yoder", respiration = "PMC", pmn = "PMN"),
         Curve = factor(Curve, levels = c("20% Row Crop", "57% Row Crop",
                                          "75% Row Crop", "80% Row Crop"))) %>%
  mutate(Mean_Value = round(Mean_Value, 1), Score = round(Score)) %>%
  pivot_wider(names_from = Curve, values_from = Score) %>%
  dplyr::rename(`Mean Row Crop Value` = Mean_Value)

table6

write.csv(table6, 'final/outputs/table6_design_scores.csv', row.names = FALSE)

# Figure S2. Mapped vs laboratory-measured particle size class

library(ggalluvial)

dat_db # REDACTED due to privacy concerns

dat_db <- dat_db %>%
  mutate(new_psc = if_else(is.na(particlesizeclass_sieve), particlesizeclass_hyd, particlesizeclass_sieve)) %>%
  filter(depthclass == '0-15') %>%
  dplyr::select(slocid, new_psc)

mapped #redacted
legend #redacted

dat_psc <- dat_r %>%
  left_join(mapped, by = 'slocid') %>%
  left_join(dat_db, by = 'slocid') %>%
  mutate(shg = case_when(substr(rast_code, 1, 1) == 1 ~ 'clayey',
                         substr(rast_code, 1, 1) == 2 ~ 'fine-silty',
                         substr(rast_code, 1, 1) == 3 ~ 'fine-loamy',
                         substr(rast_code, 1, 1) == 4 ~ 'coarse-silty',
                         substr(rast_code, 1, 1) == 5 ~ 'coarse-loamy',
                         substr(rast_code, 1, 1) == 6 ~ 'sandy',
                         .default = "other"),
         psc_grouped = factor(if_else(new_psc %in% c("fine-loamy", "fine-silty", "clayey", "coarse-loamy", "coarse-silty", "sandy"), new_psc, "other")),
         shg_grouped = factor(if_else(shg %in% c("fine-loamy", "fine-silty", "clayey", "coarse-loamy", "coarse-silty", "sandy"), shg, "other")))

p_psc <- dat_psc %>%
  ggplot(aes(axis1 = shg_grouped, axis2 = psc_grouped)) +
  geom_alluvium(aes(fill = shg_grouped), stat = 'alluvium', alpha = 0.7) +
  geom_stratum(aes(fill = after_stat(stratum)), stat = 'stratum') +
  annotate("text", x = 1, y = -10, label = "Mapped", size = 4) +
  annotate("text", x = 2, y = -10, label = "Lab", size = 4) +
  theme_minimal() +
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text.x = element_blank()) +
  ylab('') + xlab('') +
  scale_fill_manual(name = 'Particle size class',
                    values = c(legend$hexcode[2], legend$hexcode[6], legend$hexcode[4],
                               legend$hexcode[3], legend$hexcode[1], legend$hexcode[7], legend$hexcode[5]))

p_psc

ggsave('final/figs/figS2_psc_mapped_vs_lab.png', p_psc, width = 7, height = 5, dpi = 300)

# percentage of samples whose mapped class matched the lab class
dat_psc %>%
  filter(!is.na(shg), !is.na(new_psc)) %>%
  count(psc_match = if_else(shg == new_psc, "same", "different")) %>%
  mutate(pct = n / sum(n))

# Figure 5 (right). Texture triangle with regions A and B

library(ggtern)

data(USDA)

USDA <- USDA %>%
  mutate(tex_group = case_when(
    Label %in% c("Sand", "Loamy Sand", "Sandy Loam") ~ "S/LS/SL",
    Label %in% c("Sandy Clay Loam", "Loam", "Silt Loam", "Silt") ~ "SCL/L/SiL/Si",
    Label %in% c("Heavy Clay", "Clay", "Silty Clay", "Sandy Clay", "Clay Loam", "Silty Clay Loam") ~ "HC/C/SiC/SC/CL/SiCL",
    .default = NA_character_))

# stop if this ggtern version uses different class labels
if (any(is.na(USDA$tex_group))) {
  stop("Unassigned USDA texture labels: ",
       paste(unique(USDA$Label[is.na(USDA$tex_group)]), collapse = ", "))
}

pal_groups <- c("S/LS/SL" = "#F28E2B",
                "SCL/L/SiL/Si" = "#EDC948",
                "HC/C/SiC/SC/CL/SiCL" = "#76B7B2")

# label only loam, clay loam and sandy loam, at the tile midpoints
USDA_lab <- USDA %>%
  group_by(Label) %>%
  summarise(Sand = mean(Sand), Silt = mean(Silt), Clay = mean(Clay), .groups = "drop") %>%
  filter(Label %in% c("Loam", "Clay Loam", "Sandy Loam")) %>%
  mutate(short = case_when(Label == "Loam" ~ "L",
                           Label == "Clay Loam" ~ "CL",
                           Label == "Sandy Loam" ~ "SL")) %>%
  as.data.frame()

regions_poly <- bind_rows(
  data.frame(region = "A", Sand = c(49, 55, 55, 49), Clay = c(12, 12, 18, 18)),
  data.frame(region = "B",
             Sand = c(min(dat_regB$newsand_0_15), max(dat_regB$newsand_0_15),
                      max(dat_regB$newsand_0_15), min(dat_regB$newsand_0_15)),
             Clay = c(25, 25, 31, 31))
) %>%
  mutate(Silt = 100 - Sand - Clay)

regions_centers <- regions_poly %>%
  group_by(region) %>%
  summarise(Sand = mean(Sand), Clay = mean(Clay), Silt = mean(Silt), .groups = "drop") %>%
  as.data.frame()

p_tern <- ggtern(data = USDA, aes(Sand, Clay, Silt)) +
  geom_polygon(aes(group = Label, fill = tex_group), alpha = 0.75, linewidth = 0.5, color = "black") +
  scale_fill_manual(values = pal_groups, guide = "none") +
  geom_mask() +
  theme_bw() +
  geom_point(data = dat_crop %>% mutate(Sand = newsand_0_15, Clay = clay_0_15, Silt = 100 - Sand - Clay),
             aes(Sand, Clay, Silt), inherit.aes = FALSE, size = 2, color = "black", alpha = 0.2) +
  geom_polygon(data = regions_poly, aes(Sand, Clay, Silt, group = region), inherit.aes = FALSE,
               fill = "grey10", alpha = 0.18, color = "grey10", linewidth = 0.6) +
  geom_point(data = dat_regions %>% mutate(Sand = newsand_0_15, Clay = clay_0_15, Silt = 100 - Sand - Clay),
             aes(Sand, Clay, Silt), inherit.aes = FALSE, size = 2, alpha = 0.85) +
  geom_text(data = regions_centers, aes(Sand, Clay, Silt, label = region), inherit.aes = FALSE, position = "identity",
            size = 6, fontface = "bold", color = "white") +
  geom_text(data = USDA_lab, aes(Sand, Clay, Silt, label = short), inherit.aes = FALSE, position = "identity",
            fontface = "bold", color = "black", size = 5) +
  theme_showarrows() +
  weight_percent() +
  theme_legend_position("topleft") +
  labs(title = "", yarrow = "Clay (%)", zarrow = "Silt (%)", xarrow = "Sand (%)") +
  theme_hidearrows() +
  theme_clockwise() +
  theme(tern.panel.grid.major = element_blank(),
        tern.panel.grid.minor = element_blank(),
        tern.axis.title = element_text(size = 16))

p_tern

ggsave('final/figs/fig5_right_texture_triangle.png', p_tern, width = 5, height = 5, dpi = 300)
