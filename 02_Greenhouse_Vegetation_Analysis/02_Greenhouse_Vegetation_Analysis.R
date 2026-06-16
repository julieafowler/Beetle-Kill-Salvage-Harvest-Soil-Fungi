#### SFSP - Vegetation Analysis - Greenhouse ####

setwd()

#### Emergence by Treatment ####

# ============================================================
# Binomial GLM: Seedling Emergence by Land Use History
# ============================================================
# Analyzes emergence proportions (Number_Emerged / Number_Total_Containers)
# across six treatment groups using a binomial GLM with logit link.
# Includes planned contrasts, back-transformed emmeans, and visualization.

## Made with the help of Claude Sonnet 4.6
# ============================================================

library(tidyverse)
library(emmeans)

setwd()

# ── 1. Load data ─────────────────────────────────────────────

dat <- read_tsv("Emergence_Barplot.txt")

# Set factor level order (logical: Old_Growth → Pre_MPB → Post_MPB → Recent_Cut)
dat <- dat |>
  mutate(
    Treatment_LandUse = factor(
      Treatment_LandUse,
      levels = c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB",
                 "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut")
    ),
    Salvage_Harvest_Status = factor(
      Salvage_Harvest_Status,
      levels = c("Non_Salvage_Harvested", "Salvage_Harvested")
    )
  )


# ── 2. Fit the saturated binomial GLM ────────────────────────
# Treatment_LandUse is the only usable predictor given n=6 groups.
# Salvage_Harvested is perfectly collinear with Treatment_LandUse,
# so it is used for contrast definitions only (see Section 4).

m <- glm(
  cbind(Number_Emerged, Number_Total_Containers - Number_Emerged) ~ Treatment_LandUse,
  data   = dat,
  family = binomial(link = "logit")
)

summary(m)
# Note: with 6 groups and one obs per group, the model is saturated —
# it perfectly reproduces observed proportions. Formal omnibus tests
# are uninformative; inference proceeds through planned contrasts below.


# ── 3. Check for overdispersion ──────────────────────────────
# Pearson dispersion statistic: values >> 1 suggest overdispersion.
# With a saturated model this will equal 0 (perfect fit), but check
# after any model simplification.

dispersion <- sum(residuals(m, type = "pearson")^2) / m$df.residual
cat("Pearson dispersion statistic:", round(dispersion, 3), "\n")
# If > ~1.5 after adding structure, refit with family = quasibinomial


# ── 4. Estimated marginal means (back-transformed to probability scale) ──

em <- emmeans(m, ~ Treatment_LandUse, type = "response")
# type = "response" back-transforms from log-odds to probability scale
print(em)


# ── 5. Pairwise comparisons (exploratory / supplementary) ────
# All pairwise differences with Tukey correction.
# More conservative than planned contrasts; treat as supplementary.

pairwise_results <- pairs(em, adjust = "tukey")
summary(pairwise_results, infer = TRUE)


# ── 6. Visualization ─────────────────────────────────────────

# --- 6a. Dot-and-whisker plot: model-estimated emergence probabilities ---

em_df <- as.data.frame(em) |>
  rename(prob = prob, LCL = asymp.LCL, UCL = asymp.UCL)

# Define fill color by salvage status for visual grouping
em_df <- em_df |>
  mutate(
    Salvage = ifelse(
      Treatment_LandUse %in% c("1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut"),
      "Salvage Harvested", "Not Salvage Harvested"
    )
  )

p1 <- ggplot(em_df,
             aes(x = Treatment_LandUse, y = prob,
                 ymin = LCL, ymax = UCL, color = Salvage)) +
  geom_pointrange(size = 0.8, linewidth = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_color_manual(values = c("Not Salvage Harvested" = "#4E79A7",
                                "Salvage Harvested"     = "#F28E2B")) +
  labs(
    title    = "Estimated Seedling Emergence Probability by Land Use History",
    subtitle = "Binomial GLM estimates with 95% confidence intervals",
    x        = "Land Use History",
    y        = "Emergence Probability",
    color    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "top",
    axis.text.x       = element_text(angle = 30, hjust = 1),
    panel.grid.minor  = element_blank()
  )

print(p1)
ggsave("emergence_emmeans.pdf", p1, width = 7, height = 5, dpi = 300)


# --- 6b. Observed proportions barplot with model CI overlay ---
# Combines raw data (bars) with GLM uncertainty (error bars) in one plot.

obs_df <- dat |>
  select(Treatment_LandUse, Percent_Emerged, Salvage_Harvest_Status) |>
  mutate(
    Salvage = ifelse(Salvage_Harvested == "Salvage_Harvested",
                     "Salvage Harvested", "Not Salvage Harvested")
  )

p2 <- ggplot() +
  geom_col(data = obs_df,
           aes(x = Treatment_LandUse, y = Percent_Emerged, fill = Salvage),
           alpha = 0.6, width = 0.6) +
  geom_pointrange(data = em_df,
                  aes(x = Treatment_LandUse, y = prob, ymin = LCL, ymax = UCL),
                  color = "black", size = 0.5, linewidth = 0.7) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Not Salvage Harvested" = "#4E79A7",
                               "Salvage Harvested"     = "#F28E2B")) +
  labs(
    title    = "Seedling Emergence: Observed Proportions with Model Estimates",
    subtitle = "Bars = observed; points + whiskers = GLM-estimated probability ± 95% CI",
    x        = "Land Use History",
    y        = "Emergence",
    fill     = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "top",
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p2)
ggsave("emergence_observed_vs_model.pdf", p2, width = 7, height = 5, dpi = 300)



# ============================================================
# Binomial GLM: Seedling Emergence — Salvage vs. Non-Salvage
# ============================================================
# Collapses the six treatment groups into two categories by
# pooling raw container counts, then fits a binomial GLM to
# test whether salvage harvest status predicts emergence.
# Intended as a targeted complement to the six-group analysis.

## Made with the help of Claude Sonnet 4.6
# ============================================================

library(tidyverse)
library(emmeans)

setwd()

# ── 1. Load and collapse data ─────────────────────────────────
# Pool raw counts across treatments within each salvage category.
# This preserves the full n (~91 containers per group) rather than
# averaging proportions, which would give only n=2 data points.

dat <- read_tsv("Emergence_Barplot.txt") |>
  mutate(
    Salvage_Harvest_Status = factor(
      Salvage_Harvest_Status,
      levels = c("Non_Salvage_Harvested", "Salvage_Harvested")
    )
  )

dat_collapsed <- dat |>
  group_by(Salvage_Harvest_Status) |>
  summarise(
    Total   = sum(Number_Total_Containers),
    Emerged = sum(Number_Emerged),
    .groups = "drop"
  )

print(dat_collapsed)
# Non_Salvage_Harvested: n containers, n emerged
# Salvage_Harvested:     n containers, n emerged


# ── 2. Fit the binomial GLM ───────────────────────────────────

m_col <- glm(
  cbind(Emerged, Total - Emerged) ~ Salvage_Harvest_Status,
  data   = dat_collapsed,
  family = binomial(link = "logit")
)

summary(m_col)


# ── 3. Check overdispersion ───────────────────────────────────
# Now meaningful — model is no longer saturated (1 df residual).

dispersion_col <- sum(residuals(m_col, type = "pearson")^2) / m_col$df.residual
cat("Pearson dispersion statistic:", round(dispersion_col, 3), "\n")
# If > ~1.5, refit with family = quasibinomial


# ── 4. Estimated marginal means (probability scale) ───────────

em_col <- emmeans(m_col, ~ Salvage_Harvest_Status, type = "response")
print(em_col)


# ── 5. Contrast: salvage vs. non-salvage ─────────────────────

contrast_col <- contrast(em_col, list(
  salvage_vs_not = c(-1, 1)   # Salvage_Harvested minus Non_Salvage_Harvested
), adjust = "none")

summary(contrast_col, infer = TRUE)


# ── 6. Visualization ──────────────────────────────────────────

em_col_df <- as.data.frame(em_col) |>
  rename(LCL = asymp.LCL, UCL = asymp.UCL) |>
  mutate(
    label = ifelse(Salvage_Harvest_Status == "Salvage_Harvested",
                   "Salvage\nHarvested", "Not Salvage\nHarvested")
  )

# --- 6a. Dot-and-whisker ---

p_col1 <- ggplot(em_col_df,
                 aes(x = Salvage_Harvest_Status, y = prob,
                     ymin = LCL, ymax = UCL,
                     color = Salvage_Harvest_Status)) +
  geom_pointrange(size = 1, linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "gray50", linewidth = 0.4) +
  scale_x_discrete(labels = c("Non_Salvage_Harvested" = "Not Salvage\nHarvested",
                              "Salvage_Harvested"     = "Salvage\nHarvested")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_color_manual(values = c("Non_Salvage_Harvested" = "#4E79A7",
                                "Salvage_Harvested"     = "#F28E2B"),
                     guide = "none") +
  labs(
    title    = "Estimated Emergence Probability by Salvage Harvest Status",
    subtitle = "Binomial GLM estimates with 95% confidence intervals\n(pooled container counts within each category)",
    x        = NULL,
    y        = "Emergence Probability"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_col1)
ggsave("emergence_collapsed_emmeans.pdf", p_col1, width = 2, height = 5, dpi = 300)


# --- 6b. Barplot with CI overlay ---

obs_col_df <- dat_collapsed |>
  mutate(
    Percent_Emerged   = Emerged / Total,
    Salvage_Harvest_Status = factor(Salvage_Harvest_Status,
                               levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
  )

p_col2 <- ggplot() +
  geom_col(data = obs_col_df,
           aes(x = Salvage_Harvest_Status, y = Percent_Emerged,
               fill = Salvage_Harvest_Status),
           alpha = 0.6, width = 0.5) +
  geom_pointrange(data = em_col_df,
                  aes(x = Salvage_Harvest_Status, y = prob,
                      ymin = LCL, ymax = UCL),
                  color = "black", size = 0.6, linewidth = 0.8) +
  scale_x_discrete(labels = c("Non_Salvage_Harvested" = "Not Salvage\nHarvested",
                              "Salvage_Harvested"     = "Salvage\nHarvested")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Non_Salvage_Harvested" = "#4E79A7",
                               "Salvage_Harvested"     = "#F28E2B"),
                    guide = "none") +
  labs(
    title    = "Seedling Emergence by Salvage Harvest Status",
    subtitle = "Bars = observed (pooled); points + whiskers = GLM estimate ± 95% CI",
    x        = NULL,
    y        = "Emergence"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_col2)
ggsave("emergence_collapsed_observed_vs_model.pdf", p_col2, width = 4, height = 5, dpi = 300)


# ── 7. Summary table ──────────────────────────────────────────

contrast_table_col <- summary(contrast_col, infer = TRUE) |>
  as.data.frame() |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

print(contrast_table_col)



#### Vegetation Boxplots ####

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(forcats)

# Read the .txt file
data <- read_tsv("Vegetation_Boxplots.txt")  # Replace with your actual file path

# Ensure column names are proper
colnames(data)[1:4] <- c("SampleID", "Treatment", "Moisture", "Salvage_Harvest_Status")

# Define treatment levels and colors
treatment_levels <- c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut")

# List of response variables to plot
response_vars <- c("Height", "MeanCanWidth", "RootDryMass", "Aboveground_Biomass", "RoottoShootRatio")

data_filtered$Salvage_Harvest_Status = factor(data_filtered$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
treatment_levels_SH <- c("Non_Salvage_Harvested", "Salvage_Harvested")

data_filtered$Treatment = factor(data_filtered$Treatment, levels = c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut"))
treatment_levels <- c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut")


# Function to generate each plot
generate_plot <- function(response) {
  p <- ggplot(data_filtered, aes(x = Salvage_Harvest_Status, y = .data[[response]], fill = Salvage_Harvest_Status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Salvage_Harvest_Status), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = response, x = "Treatment", y = response) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add significance test between treatments
  comparisons <- combn(treatment_levels_SH, 2, simplify = FALSE)
  p <- p + stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    aes(group = Treatment),
    bracket.size = 0.5,
    size = 3
  )
  
  return(p)
}

# Generate all plots
plot_list <- lapply(response_vars, generate_plot)

# Arrange in grid with 3 on top row, 2 on second row
final_plot <- ggarrange(
  plotlist = plot_list,
  ncol = 5,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

# Print final plot
print(final_plot)

ggsave("boxplot_vegetation_wtests_SH_UPDATED117125.pdf", final_plot, width = 8, height = 5, units = "in")




#### Root Tips - Raw Counts and Metric (Per Dry Root Mass) ####

RootTips_RawandMetric <- read_tsv("RootTipStaining_Results_forBoxplot.txt")

RootTips_RawandMetric$Treatment_LandUseHistory <- factor(RootTips_RawandMetric$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
RootTips_RawandMetric$Salvage_Harvest_Status_Status <- factor(RootTips_RawandMetric$Salvage_Harvest_Status_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_RawandMetric$Salvage_Harvest_Status_Status = factor(RootTips_RawandMetric$Salvage_Harvest_Status_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
treatment_levels_SH <- c("Non_Salvage_Harvested", "Salvage_Harvested")
comparisons <- combn(treatment_levels_SH, 2, simplify = FALSE)


RootTips_Boxplot_RawCounts <- ggplot(RootTips_RawandMetric, aes(x = Salvage_Harvest_Status_Status, y = Colonized_RootTips_Counts)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons, 
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) + 
  labs(title = "ECM Colonized Root Tips - Raw Counts",
       y = "Colonized Root Tips (Counts)") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_RawCounts

ggsave("RootTips_Boxplot_RawCounts_bySH.pdf", RootTips_Boxplot_RawCounts, width = 2, height = 6, units = "in")


RootTips_Boxplot_RawCounts_notests <- ggplot(RootTips_RawandMetric, aes(x = Treatment_LandUseHistory, y = Colonized_RootTips_Counts)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) +
  labs(title = "ECM Colonized Root Tips - Raw Counts",
       y = "Colonized Root Tips (Counts)") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_RawCounts_notests

ggsave("RootTips_Boxplot_RawCounts_notests.pdf", RootTips_Boxplot_RawCounts_notests, width = 4, height = 6, units = "in")



treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_RawandMetric$Salvage_Harvest_Status_Status = factor(RootTips_RawandMetric$Salvage_Harvest_Status_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
treatment_levels_SH <- c("Non_Salvage_Harvested", "Salvage_Harvested")
comparisons <- combn(treatment_levels_SH, 2, simplify = FALSE)

RootTips_Boxplot_Metric <- ggplot(RootTips_RawandMetric_filtered, aes(x = Salvage_Harvest_Status_Status, y = Metric_ColonizedTipsper_RootDryMass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) + 
  labs(title = "ECM Colonized Root Tips - Colonized Root Tips per Average Dry Root Mass",
       y = "Colonized Root Tips per Average Dry Root Mass") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_Metric

ggsave("RootTips_Boxplot_Metric_bySH.pdf", RootTips_Boxplot_Metric, width = 2, height = 6, units = "in")


RootTips_Boxplot_Metric_notests <- ggplot(RootTips_RawandMetric, aes(x = Treatment_LandUseHistory, y = Metric_ColonizedTipsper_RootDryMass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) +
  labs(title = "ECM Colonized Root Tips - Colonized Root Tips per Average Dry Root Mass",
       y = "Colonized Root Tips per Average Dry Root Mass") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_Metric_notests

ggsave("RootTips_Boxplot_Metric_notests.pdf", RootTips_Boxplot_Metric_notests, width = 4, height = 6, units = "in")




treatment_levels <- c("Salvage_Harvested", "Non_Salvage_Harvested")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_Boxplot_Metric_SHComparisons <- ggplot(RootTips_RawandMetric, aes(x = Salvage_Harvest_Status_Status, y = Metric_ColonizedTipsper_RootDryMass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) + 
  labs(title = "ECM Colonized Root Tips - Colonized Root Tips per Average Dry Root Mass",
       y = "Colonized Root Tips per Average Dry Root Mass") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_Metric_SHComparisons

ggsave("RootTips_Boxplot_Metric_SHComparisons.pdf", RootTips_Boxplot_Metric_SHComparisons, width = 14, height = 10, units = "in")



treatment_levels <- c("Salvage_Harvested", "Non_Salvage_Harvested")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_Boxplot_RawCounts_SHComparisons <- ggplot(RootTips_RawandMetric, aes(x = Salvage_Harvest_Status_Status, y = Colonized_RootTips_Counts)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.75) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test", ## Non-parametric
                     p.adjust.method = "BH", ## Adjustment for multiple tests
                     label = "p.signif",
                     hide.ns = TRUE) + 
  labs(title = "ECM Colonized Root Tips - Colonized Root Tips per Average Dry Root Mass",
       y = "Colonized Root Tips per Average Dry Root Mass") +
  theme_minimal() 
#theme(legend.position = "none")

RootTips_Boxplot_RawCounts_SHComparisons

ggsave("RootTips_Boxplot_RawCounts_SHComparisons.pdf", RootTips_Boxplot_RawCounts_SHComparisons, width = 14, height = 10, units = "in")



