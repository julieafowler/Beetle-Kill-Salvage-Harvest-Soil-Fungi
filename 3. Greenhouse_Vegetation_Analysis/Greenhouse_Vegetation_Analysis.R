#### SFSP - Vegetation Analysis - Greenhouse ####

setwd()

#### Emergence by Treatment ####

# Load required packages
library(ggplot2)
library(dplyr)
library(readr)

# Load the data (adjust file name/path as needed)
df <- read_tsv("Emergence_Barplot.txt")

# Set factors and desired order
df$Treatment_LandUse <- factor(df$Treatment_LandUse, levels = c(
  "Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut"
))


# Plot
ggplot(df_filtered, aes(x = Treatment_LandUse, y = Percent_Emerged)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Land Use Treatment",
    y = "Percent Emerged",
    title = "Percent Emerged by Treatment"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


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
colnames(data)[1:4] <- c("SampleID", "Treatment", "Moisture", "MPB")

# Define treatment levels and colors
treatment_levels <- c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut")

# List of response variables to plot
response_vars <- c("Height", "MeanCanWidth", "RootDryMass", "Aboveground_Biomass", "RoottoShootRatio")

data_filtered$MPB = factor(data_filtered$MPB, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")

data_filtered$Treatment = factor(data_filtered$Treatment, levels = c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut"))
treatment_levels <- c("Old_Growth", "1st_Pre_MPB", "2nd_Pre_MPB", "1st_Post_MPB", "2nd_Post_MPB", "Recent_Cut")


# Function to generate each plot
generate_plot <- function(response) {
  p <- ggplot(data_filtered, aes(x = MPB, y = .data[[response]], fill = MPB)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = MPB), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = response, x = "Treatment", y = response) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add significance test between treatments
  comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)
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

ggsave("boxplot_vegetation_wtests_MPB_UPDATED117125.pdf", final_plot, width = 8, height = 5, units = "in")




#### Root Tips - Raw Counts and Metric (Per Dry Root Mass) ####

RootTips_RawandMetric <- read_tsv("RootTipStaining_Results_forBoxplot.txt")

RootTips_RawandMetric$Treatment_LandUseHistory <- factor(RootTips_RawandMetric$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
RootTips_RawandMetric$MPB_Status <- factor(RootTips_RawandMetric$MPB_Status, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_RawandMetric$MPB_Status = factor(RootTips_RawandMetric$MPB_Status, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)


RootTips_Boxplot_RawCounts <- ggplot(RootTips_RawandMetric, aes(x = MPB_Status, y = Colonized_RootTips_Counts)) +
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

ggsave("RootTips_Boxplot_RawCounts_byMPB.pdf", RootTips_Boxplot_RawCounts, width = 2, height = 6, units = "in")


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

RootTips_RawandMetric$MPB_Status = factor(RootTips_RawandMetric$MPB_Status, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

RootTips_Boxplot_Metric <- ggplot(RootTips_RawandMetric_filtered, aes(x = MPB_Status, y = Metric_ColonizedTipsper_RootDryMass)) +
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

ggsave("RootTips_Boxplot_Metric_byMPB.pdf", RootTips_Boxplot_Metric, width = 2, height = 6, units = "in")


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




treatment_levels <- c("MPB", "Non_MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_Boxplot_Metric_MPBComparisons <- ggplot(RootTips_RawandMetric, aes(x = MPB_Status, y = Metric_ColonizedTipsper_RootDryMass)) +
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

RootTips_Boxplot_Metric_MPBComparisons

ggsave("RootTips_Boxplot_Metric_MPBComparisons.pdf", RootTips_Boxplot_Metric_MPBComparisons, width = 14, height = 10, units = "in")



treatment_levels <- c("MPB", "Non_MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

RootTips_Boxplot_RawCounts_MPBComparisons <- ggplot(RootTips_RawandMetric, aes(x = MPB_Status, y = Colonized_RootTips_Counts)) +
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

RootTips_Boxplot_RawCounts_MPBComparisons

ggsave("RootTips_Boxplot_RawCounts_MPBComparisons.pdf", RootTips_Boxplot_RawCounts_MPBComparisons, width = 14, height = 10, units = "in")



