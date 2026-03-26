#### Bioassays - Chemistry Boxplots & Texture ####

#### Textural Triangle ####

## Tutorial from: https://saryace.github.io/flipbook_soiltexture_en/#36, https://stackoverflow.com/questions/29136168/add-texture-classes-to-soil-classification-triangle-via-ggplot2 
## Confirmed with: https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/?cid=nrcs142p2_054167

setwd()

library(ggtern)
data(USDA)
head(USDA, 10)


library(dplyr)
USDA_text <- USDA %>% group_by(Label) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
USDA_text


soil_data <- read.delim("SoilTexture_forR_TernaryPlot.txt", header = T, check.names=FALSE)
summary(soil_data)


# Set factors and desired order
soil_data$Treatment_LandUseHistory <- factor(soil_data$Treatment_LandUseHistory, levels = c(
  "Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut"
))

Soil_Textural_Triangle <- ggplot(data = USDA, aes(
  y = Clay,
  x = Sand,
  z = Silt
)) +
  coord_tern(L = "x", T = "y", R = "z") +
  geom_polygon(
    aes(fill = Label),
    alpha = 0.0,
    size = 0.5,
    color = "black"
  ) +
  geom_text(data = USDA_text,
            aes(label = Label),
            color = 'black',
            size = 3) +
  theme_showarrows() +
  theme_clockwise() +
  guides(fill=FALSE, color=FALSE) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(size= 10)) +
  theme(axis.title = element_text(size = 15)) +
  theme_nomask() +
  theme(legend.position="right") 


Soil_Textural_Triangle <- Soil_Textural_Triangle + 
  geom_point(data = soil_data, size=5,
             aes(x = Sand, y = Clay, z = Silt)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + 
  theme(legend.position = "right")


Soil_Textural_Triangle


#### Chemistry Boxplots ####

setwd()

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(forcats)

# Read the .txt file
data <- read_tsv("Bioassays_Chemistry_Boxplots.txt")

# Ensure column names are proper
colnames(data)[1:6] <- c("SampleID", "RepNo", "Treatment", "Moisture", "Tree_or_NoTree", "MPB_Impacted")

# Define treatment levels and colors
treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")

# List of response variables to plot
response_vars <- c("pH", "Total_C_percent", "Total_N_percent", "Na_waterextract_mgperL", 
                   "NH4_waterextract_mgperL", "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                   "Ca_waterextract_mgperL", "Cl_waterextract_mgperL", "NO3_waterextract_mgperL",
                   "PO4_waterextract_mgperL", "SO4_waterextract_mgperL")


data$MPB_Impacted = factor(data$MPB_Impacted, levels = c("Not_MPB_Impacted", "MPB_Impacted"))
treatment_levels_MPB <- c("Not_MPB_Impacted", "MPB_Impacted")



# Function to generate each plot
generate_plot <- function(response) {
  p <- ggplot(data, aes(x = MPB_Impacted, y = .data[[response]], fill = MPB_Impacted)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(y = .data[[response]], color = MPB_Impacted), 
                position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8) +
    theme_minimal(base_size = 12) +
    labs(title = response, x = "Treatment", y = response) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p <- p + stat_compare_means(
    comparisons = combn(treatment_levels_MPB, 2, simplify = FALSE),
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3,
    aes(group = MPB_Impacted) 
  )
  
  return(p)
}


# Generate all plots
plot_list <- lapply(response_vars, generate_plot)

# Arrange in grid with 4 on each row
final_plot <- ggarrange(
  plotlist = plot_list,
  ncol = 4,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom"
)

# Print final plot
print(final_plot)

ggsave("boxplot_chemistry_stats_summary_MPBComparisons_withpvaluecorrections.pdf", final_plot, width = 7, height = 15, units = "in")







