#### SFSP - Vegetation Analysis - Field ####

library(tidyverse)
library(readr)
library(ggpubr)

setwd()


#### Regen Data - LPP Only - Regeneration Data by HT class ####

df <- read.delim("FieldVegetation_LPP_Simplified_RegenData.txt")

long_df <- df %>%
  pivot_longer(
    cols = starts_with("LPP"),
    names_to = "ResponseVar",
    values_to = "Value"
  )

long_df$Treatment_Field <- factor(long_df$Treatment_Field, levels = c(
  "PreMPB", "EarlyMPB", "LateMPB", "RecentCut"
))

long_df$ResponseVar <- factor(long_df$ResponseVar, levels = c(
  "LPP_HTClass1", "LPP_HTClass2", "LPP_HtClass3"   
))

field_LPP <- ggplot(long_df, aes(x = Treatment_Field, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +  # avoids double-plotting outliers
  geom_jitter(aes(shape = as.factor(Transect)), 
              width = 0.2, alpha = 0.7, size = 3) +  # jitter points by unit/transect
  facet_wrap(~ ResponseVar, scales = "free_y") +  # one plot per response variable
  theme_bw() +
  labs(shape = "Transect") +
  theme(strip.text = element_text(face = "bold"))

field_LPP

treatment_levels <- c("PreMPB", "EarlyMPB", "LateMPB", "RecentCut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)
field_LPP <- field_LPP + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment_Field),
  bracket.size = 0.5,
  size = 3
)

field_LPP

ggsave("Field_LPP_Boxplots_Regen.pdf", field_LPP, width = 10, height = 7, units = "in")



#### Regen Data - All Tree Types ####

df_all <- read.delim("FieldVegetation_ALLTREES_Simplified_RegenData.txt")

long_df_all <- df_all %>%
  pivot_longer(
    cols = starts_with("Trees"),
    names_to = "ResponseVar",
    values_to = "Value"
  )

long_df_all$Treatment_Field <- factor(long_df_all$Treatment_Field, levels = c(
  "PreMPB", "EarlyMPB", "LateMPB", "RecentCut"
))

long_df_all$ResponseVar <- factor(long_df_all$ResponseVar, levels = c(
  "Trees_LPP_HTClass1", "Trees_LPP_HTClass2", "Trees_LPP_HtClass3", "Trees_ASP_HTClass1", "Trees_ASP_HTClass2", "Trees_ASP_HtClass3",	"Trees_SAF_HTClass1",	"Trees_SAF_HTClass2",	"Trees_SAF_HtClass3",	"Trees_ES_HTClass1", "Trees_ES_HTClass2", "Trees_ES_HtClass3"   
))

field_LPP_all <- ggplot(long_df_all, aes(x = Treatment_Field, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +  # avoids double-plotting outliers
  geom_jitter(aes(shape = as.factor(Transect)), 
              width = 0.2, alpha = 0.7, size = 3) +  # jitter points by unit/transect
  facet_wrap(~ ResponseVar, scales = "free_y", ncol = 3) +  # one plot per response variable
  theme_bw() +
  labs(shape = "Transect") +
  theme(strip.text = element_text(face = "bold"))

field_LPP_all

treatment_levels <- c("PreMPB", "EarlyMPB", "LateMPB", "RecentCut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)
field_LPP_all <- field_LPP_all + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment_Field),
  bracket.size = 1,
  size = 5
)

field_LPP_all

ggsave("Field_ALLTREES_Boxplots_Regen.pdf", field_LPP_all, width = 20, height = 20, units = "in")



#### Substrate Ground Cover ####

library(tidyverse)
library(readxl)
library(ggpubr)

df_substrate <- read_excel("Substrate_GroundCover.xlsx")

long_df_substrate <- df_substrate %>%
  pivot_longer(
    cols = c("Litter Duff", "Soil Gravel", "Rock (>1cm)", "1000 hr fuel", 
             "Moss Lichen", "Woody Basal Stump", "Herb Basal"),
    names_to = "ResponseVar",
    values_to = "Value"
  )

long_df_substrate$Age <- factor(long_df_substrate$Age, levels = c(
  "PMPB", "EMPB", "LMPB", "RC"
))

long_df_substrate$ResponseVar <- factor(long_df_substrate$ResponseVar, levels = c(
  "Litter Duff", "Soil Gravel", "Rock (>1cm)", "1000 hr fuel", 
  "Moss Lichen", "Woody Basal Stump", "Herb Basal"
))

ground_cover_plot <- ggplot(long_df_substrate, aes(x = Age, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray", color = "black") +  # avoids double-plotting outliers
  geom_jitter(aes(shape = as.factor(Trans)), 
              width = 0.2, alpha = 0.7, size = 2) +  # jitter points by transect
  facet_wrap(~ ResponseVar, scales = "free_y", ncol = 4) +  # one plot per response variable
  theme_bw() +
  labs(shape = "Transect", 
       x = "Age", 
       y = "Cover (%)") +
  theme(strip.text = element_text(face = "bold"))

ground_cover_plot

age_levels <- c("PMPB", "EMPB", "LMPB", "RC")
comparisons <- combn(age_levels, 2, simplify = FALSE)

ground_cover_plot <- ground_cover_plot + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Age),
  bracket.size = 1,
  size = 4
)

ground_cover_plot

ggsave("GroundCover_Boxplots.pdf", ground_cover_plot, width = 20, height = 10, units = "in")



#### Understory Cover ####

library(tidyverse)
library(readxl)
library(ggpubr)

df_GFS <- read_excel("GFS_UnderstoryCover.xlsx")

df_GFS <- df_GFS %>%
  mutate(
    Graminoids = as.numeric(Graminoids),
    Forbs = as.numeric(Forbs),
    Shrubs = as.numeric(Shrubs)
  )

long_df_GFS <- df_GFS %>%
  pivot_longer(
    cols = c(Graminoids, Forbs, Shrubs),
    names_to = "ResponseVar",
    values_to = "Value"
  )

long_df_GFS$Value <- as.numeric(long_df_GFS$Value)

long_df_GFS$Age <- factor(long_df_GFS$Age, levels = c(
  "PMPB", "EMPB", "LMPB", "RC"
))

long_df_GFS$ResponseVar <- factor(long_df_GFS$ResponseVar, levels = c(
  "Graminoids", "Forbs", "Shrubs"
))

cat("Structure of long_df_GFS:\n")
str(long_df_GFS)
cat("\nValue column class:", class(long_df_GFS$Value), "\n")
cat("First few values:", head(long_df_GFS$Value, 10), "\n")

understory_cover_plot <- ggplot(long_df_GFS, aes(x = Age, y = Value, group = Age)) +
  geom_boxplot(fill = "lightgray", color = "black", outlier.shape = NA) +  
  geom_jitter(aes(shape = as.factor(Trans)), 
              width = 0.2, alpha = 0.7, size = 2) +  
  facet_wrap(~ ResponseVar, scales = "free_y", ncol = 3) +  
  theme_bw() +
  labs(shape = "Transect", 
       x = "Age", 
       y = "Count") +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

understory_cover_plot

age_levels <- c("PMPB", "EMPB", "LMPB", "RC")
comparisons <- combn(age_levels, 2, simplify = FALSE)

understory_cover_plot <- understory_cover_plot + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  bracket.size = 0.5,
  size = 4,
  step.increase = 0.08
)

understory_cover_plot

ggsave("UnderstoryCover_Boxplots.pdf", understory_cover_plot, width = 15, height = 8, units = "in")



#### Over/Under - Overstory Data from 10 m radius plots ####

df <- read_tsv("FieldVegetation_LPP_Simplified_OverUnderData.txt")

long_df <- df %>%
  pivot_longer(
    cols = ends_with("_tperha"),
    names_to = "ResponseVar",
    values_to = "Value"
  )

long_df$Treatment_Field <- factor(long_df$Treatment_Field, levels = c(
  "PreMPB", "EarlyMPB", "LateMPB", "RecentCut"
))

long_df$ResponseVar <- factor(long_df$ResponseVar, levels = c(
  "LPP_lessthan10cm_saplings_tperha", "Overstory_Trees_tperha", "Overstory_Stump_tperha", 
  "Overstory_Sum_tperha", "AllSum_tperha"   
))

field_LPP <- ggplot(long_df, aes(x = Treatment_Field, y = Value)) +
  geom_boxplot(fill = "lightgray", color = "black", outlier.shape = NA) +  # avoids double-plotting outliers
  geom_jitter(aes(shape = as.factor(Transect)), 
              width = 0.2, alpha = 0.7) +  # jitter points by unit/transect
  facet_wrap(~ ResponseVar, scales = "free_y") +  # one plot per response variable
  theme_bw() +
  labs(shape = "Transect") +
  theme(strip.text = element_text(face = "bold"))

treatment_levels <- c("PreMPB", "EarlyMPB", "LateMPB", "RecentCut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)
field_LPP <- field_LPP + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment_Field),
  bracket.size = 0.5,
  size = 3
)

field_LPP

ggsave("Field_LPP_Boxplots_OverUnder_wtests.pdf", field_LPP, width = 14, height = 7, units = "in")

