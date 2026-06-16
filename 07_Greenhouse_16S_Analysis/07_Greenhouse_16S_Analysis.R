#### Greenhouse Bioassays - 16S Data ####

setwd()

library(phyloseq) ## if needed to download, will have to get from Bioconductor 
library(ggplot2)
library(grid)
library(plyr)
library(vegan)
library(Hmisc)
library(reshape2)
library(ggpubr)
library(skimr)
library(ggthemr) 


#### Rarefying ####

setwd()
Prokaryote_rhizo <- read.delim("16S_FeatureTable_noMitochondriaChloroplastUnassignedEukaryota_Rhizosphere.txt",header=T,row.names=1, check.names=FALSE)
Prokaryote_rhizo_t <- t(Prokaryote_rhizo)
sample_sums <- rowSums(Prokaryote_rhizo_t)
summary(sample_sums)
rarecurve(Prokaryote_rhizo_t, step = 500, label = FALSE)
which(sample_sums < 18000)
sum(sample_sums < 18000)
Prokaryote_rhizo_t_filtered <- Prokaryote_rhizo_t[rowSums(Prokaryote_rhizo_t) >= 18000, ]
Prokaryote_rhizo_t_rarefied_18000 <- rrarefy(Prokaryote_rhizo_t_filtered, sample = 18000)
Prokaryote_rhizo_t_rarefied_18000 <- t(Prokaryote_rhizo_t_rarefied_18000)



#### Rhizosphere ####

setwd()

## Use rarefied counts
otus <- Prokaryote_rhizo_t_rarefied_18000
otus <- as.data.frame(otus)
sample_names <- colnames(otus)
map_file <- read.delim("16S_metadata_Rhizosphere.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)

all.equal(names(otus),row.names(map_file))

## Ordering variable for ggplot
map_file$Treatment_LandUseHistory <- factor(map_file$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
map_file$Salvage_Harvest_Status <- factor(map_file$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 



#### Alpha Diversity ####

library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)

asv_table <- t(otus)
metadata <- map_file

# Ensure sample IDs match
metadata$SampleID <- rownames(metadata)
asv_table <- as.data.frame(asv_table)
asv_table$SampleID <- rownames(asv_table)

# Merge ASV table and metadata
merged <- inner_join(asv_table, metadata, by = "SampleID")

# Step 1: Save the ASV count columns
counts <- merged %>%
  select(where(is.numeric)) %>% select(-Replicate)

# Step 2: Save metadata columns
meta <- merged %>%
  select(where(~!is.numeric(.)))

# Calculate alpha diversity metrics
alpha_div <- data.frame(
  Observed = rowSums(counts > 0),
  Shannon  = diversity(counts, index = "shannon")
)

# Add metadata back
alpha_div$SampleID <- merged$SampleID
alpha_div <- left_join(alpha_div, metadata, by = "SampleID")


## Details for making figures

alpha_div$Treatment_LandUseHistory <- factor(alpha_div$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
alpha_div$Salvage_Harvest_Status <- factor(alpha_div$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

alpha_div$Treatment_All <- alpha_div$Treatment_LandUseHistory
unique(alpha_div$Treatment_All)

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

alpha_div <- alpha_div %>%
  relocate(Observed, Shannon, .after = last_col())

write.csv(alpha_div, "16S_alpha_div_rhizosphere_ObservedShannons.csv", row.names = FALSE)


## Figures 

## Observed ASVs

AD_Rhizosphere_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Observed)) +
  geom_boxplot(aes(fill = Treatment_All), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_All), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Observed Species by Treatment",
       y = "Observed Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = Treatment_LandUseHistory),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_ObservedSpecies_Rarefied

ggsave("AD_Rhizosphere_ObservedSpecies_Rarefied_notests.pdf", AD_Rhizosphere_ObservedSpecies_Rarefied, width = 3, height = 5, units = "in")



## Shannons 

AD_Rhizosphere_Shannons_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Shannon)) +
  geom_boxplot(aes(fill = Treatment_All), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_All), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Shannon Diversity by Treatment",
       y = "Shannon Index") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = Treatment_LandUseHistory),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_Shannons_Rarefied

ggsave("AD_Rhizosphere_Shannons_Rarefied_notests.pdf", AD_Rhizosphere_Shannons_Rarefied, width = 3, height = 5, units = "in")

## Final plots together 

library(ggpubr)

final_plot <- ggarrange(
  AD_Rhizosphere_ObservedSpecies_Rarefied, AD_Rhizosphere_Shannons_Rarefied,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)

# Print final plot
print(final_plot)

ggsave("alphadiversity_boxplots_Rhizosphere_rarefied_tests.pdf", final_plot, width = 14, height = 5, units = "in")





















