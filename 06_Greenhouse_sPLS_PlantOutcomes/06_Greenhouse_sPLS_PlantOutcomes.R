#### sPLS with Plant Outcomes - ITS ####

## Code from Laura Moore, Wrighton Lab, https://github.com/Lmoore45/Grazing_Microbiome_2025/tree/main 
## VIP code from someone else, https://github.com/Lmoore45/Grazing_Microbiome_2025/blob/main/6.0_sPLS_VIP.R
## Preprocessing guidelines from the mixOmics website: https://mixomics.org/mixmc/mixmc-preprocessing/


## Looking at “Which small set of microbes best predicts height/biomass/other plant outcomes?”
## Doing a separate run for each of my five outcome variables of interest because of inconsistent NAs 
## Next, following Laura Moore's example, where then I collapse these top ASVs/features (positive relationships and negative relationships seperately) into sum and mean abundances per samples of these top ASVs
## These are then used as predictors in my Bayesian models

set.seed(1234)

#### Rhizosphere ####

setwd()

#### Rhizosphere - Top Feature Selection ####

#### sPLS - Top Feature/ASV Selection - Rhizosphere ####

#### sPLS - Top Feature/ASV Selection - Rhizosphere - Height ####


##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_Rhizosphere.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, Height) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER
# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable
hist(metadata_filtered$Height)
# print(shapiro.test(metadata_filtered$Height))
# response_variable <- log(metadata_filtered$Height)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Height", xlab = "Log(Height)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$Height

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 2, keepX = c(30,10)
## Run2: ncomp = 2, keepX = c(40,10)
## Run3: ncomp = 2, keepX = c(40,10)

ncomp <- 2
keepX <- c(40,10)


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_rhizo_height <- ggplot(df_vip[1:36,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 36 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_rhizo_height

ggsave("TopASVs_rhizo_height_UPDATED112125.pdf", TopASVs_rhizo_height, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
height_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_phylum_plot)

ggsave("curveplot_rhizo_height_phyla_UPDATED112125.pdf", height_phylum_plot, width = 8, height = 7, units = "in")

height_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_genus_plot)

ggsave("curveplot_rhizo_height_genus_UPDATED112125.pdf", height_genus_plot, width = 8, height = 7, units = "in")


height_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_species_plot)

ggsave("curveplot_rhizo_height_species_UPDATED112125.pdf", height_species_plot, width = 8, height = 7, units = "in")


height_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_guild_plot)

ggsave("curveplot_rhizo_height_guild_UPDATED112125.pdf", height_guild_plot, width = 8, height = 7, units = "in")




##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%
  

#transposed_asv_height <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, Height), by = "Sample")
transposed_asv_height <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, Height), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_height %>%
  pivot_longer(cols = -c(Sample, Height), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(Height ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Positive"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_rhizo_height_negativeassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_rhizo_height_negativeassociations_UPDATED112125.csv")






#### sPLS - Top Feature/ASV Selection - Rhizosphere - MeanCanWidth ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_Rhizosphere.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, MeanCanWidth) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)


# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable
hist(metadata_filtered$MeanCanWidth)
# print(shapiro.test(metadata_filtered$MeanCanWidth))
# response_variable <- log(metadata_filtered$MeanCanWidth)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine MeanCanWidth", xlab = "Log(MeanCanWidth)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$MeanCanWidth

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 10
## Run2: ncomp = 1, keepX = 20
## Run3: ncomp = 1, keepX = 10


ncomp <- 1
keepX <- 10


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_rhizo_meancanwidth <- ggplot(df_vip[1:9,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 9 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_rhizo_meancanwidth

ggsave("TopASVs_rhizo_meancanwidth_UPDATED112125.pdf", TopASVs_rhizo_meancanwidth, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
meancanwidth_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_phylum_plot)

ggsave("curveplot_rhizo_meancanwidth_phyla_UPDATED112125.pdf", meancanwidth_phylum_plot, width = 8, height = 7, units = "in")

meancanwidth_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_genus_plot)

ggsave("curveplot_rhizo_meancanwidth_genus_UPDATED112125.pdf", meancanwidth_genus_plot, width = 8, height = 7, units = "in")


meancanwidth_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_species_plot)

ggsave("curveplot_rhizo_meancanwidth_species_UPDATED112125.pdf", meancanwidth_species_plot, width = 8, height = 7, units = "in")


meancanwidth_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_guild_plot)

ggsave("curveplot_rhizo_meancanwidth_guild_UPDATED112125.pdf", meancanwidth_guild_plot, width = 8, height = 7, units = "in")




##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_meancanwidth <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, MeanCanWidth), by = "Sample")
transposed_asv_meancanwidth <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, MeanCanWidth), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_meancanwidth %>%
  pivot_longer(cols = -c(Sample, MeanCanWidth), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(MeanCanWidth ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_rhizo_meancanwidth_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_rhizo_meancanwidth_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - Rhizosphere - RootDryMass ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_Rhizosphere.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, RootDryMass) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$RootDryMass)
# print(shapiro.test(metadata_filtered$RootDryMass))
# response_variable <- log(metadata_filtered$RootDryMass)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Root Dry Mass", xlab = "Log(RootDryMass)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$RootDryMass

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 100
## Run2: ncomp = 1, keepX = 90
## Run3: ncomp = 1, keepX = 90

ncomp <- 1
keepX <- 90


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_rhizo_rootdrymass <- ggplot(df_vip[1:51,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 51 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_rhizo_rootdrymass

ggsave("TopASVs_rhizo_rootdrymass_UPDATED112125.pdf", TopASVs_rhizo_rootdrymass, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
rootdrymass_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_phylum_plot)

ggsave("curveplot_rhizo_rootdrymass_phyla_UPDATED112125.pdf", rootdrymass_phylum_plot, width = 8, height = 7, units = "in")

rootdrymass_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_genus_plot)

ggsave("curveplot_rhizo_rootdrymass_genus_UPDATED112125.pdf", rootdrymass_genus_plot, width = 8, height = 7, units = "in")



rootdrymass_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_species_plot)

ggsave("curveplot_rhizo_rootdrymass_species_UPDATED112125.pdf", rootdrymass_species_plot, width = 8, height = 7, units = "in")


rootdrymass_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_guild_plot)

ggsave("curveplot_rhizo_rootdrymass_guild_UPDATED112125.pdf", rootdrymass_guild_plot, width = 8, height = 7, units = "in")




##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_rootdrymass <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, RootDryMass), by = "Sample")
transposed_asv_rootdrymass <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, RootDryMass), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_rootdrymass %>%
  pivot_longer(cols = -c(Sample, RootDryMass), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(RootDryMass ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_rhizo_rootdrymass_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_rhizo_rootdrymass_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - Rhizosphere - Aboveground_Biomass ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_Rhizosphere.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, Aboveground_Biomass) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$Aboveground_Biomass)
# print(shapiro.test(metadata_filtered$Aboveground_Biomass))
# response_variable <- log(metadata_filtered$Aboveground_Biomass)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Aboveground_Biomass", xlab = "Log(Aboveground_Biomass)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$Aboveground_Biomass

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 40
## Run2: ncomp = 1, keepX = 30
## Run3: ncomp = 1, keepX = 30

ncomp <- 1
keepX <- 30


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_rhizo_aboveground_biomass <- ggplot(df_vip[1:21,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 21 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_rhizo_aboveground_biomass

ggsave("TopASVs_rhizo_aboveground_biomass_UPDATED112125.pdf", TopASVs_rhizo_aboveground_biomass, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
aboveground_biomass_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_phylum_plot)

ggsave("curveplot_rhizo_aboveground_biomass_phyla_UPDATED112125.pdf", aboveground_biomass_phylum_plot, width = 8, height = 7, units = "in")

aboveground_biomass_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_genus_plot)

ggsave("curveplot_rhizo_aboveground_biomass_genus_UPDATED112125.pdf", aboveground_biomass_genus_plot, width = 8, height = 7, units = "in")


aboveground_biomass_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_species_plot)

ggsave("curveplot_rhizo_aboveground_biomass_species_UPDATED112125.pdf", aboveground_biomass_species_plot, width = 8, height = 7, units = "in")


aboveground_biomass_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_guild_plot)

ggsave("curveplot_rhizo_aboveground_biomass_guild_UPDATED112125.pdf", aboveground_biomass_guild_plot, width = 8, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_aboveground_biomass <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, Aboveground_Biomass), by = "Sample")
transposed_asv_aboveground_biomass <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, Aboveground_Biomass), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_aboveground_biomass %>%
  pivot_longer(cols = -c(Sample, Aboveground_Biomass), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(Aboveground_Biomass ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_rhizo_aboveground_biomass_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_rhizo_aboveground_biomass_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - Rhizosphere - RoottoShootRatio ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_Rhizosphere.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, RoottoShootRatio) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$RoottoShootRatio)
# print(shapiro.test(metadata_filtered$RoottoShootRatio))
# response_variable <- log(metadata_filtered$RoottoShootRatio)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine RoottoShootRatio", xlab = "Log(RoottoShootRatio)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$RoottoShootRatio

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 100
## Run2: ncomp = 1, keepX = 120
## Run3: ncomp = 1, keepX = 90

ncomp <- 1
keepX <- 100


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_rhizo_roottoshootratio <- ggplot(df_vip[1:70,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 70 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_rhizo_roottoshootratio

ggsave("TopASVs_rhizo_roottoshootratio_UPDATED112125.pdf", TopASVs_rhizo_roottoshootratio, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
roottoshootratio_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_phylum_plot)

ggsave("curveplot_rhizo_roottoshootratio_phyla_UPDATED112125.pdf", roottoshootratio_phylum_plot, width = 8, height = 7, units = "in")

roottoshootratio_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_genus_plot)

ggsave("curveplot_rhizo_roottoshootratio_genus_UPDATED112125.pdf", roottoshootratio_genus_plot, width = 11, height = 7, units = "in")


roottoshootratio_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_species_plot)

ggsave("curveplot_rhizo_roottoshootratio_species_UPDATED112125.pdf", roottoshootratio_species_plot, width = 11, height = 7, units = "in")


roottoshootratio_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_guild_plot)

ggsave("curveplot_rhizo_roottoshootratio_guild_UPDATED112125.pdf", roottoshootratio_guild_plot, width = 11, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_roottoshootratio <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, RoottoShootRatio), by = "Sample")
transposed_asv_roottoshootratio <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, RoottoShootRatio), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_roottoshootratio %>%
  pivot_longer(cols = -c(Sample, RoottoShootRatio), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(RoottoShootRatio ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_rhizo_roottoshootratio_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_rhizo_roottoshootratio_positiveassociations_UPDATED112125.csv")




#### RootTips ####

set.seed(1234)

setwd()


#### RootTips - Top Feature Selection ####

#### sPLS - Top Feature/ASV Selection - RootTips ####

#### sPLS - Top Feature/ASV Selection - RootTips - Height ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_RootTipsOnly.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, Height) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$Height)
# print(shapiro.test(metadata_filtered$Height))
# response_variable <- log(metadata_filtered$Height)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Height", xlab = "Log(Height)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$Height

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 10
## Run2: ncomp = 1, keepX = 10
## Run3: ncomp = 1, keepX = 10

ncomp <- 1
keepX <- 10


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 5, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_roottips_height <- ggplot(df_vip[1:9,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 9 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_roottips_height

ggsave("TopASVs_roottips_height_UPDATED112125.pdf", TopASVs_roottips_height, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
height_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_phylum_plot)

ggsave("curveplot_roottips_height_phyla_UPDATED112125.pdf", height_phylum_plot, width = 8, height = 7, units = "in")

height_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_genus_plot)

ggsave("curveplot_roottips_height_genus_UPDATED112125.pdf", height_genus_plot, width = 8, height = 7, units = "in")


height_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_species_plot)

ggsave("curveplot_roottips_height_species_UPDATED112125.pdf", height_species_plot, width = 8, height = 7, units = "in")


height_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Height (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(height_guild_plot)

ggsave("curveplot_roottips_height_guild_UPDATED112125.pdf", height_guild_plot, width = 8, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_height <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, Height), by = "Sample")
transposed_asv_height <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, Height), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_height %>%
  pivot_longer(cols = -c(Sample, Height), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(Height ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Positive"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_roottips_height_negativeassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_roottips_height_negativeassociations_UPDATED112125.csv")






#### sPLS - Top Feature/ASV Selection - RootTips - MeanCanWidth ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_RootTipsOnly.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, MeanCanWidth) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$MeanCanWidth)
# print(shapiro.test(metadata_filtered$MeanCanWidth))
# response_variable <- log(metadata_filtered$MeanCanWidth)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine MeanCanWidth", xlab = "Log(MeanCanWidth)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$MeanCanWidth

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 10
## Run2: ncomp = 1, keepX = 10
## Run3: ncomp = 1, keepX = 10

ncomp <- 1
keepX <- 10


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_roottips_meancanwidth <- ggplot(df_vip[1:9,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 9 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_roottips_meancanwidth

ggsave("TopASVs_roottips_meancanwidth_UPDATED112125.pdf", TopASVs_roottips_meancanwidth, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
meancanwidth_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_phylum_plot)

ggsave("curveplot_roottips_meancanwidth_phyla_UPDATED112125.pdf", meancanwidth_phylum_plot, width = 8, height = 7, units = "in")

meancanwidth_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_genus_plot)

ggsave("curveplot_roottips_meancanwidth_genus_UPDATED112125.pdf", meancanwidth_genus_plot, width = 8, height = 7, units = "in")


meancanwidth_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_species_plot)

ggsave("curveplot_roottips_meancanwidth_species_UPDATED112125.pdf", meancanwidth_species_plot, width = 8, height = 7, units = "in")


meancanwidth_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine MeanCanWidth (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(meancanwidth_guild_plot)

ggsave("curveplot_roottips_meancanwidth_guild_UPDATED112125.pdf", meancanwidth_guild_plot, width = 8, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_meancanwidth <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, MeanCanWidth), by = "Sample")
transposed_asv_meancanwidth <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, MeanCanWidth), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_meancanwidth %>%
  pivot_longer(cols = -c(Sample, MeanCanWidth), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(MeanCanWidth ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_roottips_meancanwidth_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_roottips_meancanwidth_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - RootTips - RootDryMass ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_RootTipsOnly.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, RootDryMass) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$RootDryMass)
# print(shapiro.test(metadata_filtered$RootDryMass))
# response_variable <- log(metadata_filtered$RootDryMass)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Root Dry Mass", xlab = "Log(RootDryMass)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$RootDryMass

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 150
## Run2: ncomp = 1, keepX = 10
## Run3: ncomp = 1, keepX = 150

ncomp <- 1
keepX <- 150


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_roottips_rootdrymass <- ggplot(df_vip[1:69,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 69 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_roottips_rootdrymass

ggsave("TopASVs_roottips_rootdrymass_UPDATED112125.pdf", TopASVs_roottips_rootdrymass, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
rootdrymass_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_phylum_plot)

ggsave("curveplot_roottips_rootdrymass_phyla_UPDATED112125.pdf", rootdrymass_phylum_plot, width = 8, height = 7, units = "in")

rootdrymass_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_genus_plot)

ggsave("curveplot_roottips_rootdrymass_genus_UPDATED112125.pdf", rootdrymass_genus_plot, width = 8, height = 7, units = "in")


rootdrymass_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_species_plot)

ggsave("curveplot_roottips_rootdrymass_species_UPDATED112125.pdf", rootdrymass_species_plot, width = 8, height = 7, units = "in")


rootdrymass_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Root Dry Mass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(rootdrymass_guild_plot)

ggsave("curveplot_roottips_rootdrymass_guild_UPDATED112125.pdf", rootdrymass_guild_plot, width = 8, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_rootdrymass <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, RootDryMass), by = "Sample")
transposed_asv_rootdrymass <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, RootDryMass), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_rootdrymass %>%
  pivot_longer(cols = -c(Sample, RootDryMass), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(RootDryMass ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_roottips_rootdrymass_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_roottips_rootdrymass_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - RootTips - Aboveground_Biomass ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_RootTipsOnly.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, Aboveground_Biomass) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$Aboveground_Biomass)
# print(shapiro.test(metadata_filtered$Aboveground_Biomass))
# response_variable <- log(metadata_filtered$Aboveground_Biomass)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine Aboveground_Biomass", xlab = "Log(Aboveground_Biomass)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$Aboveground_Biomass

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 10
## Run2: ncomp = 1, keepX = 10
## Run3: ncomp = 1, keepX = 10

ncomp <- 1
keepX <- 10


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_roottips_aboveground_biomass <- ggplot(df_vip[1:6,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 6 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_roottips_aboveground_biomass

ggsave("TopASVs_roottips_aboveground_biomass_UPDATED112125.pdf", TopASVs_roottips_aboveground_biomass, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
aboveground_biomass_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_phylum_plot)

ggsave("curveplot_roottips_aboveground_biomass_phyla_UPDATED112125.pdf", aboveground_biomass_phylum_plot, width = 8, height = 7, units = "in")

aboveground_biomass_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_genus_plot)

ggsave("curveplot_roottips_aboveground_biomass_genus_UPDATED112125.pdf", aboveground_biomass_genus_plot, width = 8, height = 7, units = "in")


aboveground_biomass_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_species_plot)

ggsave("curveplot_roottips_aboveground_biomass_species_UPDATED112125.pdf", aboveground_biomass_species_plot, width = 8, height = 7, units = "in")


aboveground_biomass_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine Aboveground_Biomass (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(aboveground_biomass_guild_plot)

ggsave("curveplot_roottips_aboveground_biomass_guild_UPDATED112125.pdf", aboveground_biomass_guild_plot, width = 8, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_aboveground_biomass <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, Aboveground_Biomass), by = "Sample")
transposed_asv_aboveground_biomass <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, Aboveground_Biomass), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_aboveground_biomass %>%
  pivot_longer(cols = -c(Sample, Aboveground_Biomass), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(Aboveground_Biomass ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_roottips_aboveground_biomass_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_roottips_aboveground_biomass_positiveassociations_UPDATED112125.csv")




#### sPLS - Top Feature/ASV Selection - RootTips - RoottoShootRatio ####

##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script from Laura Mason to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt")
asv_ITS <- read.delim("ITS_FeatureTable_RootTipsOnly.txt")


##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxonomy) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(SampleID, RoottoShootRatio) %>% drop_na() %>% filter(SampleID %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
# asv_matrix[asv_matrix == 0] <- 0.0001
# asv_clr <- clr(asv_matrix)



# STEP 1: OFFSET
asv_matrix <- asv_matrix+1
sum(which(asv_matrix == 0))
dim(asv_matrix)

# STEP 2: PRE-FILTER

# remove low count OTUs
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent=0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(asv_matrix, percent=0.01)
data.filter <- result.filter$data.filter

# check the number of variables kept after filtering
# but from now on work with 'data.filter'
length(result.filter$keep.otu)

lib.size <- apply(data.filter, 1, sum) # determine total count for each sample
barplot(lib.size) # and plot as bar plot

## Setting a maximum library size to exclude the couple that are way larger than others if needed
# maximum.lib.size <- 60000
# 
# data.filter <- data.filter[-which(lib.size > maximum.lib.size),]

## CLR

asv_matrix <- as.matrix(data.filter)
asv_clr <- clr(data.filter)


# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(SampleID, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable

hist(metadata_filtered$RoottoShootRatio)
# print(shapiro.test(metadata_filtered$RoottoShootRatio))
# response_variable <- log(metadata_filtered$RoottoShootRatio)
# hist(response_variable, main = "Histogram of Log-Transformed Lodgepole Pine RoottoShootRatio", xlab = "Log(RoottoShootRatio)", breaks = 20)
# print(shapiro.test(response_variable))

## Not log transforming 
response_variable <- metadata_filtered$RoottoShootRatio

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation resultsx
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
# ncomp <- 1
# keepX <- rep(50, ncomp)

## Tuning model
grid_keepX <- seq(10, 150, by = 10)

tune_res <- tune.spls(
  X = ASV_table, Y = response_variable,
  ncomp = 5,
  test.keepX = grid_keepX,
  validation = "Mfold", 
  folds = 5, nrepeat = 10,
  measure = "MSE",
  progressBar = TRUE
)

plot(tune_res)

best_ncomp  <- tune_res$choice.ncomp$ncomp
best_keepX <- tune_res$choice.keepX[1:best_ncomp]  # vector per component
best_ncomp; best_keepX

## How to pick: Run tuning a few times and choose the simplest common answer:
## Pick the most frequent ncomp across runs.
## For keepX per component, take the median across runs and round.
## (If two options are basically tied, use the 1-SE rule: pick the smaller model—fewer components / smaller keepX.)

## Run1: ncomp = 1, keepX = 10
## Run2: ncomp = 1, keepX = 10
## Run3: ncomp = 1, keepX = 10

ncomp <- 1
keepX <- 10


spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)
perf_result$measures$Q2.total$summary ## Maximize Q2, minimize MSEP
perf_result$measures$MSEP$summary


##### Step 9: Predicted vs Observed 
# Evaluate model performance visually and quantitatively

# GOOD observed vs predicted: points near the 1:1 line; narrow, symmetric scatter.
# BAD observed vs predicted: systematic bias (above/below the line), curvature, very wide scatter or clear subgroups.

predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

## Also check Residuals vs. Fitted 

# GOOD residuals vs fitted: shapeless cloud around 0, roughly constant vertical spread.
# BAD residuals vs fitted: curves (model misspecification), fan (heteroscedasticity), bands/clusters (batch effects), big outliers.

fitted <- drop(predict(spls_result, ASV_table, ncomp = ncomp)$predict[,1,ncomp])
plot(fitted, response_variable - fitted, xlab = "Fitted", ylab = "Residuals"); abline(h = 0, lty = 2)


##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

TopASVs_roottips_roottoshootratio <- ggplot(df_vip[1:10,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

TopASVs_roottips_roottoshootratio

ggsave("TopASVs_roottips_roottoshootratio_UPDATED112125.pdf", TopASVs_roottips_roottoshootratio, width = 8, height = 7, units = "in")



##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SH", "Guild"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
roottoshootratio_phylum_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_phylum_plot)

ggsave("curveplot_roottips_roottoshootratio_phyla_UPDATED112125.pdf", roottoshootratio_phylum_plot, width = 8, height = 7, units = "in")

roottoshootratio_genus_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Genus), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Genus), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Genus), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Genus", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_genus_plot)

ggsave("curveplot_roottips_roottoshootratio_genus_UPDATED112125.pdf", roottoshootratio_genus_plot, width = 11, height = 7, units = "in")


roottoshootratio_species_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Species), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Species), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Species), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Species", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_species_plot)

ggsave("curveplot_roottips_roottoshootratio_species_UPDATED112125.pdf", roottoshootratio_species_plot, width = 11, height = 7, units = "in")


roottoshootratio_guild_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Guild), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Guild), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Guild), shape = 21, size = 4, color = "black") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Guild", x = "VIP Score", y = "Correlation to Lodgepole Pine RoottoShootRatio (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(roottoshootratio_guild_plot)

ggsave("curveplot_roottips_roottoshootratio_guild_UPDATED112125.pdf", roottoshootratio_guild_plot, width = 11, height = 7, units = "in")



##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-c(ASV, Taxonomy), ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

filtered_asv <- filtered_asv[, -2]

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) #%>%


#transposed_asv_roottoshootratio <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, RoottoShootRatio), by = "Sample")
transposed_asv_roottoshootratio <- transposed_asv %>% left_join(metadata %>% dplyr::select(SampleID, RoottoShootRatio), by = c("Sample" = "SampleID"))


# Identify slope direction per ASV
asv_relationships <- transposed_asv_roottoshootratio %>%
  pivot_longer(cols = -c(Sample, RoottoShootRatio), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(RoottoShootRatio ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))


## Save the features positively associated with outcome of interest 
## Do this for both positive & negative associations, but separately 
plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
write_csv(asv_sums, "ITS_summed_relab_vip_roottips_roottoshootratio_positiveassociations_UPDATED112125.csv")
write_csv(asv_means, "ITS_mean_relab_vip_roottips_roottoshootratio_positiveassociations_UPDATED112125.csv")





