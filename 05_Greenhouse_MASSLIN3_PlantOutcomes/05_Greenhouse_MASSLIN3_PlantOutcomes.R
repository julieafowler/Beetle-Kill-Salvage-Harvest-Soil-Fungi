#### MASSLIN3 with Plant Outcomes ####

## Testing only associations with continuous plant output variables
## Basically, "I know treatments impact our plant variables, how do microbial features align around these plant changes?"
## Includes production of volcano plots 


#### Rhizosphere ####

library(maaslin3)
library(EnhancedVolcano)
library(dplyr)
library(readr)

setwd()

#### ASVs - Rhizosphere ####

# --- Load ---
feat <- read.delim("ITS_FeatureTable_Rhizosphere.txt", row.names = 1, check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt", row.names = 1, check.names = FALSE)


# If feature table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(feat))) stop("Sample IDs don't line up")
feat <- t(feat)  # now rows = samples, cols = features

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(feat)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - feat: samples (rows) × features (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(feat))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Height_volcano_abundance

ggsave("ASVs_Height_volcano_abundance.pdf", ASVs_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Height_volcano_prevalence

ggsave("ASVs_Height_volcano_prevalence.pdf", ASVs_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_MeanCanWidth_volcano_abundance

ggsave("ASVs_MeanCanWidth_volcano_abundance.pdf", ASVs_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_MeanCanWidth_volcano_prevalence

ggsave("ASVs_MeanCanWidth_volcano_prevalence.pdf", ASVs_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RootDryMass_volcano_abundance

ggsave("ASVs_RootDryMass_volcano_abundanc.pdf", ASVs_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RootDryMass_volcano_prevalence

ggsave("ASVs_RootDryMass_volcano_prevalence.pdf", ASVs_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Aboveground_Biomass_volcano_abundance

ggsave("ASVs_Aboveground_Biomass_volcano_abundance.pdf", ASVs_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Aboveground_Biomass_volcano_prevalence

ggsave("ASVs_Aboveground_Biomass_volcano_prevalence_111825.pdf", ASVs_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RoottoShootRatio_volcano_abundance

ggsave("ASVs_RoottoShootRatio_volcano_abundance.pdf", ASVs_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RoottoShootRatio_volcano_prevalence

ggsave("ASVs_RoottoShootRatio_volcano_prevalence.pdf", ASVs_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")




#### Guilds - Rhizosphere ####

setwd()

library(maaslin3)

# --- Load ---
guild <- read.delim("ITS_GuildTable_Rhizosphere.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt", row.names = 1, check.names = FALSE)


sample_cols <- names(guild)[-1]
guild[ sample_cols ] <- lapply(guild[ sample_cols ], function(x) as.numeric(as.character(x)))
guild_sum <- rowsum(as.matrix(guild[ sample_cols ]), group = guild$Guild, na.rm = TRUE)
guild_sum_df <- data.frame(Guild = rownames(guild_sum), guild_sum, row.names = NULL)
guild <- guild_sum_df
rownames(guild) <- guild[, 1]
guild <- guild[, -1]


# If your guild table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(guild))) stop("Sample IDs don't line up")
guild <- t(guild)  # now rows = samples, cols = guilds

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(guild)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - guild: samples (rows) × guilds (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(guild))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_Height_volcano_abundance

ggsave("guild_Height_volcano_abundance.pdf", guild_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

guild_Height_volcano_prevalence

ggsave("guild_Height_volcano_prevalence.pdf", guild_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_MeanCanWidth_volcano_abundance

ggsave("guild_MeanCanWidth_volcano_abundance.pdf", guild_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

guild_MeanCanWidth_volcano_prevalence

ggsave("guild_MeanCanWidth_volcano_prevalence_UPDATED111824.pdf", guild_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_RootDryMass_volcano_abundance

ggsave("guild_RootDryMass_volcano_abundance.pdf", guild_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

guild_RootDryMass_volcano_prevalence

ggsave("guild_RootDryMass_volcano_prevalence_111825.pdf", guild_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_Aboveground_Biomass_volcano_abundance

ggsave("guild_Aboveground_Biomass_volcano_abundance.pdf", guild_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

guild_Aboveground_Biomass_volcano_prevalence

ggsave("guild_Aboveground_Biomass_volcano_prevalence_111825.pdf", guild_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_RoottoShootRatio_volcano_abundance

ggsave("guild_RoottoShootRatio_volcano_abundance_111825.pdf", guild_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

guild_RoottoShootRatio_volcano_prevalence

ggsave("guild_RoottoShootRatio_volcano_prevalence_111825.pdf", guild_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")




#### Species - Rhizosphere ####

setwd()

library(maaslin3)

# --- Load ---
species <- read.delim("ITS_SpeciesTable_Rhizosphere.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt", row.names = 1, check.names = FALSE)




sample_cols <- names(species)[-1]
species[ sample_cols ] <- lapply(species[ sample_cols ], function(x) as.numeric(as.character(x)))
species_sum <- rowsum(as.matrix(species[ sample_cols ]), group = species$Species, na.rm = TRUE)
species_sum_df <- data.frame(Species = rownames(species_sum), species_sum, row.names = NULL)
species <- species_sum_df
rownames(species) <- species[, 1]
species <- species[, -1]


# If your species table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(species))) stop("Sample IDs don't line up")
species <- t(species)  # now rows = samples, cols = species

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(species)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - species: samples (rows) × speciess (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(species))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)

## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_Height_volcano_abundance

ggsave("species_Height_volcano_abundance_111825.pdf", species_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_Height_volcano_prevalence

ggsave("species_Height_volcano_prevalence_111825.pdf", species_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_MeanCanWidth_volcano_abundance

ggsave("species_MeanCanWidth_volcano_abundance_111825.pdf", species_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_MeanCanWidth_volcano_prevalence

ggsave("species_MeanCanWidth_volcano_prevalence_111825.pdf", species_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_RootDryMass_volcano_abundance

ggsave("species_RootDryMass_volcano_abundance_111825.pdf", species_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_RootDryMass_volcano_prevalence

ggsave("species_RootDryMass_volcano_prevalence_111825.pdf", species_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_Aboveground_Biomass_volcano_abundance

ggsave("species_Aboveground_Biomass_volcano_abundance_111825.pdf", species_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_Aboveground_Biomass_volcano_prevalence

ggsave("species_Aboveground_Biomass_volcano_prevalence_111825.pdf", species_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_RoottoShootRatio_volcano_abundance

ggsave("species_RoottoShootRatio_volcano_abundance.pdf", species_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_RoottoShootRatio_volcano_prevalence

ggsave("species_RoottoShootRatio_volcano_prevalence_111825.pdf", species_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")



#### Genus - Rhizosphere ####

setwd()

library(maaslin3)

# --- Load ---
genus <- read.delim("ITS_GenusTable_Rhizosphere.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt", row.names = 1, check.names = FALSE)




sample_cols <- names(genus)[-1]
genus[ sample_cols ] <- lapply(genus[ sample_cols ], function(x) as.numeric(as.character(x)))
genus_sum <- rowsum(as.matrix(genus[ sample_cols ]), group = genus$Genus, na.rm = TRUE)
genus_sum_df <- data.frame(Genus = rownames(genus_sum), genus_sum, row.names = NULL)
genus <- genus_sum_df
rownames(genus) <- genus[, 1]
genus <- genus[, -1]


# If your genus table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(genus))) stop("Sample IDs don't line up")
genus <- t(genus)  # now rows = samples, cols = genuss

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(genus)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - genus: samples (rows) × genuss (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(genus))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)

## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_Height_volcano_abundance

ggsave("genus_Height_volcano_abundance.pdf", genus_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_Height_volcano_prevalence

ggsave("genus_Height_volcano_prevalence.pdf", genus_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_MeanCanWidth_volcano_abundance

ggsave("genus_MeanCanWidth_volcano_abundance_111825.pdf", genus_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_MeanCanWidth_volcano_prevalence

ggsave("genus_MeanCanWidth_volcano_prevalence.pdf", genus_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_RootDryMass_volcano_abundance

ggsave("genus_RootDryMass_volcano_abundance.pdf", genus_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_RootDryMass_volcano_prevalence

ggsave("genus_RootDryMass_volcano_prevalence.pdf", genus_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")




## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_Aboveground_Biomass_volcano_abundance

ggsave("genus_Aboveground_Biomass_volcano_abundance.pdf", genus_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_Aboveground_Biomass_volcano_prevalence

ggsave("genus_Aboveground_Biomass_volcano_prevalence.pdf", genus_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")




## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_RoottoShootRatio_volcano_abundance

ggsave("genus_RoottoShootRatio_volcano_abundance_111825.pdf", genus_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_RoottoShootRatio_volcano_prevalence

ggsave("genus_RoottoShootRatio_volcano_prevalence_111825.pdf", genus_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")




#### EMF Only - ASVs - Rhizosphere ####

# --- Load ---
feat <- read.delim("ITS_FeatureTable_Rhizosphere_EMFOnly.txt", row.names = 1, check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt", row.names = 1, check.names = FALSE)




# If feature table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(feat))) stop("Sample IDs don't line up")
feat <- t(feat)  # now rows = samples, cols = features

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(feat)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - feat: samples (rows) × features (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(feat))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Height_EMFONLY_RHIZO",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Height_EMFONLY_RHIZO/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Height_volcano_abundance

ggsave("ASVs_Height_volcano_abundance_EMFONLY_RHIZO.pdf", ASVs_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Height_volcano_prevalence

ggsave("ASVs_Height_volcano_prevalence_EMFONLY_RHIZO.pdf", ASVs_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_MeanCanWidth_EMFONLY_RHIZO",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_MeanCanWidth_EMFONLY_RHIZO/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_MeanCanWidth_volcano_abundance

ggsave("ASVs_MeanCanWidth_volcano_abundance_EMFONLY_RHIZO.pdf", ASVs_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_MeanCanWidth_volcano_prevalence

ggsave("ASVs_MeanCanWidth_volcano_prevalence_EMFONLY_RHIZO.pdf", ASVs_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RootDryMass_EMFONLY_RHIZO",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RootDryMass_EMFONLY_RHIZO/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RootDryMass_volcano_abundance

ggsave("ASVs_RootDryMass_volcano_abundanc_EMFONLY_RHIZO.pdf", ASVs_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RootDryMass_volcano_prevalence

ggsave("ASVs_RootDryMass_volcano_prevalence_EMFONLY_RHIZO.pdf", ASVs_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Aboveground_Biomass_EMFONLY_RHIZO",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Aboveground_Biomass_EMFONLY_RHIZO/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Aboveground_Biomass_volcano_abundance

ggsave("ASVs_Aboveground_Biomass_volcano_abundance_EMFONLY_RHIZO.pdf", ASVs_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Aboveground_Biomass_volcano_prevalence

ggsave("ASVs_Aboveground_Biomass_volcano_prevalence_111825_EMFONLY_RHIZO.pdf", ASVs_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RoottoShootRatio_EMFONLY_RHIZO",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RoottoShootRatio_EMFONLY_RHIZO/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RoottoShootRatio_volcano_abundance

ggsave("ASVs_RoottoShootRatio_volcano_abundance_EMFONLY_RHIZO.pdf", ASVs_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RoottoShootRatio_volcano_prevalence

ggsave("ASVs_RoottoShootRatio_volcano_prevalence_EMFONLY_RHIZO.pdf", ASVs_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")





#### Root Tips ####

library(maaslin3)
library(EnhancedVolcano)
library(dplyr)
library(readr)

setwd()

#### ASVs - Root Tips ####

# --- Load ---
feat <- read.delim("ITS_FeatureTable_RootTipsOnly.txt", row.names = 1, check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt", row.names = 1, check.names = FALSE)

# If feature table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(feat))) stop("Sample IDs don't line up")
feat <- t(feat)  # now rows = samples, cols = features

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
#meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(feat)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - feat: samples (rows) × features (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(feat))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Height_volcano_abundance

ggsave("ASVs_Height_volcano_abundance_RootTips.pdf", ASVs_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Height_volcano_prevalence

ggsave("ASVs_Height_volcano_prevalence_RootTips.pdf", ASVs_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_MeanCanWidth_volcano_abundance

ggsave("ASVs_MeanCanWidth_volcano_abundance_RootTips.pdf", ASVs_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_MeanCanWidth_volcano_prevalence

ggsave("ASVs_MeanCanWidth_volcano_prevalence_RootTips.pdf", ASVs_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RootDryMass_volcano_abundance

ggsave("ASVs_RootDryMass_volcano_abundance_RootTips.pdf", ASVs_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RootDryMass_volcano_prevalence

ggsave("ASVs_RootDryMass_volcano_prevalence_RootTips.pdf", ASVs_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Aboveground_Biomass_volcano_abundance

ggsave("ASVs_Aboveground_Biomass_volcano_abundance_RootTips.pdf", ASVs_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Aboveground_Biomass_volcano_prevalence

ggsave("ASVs_Aboveground_Biomass_volcano_prevalence_RootTips.pdf", ASVs_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RoottoShootRatio_volcano_abundance

ggsave("ASVs_RoottoShootRatio_volcano_abundance_RootTips.pdf", ASVs_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RoottoShootRatio_volcano_prevalence

ggsave("ASVs_RoottoShootRatio_volcano_prevalence_RootTips.pdf", ASVs_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")




#### Guilds - Root Tips ####

setwd()

library(maaslin3)

# --- Load ---
guild <- read.delim("ITS_GuildTable_RootTipsOnly.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt", row.names = 1, check.names = FALSE)
meta <- meta[, c(2, 1, 3:ncol(meta))]

sample_cols <- names(guild)[-1]
guild[ sample_cols ] <- lapply(guild[ sample_cols ], function(x) as.numeric(as.character(x)))
guild_sum <- rowsum(as.matrix(guild[ sample_cols ]), group = guild$Guild, na.rm = TRUE)
guild_sum_df <- data.frame(Guild = rownames(guild_sum), guild_sum, row.names = NULL)
guild <- guild_sum_df
rownames(guild) <- guild[, 1]
guild <- guild[, -1]


# If your guild table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(guild))) stop("Sample IDs don't line up")
guild <- t(guild)  # now rows = samples, cols = guilds

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
#meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(guild)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
#meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - guild: samples (rows) × guilds (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(guild))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_Height_volcano_abundance

ggsave("guild_Height_volcano_abundance_RootTips.pdf", guild_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = FALSE
)

guild_Height_volcano_prevalence

ggsave("guild_Height_volcano_prevalence_RootTips.pdf", guild_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_MeanCanWidth_volcano_abundance

ggsave("guild_MeanCanWidth_volcano_abundance_RootTips.pdf", guild_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = FALSE
)

guild_MeanCanWidth_volcano_prevalence

ggsave("guild_MeanCanWidth_volcano_prevalence_RootTips.pdf", guild_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_RootDryMass_volcano_abundance

ggsave("guild_RootDryMass_volcano_abundance_RootTips.pdf", guild_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = FALSE
)

guild_RootDryMass_volcano_prevalence

ggsave("guild_RootDryMass_volcano_prevalence_RootTips.pdf", guild_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_Aboveground_Biomass_volcano_abundance

ggsave("guild_Aboveground_Biomass_volcano_abundance_RootTips.pdf", guild_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = FALSE
)

guild_Aboveground_Biomass_volcano_prevalence

ggsave("guild_Aboveground_Biomass_volcano_prevalence_RootTips.pdf", guild_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  selectLab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_RoottoShootRatio_volcano_abundance

ggsave("guild_RoottoShootRatio_volcano_abundance_RootTips.pdf", guild_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  selectLab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = FALSE
)

guild_RoottoShootRatio_volcano_prevalence

ggsave("guild_RoottoShootRatio_volcano_prevalence_RootTips.pdf", guild_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")




#### Species - Root Tips ####

setwd()

library(maaslin3)

# --- Load ---
species <- read.delim("ITS_SpeciesTable_RootTipsOnly.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt", row.names = 1, check.names = FALSE)

sample_cols <- names(species)[-1]
species[ sample_cols ] <- lapply(species[ sample_cols ], function(x) as.numeric(as.character(x)))
species_sum <- rowsum(as.matrix(species[ sample_cols ]), group = species$Species, na.rm = TRUE)
species_sum_df <- data.frame(Species = rownames(species_sum), species_sum, row.names = NULL)
species <- species_sum_df
rownames(species) <- species[, 1]
species <- species[, -1]


# If your species table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(species))) stop("Sample IDs don't line up")
species <- t(species)  # now rows = samples, cols = speciess

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
#meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(species)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
#meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - species: samples (rows) × speciess (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(species))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)

## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_Height_volcano_abundance

ggsave("species_Height_volcano_abundance_RootTips.pdf", species_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_Height_volcano_prevalence

ggsave("species_Height_volcano_prevalence_RootTips.pdf", species_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_MeanCanWidth_volcano_abundance

ggsave("species_MeanCanWidth_volcano_abundance_RootTips.pdf", species_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_MeanCanWidth_volcano_prevalence

ggsave("species_MeanCanWidth_volcano_prevalence_RootTips.pdf", species_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_RootDryMass_volcano_abundance

ggsave("species_RootDryMass_volcano_abundance_RootTips.pdf", species_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_RootDryMass_volcano_prevalence

ggsave("species_RootDryMass_volcano_prevalence_RootTips.pdf", species_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_Aboveground_Biomass_volcano_abundance

ggsave("species_Aboveground_Biomass_volcano_abundance_RootTips.pdf", species_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_Aboveground_Biomass_volcano_prevalence

ggsave("species_Aboveground_Biomass_volcano_prevalence_RootTips.pdf", species_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(species))
maaslin3(
  input_data      = species[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_species_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_species_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
species_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

species_RoottoShootRatio_volcano_abundance

ggsave("species_RoottoShootRatio_volcano_abundance_RootTips.pdf", species_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

species_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

species_RoottoShootRatio_volcano_prevalence

ggsave("species_RoottoShootRatio_volcano_prevalence_RootTips.pdf", species_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")







#### Genus - Root Tips ####

setwd()

library(maaslin3)

# --- Load ---
genus <- read.delim("ITS_GenusTable_RootTipsOnly.txt", check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt", row.names = 1, check.names = FALSE)

sample_cols <- names(genus)[-1]
genus[ sample_cols ] <- lapply(genus[ sample_cols ], function(x) as.numeric(as.character(x)))
genus_sum <- rowsum(as.matrix(genus[ sample_cols ]), group = genus$Genus, na.rm = TRUE)
genus_sum_df <- data.frame(Genus = rownames(genus_sum), genus_sum, row.names = NULL)
genus <- genus_sum_df
rownames(genus) <- genus[, 1]
genus <- genus[, -1]


# If your genus table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(genus))) stop("Sample IDs don't line up")
genus <- t(genus)  # now rows = samples, cols = genuss

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
#meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(genus)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")

meta$Salvage_Harvest_Status <- relevel(factor(meta$Salvage_Harvest_Status), ref = "Non_Salvage_Harvested")
#meta$Logging <- relevel(factor(meta$Logging), ref = "Non_Logged")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - genus: samples (rows) × genuss (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(genus))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_Height",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)

## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_Height/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_Height_volcano_abundance

ggsave("genus_Height_volcano_abundance_RootTips.pdf", genus_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_Height_volcano_prevalence

ggsave("genus_Height_volcano_prevalence_RootTips.pdf", genus_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_MeanCanWidth",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_MeanCanWidth/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_MeanCanWidth_volcano_abundance

ggsave("genus_MeanCanWidth_volcano_abundance_RootTips.pdf", genus_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_MeanCanWidth_volcano_prevalence

ggsave("genus_MeanCanWidth_volcano_prevalence_RootTips.pdf", genus_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_RootDryMass",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_RootDryMass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_RootDryMass_volcano_abundance

ggsave("genus_RootDryMass_volcano_abundance_RootTips.pdf", genus_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_RootDryMass_volcano_prevalence

ggsave("genus_RootDryMass_volcano_prevalence_RootTips.pdf", genus_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")




## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_Aboveground_Biomass",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_Aboveground_Biomass/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_Aboveground_Biomass_volcano_abundance

ggsave("genus_Aboveground_Biomass_volcano_abundance_RootTips.pdf", genus_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_Aboveground_Biomass_volcano_prevalence

ggsave("genus_Aboveground_Biomass_volcano_prevalence_RootTips.pdf", genus_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")




## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(genus))
maaslin3(
  input_data      = genus[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_genus_RoottoShootRatio",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_genus_RoottoShootRatio/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
genus_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

genus_RoottoShootRatio_volcano_abundance

ggsave("genus_RoottoShootRatio_volcano_abundance_RootTips.pdf", genus_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

genus_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

genus_RoottoShootRatio_volcano_prevalence

ggsave("genus_RoottoShootRatio_volcano_prevalence_RootTips.pdf", genus_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")





#### EMF Only - ASVs - Root Tips ####

# --- Load ---
feat <- read.delim("ITS_FeatureTable_RootTipsOnly_EMFOnly.txt", row.names = 1, check.names = FALSE)  # ASVs x samples OK
meta <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt", row.names = 1, check.names = FALSE)

# If feature table is ASVs as rows, flip to samples-as-rows (optional clarity)
if (!all(rownames(meta) %in% colnames(feat))) stop("Sample IDs don't line up")
feat <- t(feat)  # now rows = samples, cols = features

# Outcomes and covariates
outs <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio")
meta[outs] <- lapply(meta[outs], function(x) as.numeric(as.character(x)))

meta$Treatment_LandUseHistory <- factor(meta$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB","1stPostMPB", "2ndPostMPB","Recent_Cut"))
meta$Salvage_Harvest_Status  <- factor(meta$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested"))
#meta$Logging  <- factor(meta$Logging, levels = c("Non_Logged", "1980s", "1990s", "2000s", "2010s", "2018_orLater"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(feat)  # remove if already supplied relative abundances



meta$Treatment_LandUseHistory <- relevel(factor(meta$Treatment_LandUseHistory), ref = "Old_Forest")


meta$Height <- as.numeric(as.character(meta$Height))
meta$MeanCanWidth <- as.numeric(as.character(meta$MeanCanWidth))
meta$RootDryMass <- as.numeric(as.character(meta$RootDryMass))
meta$Aboveground_Biomass <- as.numeric(as.character(meta$Aboveground_Biomass))
meta$RoottoShootRatio <- as.numeric(as.character(meta$RoottoShootRatio))



# Assumes:
# - feat: samples (rows) × features (cols)
# - meta: samples (rows) × metadata (cols)
# - meta has Treatment_LandUseHistory, and ReadDepth
#   (If need ReadDepth from counts: meta$ReadDepth <- rowSums(feat))


## 1) Height
this_meta <- subset(meta, !is.na(Height))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Height_EMFONLY_ROOTTIPS",
  formula         = ~ Height + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Height_EMFONLY_ROOTTIPS/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Height"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Height_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Height_volcano_abundance

ggsave("ASVs_Height_volcano_abundance_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_Height_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Height_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Height_volcano_prevalence

ggsave("ASVs_Height_volcano_prevalence_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_Height_volcano_prevalence, width = 8, height = 7, units = "in")



## 2) MeanCanWidth
this_meta <- subset(meta, !is.na(MeanCanWidth))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_MeanCanWidth_EMFONLY_ROOTTIPS",
  formula         = ~ MeanCanWidth + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_MeanCanWidth_EMFONLY_ROOTTIPS/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MeanCanWidth"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_MeanCanWidth_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_MeanCanWidth_volcano_abundance

ggsave("ASVs_MeanCanWidth_volcano_abundance_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_MeanCanWidth_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_MeanCanWidth_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_MeanCanWidth_volcano_prevalence

ggsave("ASVs_MeanCanWidth_volcano_prevalence_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_MeanCanWidth_volcano_prevalence, width = 8, height = 7, units = "in")



## 3) RootDryMass
this_meta <- subset(meta, !is.na(RootDryMass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RootDryMass_EMFONLY_ROOTTIPS",
  formula         = ~ RootDryMass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RootDryMass_EMFONLY_ROOTTIPS/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RootDryMass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RootDryMass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RootDryMass_volcano_abundance

ggsave("ASVs_RootDryMass_volcano_abundance_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_RootDryMass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RootDryMass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RootDryMass_volcano_prevalence

ggsave("ASVs_RootDryMass_volcano_prevalence_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_RootDryMass_volcano_prevalence, width = 8, height = 7, units = "in")



## 4) Aboveground_Biomass
this_meta <- subset(meta, !is.na(Aboveground_Biomass))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_Aboveground_Biomass_EMFONLY_ROOTTIPS",
  formula         = ~ Aboveground_Biomass + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_Aboveground_Biomass_EMFONLY_ROOTTIPS/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "Aboveground_Biomass"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_Aboveground_Biomass_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_Aboveground_Biomass_volcano_abundance

ggsave("ASVs_Aboveground_Biomass_volcano_abundance_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_Aboveground_Biomass_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_Aboveground_Biomass_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_Aboveground_Biomass_volcano_prevalence

ggsave("ASVs_Aboveground_Biomass_volcano_prevalence_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_Aboveground_Biomass_volcano_prevalence, width = 8, height = 7, units = "in")



## 5) RoottoShootRatio
this_meta <- subset(meta, !is.na(RoottoShootRatio))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_RoottoShootRatio_EMFONLY_ROOTTIPS",
  formula         = ~ RoottoShootRatio + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_RoottoShootRatio_EMFONLY_ROOTTIPS/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "RoottoShootRatio"
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,                 # abundance coef is already log2-scale
    q_plot = qval_individual
  )

# Optional: set rownames to labels EnhancedVolcano will use
tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_RoottoShootRatio_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,                  # q-value cutoff want to highlight
  FCcutoff = log2(1.5),            # ~1.5x per SD threshold
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_RoottoShootRatio_volcano_abundance

ggsave("ASVs_RoottoShootRatio_volcano_abundance_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_RoottoShootRatio_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),   # log-odds -> log2(odds ratio)
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_RoottoShootRatio_volcano_prevalence <- EnhancedVolcano(
  tt_prev,
  lab = rownames(tt_prev),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05,
  FCcutoff = log2(1.5),
  title = paste0("Volcano: ", metadatum, " (prevalence)"),
  subtitle = "Effect = log2(odds ratio) per 1 SD in outcome",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6, labSize = 3.0, drawConnectors = TRUE
)

ASVs_RoottoShootRatio_volcano_prevalence

ggsave("ASVs_RoottoShootRatio_volcano_prevalence_RootTips_EMFONLY_ROOTTIPS.pdf", ASVs_RoottoShootRatio_volcano_prevalence, width = 8, height = 7, units = "in")







