#### Greenhouse Bioassays - ITS Data ####

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

library(vegan)

setwd()
ITS_rhizo <- read.delim("ITS_FeatureTable_Rhizosphere.txt",header=T,row.names=1, check.names=FALSE)
ITS_rhizo_t <- t(ITS_rhizo)
sample_sums <- rowSums(ITS_rhizo_t)
summary(sample_sums)
rarecurve(ITS_rhizo_t, step = 500, label = FALSE)
which(sample_sums < 28000)
sum(sample_sums < 28000)
ITS_rhizo_t_filtered <- ITS_rhizo_t[rowSums(ITS_rhizo_t) >= 28000, ]
ITS_rhizo_t_rarefied_28000 <- rrarefy(ITS_rhizo_t_filtered, sample = 28000)
ITS_rhizo_t_rarefied_28000 <- t(ITS_rhizo_t_rarefied_28000)


setwd()
ITS_roottips <- read.delim("ITS_FeatureTable_RootTipsOnly.txt",header=T,row.names=1, check.names=FALSE)
ITS_roottips_t <- t(ITS_roottips)
sample_sums <- rowSums(ITS_roottips_t)
summary(sample_sums)
rarecurve(ITS_roottips_t, step = 500, label = FALSE)
which(sample_sums < 22000)
sum(sample_sums < 22000)
ITS_roottips_t_filtered <- ITS_roottips_t[rowSums(ITS_roottips_t) >= 22000, ]
ITS_roottips_t_rarefied_22000 <- rrarefy(ITS_roottips_t_filtered, sample = 22000)
ITS_roottips_t_rarefied_22000 <- t(ITS_roottips_t_rarefied_22000)



#### **Rhizosphere** ####

setwd()

## Use rarefied counts
otus <- ITS_rhizo_t_rarefied_28000
otus <- as.data.frame(otus)
sample_names <- colnames(otus)
# Filter metadata to only samples present in rarefied count data
map_file <- read.delim("ITS_metadata_Rhizosphere.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)

all.equal(names(otus),row.names(map_file)) 

map_file$Treatment_LandUseHistory <- factor(map_file$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
map_file$MPB_Impacted <- factor(map_file$MPB_Impacted, levels = c("Non_MPB", "MPB")) 


otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("ITS_taxonomy_Rhizosphere.txt",header=T,row.names=1) 
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat))
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) 
all.equal(row.names(map_file),sample_names(physeq)) 
sampledata<-sample_data(map_file)
mgd<-merge_phyloseq(physeq,sampledata)



#### NMDS ####

mgd_ge5K<-mgd
mgd_ge5K_relabund<-transform_sample_counts(mgd_ge5K,function(x)x/sum(x)) 
mgd_ge5K_relabund.bray<-distance(mgd_ge5K_relabund,"bray")
mgd_ge5K_relabund.bray.nmds<-ordinate(mgd_ge5K_relabund,"NMDS",mgd_ge5K_relabund.bray)
mgd_ge5K_relabund.bray

mgd_ge5K_relabund.bray.nmds$stress

mgd_relabund_map=as(sample_data(mgd_ge5K_relabund),"data.frame")
sample_tab<-mgd_relabund_map
head(sample_tab)

sample_tab$NMDS1<-mgd_ge5K_relabund.bray.nmds$points[,1]
sample_tab$NMDS2<-mgd_ge5K_relabund.bray.nmds$points[,2]

plot(mgd_ge5K_relabund.bray.nmds)

NMDS <- ggplot(sample_tab) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=Treatment_LandUseHistory), size=4)+
  theme(text=element_text(size = 20)) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = MPB_Impacted, color=MPB_Impacted),linetype = 2) +
  #ylim(-0.4, 0.4) + xlim(-0.6,0.6) 
  theme(legend.position = "none")

NMDS

ggsave("NMDS_Rhizosphere_Rarefied_nokey.pdf", NMDS, width = 14, height = 10, units = "in")



#### Beta Dispersion #### 

# Load necessary packages
library(vegan)
library(ggplot2)
library(dplyr)

mod <- betadisper(mgd_ge5K_relabund.bray, sampledata$MPB_Impacted, bias.adjust = TRUE)
plot(mod)

beta_df <- data.frame(
  Group = mod$group,
  DistanceToCentroid = mod$distances
)

set.seed(123)
perm_test <- permutest(mod, permutations = 999)

betadispersion_boxplot <- ggplot(beta_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    y = "Distance to Group Centroid",
    x = "Group",
    title = "Beta Dispersion by Group (Bray-Curtis)"
  ) +
  theme(legend.position="none")

betadispersion_boxplot

print(perm_test)

ggsave("NMDS_Rhizosphere_Rarefied_BetaDispersion_MPBImpacted.pdf", betadispersion_boxplot, width = 2, height = 5, units = "in")



#### PERMANOVA ####

adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ MPB_Impacted, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory*MPB_Impacted, data = sample_tab, permutations = 999, method = "bray")


## P-value adjustments:
library(stats)
p = c(0.001, 0.001) ## Change with above findings 
p.adjusted = p.adjust(p, method = "BH")
p.adjusted

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)


#### NMDS w/ Envfit for FUNGuild ####


## Making a spreadsheet with my OTU table where I kept the same layout and metadata as map_file 
## but with each column representing an OTU, but OTU code replaced with guild of that OTU
## Now going to sum all columns of the same guild

## Then need to normalize each row to one 
## Then move forward with envfit 


## First step - OTU table with a column added for FUNGuild guild to replace the OTU code 
envfit_FUNGuild <- read.delim("ITS_FeatureTable_FUNGuild_Envfit_Rhizosphere.txt",header = T)
# Make sure only those samples that made it past rarefying are included 
sample_names <- colnames(otus)
# Extract the first column name of envfit_FUNGuild to keep it 
guild_column <- colnames(envfit_FUNGuild)[1]
# Find which of envfit_FUNGuild 's columns (after the first) match rarefied otus's column names
matching_samples <- intersect(sample_names, colnames(envfit_FUNGuild))
# Subset envfit_FUNGuild to keep the first column and matching sample columns
envfit_FUNGuild <- envfit_FUNGuild[, c(guild_column, matching_samples)]

## Collapse all rows of the same guild into one for each guild
envfit_FUNGuild_agg <- aggregate(. ~  Guild, data = envfit_FUNGuild, sum)

## Transpose the table so that each column is a guild
envfit_FUNGuild_pivot <- t(envfit_FUNGuild_agg)


library(tibble)

# Make the first row the column names
colnames(envfit_FUNGuild_pivot) <- as.character(envfit_FUNGuild_pivot[1, ])
envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[-1, ]

all.equal(row.names(envfit_FUNGuild_pivot),row.names(map_file))
identical(row.names(envfit_FUNGuild_pivot), row.names(map_file)) 

if ("Unknown" %in% colnames(envfit_FUNGuild_pivot)) {
  envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[, colnames(envfit_FUNGuild_pivot) != "Unknown"]
}

envfit_FUNGuild_pivot <- as.data.frame(envfit_FUNGuild_pivot)
# Correct skew & calculate z scores through scaling 
envfit_FUNGuild_pivot[] <- lapply(envfit_FUNGuild_pivot, as.numeric)  # ensure numeric
str(envfit_FUNGuild_pivot)

# Get rid of columns where there are only 0s (happens sometimes when dropping samples due to rarefying)
envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[, colSums(envfit_FUNGuild_pivot) != 0]
# Log-transform to handle skew and zeros
envfit_FUNGuild_pivot_log <- log1p(envfit_FUNGuild_pivot)
# Standardize to mean 0, SD 1 for comparability across guilds
envfit_FUNGuild_pivot_log_scaled <- scale(envfit_FUNGuild_pivot_log)

envfit_FUNGuild_pivot_log_scaled <- as.data.frame(envfit_FUNGuild_pivot_log_scaled)

nmds_samples <- rownames(mgd_ge5K_relabund.bray.nmds$points)
envfit_samples <- rownames(envfit_FUNGuild_pivot_log_scaled)
identical(nmds_samples, envfit_samples)


options(ggrepel.max.overlaps = Inf)

ord_FUNGuild <- mgd_ge5K_relabund.bray.nmds
ord_scores <- scores(ord_FUNGuild) 
fit_FUNGuild <- envfit(ord_scores, envfit_FUNGuild_pivot_log_scaled, permutations = 10000, na.rm=TRUE)
fit_FUNGuild

p.adjust.envfit <- function (x, method = 'bonferroni', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}

fit_adjust <- p.adjust.envfit(fit_FUNGuild)
fit_adjust

fit1 <- as.data.frame(scores(fit_adjust, display = "vectors"))
fit1 <- cbind(fit1, env.variables = rownames(fit1))
fit1 <- cbind(fit1, pval = fit_adjust$vectors$pvals) 
env.scores.fit <- cbind(fit1, r = fit_adjust$vectors$r)
sig.env.scrs <- subset(env.scores.fit, pval<=0.05)

library(ggrepel)

NMDS_Envfit_FUNGuild <- NMDS +
  geom_segment(data=sig.env.scrs,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow=arrow(length=unit(0.5,'cm')),color='bisque4',inherit.aes=FALSE) + 
  coord_fixed() +
  geom_text_repel(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  theme(legend.position = "none") 

NMDS_Envfit_FUNGuild

ggsave("NMDS_Envfit_FUNGuild_Rhizosphere_Rarefied.pdf", NMDS_Envfit_FUNGuild, width = 7, height = 5, units = "in")


#### NMDS w/ Envfit for Plant Outcomes ####

## Envfit with how patchy and frequent my NAs are in the plant outcome data make envfit difficult
## The way to do this may be iterating envfit, one variable at a time 
## Best option to maximize sample use for each outcome metric.

library(dplyr)

envfit_PlantOutcomes <- read.delim("ITS_metadata_Rhizosphere_PlantOutcomes.txt",header = T)

rownames(envfit_PlantOutcomes) <- envfit_PlantOutcomes$SampleID
outcome_vars <- c("Height","MeanCanWidth","RootDryMass","Aboveground_Biomass","RoottoShootRatio","RootTips_Colonized_RawCounts","RootTips_Metric_ColonizedTipsper_RootDryMass")

# Remove metadata rows that were removed by rarefying in NMDS
envfit_PlantOutcomes_aligned <- envfit_PlantOutcomes[rownames(envfit_PlantOutcomes) %in% rownames(ord_scores), ]
identical(rownames(envfit_PlantOutcomes_aligned), rownames(ord_scores))


run_envfit_iterative <- function(ord_scores, envfit_PlantOutcomes_aligned, outcome_vars) {
  results <- list()
  vector_list <- list()
  
  for (var in outcome_vars) {
    cat("Now running:", var, "\n")
    
    if (!var %in% colnames(envfit_PlantOutcomes_aligned)) {
      warning(paste("Column not found in envfit_PlantOutcomes_aligned:", var))
      next
    }
    
    vec <- envfit_PlantOutcomes_aligned[[var]]
    keep <- !is.na(vec)
    
    if (sum(keep) >= 3) {
      ord_sub <- ord_scores[keep, , drop = FALSE]
      
      vec_transformed <- scale(log1p(vec[keep]))
      env_sub <- data.frame(var = vec_transformed)
      rownames(env_sub) <- rownames(ord_sub)
      
      fit <- envfit(ord_sub, env_sub, permutations = 999)
      
      if (is.null(fit$vectors) || is.null(scores(fit, display = "vectors"))) {
        warning(paste("envfit failed for", var))
        next
      }
      
      fit_scores <- scores(fit, display = "vectors")
      
      if (length(fit_scores) < 2) {
        warning(paste("fit_scores too short for", var))
        next
      }
      
      vec_coords <- data.frame(
        var = var,
        Dim1 = fit_scores[1],
        Dim2 = fit_scores[2],
        r2 = unname(fit$vectors$r),
        pval = unname(fit$vectors$pvals),
        n_samples = sum(keep),
        stringsAsFactors = FALSE
      )
      
      results[[var]] <- fit
      vector_list[[length(vector_list) + 1]] <- vec_coords
    } else {
      warning(paste("Too few non-NA values for", var))
    }
  }
  
  if (length(vector_list) == 0) {
    warning("No successful envfit results. Returning empty output.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors <- do.call(rbind, vector_list)
  cat("Number of successful vectors:", nrow(vectors), "\n")
  print(vectors)
  
  if (!"pval" %in% colnames(vectors)) {
    warning("p-values missing in vectors. Returning empty.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors$pval_adj <- p.adjust(vectors$pval, method = "BH")
  print("Adjusted p-values:")
  print(vectors$pval_adj)
  
  vectors_sig <- vectors %>% filter(pval_adj <= 0.05)
  
  return(list(fits = results, vectors = vectors_sig, summary = vectors))
}


fit_result <- run_envfit_iterative(ord_scores, envfit_PlantOutcomes_aligned, outcome_vars)
str(fit_result$summary)
fit_result$summary


library(ggplot2)
library(ggrepel)

sig.env.scrs <- fit_result$vectors
sig.env.scrs$NMDS1 <- sig.env.scrs$Dim1
sig.env.scrs$NMDS2 <- sig.env.scrs$Dim2
sig.env.scrs$env.variables <- sig.env.scrs$var

NMDS


NMMS_Rhizo_Plant_Significant <- NMDS +
  geom_segment(data = sig.env.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'bisque4',
               inherit.aes = FALSE) +
  geom_text_repel(data = sig.env.scrs, 
                  aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  coord_fixed()

NMMS_Rhizo_Plant_Significant

ggsave("NMMS_Rhizo_Plant_Significant_wRootTips.pdf", NMMS_Rhizo_Plant_Significant, width = 8, height = 6, units = "in")




## Plotting all vectors regardless of significance

all.env.scrs <- fit_result$summary
all.env.scrs$NMDS1 <- all.env.scrs$Dim1
all.env.scrs$NMDS2 <- all.env.scrs$Dim2
all.env.scrs$env.variables <- all.env.scrs$var

NMDS_Rhizo_Plant_AllVectors_RegardlessofSig <- NMDS +
  geom_segment(data = all.env.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'darkgray',
               inherit.aes = FALSE) +
  geom_text_repel(data = all.env.scrs,
                  aes(x = NMDS1, y = NMDS2, label = env.variables),
                  colour = "darkgray", fontface = "italic") +
  coord_fixed()

NMDS_Rhizo_Plant_AllVectors_RegardlessofSig

ggsave("NMDS_Rhizo_Plant_AllVectors_RegardlessofSig.pdf", NMDS_Rhizo_Plant_AllVectors_RegardlessofSig, width = 14, height = 10, units = "in")



#### NMDS w/ Envfit for Soil Chemistry #### 

## Envfit with how patchy and frequent my NAs are in the chemistry data make envfit difficult
## The way to do this may be iterating envfit, one variable at a time 
## Best option to maximize sample use for each outcome metric.

library(dplyr)

envfit_ChemistryEnvfit <- read.delim("ITS_metadata_Rhizosphere_ChemistryEnvfit.txt",header = T)

# Set rownames to sample IDs
rownames(envfit_ChemistryEnvfit) <- envfit_ChemistryEnvfit$SampleID

# Select only the numerical outcome columns
outcome_vars <- c("pH", "Total_C_percent", "Total_N_percent", "Na_waterextract_mgperL",
                  "NH4_waterextract_mgperL", "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                  "Ca_waterextract_mgperL", "Cl_waterextract_mgperL", "NO3_waterextract_mgperL",
                  "PO4_waterextract_mgperL", "SO4_waterextract_mgperL")

# Remove metadata rows that were removed by rarefying in NMDS
envfit_ChemistryEnvfit_aligned <- envfit_ChemistryEnvfit[rownames(envfit_ChemistryEnvfit) %in% rownames(ord_scores), ]
identical(rownames(envfit_ChemistryEnvfit_aligned), rownames(ord_scores))


run_envfit_iterative <- function(ord_scores, envfit_ChemistryEnvfit_aligned, outcome_vars) {
  results <- list()
  vector_list <- list()
  
  for (var in outcome_vars) {
    cat("Now running:", var, "\n")
    
    if (!var %in% colnames(envfit_ChemistryEnvfit_aligned)) {
      warning(paste("Column not found in envfit_ChemistryEnvfit_aligned:", var))
      next
    }
    
    vec <- envfit_ChemistryEnvfit_aligned[[var]]
    keep <- !is.na(vec)
    
    if (sum(keep) >= 3) {
      ord_sub <- ord_scores[keep, , drop = FALSE]
      
      vec_transformed <- scale(log1p(vec[keep]))
      env_sub <- data.frame(var = vec_transformed)
      rownames(env_sub) <- rownames(ord_sub)
      
      fit <- envfit(ord_sub, env_sub, permutations = 999)
      
      if (is.null(fit$vectors) || is.null(scores(fit, display = "vectors"))) {
        warning(paste("envfit failed for", var))
        next
      }
      
      fit_scores <- scores(fit, display = "vectors")
      
      if (length(fit_scores) < 2) {
        warning(paste("fit_scores too short for", var))
        next
      }
      
      vec_coords <- data.frame(
        var = var,
        Dim1 = fit_scores[1],
        Dim2 = fit_scores[2],
        r2 = unname(fit$vectors$r),
        pval = unname(fit$vectors$pvals),
        n_samples = sum(keep),
        stringsAsFactors = FALSE
      )
      
      results[[var]] <- fit
      vector_list[[length(vector_list) + 1]] <- vec_coords
    } else {
      warning(paste("Too few non-NA values for", var))
    }
  }
  
  if (length(vector_list) == 0) {
    warning("No successful envfit results. Returning empty output.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors <- do.call(rbind, vector_list)
  cat("Number of successful vectors:", nrow(vectors), "\n")
  print(vectors)
  
  if (!"pval" %in% colnames(vectors)) {
    warning("p-values missing in vectors. Returning empty.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors$pval_adj <- p.adjust(vectors$pval, method = "BH")
  print("Adjusted p-values:")
  print(vectors$pval_adj)
  
  vectors_sig <- vectors %>% filter(pval_adj <= 0.05)
  
  return(list(fits = results, vectors = vectors_sig, summary = vectors))
}


fit_result <- run_envfit_iterative(ord_scores, envfit_ChemistryEnvfit_aligned, outcome_vars)
str(fit_result$summary)
fit_result$summary


library(ggplot2)
library(ggrepel)

sig.env.scrs <- fit_result$vectors
sig.env.scrs$NMDS1 <- sig.env.scrs$Dim1
sig.env.scrs$NMDS2 <- sig.env.scrs$Dim2
sig.env.scrs$env.variables <- sig.env.scrs$var

NMDS


NMDS_Rhizo_Chemistry_Significant <- NMDS +
  geom_segment(data = sig.env.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'bisque4',
               inherit.aes = FALSE) +
  geom_text_repel(data = sig.env.scrs, 
                  aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  coord_fixed()

NMDS_Rhizo_Chemistry_Significant

ggsave("NMDS_Rhizo_Chemistry_Significant.pdf", NMDS_Rhizo_Chemistry_Significant, width = 14, height = 10, units = "in")



## Plotting all vectors regardless of significance

all.env.scrs <- fit_result$summary
all.env.scrs$NMDS1 <- all.env.scrs$Dim1
all.env.scrs$NMDS2 <- all.env.scrs$Dim2
all.env.scrs$env.variables <- all.env.scrs$var

NMDS_Rhizo_Chemistry_AllVectors_RegardlessofSig <- NMDS +
  geom_segment(data = all.env.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'darkgray',
               inherit.aes = FALSE) +
  geom_text_repel(data = all.env.scrs,
                  aes(x = NMDS1, y = NMDS2, label = env.variables),
                  colour = "darkgray", fontface = "italic") +
  coord_fixed()

NMDS_Rhizo_Chemistry_AllVectors_RegardlessofSig

ggsave("NMDS_Rhizo_Chemistry_AllVectors_RegardlessofSig.pdf", NMDS_Rhizo_Chemistry_AllVectors_RegardlessofSig, width = 14, height = 10, units = "in")





#### EMF RA Boxplot ####

#### EMF RA by 6 treatments ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         EMF_RA = Ectomycorrhizal / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

EMF_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = EMF_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere EMF RA by Treatment",
       y = "EMF RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

EMF_RA_Greenhouse

ggsave("Greenhouse_EMF_RA.pdf", EMF_RA_Greenhouse, width = 3, height = 5, units = "in")

EMF_RA_Greenhouse <- EMF_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

EMF_RA_Greenhouse

ggsave("Greenhouse_EMF_RA_wtests.pdf", EMF_RA_Greenhouse, width = 5, height = 5, units = "in")


#### EMF RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         EMF_RA = Ectomycorrhizal / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

EMF_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = EMF_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere EMF RA by MPB Impact",
       y = "EMF RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

EMF_RA_Greenhouse_byMPB

ggsave("Greenhouse_EMF_RA_byMPB.pdf", EMF_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

EMF_RA_Greenhouse_byMPB <- EMF_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

EMF_RA_Greenhouse_byMPB

ggsave("Greenhouse_EMF_RA_byMPB_wtests.pdf", EMF_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "Rhizo_EMF_RA.csv", row.names = FALSE)



#### Plant Pathogen RA by 6 treatments ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         PlantPathogen_RA = Plant.Pathogen / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

PlantPathogen_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = PlantPathogen_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere PlantPathogen RA by Treatment",
       y = "PlantPathogen RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

PlantPathogen_RA_Greenhouse

ggsave("Greenhouse_PlantPathogen_RA.pdf", PlantPathogen_RA_Greenhouse, width = 3, height = 5, units = "in")

PlantPathogen_RA_Greenhouse <- PlantPathogen_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

PlantPathogen_RA_Greenhouse

ggsave("Greenhouse_PlantPathogen_RA_wtests.pdf", PlantPathogen_RA_Greenhouse, width = 5, height = 5, units = "in")


#### Plant Pathogen RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         PlantPathogen_RA = Plant.Pathogen / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

PlantPathogen_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = PlantPathogen_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere PlantPathogen RA by MPB Impact",
       y = "PlantPathogen RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

PlantPathogen_RA_Greenhouse_byMPB

ggsave("Greenhouse_PlantPathogen_RA_byMPB.pdf", PlantPathogen_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

PlantPathogen_RA_Greenhouse_byMPB <- PlantPathogen_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

PlantPathogen_RA_Greenhouse_byMPB

ggsave("Greenhouse_PlantPathogen_RA_byMPB_wtests.pdf", PlantPathogen_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "Rhizo_PlantPathogen_RA.csv", row.names = FALSE)


#### Fungal Parasite RA by 6 treatments ####

setwd() 

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         FungalParasite_RA = Fungal.Parasite / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

FungalParasite_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = FungalParasite_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere FungalParasite RA by Treatment",
       y = "FungalParasite RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

FungalParasite_RA_Greenhouse

ggsave("Greenhouse_FungalParasite_RA.pdf", FungalParasite_RA_Greenhouse, width = 3, height = 5, units = "in")

FungalParasite_RA_Greenhouse <- FungalParasite_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

FungalParasite_RA_Greenhouse

ggsave("Greenhouse_FungalParasite_RA_wtests.pdf", FungalParasite_RA_Greenhouse, width = 5, height = 5, units = "in")


#### Fungal Parasite RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")

library(ggplot2)
library(dplyr)

# Calculate relative abundance for the guild
fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         FungalParasite_RA = Fungal.Parasite / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

FungalParasite_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = FungalParasite_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere FungalParasite RA by MPB Impact",
       y = "FungalParasite RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

FungalParasite_RA_Greenhouse_byMPB

ggsave("Greenhouse_FungalParasite_RA_byMPB.pdf", FungalParasite_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

FungalParasite_RA_Greenhouse_byMPB <- FungalParasite_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

FungalParasite_RA_Greenhouse_byMPB

ggsave("Greenhouse_FungalParasite_RA_byMPB_wtests.pdf", FungalParasite_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "Rhizo_FungalParasite_RA.csv", row.names = FALSE)




#### Alpha Diversity ####

library(vegan)
library(ggplot2)
library(dplyr)
library(ggpubr)

asv_table <- t(otus)   # samples x ASVs
metadata <- map_file    # samples x variables

# Ensure sample IDs match
metadata$SampleID <- rownames(metadata)
asv_table <- as.data.frame(asv_table)
asv_table$SampleID <- rownames(asv_table)

# Merge ASV table and metadata
merged <- inner_join(asv_table, metadata, by = "SampleID")

counts <- merged %>%
  select(where(is.numeric)) %>% select(-Replicate)

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


alpha_div$Treatment_LandUseHistory <- factor(alpha_div$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
alpha_div$MPB_Impacted <- factor(alpha_div$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

alpha_div$Treatment_All <- Treatment_LandUseHistory
unique(alpha_div$Treatment_All)

alpha_div <- alpha_div %>%
  relocate(Observed, Shannon, .after = last_col())

write.csv(alpha_div, "ITS_alpha_div_rhizosphere_ObservedShannons.csv", row.names = FALSE)


## Figures 

## Observed ASVs

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_Rhizosphere_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Observed)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
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

ggsave("AD_Rhizosphere_ObservedSpecies_Rarefied_wtests.pdf", AD_Rhizosphere_ObservedSpecies_Rarefied, width = 3, height = 5, units = "in")


treatment_levels <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_Rhizosphere_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = MPB_Impacted, y = Observed)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Observed Species by MPB Impact",
       y = "Observed Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impacted),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_ObservedSpecies_Rarefied

ggsave("AD_Rhizosphere_ObservedSpecies_Rarefied_byMPB_wtests.pdf", AD_Rhizosphere_ObservedSpecies_Rarefied, width = 2, height = 5, units = "in")



## Shannons 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_Rhizosphere_Shannons_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Shannon)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
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

ggsave("AD_Rhizosphere_Shannons_Rarefied_wtests.pdf", AD_Rhizosphere_Shannons_Rarefied, width = 3, height = 5, units = "in")


treatment_levels <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_Rhizosphere_Shannons_Rarefied <- ggplot(alpha_div, aes(x = MPB_Impacted, y = Shannon)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Shannon Diversity by MPB Impact",
       y = "Shannon Index") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impacted),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_Shannons_Rarefied

ggsave("AD_Rhizosphere_Shannons_Rarefied_byMPB_wtests.pdf", AD_Rhizosphere_Shannons_Rarefied, width = 2, height = 5, units = "in")


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



#### FUNGuild Barplot ####

library(dplyr)
library(tidyr)
library(ggplot2)

## FUNGuild: Ectomycorrhizal fungi by soil depth and treatment

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_Rhizosphere.txt")


group_var <- "Treatment_LandUseHistory"
facet_row <- "."
facet_col <- "."
guild_start_col <- 9

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Guild",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$MPB_Impacted <- factor(df_long$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

unique(df_long$Guild)
guild_levels <- c(
  "Animal.Parasite", "Animal.Pathogen", "Arbuscular.Mycorrhizal",
  "Dung.Saprotroph", "Ectomycorrhizal", "Endophyte",          
  "Epiphyte", "Fungal.Parasite", "Lichen.Parasite",     
  "Plant.Pathogen", "Plant.Saprotroph", "Soil.Saprotroph",    
  "Undefined.Saprotroph", "Wood.Saprotroph", "Unknown"   
)
df_long$Guild <- factor(df_long$Guild, levels = guild_levels)


# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Guild")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()

final_plot <- ggplot(df_relabund, aes_string(x = group_var, y = "RelAbund", fill = "Guild")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_Rhizosphere.pdf", final_plot, width = 14, height = 10, units = "in")







#### **Root Tips** ####

setwd()

otus <- ITS_roottips_t_rarefied_22000
otus <- as.data.frame(otus)
sample_names <- colnames(otus)
map_file <- read.delim("ITS_metadata_RootTips.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)

all.equal(names(otus),row.names(map_file))

map_file$Treatment_LandUseHistory <- factor(map_file$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
map_file$MPB_Impacted <- factor(map_file$MPB_Impacted, levels = c("Non_MPB", "MPB")) 



otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("ITS_taxonomy_RootTips.txt",header=T,row.names=1)
taxmat<-as.matrix(taxa)
all.equal(row.names(taxmat),row.names(otumat))
TAX = tax_table(taxmat)
physeq<-phyloseq(OTU,TAX) 
all.equal(row.names(map_file),sample_names(physeq))
sampledata<-sample_data(map_file)
mgd<-merge_phyloseq(physeq,sampledata)



#### NMDS ####

mgd_ge5K<-mgd
mgd_ge5K_relabund<-transform_sample_counts(mgd_ge5K,function(x)x/sum(x))
mgd_ge5K_relabund.bray<-distance(mgd_ge5K_relabund,"bray") 
mgd_ge5K_relabund.bray.nmds<-ordinate(mgd_ge5K_relabund,"NMDS",mgd_ge5K_relabund.bray)
mgd_ge5K_relabund.bray

mgd_ge5K_relabund.bray.nmds$stress

mgd_relabund_map=as(sample_data(mgd_ge5K_relabund),"data.frame") 
sample_tab<-mgd_relabund_map
head(sample_tab)

sample_tab$NMDS1<-mgd_ge5K_relabund.bray.nmds$points[,1]
sample_tab$NMDS2<-mgd_ge5K_relabund.bray.nmds$points[,2]

plot(mgd_ge5K_relabund.bray.nmds)

NMDS <- ggplot(sample_tab) +
  geom_point(aes(x=NMDS1, y=NMDS2, color=Treatment_LandUseHistory), size=4)+
  theme(text=element_text(size = 20)) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = MPB_Impacted, color=MPB_Impacted),linetype = 2) 

NMDS

ggsave("NMDS_RootTips_Rarefied_key.pdf", NMDS, width = 14, height = 10, units = "in")



#### Beta Dispersion #### 

## For groups with only 2 factors (MPB vs. MPB_Impacted)

# Load necessary packages
library(vegan)
library(ggplot2)
library(dplyr)

mod <- betadisper(mgd_ge5K_relabund.bray, sampledata$MPB_Impacted, bias.adjust = TRUE)
plot(mod)

beta_df <- data.frame(
  Group = mod$group,
  DistanceToCentroid = mod$distances
)

set.seed(123)
perm_test <- permutest(mod, permutations = 999)

betadispersion_boxplot <- ggplot(beta_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    y = "Distance to Group Centroid",
    x = "Group",
    title = "Beta Dispersion by Group (Bray-Curtis)"
  ) +
  theme(legend.position="none")

betadispersion_boxplot

print(perm_test)

ggsave("NMDS_RootTips_Rarefied_BetaDispersion_MPBImpacted.pdf", betadispersion_boxplot, width = 2, height = 5, units = "in")



#### PERMANOVA ####

adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ MPB_Impacted, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory*MPB_Impacted, data = sample_tab, permutations = 999, method = "bray")


## P-value adjustments:
library(stats)
p = c(0.001, 0.006) ## Change with above findings 
p.adjusted = p.adjust(p, method = "BH")
p.adjusted

# Adjusted: 0.001666667 0.001666667 0.002500000 0.001666667 0.003000000


library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)


#### NMDS w/ Envfit for FUNGuild ####

envfit_FUNGuild <- read.delim("ITS_FeatureTable_FUNGuild_Envfit_RootTips.txt",header = T)
sample_names <- colnames(otus)
guild_column <- colnames(envfit_FUNGuild)[1]
matching_samples <- intersect(sample_names, colnames(envfit_FUNGuild))
envfit_FUNGuild <- envfit_FUNGuild[, c(guild_column, matching_samples)]

envfit_FUNGuild_agg <- aggregate(. ~  Guild, data = envfit_FUNGuild, sum)

envfit_FUNGuild_pivot <- t(envfit_FUNGuild_agg)


library(tibble)

colnames(envfit_FUNGuild_pivot) <- as.character(envfit_FUNGuild_pivot[1, ])
envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[-1, ]

all.equal(row.names(envfit_FUNGuild_pivot),row.names(map_file))
identical(row.names(envfit_FUNGuild_pivot), row.names(map_file))

if ("Unknown" %in% colnames(envfit_FUNGuild_pivot)) {
  envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[, colnames(envfit_FUNGuild_pivot) != "Unknown"]
}

envfit_FUNGuild_pivot <- as.data.frame(envfit_FUNGuild_pivot)
envfit_FUNGuild_pivot[] <- lapply(envfit_FUNGuild_pivot, as.numeric)
str(envfit_FUNGuild_pivot)

envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[, colSums(envfit_FUNGuild_pivot) != 0]
envfit_FUNGuild_pivot_log <- log1p(envfit_FUNGuild_pivot)
envfit_FUNGuild_pivot_log_scaled <- scale(envfit_FUNGuild_pivot_log)

envfit_FUNGuild_pivot_log_scaled <- as.data.frame(envfit_FUNGuild_pivot_log_scaled)

nmds_samples <- rownames(mgd_ge5K_relabund.bray.nmds$points)
envfit_samples <- rownames(envfit_FUNGuild_pivot_log_scaled)
identical(nmds_samples, envfit_samples)



options(ggrepel.max.overlaps = Inf)


ord_FUNGuild <- mgd_ge5K_relabund.bray.nmds
ord_scores <- scores(ord_FUNGuild) 
fit_FUNGuild <- envfit(ord_scores, envfit_FUNGuild_pivot_log_scaled, permutations = 10000, na.rm=TRUE)
fit_FUNGuild

p.adjust.envfit <- function (x, method = 'bonferroni', n)
{
  x.new <- x
  if (!is.null (x$vectors)) pval.vectors <- x$vectors$pvals else pval.vectors <- NULL
  if (!is.null (x$factors)) pval.factors <- x$factors$pvals else pval.factors <- NULL
  if (missing (n)) n <- length (pval.vectors) + length (pval.factors)
  if (!is.null (x$vectors)) x.new$vectors$pvals <- p.adjust (x$vectors$pvals, method = method, n = n)
  if (!is.null (x$factors)) x.new$factors$pvals <- p.adjust (x$factors$pvals, method = method, n = n)
  cat ('Adjustment of significance by', method, 'method')
  return (x.new)
}

fit_adjust <- p.adjust.envfit(fit_FUNGuild)
fit_adjust

fit1 <- as.data.frame(scores(fit_adjust, display = "vectors"))
fit1 <- cbind(fit1, env.variables = rownames(fit1))
fit1 <- cbind(fit1, pval = fit_adjust$vectors$pvals)
env.scores.fit <- cbind(fit1, r = fit_adjust$vectors$r)
sig.env.scrs <- subset(env.scores.fit, pval<=0.05)

library(ggrepel)

## Plot 
NMDS_Envfit_FUNGuild <- NMDS +
  geom_segment(data=sig.env.scrs,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow=arrow(length=unit(0.5,'cm')),color='bisque4',inherit.aes=FALSE) + 
  coord_fixed() +
  geom_text_repel(data = sig.env.scrs, aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  theme(legend.position = "none") 
  # ylim(-0.5, 0.7) + xlim(-0.8,0.8)

NMDS_Envfit_FUNGuild

ggsave("NMDS_Envfit_FUNGuild_RootTips_Rarefied.pdf", NMDS_Envfit_FUNGuild, width = 7, height = 5, units = "in")


#### NMDS w/ Envfit for Plant Outcomes ####

## Envfit with how patchy and frequent my NAs are in the plant outcome data make envfit difficult
## The way to do this may be iterating envfit, one variable at a time 
## Best option to maximize sample use for each outcome metric.

library(dplyr)

envfit_PlantOutcomes <- read.delim("ITS_metadata_RootTips_PlantOutcomes.txt",header = T)

rownames(envfit_PlantOutcomes) <- envfit_PlantOutcomes$SampleID

outcome_vars <- c("Height", "MeanCanWidth",	"RootDryMass", "Aboveground_Biomass",	"RoottoShootRatio")

envfit_PlantOutcomes_aligned <- envfit_PlantOutcomes[rownames(envfit_PlantOutcomes) %in% rownames(ord_scores), ]
identical(rownames(envfit_PlantOutcomes_aligned), rownames(ord_scores))


run_envfit_iterative <- function(ord_scores, envfit_PlantOutcomes_aligned, outcome_vars) {
  results <- list()
  vector_list <- list()
  
  for (var in outcome_vars) {
    cat("Now running:", var, "\n")
    
    if (!var %in% colnames(envfit_PlantOutcomes_aligned)) {
      warning(paste("Column not found in envfit_PlantOutcomes_aligned:", var))
      next
    }
    
    vec <- envfit_PlantOutcomes_aligned[[var]]
    keep <- !is.na(vec)
    
    if (sum(keep) >= 3) {
      ord_sub <- ord_scores[keep, , drop = FALSE]
      
      vec_transformed <- scale(log1p(vec[keep]))
      env_sub <- data.frame(var = vec_transformed)
      rownames(env_sub) <- rownames(ord_sub)
      
      fit <- envfit(ord_sub, env_sub, permutations = 999)
      
      if (is.null(fit$vectors) || is.null(scores(fit, display = "vectors"))) {
        warning(paste("envfit failed for", var))
        next
      }
      
      fit_scores <- scores(fit, display = "vectors")
      
      if (length(fit_scores) < 2) {
        warning(paste("fit_scores too short for", var))
        next
      }
      
      vec_coords <- data.frame(
        var = var,
        Dim1 = fit_scores[1],
        Dim2 = fit_scores[2],
        r2 = unname(fit$vectors$r),
        pval = unname(fit$vectors$pvals),
        n_samples = sum(keep),
        stringsAsFactors = FALSE
      )
      
      results[[var]] <- fit
      vector_list[[length(vector_list) + 1]] <- vec_coords
    } else {
      warning(paste("Too few non-NA values for", var))
    }
  }
  
  if (length(vector_list) == 0) {
    warning("No successful envfit results. Returning empty output.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors <- do.call(rbind, vector_list)
  cat("Number of successful vectors:", nrow(vectors), "\n")
  print(vectors)
  
  if (!"pval" %in% colnames(vectors)) {
    warning("p-values missing in vectors. Returning empty.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors$pval_adj <- p.adjust(vectors$pval, method = "BH")
  print("Adjusted p-values:")
  print(vectors$pval_adj)
  
  vectors_sig <- vectors %>% filter(pval_adj <= 0.05)
  
  return(list(fits = results, vectors = vectors_sig, summary = vectors))
}


fit_result <- run_envfit_iterative(ord_scores, envfit_PlantOutcomes_aligned, outcome_vars)
str(fit_result$summary)
fit_result$summary


library(ggplot2)
library(ggrepel)

sig.env.scrs <- fit_result$vectors
sig.env.scrs$NMDS1 <- sig.env.scrs$Dim1
sig.env.scrs$NMDS2 <- sig.env.scrs$Dim2
sig.env.scrs$env.variables <- sig.env.scrs$var

NMDS

NMMS_RootTips_Plant_Significant <- NMDS +
  geom_segment(data = sig.env.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'bisque4',
               inherit.aes = FALSE) +
  geom_text_repel(data = sig.env.scrs, 
                  aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  coord_fixed() 
  #theme(legend.position = "none") 

NMMS_RootTips_Plant_Significant

ggsave("NMMS_RootTips_Plant_Significant_wRootTips.pdf", NMMS_RootTips_Plant_Significant, width = 8, height = 6, units = "in")




## Plotting all vectors regardless of significance

all.env.scrs <- fit_result$summary
all.env.scrs$NMDS1 <- all.env.scrs$Dim1
all.env.scrs$NMDS2 <- all.env.scrs$Dim2
all.env.scrs$env.variables <- all.env.scrs$var

NMDS_RootTips_Plant_AllVectors_RegardlessofSig <- NMDS +
  geom_segment(data = all.env.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'darkgray',
               inherit.aes = FALSE) +
  geom_text_repel(data = all.env.scrs,
                  aes(x = NMDS1, y = NMDS2, label = env.variables),
                  colour = "darkgray", fontface = "italic") +
  coord_fixed()

NMDS_RootTips_Plant_AllVectors_RegardlessofSig

ggsave("NMDS_RootTips_Plant_AllVectors_RegardlessofSig.pdf", NMDS_RootTips_Plant_AllVectors_RegardlessofSig, width = 14, height = 10, units = "in")





#### NMDS w/ Envfit for Soil Chemistry #### 

## Envfit with how patchy and frequent my NAs are in the chemistry data make envfit difficult
## The way to do this may be iterating envfit, one variable at a time 
## Best option to maximize sample use for each outcome metric.

library(dplyr)
 
envfit_ChemistryEnvfit <- read.delim("ITS_metadata_RootTips_ChemistryEnvfit.txt",header = T)

rownames(envfit_ChemistryEnvfit) <- envfit_ChemistryEnvfit$SampleID

outcome_vars <- c("pH", "Total_C_percent", "Total_N_percent", "Na_waterextract_mgperL",
                  "NH4_waterextract_mgperL", "K_waterextract_mgperL", "Mg_waterextract_mgperL",
                  "Ca_waterextract_mgperL", "Cl_waterextract_mgperL", "NO3_waterextract_mgperL",
                  "PO4_waterextract_mgperL", "SO4_waterextract_mgperL")

envfit_ChemistryEnvfit_aligned <- envfit_ChemistryEnvfit[rownames(envfit_ChemistryEnvfit) %in% rownames(ord_scores), ]
identical(rownames(envfit_ChemistryEnvfit_aligned), rownames(ord_scores))


run_envfit_iterative <- function(ord_scores, envfit_ChemistryEnvfit_aligned, outcome_vars) {
  results <- list()
  vector_list <- list()
  
  for (var in outcome_vars) {
    cat("Now running:", var, "\n")
    
    if (!var %in% colnames(envfit_ChemistryEnvfit_aligned)) {
      warning(paste("Column not found in envfit_ChemistryEnvfit_aligned:", var))
      next
    }
    
    vec <- envfit_ChemistryEnvfit_aligned[[var]]
    keep <- !is.na(vec)
    
    if (sum(keep) >= 3) {
      ord_sub <- ord_scores[keep, , drop = FALSE]
      
      vec_transformed <- scale(log1p(vec[keep]))
      env_sub <- data.frame(var = vec_transformed)
      rownames(env_sub) <- rownames(ord_sub)
      
      fit <- envfit(ord_sub, env_sub, permutations = 999)
      
      if (is.null(fit$vectors) || is.null(scores(fit, display = "vectors"))) {
        warning(paste("envfit failed for", var))
        next
      }
      
      fit_scores <- scores(fit, display = "vectors")
      
      if (length(fit_scores) < 2) {
        warning(paste("fit_scores too short for", var))
        next
      }
      
      vec_coords <- data.frame(
        var = var,
        Dim1 = fit_scores[1],
        Dim2 = fit_scores[2],
        r2 = unname(fit$vectors$r),
        pval = unname(fit$vectors$pvals),
        n_samples = sum(keep),
        stringsAsFactors = FALSE
      )
      
      results[[var]] <- fit
      vector_list[[length(vector_list) + 1]] <- vec_coords
    } else {
      warning(paste("Too few non-NA values for", var))
    }
  }
  
  if (length(vector_list) == 0) {
    warning("No successful envfit results. Returning empty output.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors <- do.call(rbind, vector_list)
  cat("Number of successful vectors:", nrow(vectors), "\n")
  print(vectors)
  
  if (!"pval" %in% colnames(vectors)) {
    warning("p-values missing in vectors. Returning empty.")
    return(list(fits = results, vectors = data.frame(), summary = data.frame()))
  }
  
  vectors$pval_adj <- p.adjust(vectors$pval, method = "BH")
  print("Adjusted p-values:")
  print(vectors$pval_adj)
  
  vectors_sig <- vectors %>% filter(pval_adj <= 0.05)
  
  return(list(fits = results, vectors = vectors_sig, summary = vectors))
}


fit_result <- run_envfit_iterative(ord_scores, envfit_ChemistryEnvfit_aligned, outcome_vars)
str(fit_result$summary)
fit_result$summary


library(ggplot2)
library(ggrepel)

sig.env.scrs <- fit_result$vectors
sig.env.scrs$NMDS1 <- sig.env.scrs$Dim1
sig.env.scrs$NMDS2 <- sig.env.scrs$Dim2
sig.env.scrs$env.variables <- sig.env.scrs$var

NMDS


NMDS_RootTips_Chemistry_Significant <- NMDS +
  geom_segment(data = sig.env.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'bisque4',
               inherit.aes = FALSE) +
  geom_text_repel(data = sig.env.scrs, 
                  aes(x = NMDS1, y = NMDS2, label = env.variables), 
                  colour = "bisque4", fontface = "bold") +
  coord_fixed()

NMDS_RootTips_Chemistry_Significant

ggsave("NMDS_RootTips_Chemistry_Significant.pdf", NMDS_RootTips_Chemistry_Significant, width = 14, height = 10, units = "in")



## Plotting all vectors regardless of significance

all.env.scrs <- fit_result$summary
all.env.scrs$NMDS1 <- all.env.scrs$Dim1
all.env.scrs$NMDS2 <- all.env.scrs$Dim2
all.env.scrs$env.variables <- all.env.scrs$var

NMDS_RootTips_Chemistry_AllVectors_RegardlessofSig <- NMDS +
  geom_segment(data = all.env.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, 'cm')),
               color = 'darkgray',
               inherit.aes = FALSE) +
  geom_text_repel(data = all.env.scrs,
                  aes(x = NMDS1, y = NMDS2, label = env.variables),
                  colour = "darkgray", fontface = "italic") +
  coord_fixed()

NMDS_RootTips_Chemistry_AllVectors_RegardlessofSig

ggsave("NMDS_RootTips_Chemistry_AllVectors_RegardlessofSig.pdf", NMDS_RootTips_Chemistry_AllVectors_RegardlessofSig, width = 14, height = 10, units = "in")




#### EMF RA Boxplot ####

#### EMF RA by 6 treatments ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         EMF_RA = Ectomycorrhizal / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

EMF_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = EMF_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips EMF RA by Treatment",
       y = "EMF RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

EMF_RA_Greenhouse

ggsave("Greenhouse_EMF_RA.pdf", EMF_RA_Greenhouse, width = 3, height = 5, units = "in")

EMF_RA_Greenhouse <- EMF_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

EMF_RA_Greenhouse

ggsave("Greenhouse_EMF_RA_wtests.pdf", EMF_RA_Greenhouse, width = 5, height = 5, units = "in")


#### EMF RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         EMF_RA = Ectomycorrhizal / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

EMF_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = EMF_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips EMF RA by MPB Impact",
       y = "EMF RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

EMF_RA_Greenhouse_byMPB

ggsave("Greenhouse_EMF_RA_byMPB.pdf", EMF_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

EMF_RA_Greenhouse_byMPB <- EMF_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

EMF_RA_Greenhouse_byMPB

ggsave("Greenhouse_EMF_RA_byMPB_wtests.pdf", EMF_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "RootTips_EMF_RA.csv", row.names = FALSE)

#### Plant Pathogen RA by 6 treatments ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         PlantPathogen_RA = Plant.Pathogen / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

PlantPathogen_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = PlantPathogen_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips PlantPathogen RA by Treatment",
       y = "PlantPathogen RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

PlantPathogen_RA_Greenhouse

ggsave("Greenhouse_PlantPathogen_RA.pdf", PlantPathogen_RA_Greenhouse, width = 3, height = 5, units = "in")

PlantPathogen_RA_Greenhouse <- PlantPathogen_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

PlantPathogen_RA_Greenhouse

ggsave("Greenhouse_PlantPathogen_RA_wtests.pdf", PlantPathogen_RA_Greenhouse, width = 5, height = 5, units = "in")


#### Plant Pathogen RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         PlantPathogen_RA = Plant.Pathogen / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

PlantPathogen_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = PlantPathogen_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips PlantPathogen RA by MPB Impact",
       y = "PlantPathogen RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

PlantPathogen_RA_Greenhouse_byMPB

ggsave("Greenhouse_PlantPathogen_RA_byMPB.pdf", PlantPathogen_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

PlantPathogen_RA_Greenhouse_byMPB <- PlantPathogen_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

PlantPathogen_RA_Greenhouse_byMPB

ggsave("Greenhouse_PlantPathogen_RA_byMPB_wtests.pdf", PlantPathogen_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "RootTips_PlantPathogen_RA.csv", row.names = FALSE)



#### Fungal Parasite RA by 6 treatments ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         FungalParasite_RA = Fungal.Parasite / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

FungalParasite_RA_Greenhouse <- ggplot(fungal_guild_df, aes(x = Treatment_LandUseHistory, y = FungalParasite_RA)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips FungalParasite RA by Treatment",
       y = "FungalParasite RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

FungalParasite_RA_Greenhouse

ggsave("Greenhouse_FungalParasite_RA.pdf", FungalParasite_RA_Greenhouse, width = 3, height = 5, units = "in")

FungalParasite_RA_Greenhouse <- FungalParasite_RA_Greenhouse + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

FungalParasite_RA_Greenhouse

ggsave("Greenhouse_FungalParasite_RA_wtests.pdf", FungalParasite_RA_Greenhouse, width = 5, height = 5, units = "in")


#### Fungal Parasite RA by MPB Impact ####

setwd()

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         FungalParasite_RA = Fungal.Parasite / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Treatment_LandUseHistory <- factor(fungal_guild_df$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
fungal_guild_df$MPB_Impacted <- factor(fungal_guild_df$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

df$MPB_Impacted = factor(df$MPB_Impacted, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

FungalParasite_RA_Greenhouse_byMPB <- ggplot(fungal_guild_df, aes(x = MPB_Impacted, y = FungalParasite_RA)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips FungalParasite RA by MPB Impact",
       y = "FungalParasite RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle =  45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

FungalParasite_RA_Greenhouse_byMPB

ggsave("Greenhouse_FungalParasite_RA_byMPB.pdf", FungalParasite_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

FungalParasite_RA_Greenhouse_byMPB <- FungalParasite_RA_Greenhouse_byMPB + stat_compare_means(
  comparisons = comparisons,
  method = "wilcox.test",  # Non-parametric pairwise comparison
  p.adjust.method = "BH",
  label = "p.signif",
  hide.ns = TRUE,
  aes(group = Treatment),
  bracket.size = 0.5,
  size = 3)

FungalParasite_RA_Greenhouse_byMPB

ggsave("Greenhouse_FungalParasite_RA_byMPB_wtests.pdf", FungalParasite_RA_Greenhouse_byMPB, width = 2, height = 5, units = "in")

write.csv(fungal_guild_df, "RootTips_FungalParasite_RA.csv", row.names = FALSE)




#### Alpha Diversity ####

library(vegan) 
library(ggplot2) 
library(dplyr) 
library(ggpubr) 

asv_table <- t(otus)  
metadata <- map_file 


metadata$SampleID <- rownames(metadata)
asv_table <- as.data.frame(asv_table)
asv_table$SampleID <- rownames(asv_table)

merged <- inner_join(asv_table, metadata, by = "SampleID")

counts <- merged %>%
  select(where(is.numeric)) %>% select(-Replicate)

meta <- merged %>%
  select(where(~!is.numeric(.)))

alpha_div <- data.frame(
  Observed = rowSums(counts > 0),
  Shannon  = diversity(counts, index = "shannon")
)

alpha_div$SampleID <- merged$SampleID
alpha_div <- left_join(alpha_div, metadata, by = "SampleID")


alpha_div$Treatment_LandUseHistory <- factor(alpha_div$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
alpha_div$MPB_Impacted <- factor(alpha_div$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

alpha_div$Treatment_All <- Treatment_LandUseHistory
unique(alpha_div$Treatment_All)

alpha_div <- alpha_div %>%
  relocate(Observed, Shannon, .after = last_col())

write.csv(alpha_div, "ITS_alpha_div_RootTips_ObservedShannons.csv", row.names = FALSE)


## Figures 

## Observed ASVs

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_RootTips_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Observed)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips Observed Species by Treatment",
       y = "Observed Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

AD_RootTips_ObservedSpecies_Rarefied

ggsave("AD_RootTips_ObservedSpecies_Rarefied.pdf", AD_RootTips_ObservedSpecies_Rarefied, width = 3, height = 5, units = "in")


treatment_levels <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_RootTips_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = MPB_Impacted, y = Observed)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips Observed Species by MPB Impact",
       y = "Observed Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impacted),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_RootTips_ObservedSpecies_Rarefied

ggsave("AD_RootTips_ObservedSpecies_Rarefied_byMPB_wtests.pdf", AD_RootTips_ObservedSpecies_Rarefied, width = 2, height = 5, units = "in")



## Shannons 

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_RootTips_Shannons_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Shannon)) +
  geom_boxplot(aes(fill = Treatment_LandUseHistory), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_LandUseHistory), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips Shannon Diversity by Treatment",
       y = "Shannon Index") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # stat_compare_means(
  #   aes(group = Treatment_LandUseHistory),
  #   comparisons = comparisons,
  #   method = "wilcox.test",  # Non-parametric pairwise comparison
  #   p.adjust.method = "BH",
  #   label = "p.signif",
  #   hide.ns = TRUE,
  #   bracket.size = 0.5,
  #   size = 3) +
  theme(legend.position = "none")

AD_RootTips_Shannons_Rarefied

ggsave("AD_RootTips_Shannons_Rarefied.pdf", AD_RootTips_Shannons_Rarefied, width = 3, height = 5, units = "in")


treatment_levels <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

AD_RootTips_Shannons_Rarefied <- ggplot(alpha_div, aes(x = MPB_Impacted, y = Shannon)) +
  geom_boxplot(aes(fill = MPB_Impacted), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impacted), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "RootTips Shannon Diversity by MPB Impact",
       y = "Shannon Index") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impacted),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_RootTips_Shannons_Rarefied

ggsave("AD_RootTips_Shannons_Rarefied_byMPB_wtests.pdf", AD_RootTips_Shannons_Rarefied, width = 2, height = 5, units = "in")


## Final plots together 

library(ggpubr)

final_plot <- ggarrange(
  AD_RootTips_ObservedSpecies_Rarefied, AD_RootTips_Shannons_Rarefied,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)

print(final_plot)

ggsave("alphadiversity_boxplots_RootTips_rarefied_tests.pdf", final_plot, width = 14, height = 5, units = "in")



#### FUNGuild Barplot ####

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.delim("envfit_FUNGuild_pivot_FeatureTable_wMetadata_ITS_RootTips.txt")


group_var <- "Treatment_LandUseHistory"
facet_row <- "." 
facet_col <- "." 
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Guild",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$MPB_Impacted <- factor(df_long$MPB_Impacted, levels = c("Non_MPB", "MPB")) 

unique(df_long$Guild)
guild_levels <- c(
  "Animal.Parasite", "Animal.Pathogen", "Arbuscular.Mycorrhizal",
  "Dung.Saprotroph", "Ectomycorrhizal", "Endophyte",          
  "Epiphyte", "Fungal.Parasite", "Lichen.Parasite",     
  "Plant.Pathogen", "Plant.Saprotroph", "Soil.Saprotroph",    
  "Undefined.Saprotroph", "Wood.Saprotroph", "Unknown"   
)
df_long$Guild <- factor(df_long$Guild, levels = guild_levels)

grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Guild")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()

final_plot <- ggplot(df_relabund, aes_string(x = group_var, y = "RelAbund", fill = "Guild")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_RootTips.pdf", final_plot, width = 14, height = 10, units = "in")








