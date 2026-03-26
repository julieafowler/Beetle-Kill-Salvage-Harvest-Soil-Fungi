#### SFSP - Fungal - ITS - Analysis ####

setwd()


#### Package Loading ####

library(phyloseq) 
library(ggplot2)
library(grid)
library(plyr)
library(vegan)
library(Hmisc)
library(reshape2)
library(ggpubr)
library(skimr)
library(ggthemr) 
library(vegan)


#### Rarefying ####

setwd()
ITS_Field_alldata <- read.delim("ITS-featuretable-updated.txt",header=T,row.names=1, check.names=FALSE)
ITS_Field_alldata_t <- t(ITS_Field_alldata)
sample_sums <- rowSums(ITS_Field_alldata_t)
summary(sample_sums)
rarecurve(ITS_Field_alldata_t, step = 500, label = FALSE)
which(sample_sums < 14000)
sum(sample_sums < 14000)
ITS_Field_alldata_t_filtered <- ITS_Field_alldata_t[rowSums(ITS_Field_alldata_t) >= 14000, ]
ITS_Field_alldata_t_rarefied_14000 <- rrarefy(ITS_Field_alldata_t_filtered, sample = 14000)
ITS_Field_alldata_t_rarefied_14000 <- t(ITS_Field_alldata_t_rarefied_14000)


#### FUNGuildR - Updating FUNGuild 09/25/25 to the new database ####

## Making a dataframe with the FUNGuild data using FUNGuildR to get the most up to date database 

install.packages("devtools")
devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)
setwd()

fung <- get_funguild_db() ## Download the most up-to-date version of FUNGuild and save
sample_fungi <- read.delim("Field_FUNGuildR_Input_Updated.txt")
fung_guilds <- funguild_assign(sample_fungi, db = fung)
write.table(fung_guilds, file = "field_fung_guilds_updated.txt", sep = "\t", row.names = FALSE)



#### Import Files and Create Object ####

setwd()

otus<-ITS_Field_alldata_t_rarefied_14000
otus <- as.data.frame(otus)
sample_names <- colnames(otus)
map_file<-read.delim("ITS-map-file.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)

all.equal(names(otus),row.names(map_file))

## Ordering variable for ggplot
map_file$Decade <- factor(map_file$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")) 
map_file$Soil_Depth <- factor(map_file$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm")) 

otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("ITS_taxonomy_updated_edited.txt",header=T,row.names=1)
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

NMDS <- ggplot(sample_tab)+
  geom_point(aes(x=NMDS1, y=NMDS2, color=Decade, shape=Soil_Depth), size=4)+
  theme(text=element_text(size = 20)) +
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = Decade, color=Decade),linetype = 2) +
  #ylim(-0.4, 0.4) + xlim(-0.6,0.6) 
  theme(legend.position = "none")
##theme(legend.position = "none")

## 10 x 14 inches 
NMDS

ggsave("ITS_NMDS_Field_Rarefied.pdf", NMDS, width = 10, height = 8, units = "in")



#### PERMANOVA ####

mgd_ge5K_relabund.bray
sample_tab

adonis2(mgd_ge5K_relabund.bray ~ Decade, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Soil_Depth, data = sample_tab, permutations = 999, method = "bray", na.action = na.omit)
adonis2(mgd_ge5K_relabund.bray ~ MPB, data = sample_tab, permutations = 999, method = "bray")

adonis2(mgd_ge5K_relabund.bray ~ Decade * Soil_Depth, data = sample_tab, permutations = 999, method = "bray", by = "margin", na.action = na.omit)

adonis2(mgd_ge5K_relabund.bray ~ MPB * Soil_Depth, data = sample_tab, permutations = 999, method = "bray", by = "margin", na.action = na.omit)

## P-value adjustments:
library(stats)
p = c(0.001, 0.001, 0.001, 0.999, 0.229) ## Change with above findings 
p.adjusted = p.adjust(p, method = "BH")
p.adjusted


library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Decade, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)





#### NMDS w/ Envfit for FUNGuild ####

envfit_FUNGuild <- read.delim("ITS-featuretable-FUNGuild-Envfit-updated.txt",header = T)
sample_names <- colnames(otus)
guild_column <- colnames(envfit_FUNGuild)[1]
matching_samples <- intersect(sample_names, colnames(envfit_FUNGuild))
envfit_FUNGuild <- envfit_FUNGuild[, c(guild_column, matching_samples)]

envfit_FUNGuild_agg <- aggregate(. ~  Guild, data = envfit_FUNGuild, sum)

envfit_FUNGuild_pivot <- t(envfit_FUNGuild_agg)


library(tibble)

# Step 1: Make the first row the column names
colnames(envfit_FUNGuild_pivot) <- as.character(envfit_FUNGuild_pivot[1, ])
envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[-1, ]

# Step 2: Check that row names match exactly between envfit_FUNGuild_pivot and metadata
all.equal(row.names(envfit_FUNGuild_pivot),row.names(map_file))
# Do the row names match AND are in the same order?
identical(row.names(envfit_FUNGuild_pivot), row.names(map_file))  # TRUE means perfectly matched and ordered

# Step 3: Remove "Unknown" column if it exists
if ("Unknown" %in% colnames(envfit_FUNGuild_pivot)) {
  envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[, colnames(envfit_FUNGuild_pivot) != "Unknown"]
}

envfit_FUNGuild_pivot <- as.data.frame(envfit_FUNGuild_pivot)
# Step 4: Correct skew & calculate z scores through scaling 
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
identical(nmds_samples, envfit_samples) ## Should be true before running 


## NA check
any(!is.finite(scores(mgd_ge5K_relabund.bray.nmds)))
which(!is.finite(scores(mgd_ge5K_relabund.bray.nmds)), arr.ind = TRUE)

any(!is.finite(as.matrix(envfit_FUNGuild_pivot_log_scaled)))
which(!is.finite(as.matrix(envfit_FUNGuild_pivot_log_scaled)), arr.ind = TRUE)




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
  theme(legend.position = "none") +
  ylim(-0.5, 0.5) + xlim(-0.4,0.8)

## 10 x 10 inches 
NMDS_Envfit_FUNGuild

ggsave("Field_NMDS_Envfit_FUNGuild_Rarefied_UPDATEDSCALEDVECTORS.pdf", NMDS_Envfit_FUNGuild, width = 10, height = 8, units = "in")




#### FUNGuild Barplot: Ectomycorrhizal fungi by soil depth and treatment ####

setwd()

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(tidyverse)

envfit_FUNGuild <- read.delim("ITS-featuretable-FUNGuild-Envfit-updated.txt", header = TRUE)
envfit_FUNGuild_agg <- aggregate(. ~ Guild, data = envfit_FUNGuild, sum)
envfit_FUNGuild_pivot <- t(envfit_FUNGuild_agg)

# Export for use in ECM RA boxplot later 
write.csv(envfit_FUNGuild_pivot, "ITS_field_envfit_FUNGuild_pivot.csv", row.names = TRUE)

colnames(envfit_FUNGuild_pivot) <- envfit_FUNGuild_pivot[1, ]
envfit_FUNGuild_pivot <- envfit_FUNGuild_pivot[-1, ]

envfit_FUNGuild_pivot <- as.data.frame(envfit_FUNGuild_pivot)
envfit_FUNGuild_pivot$SampleID <- rownames(envfit_FUNGuild_pivot)

# Load metadata
map_file <- read.delim("ITS-map-file.txt", header = TRUE, row.names = 1, check.names = FALSE)
map_file$SampleID <- rownames(map_file)

identical(envfit_FUNGuild_pivot$SampleID, map_file$SampleID)

# Merge metadata
df_joined <- merge(envfit_FUNGuild_pivot, map_file[, c("SampleID", "Decade", "Soil_Depth")], by = "SampleID")

df_long <- df_joined %>%
  pivot_longer(cols = -c(SampleID, Decade, Soil_Depth),
               names_to = "Guild",
               values_to = "Count") %>%
  mutate(Count = as.numeric(Count))

df_long_filtered <- df_long %>%
  group_by(SampleID) %>%
  filter(sum(Count, na.rm = TRUE) > 0) %>%
  ungroup()

## Calculate relative abundance and average by Decade × Soil_Depth × Guild
## If this doesn't work, it is because I've called plyr instead of dplyr 
df_rel_abund <- df_long_filtered %>%
  group_by(SampleID, Decade, Soil_Depth) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup() %>%
  group_by(Decade, Soil_Depth, Guild) %>%
  summarise(value = mean(RelAbund, na.rm = TRUE), .groups = "drop")



# Set the desired order for Soil_Depth and Site
df_rel_abund$Soil_Depth <- factor(df_rel_abund$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm"))
df_rel_abund$Decade <- factor(df_rel_abund$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts"))


# Create the plot with horizontal facets
final_plot <- ggplot(df_rel_abund, aes(x = Decade, y = value, fill = Guild)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Soil_Depth, nrow = 1) +
  labs(x = "Land Use History", y = "Relative Abundance", fill = "Fungal Guild") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("Field_FUNGuild_Barplot_Updated.pdf", final_plot, width = 14, height = 10, units = "in")



#### EMF RA Boxplot ####

df <- read.delim("ITS_field_envfit_FUNGuild_pivot_forEMFRABoxplot.txt")

df <- na.omit(df)

library(ggplot2)
library(dplyr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         EMF_RA = Ectomycorrhizal / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Decade <- factor(fungal_guild_df$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")) 
fungal_guild_df$Soil_Depth <- factor(fungal_guild_df$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm")) 
fungal_guild_df$MPB_Impact <- factor(fungal_guild_df$MPB_Impact, levels = c("Non_MPB", "MPB")) 

treatment_levels <- c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

fungal_guild_df$MPB_Impact = factor(fungal_guild_df$MPB_Impact, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

fungal_guild_df <- fungal_guild_df %>%
  filter(Soil_Depth != "OrgC")

EMF_RA_Field <- ggplot(fungal_guild_df, aes(x = MPB_Impact, y = EMF_RA)) +
  facet_grid(~Soil_Depth) +
  geom_boxplot(aes(fill = MPB_Impact), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impact), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = " EMF RA by Treatment",
       y = "EMF RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impact),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

EMF_RA_Field

ggsave("Field_EMF_RA_Boxplots_byMPB_wtests.pdf", EMF_RA_Field, width = 7, height = 3, units = "in")



#### Plant Pathogen Boxplot #### 

df <- read.delim("ITS_field_envfit_FUNGuild_pivot_forEMFRABoxplot.txt")

df <- na.omit(df)

library(ggplot2)
library(dplyr)
library(ggpubr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         PlantPath_RA = Plant.Pathogen / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Decade <- factor(fungal_guild_df$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")) 
fungal_guild_df$Soil_Depth <- factor(fungal_guild_df$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm")) 

treatment_levels <- c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

fungal_guild_df$MPB_Impact = factor(fungal_guild_df$MPB_Impact, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

fungal_guild_df <- fungal_guild_df %>%
  filter(Soil_Depth != "OrgC")

PlantPath_RA_Field <- ggplot(fungal_guild_df, aes(x = MPB_Impact, y = PlantPath_RA)) +
  facet_grid(~Soil_Depth) +
  geom_boxplot(aes(fill = MPB_Impact), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impact), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = " PlantPath RA by Treatment",
       y = "PlantPath RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impact),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

PlantPath_RA_Field

ggsave("Field_PlantPath_RA_Boxplots_MPBImpact.pdf", PlantPath_RA_Field, width = 7, height = 3, units = "in")



#### Endophyte RA Boxplots ####

df <- read.delim("ITS_field_envfit_FUNGuild_pivot_forEMFRABoxplot.txt")

df <- na.omit(df)

library(ggplot2)
library(dplyr)
library(ggpubr)

fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         Endophyte_RA = Endophyte / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Decade <- factor(fungal_guild_df$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")) 
fungal_guild_df$Soil_Depth <- factor(fungal_guild_df$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm")) 

treatment_levels <- c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

fungal_guild_df$MPB_Impact = factor(fungal_guild_df$MPB_Impact, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

fungal_guild_df <- fungal_guild_df %>%
  filter(Soil_Depth != "OrgC")

Endophyte_RA_Field <- ggplot(fungal_guild_df, aes(x = MPB_Impact, y = Endophyte_RA)) +
  facet_grid(~Soil_Depth) +
  geom_boxplot(aes(fill = MPB_Impact), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impact), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = " Endophyte RA by Treatment",
       y = "Endophyte RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impact),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

Endophyte_RA_Field

ggsave("Field_Endophyte_RA_Boxplots.pdf", Endophyte_RA_Field, width = 7, height = 5, units = "in")



#### Fungal Parasite RA Boxplots ####

df <- read.delim("ITS_field_envfit_FUNGuild_pivot_forEMFRABoxplot.txt")

df <- na.omit(df)

library(ggplot2)
library(dplyr)
library(ggpubr)


fungal_guild_df <- df %>%
  rowwise() %>%
  mutate(TotalGuildCount = sum(c_across(c(Animal.Parasite:Wood.Saprotroph)), na.rm = TRUE),
         Fungal.Parasite_RA = Fungal.Parasite / TotalGuildCount) %>%
  ungroup()

fungal_guild_df$Decade <- factor(fungal_guild_df$Decade, levels = c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")) 
fungal_guild_df$Soil_Depth <- factor(fungal_guild_df$Soil_Depth, levels = c("OrgC", "0-5cm", "5-15cm")) 

treatment_levels <- c("Old Forest", "1st PreMPB_Old", "2nd PreMPB_NEW", "Early MPB", "Late MPB", "Recent Cuts")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

fungal_guild_df$MPB_Impact = factor(fungal_guild_df$MPB_Impact, levels = c("Non_MPB", "MPB"))
treatment_levels_MPB <- c("Non_MPB", "MPB")
comparisons <- combn(treatment_levels_MPB, 2, simplify = FALSE)

fungal_guild_df <- fungal_guild_df %>%
  filter(Soil_Depth != "OrgC")

Fungal.Parasite_RA_Field <- ggplot(fungal_guild_df, aes(x = MPB_Impact, y = Fungal.Parasite_RA)) +
  facet_grid(~Soil_Depth) +
  geom_boxplot(aes(fill = MPB_Impact), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = MPB_Impact), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = " Fungal.Parasite RA by Treatment",
       y = "Fungal.Parasite RA (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = MPB_Impact),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

Fungal.Parasite_RA_Field

ggsave("Field_Fungal.Parasite_RA_Boxplots.pdf", Fungal.Parasite_RA_Field, width = 7, height = 3, units = "in")



#### RA boxplots ggarrange ####

## Combining EMF RA, Plant Path RA, and Fungal Parasite RA into a figure 

library(ggpubr)

boxplots <- ggarrange(EMF_RA_Field, PlantPath_RA_Field, Fungal.Parasite_RA_Field, ncol = 1, nrow = 3)
ggsave("EMF_Path_Parasite_Boxplots.pdf", boxplots, width = 7, height = 10, units = "in")



#### Maaslin3 ####

library(maaslin3)
library(EnhancedVolcano)
library(dplyr)
library(readr)

#### Maaslin3 - ASVs ####

# --- Load ---
feat <- read.delim("ITS-featuretable-updated.txt", row.names = 1, check.names = FALSE)
meta <- read.delim("ITS-map-file.txt", row.names = 1, check.names = FALSE)

if (!all(rownames(meta) %in% colnames(feat))) stop("Sample IDs don't line up")
feat <- t(feat)

# Outcomes and covariates
outs <- c("MPB")
meta$Soil_Depth  <- factor(meta$Soil_Depth, levels = c("OrgC","0-5cm", "5-15cm"))
meta$Decade <- factor(meta$Decade, levels = c("Old Forest", "1st PreMPB_OLD", "2nd PreMPB_NEW","Early MPB", "Late MPB","Recent Cuts"))
meta$MPB  <- factor(meta$MPB, levels = c("Non_MPB", "MPB"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(feat)

meta$Decade <- relevel(factor(meta$Decade ), ref = "Old Forest")
meta$MPB <- relevel(factor(meta$MPB), ref = "Non_MPB")


## 1) MPB
this_meta <- subset(meta, !is.na(MPB))
common <- intersect(rownames(this_meta), rownames(feat))
maaslin3(
  input_data      = feat[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_ASVs_MPB",
  formula         = ~ MPB + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaasLin3 results for this outcome
res <- read_tsv("maaslin3_ASVs_MPB/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model to visualize
metadatum <- "MPB"   # <-- edit if needed
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,
    q_plot = qval_individual
  )

tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
ASVs_MPB_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05, 
  FCcutoff = log2(1.5), 
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

ASVs_MPB_volcano_abundance

ggsave("ASVs_MPB_volcano_abundance.pdf", ASVs_MPB_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2),
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

ASVs_MPB_volcano_prevalence <- EnhancedVolcano(
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

ASVs_MPB_volcano_prevalence

ggsave("ASVs_MPB_volcano_prevalence.pdf", ASVs_MPB_volcano_prevalence, width = 8, height = 7, units = "in")



#### Maaslin3 - Guilds #### 

# --- Load ---
guild <- read.delim("ITS-featuretable-FUNGuild-Envfit-updated.txt", check.names = FALSE)
meta <- read.delim("ITS-map-file.txt", row.names = 1, check.names = FALSE)

sample_cols <- names(guild)[-1]
guild[ sample_cols ] <- lapply(guild[ sample_cols ], function(x) as.numeric(as.character(x)))
guild_sum <- rowsum(as.matrix(guild[ sample_cols ]), group = guild$Guild, na.rm = TRUE)
guild_sum_df <- data.frame(Guild = rownames(guild_sum), guild_sum, row.names = NULL)
guild <- guild_sum_df
rownames(guild) <- guild[, 1]
guild <- guild[, -1]


if (!all(rownames(meta) %in% colnames(guild))) stop("Sample IDs don't line up")
guild <- t(guild)

# Outcomes and covariates
outs <- c("MPB")
meta$Soil_Depth  <- factor(meta$Soil_Depth, levels = c("OrgC","0-5cm", "5-15cm"))
meta$Decade <- factor(meta$Decade, levels = c("Old Forest", "1st PreMPB_OLD", "2nd PreMPB_NEW","Early MPB", "Late MPB","Recent Cuts"))
meta$MPB  <- factor(meta$MPB, levels = c("Non_MPB", "MPB"))

# If started with counts, include read depth to control prevalence detection bias
meta$ReadDepth <- rowSums(guild)

meta$Decade <- relevel(factor(meta$Decade ), ref = "Old Forest")
meta$MPB <- relevel(factor(meta$MPB), ref = "Non_MPB")


## 1) MPB
this_meta <- subset(meta, !is.na(MPB))
common <- intersect(rownames(this_meta), rownames(guild))
maaslin3(
  input_data      = guild[common, , drop = FALSE],
  input_metadata  = this_meta[common, , drop = FALSE],
  output          = "maaslin3_guild_MPB",
  formula         = ~ MPB + ReadDepth,
  normalization   = "TSS", transform = "LOG", standardize = TRUE,
  min_prevalence  = 0.10, median_comparison_abundance = TRUE
)


## EnhancedVolcano plot

# 1) Load MaAsLin3 results for this outcome
res <- read_tsv("maaslin3_guild_MPB/all_results.tsv", show_col_types = FALSE)

## Feature abundance 

# 2) Keep the predictor and model you want to visualize
metadatum <- "MPB" 
df <- res %>%
  filter(metadata == metadatum, model == "abundance") %>%
  transmute(
    feature,
    log2FC = coef,
    q_plot = qval_individual
  )

tt <- as.data.frame(df)
rownames(tt) <- tt$feature

# 3) Plot
guild_MPB_volcano_abundance <-  EnhancedVolcano(
  tt,
  lab = rownames(tt),
  x = "log2FC",
  y = "q_plot",
  pCutoff = 0.05, 
  FCcutoff = log2(1.5), 
  title = paste0("Volcano: ", metadatum, " (abundance)"),
  subtitle = "Effect = log2 fold-change per 1 SD in outcome (standardize=TRUE)",
  ylab = expression(-log[10]("q")),
  pointSize = 1.6,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colAlpha = 0.9
)

guild_MPB_volcano_abundance

ggsave("guild_MPB_volcano_abundance.pdf", guild_MPB_volcano_abundance, width = 8, height = 7, units = "in")


## Feature prevalence 

df_prev <- res %>%
  filter(metadata == metadatum, model == "prevalence") %>%
  transmute(
    feature,
    log2FC = coef / log(2), 
    q_plot = qval_individual
  )
tt_prev <- as.data.frame(df_prev); rownames(tt_prev) <- tt_prev$feature

guild_MPB_volcano_prevalence <- EnhancedVolcano(
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

guild_MPB_volcano_prevalence

ggsave("guild_MPB_volcano_prevalence.pdf", guild_MPB_volcano_prevalence, width = 8, height = 7, units = "in")




