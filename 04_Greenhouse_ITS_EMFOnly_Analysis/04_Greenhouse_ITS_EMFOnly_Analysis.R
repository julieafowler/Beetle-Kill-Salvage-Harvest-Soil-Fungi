#### SFSP Greenhouse Bioassays - ITS - EMF Only Analyses ####

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


## Not rarefying for ECM only



#### Rhizosphere - EMF Community ####

#### NMDS Prep ####

## Using relative abundances for input into NMDS owing to lack of rarefying with this dataset 

setwd()

otus <- read.delim("ITS_FeatureTable_Rhizosphere_EMFOnly.txt",header=T,row.names=1, check.names=FALSE)
otus_t <- t(otus)
otus_rel_abund_t <- decostand(otus_t, method = "total")
otus_rel_abund <- t(otus_rel_abund_t)
otus <- otus_rel_abund

otus <- as.data.frame(otus)
sample_names <- colnames(otus)

# Filter metadata to only samples present in rarefied count data
map_file <- read.delim("ITS_metadata_Rhizosphere.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)

all.equal(names(otus),row.names(map_file))

map_file$Treatment_LandUseHistory <- factor(map_file$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
map_file$Salvage_Harvest_Status <- factor(map_file$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("ITS_taxonomy_Rhizosphere.txt",header=T,row.names=1)
taxmat<-as.matrix(taxa)
taxmat <- taxmat[rownames(taxmat) %in% rownames(otumat), ]
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
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = Salvage_Harvest_Status, color=Salvage_Harvest_Status),linetype = 2) 

NMDS

ggsave("NMDS_Rhizosphere_EMFOnly.pdf", NMDS, width = 14, height = 10, units = "in")


#### PERMANOVA ####

adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Salvage_Harvest_Status, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory*Salvage_Harvest_Status, data = sample_tab, permutations = 999, method = "bray")


## P-value adjustments:
library(stats)
p = c(0.001, 0.001) ## Change with above findings 
p.adjusted = p.adjust(p, method = "BH")
p.adjusted

# Adjusted: 0.001666667 0.001666667 0.002500000 0.001666667 0.003000000

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)



#### Alpha Diversity ####

library(vegan) 
library(ggplot2)
library(dplyr)
library(ggpubr)

asv_table <- t(otus)   # samples x ASVs 
metadata <- map_file    # samples x variables

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
alpha_div$Salvage_Harvest_Status <- factor(alpha_div$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

alpha_div$Treatment_All <- alpha_div$Treatment_LandUseHistory
unique(alpha_div$Treatment_All)

treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

treatment_levels <- c("Salvage_Harvested", "Non_Salvage_Harvested")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)

alpha_div <- alpha_div %>%
  relocate(Observed, Shannon, .after = last_col())

write.csv(alpha_div, "ITS_alpha_div_rhizosphere_EMFOnly_ObservedShannons.csv", row.names = FALSE)


## Figures 

## Observed ASVs

AD_Rhizosphere_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = Salvage_Harvest_Status, y = Observed)) +
  geom_boxplot(aes(fill = Salvage_Harvest_Status), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Salvage_Harvest_Status), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Observed Species by Treatment",
       y = "Observed Species") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = Salvage_Harvest_Status),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_ObservedSpecies_Rarefied

ggsave("AD_Rhizosphere_EMFOnly_ObservedSpecies_Rarefied_bySH.pdf", AD_Rhizosphere_ObservedSpecies_Rarefied, width = 2, height = 5, units = "in")



## Shannons 

AD_Rhizosphere_Shannons_Rarefied <- ggplot(alpha_div, aes(x = Salvage_Harvest_Status, y = Shannon)) +
  geom_boxplot(aes(fill = Salvage_Harvest_Status), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Salvage_Harvest_Status), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Rhizosphere Shannon Diversity by Treatment",
       y = "Shannon Index") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    aes(group = Salvage_Harvest_Status),
    comparisons = comparisons,
    method = "wilcox.test",  # Non-parametric pairwise comparison
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    bracket.size = 0.5,
    size = 3) +
  theme(legend.position = "none")

AD_Rhizosphere_Shannons_Rarefied

ggsave("AD_Rhizosphere_EMFOnly_Shannons_Rarefied_bySH.pdf", AD_Rhizosphere_Shannons_Rarefied, width = 2, height = 5, units = "in")

## Final plots together 

library(ggpubr)

final_plot <- ggarrange(
  AD_Rhizosphere_ObservedSpecies_Rarefied, AD_Rhizosphere_Shannons_Rarefied,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)

print(final_plot)

ggsave("alphadiversity_boxplots_Rhizosphere_EMFOnly_rarefied_tests_bySH.pdf", final_plot, width = 14, height = 5, units = "in")




#### Genus/Species Barplots - All Taxa ####

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.delim("ITS_GenusTable_Rhizosphere.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),          
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body))  
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_Rhizosphere_Genus_Summed.csv", row.names = TRUE)

## Prepping table - doing things in excel 
  
df <- read.delim("ITS_SpeciesTable_Rhizosphere.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),             
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body)) 
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_Rhizosphere_Species_Summed.csv", row.names = TRUE)

## Species Barplot
  
df <- read.delim("ITS_Rhizosphere_SpeciesBarplot_Summed.txt")


group_var <- "Treatment_LandUseHistory"  
facet_row <- "." 
facet_col <- "."    
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Species",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Species)
guild_levels <- c(
  "Acephala_sp","Agaricus_agrinferus","Agaricus_arvensis","Agaricus_zhangyensis",
  "Alternaria_ellipsoidea","Alternaria_multirostrata","Alternaria_nepalensis",
  "Amanita_pseudopantherina","Antrodiella_faginea","Apiotrichum_sporotrichoides",
  "Apseudocercosporella_trigonotidis","Ascochyta_medicaginicola_var._macrospora",
  "Aspergillus_sigurros","Aureobasidium_castaneae","Auricularia_subglabra",
  "Auricularia_submesenterica","Barrenia_panici","Basidioascus_undulatus",
  "Beauveria_varroae","Bergerella_atrofusca","Bovista_cretacea","Cadophora_domestica",
  "Calonarius_caesiocinctus","Calvatia_fragilis","Candida_orthopsilosis",
  "Candida_pseudoglaebosa","Candolleomyces_cacao","Ceratobasidiaceae_sp",
  "Chaetomium_cervicicola","Chaetospermum_sp","Chalara_recta","Chroogomphus_pakistanicus",
  "Chrysozyma_pseudogriseoflava","Cladophialophora_humicola","Cladosporium_aphidis",
  "Clarireedia_paspali","Clitopilus_brunneiceps","Clonostachys_divergens",
  "Coleophoma_camelliae","Coniochaeta_vineae","Coprinellus_sclerocystidiosus",
  "Coprinopsis_marcescibilis","Cortinarius_albidipes","Cortinarius_alboadustus",
  "Cortinarius_brunneoaffinis","Cortinarius_cepistipes","Cortinarius_clarisordidus",
  "Cortinarius_collangustus","Cortinarius_deceptivissimus","Cortinarius_decipiens",
  "Cortinarius_kauffmanianus","Cortinarius_majorinus","Cortinarius_miwok",
  "Cortinarius_mucosus","Cortinarius_parvannulatus","Cortinarius_paululus",
  "Cortinarius_subcarneinatus","Cortinarius_turgidipes","Cristulariella_depraedans",
  "Curvularia_siddiquii","Cutaneotrichosporon_dermatis","Cyphellophora_eucalypti",
  "Cyphellophora_vermispora","Cystobasidium_raffinophilum","Cystobasidium_terricola",
  "Darksidea_zeta","Debaryomyces_robertsiae","Didymosphaeria_variabile",
  "Diversispora_spurca","Dothiora_europaea","Emericellopsis_glabra",
  "Entoloma_aurorae.borealis","Entoloma_fuscohebes","Entoloma_tibiicystidiatum",
  "Entomocorticium_portiae","Epicoccum_thailandicum","Erysiphe_polygoni",
  "Filobasidium_magnum","Filobasidium_oeirense","Fonsecaea_nubica",
  "Fontanospora_fusiramosa","Funneliformis_caledonium","Furcasterigmium_furcatum",
  "Fusarium_algeriense","Fusarium_brevicatenulatum","Fusarium_croci","Fusarium_nelsonii",
  "Gaeumannomyces_hyphopodioides","Ganoderma_carocalcareum","Gliomastix_roseogrisea",
  "Gloeopycnis_protuberans","Gymnopilus_aurantiophyllus","Gyoerffyella_sp",
  "Hanseniaspora_nectarophila","Harzia_acremonioides","Hebeloma_erebium",
  "Hebeloma_pungens","Helotiales_sp","Humicola_zollerniae","Hyaloscyphaceae_sp",
  "Hydnum_pallidomarginatum","Hydnum_subolympicum","Hygrophorus_fuscoalboides",
  "Hypoderma_siculum","Ilyonectria_mors.panacis","Inocybe_purpureobadia",
  "Inocybe_silvae.herbaceae","Juncaceicola_padellana","Knufia_tsunedae",
  "Laccaria_fagacicola","Lachnum_impudicum","Lactarius_aurantiacopallens",
  "Lactarius_miniatosporus","Lentithecium_pseudoclioninum","Leohumicola_atra",
  "Leucosporidium_escuderoi","Linnemannia_amoeboidea","Linnemannia_camargensis",
  "Lophiostoma_rosae","Lycoperdon_echinatum","Lycoperdon_ericaeum","Lycoperdon_nigrescens",
  "Lyophyllum_shimeji","Macroconia_bulbipes","Mallocybe_sp","Mariannaea_camptospora",
  "Megalocystidium_perticatum","Metapochonia_hahajimaensis","Microbotryum_violaceum",
  "Microdochium_lycopodinum","Minimelanolocus_clavatus","Morchella_tasmanica",
  "Mortierella_alliacea","Mortierella_fluviae","Mortierella_wuyishanensis",
  "Mrakia_niccombsii","Mycena_megaspora","Mycodidymella_aesculi","Naganishia_bhutanensis",
  "Naganishia_diffluens","Naganishia_vishniacii","Neoantrodia_infirma","Neoantrodia_serialis",
  "Neoascochyta_tardicrescens","Neocucurbitaria_aetnensis","Neocucurbitaria_rhamnioides",
  "Neodidymelliopsis_negundinis","Neofabraea_salicina","Neophaeomoniella_constricta",
  "Nodulosphaeria_thalictri","Oculimacula_yallundae","Odontia_sp","Panaeolina_foenisecii",
  "Papiliotrema_nemorosa","Parabartalinia_lateralis","Paraboeremia_selaginellae",
  "Paraleptosphaeria_nitschkei","Paraphaeosphaeria_neglecta","Paraphoma_ledniceana",
  "Paraphoma_pye","Penicillium_alutaceum","Penicillium_coprobium","Penicillium_janczewskii",
  "Penicillium_madriti","Penicillium_oregonense","Penicillium_vasconiae",
  "Peziza_montirivicola","Phenoliferia_psychrophenolica","Phialocephala_amethystea",
  "Phialocephala_helenae","Phlebia_rufa","Phlebiopsis_sp","Phlegmacium_subfoetens",
  "Piloderma_bicolor","Piloderma_lanatum","Piloderma_sp","Plectosphaerella_melonis",
  "Pluteus_kovalenkoi","Polyozellus_mariae","Psathyrella_carinthiaca","Pseudeurotium_desertorum",
  "Pseudofusicoccum_olivaceum","Pseudogymnoascus_roseus","Pseudopezicula_betulae",
  "Pseudosperma_emberizanum","Pseudoteratosphaeria_perpendicularis","Pyrenochaetopsis_kuksensis",
  "Rhizopogon_ellenae","Rhizopogon_sp","Rhodosporidiobolus_azoricus","Rhodotorula_dairenensis",
  "Rhodotorula_kratochvilovae","Rickenella_indica","Russula_cremicolor","Russula_indocatillus",
  "Russula_pseudotsugarum","Sarcodon_lidongensis","Satchmopsis_metrosideri",
  "Sclerostagonospora_cycadis","Sebacina_sp","Selenophoma_mahoniae",
  "Septoriella_hibernica","Septoriella_oudemansii","Serendipita_sp",
  "Serendipita_vermifera","Serendipitaceae_sp","Simplicillium_lanosoniveum",
  "Solicoccozyma_aeria","Stagonospora_forlicesenensis","Strobilurus_conigenoides",
  "Strobilurus_pachycystidiatus","Subulicystidium_perlongisporum","Suillus_luteus",
  "Symmetrospora_proteacearum","Talaromyces_oumae.annae","Thaxterogaster_argyrionus",
  "Thaxterogaster_rhipiduranus","Thelephora_atra","Thelephora_caryophyllea",
  "Thelephora_sp","Thelephoraceae_sp","Thelonectria_blackeriella","Thyronectria_coryli",
  "Tolypocladium_geodes","Tomentella_lammiensis","Tomentella_sp","Tomentellopsis_echinospora",
  "Tomentellopsis_sp","Torula_lancangjiangensis","Trichocladium_antarcticum","Trichoderma_atlanticum",
  "Trichoderma_austrokoningii","Trichoderma_ochroleucum","Trichoderma_sinense","Tricholoma_badicephalum",
  "Tricholoma_bonii","Tricholoma_olivaceum","Tricholoma_populinum","Tricladium_marylandicum",
  "Tubaria_similis","Umbelopsis_angularis","Unknown","Vacuiphoma_bulgarica","Venturia_albae",
  "Venturia_catenospora","Venturia_mandshurica","Venturia_saliciperda","Verticillium_albo.atrum",
  "Vishniacozyma_psychrotolerans","Vuilleminia_alni","Wilcoxina_rehmii","Xenodidymella_camporesii",
  "Xenopenidiella_nigrescens"
)
df_long$Species <- factor(df_long$Species, levels = guild_levels)



# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Species")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find species that exceed 1% in any treatment combination
species_to_keep <- df_relabund %>%
  group_by(Species) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Species)

df_grouped <- df_relabund %>%
  mutate(Species_grouped = if_else(Species %in% species_to_keep,
                                   Species,
                                   "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Species_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
species_order <- df_grouped %>%
  group_by(Species_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Species_grouped = as.character(Species_grouped)) %>%
  pull(Species_grouped)

# Push Other + Unknown to bottom
species_order <- c(
  setdiff(species_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Species_grouped = factor(Species_grouped, levels = species_order))



final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Species_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_Rhizosphere_Species_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")



## Genus Barplot

df <- read.delim("ITS_Rhizosphere_GenusBarplot_Summed.txt")

group_var <- "Treatment_LandUseHistory" 
facet_row <- "." 
facet_col <- "."     
guild_start_col <- 9    

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Genus",
    values_to = "Count"
  )


df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Genus)
guild_levels <- c(
  "Acephala","Achroiostachys","Agaricus","Alternaria","Amanita","Antrodiella",
  "Apiotrichum","Apseudocercosporella","Ascochyta","Aspergillus","Aureobasidium",
  "Auricularia","Barrenia","Basidioascus","Beauveria","Bergerella","Bovista",
  "Cadophora","Calonarius","Calvatia","Candida","Candolleomyces",
  "Ceratobasidiaceae_gen_Incertae_sedis","Chaetomium","Chaetospermum","Chalara",
  "Chroogomphus","Chrysozyma","Cladophialophora","Cladosporium","Clarireedia",
  "Clitopilus","Clonostachys","Coleophoma","Coniochaeta","Coprinellus","Coprinopsis",
  "Cortinarius","Cosmospora","Cristulariella","Curvularia","Cutaneotrichosporon",
  "Cyphellophora","Cystobasidium","Darksidea","Debaryomyces","Didymosphaeria",
  "Diversispora","Dothiora","Emericellopsis","Entoloma","Entomocorticium","Epicoccum",
  "Erysiphe","Filobasidium","Fonsecaea","Fontanospora","Funneliformis",
  "Furcasterigmium","Fusarium","Gaeumannomyces","Ganoderma","Gliomastix",
  "Gloeopycnis","Gymnopilus","Gyoerffyella","Hanseniaspora","Harzia","Hebeloma",
  "Helotiales_gen_Incertae_sedis","Humicola","Hyaloscyphaceae_gen_Incertae_sedis",
  "Hydnum","Hygrophorus","Hypoderma","Ilyonectria","Inocybe","Juncaceicola","Knufia",
  "Laccaria","Lachnum","Lactarius","Lentithecium","Leohumicola","Leucosporidium",
  "Linnemannia","Lophiostoma","Lycoperdon","Lyophyllum","Macroconia","Mallocybe",
  "Mariannaea","Megalocystidium","Metapochonia","Microbotryum","Microdochium",
  "Minimelanolocus","Morchella","Mortierella","Mrakia","Mucor","Mycena",
  "Mycodidymella","Myrtapenidiella","Naganishia","Neoantrodia","Neoascochyta",
  "Neocucurbitaria","Neodidymelliopsis","Neofabraea","Neophaeomoniella","Nigrospora",
  "Nodulosphaeria","Oculimacula","Odontia","Panaeolina","Papiliotrema",
  "Parabartalinia","Paraboeremia","Paraleptosphaeria","Paraphaeosphaeria",
  "Paraphoma","Penicillium","Pezicula","Peziza","Phenoliferia","Phialocephala",
  "Phlebia","Phlebiopsis","Phlegmacium","Piloderma","Plectosphaerella","Pluteus",
  "Polyozellus","Psathyrella","Pseudeurotium","Pseudofusicoccum","Pseudogymnoascus",
  "Pseudopezicula","Pseudophaeomoniella","Pseudosperma","Pseudoteratosphaeria",
  "Pyrenochaetopsis","Rhizopogon","Rhodosporidiobolus","Rhodotorula","Rickenella",
  "Russula","Sarcodon","Satchmopsis","Sclerostagonospora","Sebacina","Selenophoma",
  "Septoriella","Serendipita","Serendipitaceae_gen_Incertae_sedis","Simplicillium",
  "Solicoccozyma","Stagonospora","Stagonosporopsis","Strobilurus","Subulicystidium",
  "Suillus","Symmetrospora","Talaromyces","Thaxterogaster","Thelephora",
  "Thelephoraceae_gen_Incertae_sedis","Thelonectria","Thyronectria","Tolypocladium",
  "Tomentella","Tomentellopsis","Torula","Trichocladium","Trichoderma","Tricholoma",
  "Tricladium","Tubaria","Umbelopsis","Unknown","Vacuiphoma","Venturia","Verticillium",
  "Vishniacozyma","Vuilleminia","Wilcoxina","Xenodidymella","Xenopenidiella"
)
df_long$Genus <- factor(df_long$Genus, levels = guild_levels)


# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Genus")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find genus that exceed 1% in any treatment combination
genus_to_keep <- df_relabund %>%
  group_by(Genus) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Genus)

df_grouped <- df_relabund %>%
  mutate(Genus_grouped = if_else(Genus %in% genus_to_keep,
                                   Genus,
                                   "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Genus_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
genus_order <- df_grouped %>%
  group_by(Genus_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Genus_grouped = as.character(Genus_grouped)) %>%
  pull(Genus_grouped)

# Push Other + Unknown to bottom
genus_order <- c(
  setdiff(genus_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Genus_grouped = factor(Genus_grouped, levels = genus_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Genus_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_Rhizosphere_Genus_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 13, height = 10, units = "in")


#### Genus/Species EMF Only Barplots ####

library(dplyr)
library(tidyr)
library(ggplot2)


## Prepping table - doing things in excel 

df <- read.delim("ITS_FeatureTable_Rhizosphere_EMFOnly_GenusTable.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body), 
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body))
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_Rhizosphere_Genus_Summed.csv", row.names = TRUE)

## Prepping genus table - doing things in excel 

df <- read.delim("ITS_FeatureTable_Rhizosphere_EMFOnly_SpeciesTable.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),  
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body)) 
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_Rhizosphere_Species_Summed.csv", row.names = TRUE)


## Species Barplot

df <- read.delim("ITS_Rhizosphere_SpeciesBarplot_EMFOnly_Summed.txt")


group_var <- "Treatment_LandUseHistory" 
facet_row <- "." 
facet_col <- "."  
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Species",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Species)
guild_levels <- c(
  "Amanita_pseudopantherina","Calonarius_caesiocinctus","Chroogomphus_pakistanicus",
  "Cortinarius_albidipes","Cortinarius_alboadustus","Cortinarius_brunneoaffinis",
  "Cortinarius_cepistipes","Cortinarius_clarisordidus","Cortinarius_collangustus",
  "Cortinarius_deceptivissimus","Cortinarius_decipiens","Cortinarius_kauffmanianus",
  "Cortinarius_majorinus","Cortinarius_miwok","Cortinarius_mucosus",
  "Cortinarius_parvannulatus","Cortinarius_paululus","Cortinarius_subcarneinatus",
  "Cortinarius_turgidipes","Hebeloma_erebium","Hebeloma_pungens",
  "Hydnum_pallidomarginatum","Hydnum_subolympicum","Hygrophorus_fuscoalboides",
  "Inocybe_purpureobadia","Inocybe_silvae.herbaceae","Laccaria_fagacicola",
  "Lactarius_aurantiacopallens","Lactarius_miniatosporus","Mallocybe_sp",
  "Phlegmacium_subfoetens","Piloderma_bicolor","Piloderma_lanatum","Piloderma_sp",
  "Polyozellus_mariae","Pseudosperma_emberizanum","Rhizopogon_ellenae",
  "Rhizopogon_sp","Russula_cremicolor","Russula_indocatillus","Russula_pseudotsugarum",
  "Sarcodon_lidongensis","Suillus_luteus","Thaxterogaster_argyrionus",
  "Thaxterogaster_rhipiduranus","Thelephora_atra","Thelephora_caryophyllea",
  "Thelephora_sp","Tomentella_lammiensis","Tomentella_sp","Tomentellopsis_echinospora",
  "Tomentellopsis_sp","Tricholoma_badicephalum","Tricholoma_bonii",
  "Tricholoma_olivaceum","Tricholoma_populinum","Unknown","Wilcoxina_rehmii"
)
df_long$Species <- factor(df_long$Species, levels = guild_levels)


# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Species")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find species that exceed 1% in any treatment combination
species_to_keep <- df_relabund %>%
  group_by(Species) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Species)

df_grouped <- df_relabund %>%
  mutate(Species_grouped = if_else(Species %in% species_to_keep,
                                   Species,
                                   "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Species_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
species_order <- df_grouped %>%
  group_by(Species_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Species_grouped = as.character(Species_grouped)) %>%
  pull(Species_grouped)

# Push Other + Unknown to bottom
species_order <- c(
  setdiff(species_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Species_grouped = factor(Species_grouped, levels = species_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Species_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("Barplot_Rhizosphere_EMFOnly_Species_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")



## Genus Barplot

df <- read.delim("ITS_Rhizosphere_GenusBarplot_EMFOnly_Summed.txt")


group_var <- "Treatment_LandUseHistory"  
facet_row <- "."  
facet_col <- "."   
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Genus",
    values_to = "Count"
  )


df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Genus)
guild_levels <- c(
  "Amanita","Calonarius","Chroogomphus","Cortinarius","Hebeloma","Hydnum",
  "Hygrophorus","Inocybe","Laccaria","Lactarius","Mallocybe","Phlegmacium",
  "Piloderma","Polyozellus","Pseudosperma","Rhizopogon","Russula","Sarcodon",
  "Suillus","Thaxterogaster","Thelephora","Tomentella","Tomentellopsis",
  "Tricholoma","Wilcoxina"
)
df_long$Genus <- factor(df_long$Genus, levels = guild_levels)


# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Genus")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find genus that exceed 1% in any treatment combination
genus_to_keep <- df_relabund %>%
  group_by(Genus) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Genus)

df_grouped <- df_relabund %>%
  mutate(Genus_grouped = if_else(Genus %in% genus_to_keep,
                                 Genus,
                                 "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Genus_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
genus_order <- df_grouped %>%
  group_by(Genus_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Genus_grouped = as.character(Genus_grouped)) %>%
  pull(Genus_grouped)

# Push Other + Unknown to bottom
genus_order <- c(
  setdiff(genus_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Genus_grouped = factor(Genus_grouped, levels = genus_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Genus_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("Barplot_Rhizosphere_EMFOnly_Genus_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")








#### Root Tips - EMF Community #### 

#### NMDS Prep ####

## Using relative abundances for input into NMDS owing to lack of rarefying with this dataset 

setwd()

otus <- read.delim("ITS_FeatureTable_RootTipsOnly_EMFOnly.txt",header=T,row.names=1, check.names=FALSE)
otus_t <- t(otus)
otus_rel_abund_t <- decostand(otus_t, method = "total")
otus_rel_abund <- t(otus_rel_abund_t)
otus <- otus_rel_abund

otus <- as.data.frame(otus)
sample_names <- colnames(otus)

map_file <- read.delim("ITS_metadata_RootTips.txt",header = T,row.names=1,check.names=FALSE)
map_file <- map_file[rownames(map_file) %in% sample_names, ]
map_file <- as.data.frame(map_file)


all.equal(names(otus),row.names(map_file))


map_file$Treatment_LandUseHistory <- factor(map_file$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
map_file$Salvage_Harvest_Status <- factor(map_file$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 



otumat<-as.matrix(otus)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
head(OTU)
taxa<-read.delim("ITS_taxonomy_RootTips.txt",header=T,row.names=1) 
taxmat<-as.matrix(taxa)
taxmat <- taxmat[rownames(taxmat) %in% rownames(otumat), ]
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
  stat_ellipse(aes(x=NMDS1, y=NMDS2, group = Salvage_Harvest_Status, color=Salvage_Harvest_Status),linetype = 2) 

NMDS

ggsave("NMDS_RootTips_EMFOnly.pdf", NMDS, width = 14, height = 10, units = "in")


#### PERMANOVA ####

adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Salvage_Harvest_Status, data = sample_tab, permutations = 999, method = "bray")
adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory*Salvage_Harvest_Status, data = sample_tab, permutations = 999, method = "bray")


## P-value adjustments:
library(stats)
p = c(0.001, 0.067) ## Change with above findings 
p.adjusted = p.adjust(p, method = "BH")
p.adjusted

# Adjusted: 0.001666667 0.001666667 0.002500000 0.001666667 0.003000000

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)
pairwise.adonis2(mgd_ge5K_relabund.bray ~ Treatment_LandUseHistory, data = sample_tab, sim.method = "bray", p.adjust.m = "BH", perm = 999)




#### Alpha Diversity ####

library(vegan)    
library(ggplot2) 
library(dplyr)  
library(ggpubr) 

asv_table <- t(otus)   # samples x ASVs 
metadata <- map_file    # samples x variables

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
alpha_div$Salvage_Harvest_Status <- factor(alpha_div$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

alpha_div$Treatment_All <- alpha_div$Treatment_LandUseHistory
unique(alpha_div$Treatment_All)


treatment_levels <- c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")
comparisons <- combn(treatment_levels, 2, simplify = FALSE)



alpha_div <- alpha_div %>%
  relocate(Observed, Shannon, .after = last_col())

write.csv(alpha_div, "ITS_alpha_div_roottips_EMFOnly_ObservedShannons.csv", row.names = FALSE)




## Figures 

## Observed ASVs

AD_RootTips_ObservedSpecies_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Observed)) +
  geom_boxplot(aes(fill = Treatment_All), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_All), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Root Tips Observed Species by Treatment",
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

AD_RootTips_ObservedSpecies_Rarefied

ggsave("AD_RootTips_EMFOnly_ObservedSpecies_6_tests.pdf", AD_RootTips_ObservedSpecies_Rarefied, width = 3, height = 5, units = "in")


## Shannons 

AD_RootTips_Shannons_Rarefied <- ggplot(alpha_div, aes(x = Treatment_LandUseHistory, y = Shannon)) +
  geom_boxplot(aes(fill = Treatment_All), outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(fill = Treatment_All), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.8, shape = 21, color = "black") +
  labs(title = "Root Tips Shannon Diversity by Treatment",
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

AD_RootTips_Shannons_Rarefied

ggsave("AD_RootTips_EMFOnly_Shannons_6_tests.pdf", AD_RootTips_Shannons_Rarefied, width = 3, height = 5, units = "in")


## Final plots together 

library(ggpubr)

final_plot <- ggarrange(
  AD_RootTips_ObservedSpecies_Rarefied, AD_RootTips_Shannons_Rarefied,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)

print(final_plot)

ggsave("alphadiversity_boxplots_RootTips_EMFOnly.pdf", final_plot, width = 14, height = 5, units = "in")







#### Genus/Species Barplots - All Taxa ####

library(dplyr)
library(tidyr)
library(ggplot2)


## Prepping table - doing things in excel 

df <- read.delim("ITS_GenusTable_RootTipsOnly.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),       
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body))  
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_RootTipsOnly_Genus_Summed.csv", row.names = TRUE)

## Prepping genus table - doing things in excel 

df <- read.delim("ITS_SpeciesTable_RootTipsOnly.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),           
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body))  
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_RootTipsOnly_Species_Summed.csv", row.names = TRUE)

## Species Barplot

df <- read.delim("ITS_RootTipsOnly_SpeciesBarplot_Summed.txt")


group_var <- "Treatment_LandUseHistory"   
facet_row <- "."      
facet_col <- "."        
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Species",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Species)
guild_levels <- c(
  "Acephala_sp",
  "Agaricus_agrinferus",
  "Agaricus_arvensis",
  "Agaricus_zhangyensis",
  "Alternaria_ellipsoidea",
  "Alternaria_multirostrata",
  "Alternaria_nepalensis",
  "Amanita_pseudopantherina",
  "Antrodiella_faginea",
  "Apiotrichum_sporotrichoides",
  "Apseudocercosporella_trigonotidis",
  "Ascochyta_medicaginicola_var._macrospora",
  "Aspergillus_sigurros",
  "Aureobasidium_castaneae",
  "Auricularia_subglabra",
  "Auricularia_submesenterica",
  "Barrenia_panici",
  "Basidioascus_undulatus",
  "Beauveria_varroae",
  "Bergerella_atrofusca",
  "Bovista_cretacea",
  "Cadophora_domestica",
  "Calonarius_caesiocinctus",
  "Calvatia_fragilis",
  "Candida_orthopsilosis",
  "Candida_pseudoglaebosa",
  "Candolleomyces_cacao",
  "Ceratobasidiaceae_sp",
  "Chaetomium_cervicicola",
  "Chaetospermum_sp",
  "Chalara_recta",
  "Chroogomphus_pakistanicus",
  "Chrysozyma_pseudogriseoflava",
  "Cladophialophora_humicola",
  "Cladosporium_aphidis",
  "Clarireedia_paspali",
  "Clitopilus_brunneiceps",
  "Clonostachys_divergens",
  "Coleophoma_camelliae",
  "Coniochaeta_vineae",
  "Coprinellus_sclerocystidiosus",
  "Coprinopsis_marcescibilis",
  "Cortinarius_albidipes",
  "Cortinarius_alboadustus",
  "Cortinarius_brunneoaffinis",
  "Cortinarius_cepistipes",
  "Cortinarius_clarisordidus",
  "Cortinarius_collangustus",
  "Cortinarius_deceptivissimus",
  "Cortinarius_decipiens",
  "Cortinarius_kauffmanianus",
  "Cortinarius_majorinus",
  "Cortinarius_miwok",
  "Cortinarius_mucosus",
  "Cortinarius_parvannulatus",
  "Cortinarius_paululus",
  "Cortinarius_subcarneinatus",
  "Cortinarius_turgidipes",
  "Cristulariella_depraedans",
  "Curvularia_siddiquii",
  "Cutaneotrichosporon_dermatis",
  "Cyphellophora_eucalypti",
  "Cyphellophora_vermispora",
  "Cystobasidium_raffinophilum",
  "Cystobasidium_terricola",
  "Darksidea_zeta",
  "Debaryomyces_robertsiae",
  "Didymosphaeria_variabile",
  "Diversispora_spurca",
  "Dothiora_europaea",
  "Emericellopsis_glabra",
  "Entoloma_aurorae.borealis",
  "Entoloma_fuscohebes",
  "Entoloma_tibiicystidiatum",
  "Entomocorticium_portiae",
  "Epicoccum_thailandicum",
  "Erysiphe_polygoni",
  "Filobasidium_magnum",
  "Filobasidium_oeirense",
  "Fonsecaea_nubica",
  "Fontanospora_fusiramosa",
  "Funneliformis_caledonium",
  "Furcasterigmium_furcatum",
  "Fusarium_algeriense",
  "Fusarium_brevicatenulatum",
  "Fusarium_croci",
  "Fusarium_nelsonii",
  "Gaeumannomyces_hyphopodioides",
  "Ganoderma_carocalcareum",
  "Gliomastix_roseogrisea",
  "Gloeopycnis_protuberans",
  "Gymnopilus_aurantiophyllus",
  "Gyoerffyella_sp",
  "Hanseniaspora_nectarophila",
  "Harzia_acremonioides",
  "Hebeloma_erebium",
  "Hebeloma_pungens",
  "Helotiales_sp",
  "Humicola_zollerniae",
  "Hyaloscyphaceae_sp",
  "Hydnum_pallidomarginatum",
  "Hydnum_subolympicum",
  "Hygrophorus_fuscoalboides",
  "Hypoderma_siculum",
  "Ilyonectria_mors.panacis",
  "Inocybe_purpureobadia",
  "Inocybe_silvae.herbaceae",
  "Juncaceicola_padellana",
  "Knufia_tsunedae",
  "Laccaria_fagacicola",
  "Lachnum_impudicum",
  "Lactarius_aurantiacopallens",
  "Lactarius_miniatosporus",
  "Lentithecium_pseudoclioninum",
  "Leohumicola_atra",
  "Leucosporidium_escuderoi",
  "Linnemannia_amoeboidea",
  "Linnemannia_camargensis",
  "Lophiostoma_rosae",
  "Lycoperdon_echinatum",
  "Lycoperdon_ericaeum",
  "Lycoperdon_nigrescens",
  "Lyophyllum_shimeji",
  "Macroconia_bulbipes",
  "Mallocybe_sp",
  "Mariannaea_camptospora",
  "Megalocystidium_perticatum",
  "Metapochonia_hahajimaensis",
  "Microbotryum_violaceum",
  "Microdochium_lycopodinum",
  "Minimelanolocus_clavatus",
  "Morchella_tasmanica",
  "Mortierella_alliacea",
  "Mortierella_fluviae",
  "Mortierella_wuyishanensis",
  "Mrakia_niccombsii",
  "Mycena_megaspora",
  "Mycodidymella_aesculi",
  "Naganishia_bhutanensis",
  "Naganishia_diffluens",
  "Naganishia_vishniacii",
  "Neoantrodia_infirma",
  "Neoantrodia_serialis",
  "Neoascochyta_tardicrescens",
  "Neocucurbitaria_aetnensis",
  "Neocucurbitaria_rhamnioides",
  "Neodidymelliopsis_negundinis",
  "Neofabraea_salicina",
  "Neophaeomoniella_constricta",
  "Nodulosphaeria_thalictri",
  "Oculimacula_yallundae",
  "Odontia_sp",
  "Panaeolina_foenisecii",
  "Papiliotrema_nemorosa",
  "Parabartalinia_lateralis",
  "Paraboeremia_selaginellae",
  "Paraleptosphaeria_nitschkei",
  "Paraphaeosphaeria_neglecta",
  "Paraphoma_ledniceana",
  "Paraphoma_pye",
  "Penicillium_alutaceum",
  "Penicillium_coprobium",
  "Penicillium_janczewskii",
  "Penicillium_madriti",
  "Penicillium_oregonense",
  "Penicillium_vasconiae",
  "Peziza_montirivicola",
  "Phenoliferia_psychrophenolica",
  "Phialocephala_amethystea",
  "Phialocephala_helenae",
  "Phlebia_rufa",
  "Phlebiopsis_sp",
  "Phlegmacium_subfoetens",
  "Piloderma_bicolor",
  "Piloderma_lanatum",
  "Piloderma_sp",
  "Plectosphaerella_melonis",
  "Pluteus_kovalenkoi",
  "Polyozellus_mariae",
  "Psathyrella_carinthiaca",
  "Pseudeurotium_desertorum",
  "Pseudofusicoccum_olivaceum",
  "Pseudogymnoascus_roseus",
  "Pseudopezicula_betulae",
  "Pseudosperma_emberizanum",
  "Pseudoteratosphaeria_perpendicularis",
  "Pyrenochaetopsis_kuksensis",
  "Rhizopogon_ellenae",
  "Rhizopogon_sp",
  "Rhodosporidiobolus_azoricus",
  "Rhodotorula_dairenensis",
  "Rhodotorula_kratochvilovae",
  "Rickenella_indica",
  "Russula_cremicolor",
  "Russula_indocatillus",
  "Russula_pseudotsugarum",
  "Sarcodon_lidongensis",
  "Satchmopsis_metrosideri",
  "Sclerostagonospora_cycadis",
  "Sebacina_sp",
  "Selenophoma_mahoniae",
  "Septoriella_hibernica",
  "Septoriella_oudemansii",
  "Serendipita_sp",
  "Serendipita_vermifera",
  "Serendipitaceae_sp",
  "Simplicillium_lanosoniveum",
  "Solicoccozyma_aeria",
  "Stagonospora_forlicesenensis",
  "Strobilurus_conigenoides",
  "Strobilurus_pachycystidiatus",
  "Subulicystidium_perlongisporum",
  "Suillus_luteus",
  "Symmetrospora_proteacearum",
  "Talaromyces_oumae.annae",
  "Thaxterogaster_argyrionus",
  "Thaxterogaster_rhipiduranus",
  "Thelephora_atra",
  "Thelephora_caryophyllea",
  "Thelephora_sp",
  "Thelephoraceae_sp",
  "Thelonectria_blackeriella",
  "Thyronectria_coryli",
  "Tolypocladium_geodes",
  "Tomentella_lammiensis",
  "Tomentella_sp",
  "Tomentellopsis_echinospora",
  "Tomentellopsis_sp",
  "Torula_lancangjiangensis",
  "Trichocladium_antarcticum",
  "Trichoderma_atlanticum",
  "Trichoderma_austrokoningii",
  "Trichoderma_ochroleucum",
  "Trichoderma_sinense",
  "Tricholoma_badicephalum",
  "Tricholoma_bonii",
  "Tricholoma_olivaceum",
  "Tricholoma_populinum",
  "Tricladium_marylandicum",
  "Tubaria_similis",
  "Umbelopsis_angularis",
  "Unknown",
  "Vacuiphoma_bulgarica",
  "Venturia_albae",
  "Venturia_catenospora",
  "Venturia_mandshurica",
  "Venturia_saliciperda",
  "Verticillium_albo.atrum",
  "Vishniacozyma_psychrotolerans",
  "Vuilleminia_alni",
  "Wilcoxina_rehmii",
  "Xenodidymella_camporesii",
  "Xenopenidiella_nigrescens"
)
df_long$Species <- factor(df_long$Species, levels = guild_levels)

# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Species")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find species that exceed 1% in any treatment combination
species_to_keep <- df_relabund %>%
  group_by(Species) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Species)

df_grouped <- df_relabund %>%
  mutate(Species_grouped = if_else(Species %in% species_to_keep,
                                   Species,
                                   "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Species_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
species_order <- df_grouped %>%
  group_by(Species_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Species_grouped = as.character(Species_grouped)) %>%
  pull(Species_grouped)

# Push Other + Unknown to bottom
species_order <- c(
  setdiff(species_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Species_grouped = factor(Species_grouped, levels = species_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Species_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_RootTipsOnly_Species_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")



## Genus Barplot

df <- read.delim("ITS_RootTipsOnly_GenusBarplot_Summed.txt")


group_var <- "Treatment_LandUseHistory"   
facet_row <- "."   
facet_col <- "."     
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Genus",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Genus)
guild_levels <- c(
  "Acephala","Achroiostachys","Agaricus","Alternaria","Amanita","Antrodiella",
  "Apiotrichum","Apseudocercosporella","Ascochyta","Aspergillus","Aureobasidium",
  "Auricularia","Barrenia","Basidioascus","Beauveria","Bergerella","Bovista",
  "Cadophora","Calonarius","Calvatia","Candida","Candolleomyces",
  "Ceratobasidiaceae_gen_Incertae_sedis","Chaetomium","Chaetospermum",
  "Chalara","Chroogomphus","Chrysozyma","Cladophialophora","Cladosporium",
  "Clarireedia","Clitopilus","Clonostachys","Coleophoma","Coniochaeta",
  "Coprinellus","Coprinopsis","Cortinarius","Cosmospora","Cristulariella",
  "Curvularia","Cutaneotrichosporon","Cyphellophora","Cystobasidium","Darksidea",
  "Debaryomyces","Didymosphaeria","Diversispora","Dothiora","Emericellopsis",
  "Entoloma","Entomocorticium","Epicoccum","Erysiphe","Filobasidium",
  "Fonsecaea","Fontanospora","Funneliformis","Furcasterigmium","Fusarium",
  "Gaeumannomyces","Ganoderma","Gliomastix","Gloeopycnis","Gymnopilus",
  "Gyoerffyella","Hanseniaspora","Harzia","Hebeloma","Helotiales_gen_Incertae_sedis",
  "Humicola","Hyaloscyphaceae_gen_Incertae_sedis","Hydnum","Hygrophorus",
  "Hypoderma","Ilyonectria","Inocybe","Juncaceicola","Knufia","Laccaria",
  "Lachnum","Lactarius","Lentithecium","Leohumicola","Leucosporidium",
  "Linnemannia","Lophiostoma","Lycoperdon","Lyophyllum","Macroconia","Mallocybe",
  "Mariannaea","Megalocystidium","Metapochonia","Microbotryum","Microdochium",
  "Minimelanolocus","Morchella","Mortierella","Mrakia","Mucor","Mycena",
  "Mycodidymella","Myrtapenidiella","Naganishia","Neoantrodia","Neoascochyta",
  "Neocucurbitaria","Neodidymelliopsis","Neofabraea","Neophaeomoniella",
  "Nigrospora","Nodulosphaeria","Oculimacula","Odontia","Panaeolina",
  "Papiliotrema","Parabartalinia","Paraboeremia","Paraleptosphaeria",
  "Paraphaeosphaeria","Paraphoma","Penicillium","Pezicula","Peziza",
  "Phenoliferia","Phialocephala","Phlebia","Phlebiopsis","Phlegmacium",
  "Piloderma","Plectosphaerella","Pluteus","Polyozellus","Psathyrella",
  "Pseudeurotium","Pseudofusicoccum","Pseudogymnoascus","Pseudopezicula",
  "Pseudophaeomoniella","Pseudosperma","Pseudoteratosphaeria","Pyrenochaetopsis",
  "Rhizopogon","Rhodosporidiobolus","Rhodotorula","Rickenella","Russula",
  "Sarcodon","Satchmopsis","Sclerostagonospora","Sebacina","Selenophoma",
  "Septoriella","Serendipita","Serendipitaceae_gen_Incertae_sedis",
  "Simplicillium","Solicoccozyma","Stagonospora","Stagonosporopsis",
  "Strobilurus","Subulicystidium","Suillus","Symmetrospora","Talaromyces",
  "Thaxterogaster","Thelephora","Thelephoraceae_gen_Incertae_sedis","Thelonectria",
  "Thyronectria","Tolypocladium","Tomentella","Tomentellopsis","Torula",
  "Trichocladium","Trichoderma","Tricholoma","Tricladium","Tubaria","Umbelopsis",
  "Unknown","Vacuiphoma","Venturia","Verticillium","Vishniacozyma","Vuilleminia",
  "Wilcoxina","Xenodidymella","Xenopenidiella"
)
df_long$Genus <- factor(df_long$Genus, levels = guild_levels)

# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

# Relative abundance calculation
df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Genus")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find genus that exceed 1% in any treatment combination
genus_to_keep <- df_relabund %>%
  group_by(Genus) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Genus)

df_grouped <- df_relabund %>%
  mutate(Genus_grouped = if_else(Genus %in% genus_to_keep,
                                 Genus,
                                 "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Genus_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
genus_order <- df_grouped %>%
  group_by(Genus_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Genus_grouped = as.character(Genus_grouped)) %>%
  pull(Genus_grouped)

# Push Other + Unknown to bottom
genus_order <- c(
  setdiff(genus_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Genus_grouped = factor(Genus_grouped, levels = genus_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Genus_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("FUNGuild_Barplot_RootTipsOnly_Genus_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")



#### Genus/Species EMF Only Barplots ####

library(dplyr)
library(tidyr)
library(ggplot2)

## Prepping species table - doing things in excel 

df <- read.delim("ITS_FeatureTable_RootTipsOnly_EMFOnly_GenusTable.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),         
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body)) 
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_RootTipsOnly_EMFOnly_Genus_Summed.csv", row.names = TRUE)

## Prepping genus table - doing things in excel 

df <- read.delim("ITS_FeatureTable_RootTipsOnly_EMFOnly_SpeciesTable.txt")
df <- t(df)

hdr <- trimws(df[1, ])

df_body <- df[-1, , drop = FALSE]

colnames(df_body) <- hdr

df_num <- matrix(
  as.numeric(df_body),       
  nrow = nrow(df_body),
  dimnames = list(rownames(df_body), colnames(df_body))  
)

col_groups <- split(seq_len(ncol(df_num)), colnames(df_num))

df_summed <- sapply(col_groups, function(idx) {
  rowSums(df_num[, idx, drop = FALSE], na.rm = TRUE)
})

write.csv(df_summed, "ITS_FeatureTable_RootTipsOnly_EMFOnly_Species_Summed.csv", row.names = TRUE)

## Species Barplot

df <- read.delim("ITS_RootTipsOnly_SpeciesBarplot_EMFOnly_Summed.txt")


group_var <- "Treatment_LandUseHistory"  
facet_row <- "."     
facet_col <- "."        
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Species",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Species)
guild_levels <- c(
  "Amanita_pseudopantherina","Calonarius_caesiocinctus","Chroogomphus_pakistanicus",
  "Cortinarius_albidipes","Cortinarius_alboadustus","Cortinarius_brunneoaffinis",
  "Cortinarius_cepistipes","Cortinarius_clarisordidus","Cortinarius_collangustus",
  "Cortinarius_deceptivissimus","Cortinarius_decipiens","Cortinarius_kauffmanianus",
  "Cortinarius_majorinus","Cortinarius_miwok","Cortinarius_mucosus",
  "Cortinarius_parvannulatus","Cortinarius_paululus","Cortinarius_subcarneinatus",
  "Cortinarius_turgidipes","Hebeloma_erebium","Hebeloma_pungens",
  "Hydnum_pallidomarginatum","Hydnum_subolympicum","Hygrophorus_fuscoalboides",
  "Inocybe_purpureobadia","Inocybe_silvae.herbaceae","Laccaria_fagacicola",
  "Lactarius_aurantiacopallens","Lactarius_miniatosporus","Mallocybe_sp",
  "Phlegmacium_subfoetens","Piloderma_bicolor","Piloderma_lanatum","Piloderma_sp",
  "Polyozellus_mariae","Pseudosperma_emberizanum","Rhizopogon_ellenae",
  "Rhizopogon_sp","Russula_cremicolor","Russula_indocatillus",
  "Russula_pseudotsugarum","Sarcodon_lidongensis","Suillus_luteus",
  "Thaxterogaster_argyrionus","Thaxterogaster_rhipiduranus","Thelephora_atra",
  "Thelephora_caryophyllea","Thelephora_sp","Tomentella_lammiensis","Tomentella_sp",
  "Tomentellopsis_echinospora","Tomentellopsis_sp","Tricholoma_badicephalum",
  "Tricholoma_bonii","Tricholoma_olivaceum","Tricholoma_populinum","Unknown",
  "Wilcoxina_rehmii"
)
df_long$Species <- factor(df_long$Species, levels = guild_levels)

# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Species")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find species that exceed 1% in any treatment combination
species_to_keep <- df_relabund %>%
  group_by(Species) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Species)

df_grouped <- df_relabund %>%
  mutate(Species_grouped = if_else(Species %in% species_to_keep,
                                   Species,
                                   "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Species_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
species_order <- df_grouped %>%
  group_by(Species_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Species_grouped = as.character(Species_grouped)) %>%
  pull(Species_grouped)

# Push Other + Unknown to bottom
species_order <- c(
  setdiff(species_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Species_grouped = factor(Species_grouped, levels = species_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Species_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("Barplot_RootTipsOnly_EMFOnly_Species_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")



## Genus Barplot

df <- read.delim("ITS_RootTipsOnly_GenusBarplot_EMFOnly_Summed.txt")

group_var <- "Treatment_LandUseHistory" 
facet_row <- "."  
facet_col <- "."    
guild_start_col <- 9        # index where relevant columns start

df_long <- df %>%
  pivot_longer(
    cols = guild_start_col:ncol(df),
    names_to = "Genus",
    values_to = "Count"
  )

df_long$Treatment_LandUseHistory <- factor(df_long$Treatment_LandUseHistory, levels = c("Old_Forest", "1stPreMPB", "2ndPreMPB", "1stPostMPB", "2ndPostMPB", "Recent_Cut")) 
df_long$Salvage_Harvest_Status <- factor(df_long$Salvage_Harvest_Status, levels = c("Non_Salvage_Harvested", "Salvage_Harvested")) 

unique(df_long$Genus)
guild_levels <- c(
  "Amanita","Calonarius","Chroogomphus","Cortinarius","Hebeloma","Hydnum",
  "Hygrophorus","Inocybe","Laccaria","Lactarius","Mallocybe","Phlegmacium",
  "Piloderma","Polyozellus","Pseudosperma","Rhizopogon","Russula","Sarcodon",
  "Suillus","Thaxterogaster","Thelephora","Tomentella","Tomentellopsis",
  "Tricholoma","Wilcoxina"
)
df_long$Genus <- factor(df_long$Genus, levels = guild_levels)

# Build grouping variable list and remove "." placeholders
grouping_vars <- c(group_var, facet_row, facet_col)
grouping_vars <- grouping_vars[grouping_vars != "."]

df_relabund <- df_long %>%
  group_by(across(all_of(c(grouping_vars, "Genus")))) %>%
  summarise(SumCount = sum(Count), .groups = "drop") %>%
  group_by(across(all_of(grouping_vars))) %>%
  mutate(RelAbund = SumCount / sum(SumCount)) %>%
  ungroup()


library(dplyr)

# Find genus that exceed 1% in any treatment combination
genus_to_keep <- df_relabund %>%
  group_by(Genus) %>%
  summarize(max_rel = max(RelAbund, na.rm = TRUE)) %>%
  filter(max_rel >= 0.01) %>%
  pull(Genus)

df_grouped <- df_relabund %>%
  mutate(Genus_grouped = if_else(Genus %in% genus_to_keep,
                                 Genus,
                                 "Other"))

df_grouped <- df_grouped %>%
  group_by(Treatment_LandUseHistory, Genus_grouped) %>%
  summarize(RelAbund = sum(RelAbund), .groups = "drop")


## Reorder to be ranked by abundance & with other/unknown at the end
genus_order <- df_grouped %>%
  group_by(Genus_grouped) %>%
  summarize(total = sum(RelAbund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(Genus_grouped = as.character(Genus_grouped)) %>%
  pull(Genus_grouped)

# Push Other + Unknown to bottom
genus_order <- c(
  setdiff(genus_order, c("Other", "Unknown")),
  "Unknown",
  "Other"
)

df_grouped <- df_grouped %>%
  mutate(Genus_grouped = factor(Genus_grouped, levels = genus_order))



# Plot using facet_grid — "." works here for no facetting
final_plot <- ggplot(df_grouped, aes_string(x = group_var, y = "RelAbund", fill = "Genus_grouped")) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(as.formula(paste(facet_row, "~", facet_col))) +
  ylab("Relative Abundance") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot

ggsave("Barplot_RootTipsOnly_EMFOnly_Genus_1percentinatleast1treatment_rearranged.pdf", final_plot, width = 10, height = 10, units = "in")











