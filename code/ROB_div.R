## ---- 1. IMPORT PACKAGE ----

#{{multivariate analysis}}
library(FactoMineR) 
library(ade4)
library(factoextra) #extract information from pca and relatives
####library(mixOmics) #discriminant analysis PLS-DA and sPLS-DA
####library(uwot) #Uniform Manifold Approximation and Projection (UMAP)
####library(phateR) #Markov Affinity-Based Graph Imputation
library(GUniFrac)

#{{hierarchical clustering with bootstrap values}}
####library(pvclust)

#{{multivariate compositional analysis}}
####library (ALDEx2)
#{{multivariate Count-Based Differential Abundance Analysis}}
library(edgeR)

#{{correlation figure}}
library(corrplot)
library(ggcorrplot)

#{{dendogram graph}}
library(ggdendro)
library(dendextend)

#{{Circular and network graph}}
####library(ggraph)
####library(igraph)

##{{Heatmap graph}}
####library(ComplexHeatmap)

#{{vendiagram analysis}}
library(ggVennDiagram)
library(UpSetR)

#{{ecological analysis; diversity indices, permanova on distance, multivariate analysis}}
library(vegan)
library(pairwiseAdonis)
library(hilldiv)

#{{normality test}}
library(nortest)

#{{data distribution}}
####library(performance)

#{{Indicator species analysis}}
library(indicspecies)

#{{Import functional database}}
library(metagMisc) #use to import FUNGuild, but many other tips as
#pairwise dissimilarity boxplots, prevalence plots, diversity profiles based on Hill numbers
library(fungaltraits) #use to import FunFun database

#{{bi-network analysis}}
library(bipartite) 

#{{all inclusive packages}}
####library(phyloseq)
####library(microeco)
####library(microbiomeExplorer)

#{{color palette}}
library(RColorBrewer)
library(viridis)
library(ggsci)

#{{figure and text manipulation extra}}
library(cowplot) 
library(ggrepel)
library(scales)
library(ggfortify)
library(ggforce)
library(ggpubr)
library(patchwork)
library(ggtext) #allow text modification, for example name in italic with *name* 
#{{Marginal distribution}}
library(ggside)
#{{convert plot in ggplot object}}
library(ggplotify)

#{{table manipulation and graph, including ggplot2}}
library(tidyverse)

## ---- 2. VARIABLE SELECTION ----
#{{.table2 files correspond to several projects, need to select ROBUST samples}}

data_info <- "metadata.txt"
data_16S <- ""
data_ITS <- "Robust_MiSeq_ITS2_ITS86F_ITS4r_20241129_126_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
data_18S <- "Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"

#{{set up parameter}}
silva_levels <- "taxa2ranks_ssu_138.1.table"

bact_taxonomic_levels <- c("domain", "phylum", "class", "order", "family", 
                           "genus", "species")

its_taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family",
                          "genus", "species")

amf_taxonomic_levels <- c("kingdom", "domain", "phylum", "class", "order", "family",
                          "genus", "species")

seed <- 1 #set up the seed for statistics


## ---- 3. IMPORT AND FORMAT METADATA ----
metadata <- data_info %>%
  read_tsv()

## ---- 4. IMPORT AND FORMAT OTU ---- 
## -------- 4.1. bacteria --------

#{{list of function to modify taxonomic ranks}}
. %>%
  read_tsv(na = c("", "NA", "*", "No_hit"),
           show_col_types = FALSE) %>%
  mutate(taxonomy_consensus = taxonomy) %>%
  separate_rows(taxonomy, sep = "[|]") -> fold_table

. %>%
  filter(taxonomy != "*" &
           ! is.na(taxonomy) &
           ! str_starts(taxonomy, "uncultured_")) ->
  filter_redundant_taxa_names

. %>%
  inner_join(silva_levels %>%
               read_tsv(show_col_types = FALSE),
             by = c("taxonomy" = "taxonomic_path")) ->
  merge_with_taxo_levels

. %>%
  distinct() -> remove_duplicated_levels

. %>%
  pivot_wider(names_from = taxonomic_rank,
              values_from = taxonomy) -> unfold_table

. %>%
  mutate(species = case_when(
    str_ends(taxonomy_consensus, "[*]") ~ NA_character_,
    is.na(genus) ~ NA_character_,
    TRUE ~ str_remove(taxonomy_consensus, str_c(".*", genus, "[|]"))
  )) -> add_species_column

bact_OTUs <- data_16S %>%
  #{{format taxonomic ranks}}
  fold_table %>%
  filter_redundant_taxa_names %>%
  merge_with_taxo_levels %>%
  remove_duplicated_levels %>%
  unfold_table %>%
  add_species_column %>%
  #{{format column names}}
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~sub("_16s", "", .)) %>%
  rename_with(~sub("_1_", "_01_", .)) %>%
  rename_with(~sub("_2_", "_02_", .)) %>%
  rename_with(~sub("_3_", "_03_", .)) %>%
  rename_with(~sub("_4_", "_04_", .)) %>%
  rename_with(~sub("_5_", "_05_", .)) %>%
  rename_with(~sub("_6_", "_06_", .)) %>%
  rename_with(~sub("_7_", "_07_", .)) %>%
  rename_with(~sub("_8_", "_08_", .)) %>%
  rename_with(~sub("_9_", "_09_", .)) %>%
  rename_with(~sub("_s", "_sk", .)) %>%
  rename_with(~sub("SU_USFLO_T_PCR2", "CONT_SU_PCR2b", .)) %>%
  rename_with(~sub("SU_USFLO_T_PCR1", "CONT_SU_PCR1b", .)) %>%
  rename_with(~sub("SU_USFLO_T_ex", "CONT_SU_T_exb", .)) %>%
  #{{select samples from SuMi project}}
  select(OTU, bact_taxonomic_levels, starts_with("CONT", ignore.case = FALSE),
         starts_with("SU", ignore.case = FALSE)) %>%
  #{{sequences with "No_hit" has to be deleted}}
  filter(grepl("^Bacteria|^Archaea",domain)) %>%
  #{{sequences affiliated to chloroplast and mitonchondria are excluded}}
  filter(order != "Chloroplast", family != "Mitochondria") %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE),
         starts_with("SU", ignore.case = FALSE)) %>%
  droplevels()

## -------- 4.2. fungi --------

ITS2_OTUs <- data_ITS %>%
  #{{import data}}
  read_tsv() %>%
  #{{format column names}}
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~sub("ROB_T_", "CONT_ROB_", .)) %>%
  rename_with(~sub("ROB_Tex", "CONT_ROB_T_ex", .)) %>%
  rename_with(~sub("_ITS", "", .)) %>%
  #{{select samples from ROBUST project}}
  select(OTU, taxonomy, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB", ignore.case = FALSE)) %>%
  #{{"Fungi" are selected}}
  filter(grepl("Fungi", taxonomy)) %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB_", ignore.case = FALSE)) %>%
  select(-CONT_ROB_PCR1_a_s) %>% #delete control due to mixing in microplate cell?
  droplevels()
  
## -------- 4.2. amf --------

amf_OTUs <- data_18S %>%
  #{{import data}}
  read_tsv() %>%
  #{{format column names}}
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~sub("_AMF_", "", .)) %>%
  rename_with(~sub("\\_AMF.*", "", .)) %>%
  rename_with(~sub("ROB_T_", "CONT_ROB_", .)) %>%
  rename_with(~sub("ROB_Tex", "CONT_ROB_T_ex", .)) %>%
  #{{select samples from ROBUST project}}
  select(OTU, taxonomy, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB", ignore.case = FALSE)) %>%
  #{{"Arbuscular mycorrhizal fungi & Mucoromycotina fine root endophytes" are selected}}
  filter(grepl("Glomeromycotina", taxonomy)|grepl("Endogonales", taxonomy)) %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB_", ignore.case = FALSE)) %>%
  droplevels()

## ---- 5. IMPORT AND FORMAT TAXONOMIC DATA ----
## -------- 5.1. bacteria --------

bact_tax <- input3b %>%
  #{{format taxonomic ranks}}
  fold_table %>%
  filter_redundant_taxa_names %>%
  merge_with_taxo_levels %>%
  remove_duplicated_levels %>%
  unfold_table %>%
  add_species_column %>%
  #{{merge taxonomy}}
  right_join(select(bact_OTUs, OTU), by = "OTU") %>%
  #{{select only OTU number and taxonomy}}
  select(OTU,bact_taxonomic_levels) %>%  
  #{{replace missing taxonomic information and homogeneize the terms}}
  mutate_if(is.character, ~sub("uncultured", "unidentified", .)) %>% 
  mutate_if(is.character, ~sub("uncultured_bacterium", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("unidentified_bacterium", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("[*]", "unidentified", .)) %>%
  mutate_if(is.character , replace_na, replace = "unidentified") %>%
  #{{replace unidentified by the higher taxonomic rank}}
  mutate(phylum = if_else(str_detect(phylum, "unidentified"), 
                          str_c("unidentified_", domain), phylum)) %>%
  mutate(class = if_else(str_detect(class, "unidentified"), 
                         str_c("unidentified_", phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "unidentified"), 
                         str_c("unidentified_", class), order)) %>%
  mutate(family = if_else(str_detect(family, "unidentified"), 
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "unidentified"), 
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "unidentified"), 
                           str_c("unidentified_", genus), species)) %>%
  mutate(phylum = if_else(str_detect(phylum, "metagenome"), 
                          str_c("unidentified_", domain), phylum)) %>%
  mutate(class = if_else(str_detect(class, "metagenome"), 
                         str_c("unidentified_", phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "metagenome"), 
                         str_c("unidentified_", class), order)) %>%
  mutate(family = if_else(str_detect(family, "metagenome"), 
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "metagenome"), 
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "metagenome"),
                           str_c("unidentified_", genus), species)) %>%
  #{{homogenize unidentified terms for taxonomic rank}}
  mutate(class = str_replace(class,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{homogenize terms for unidentified or uncharacterized species (sp)}}
  mutate(species = if_else(str_detect(species, "_sp"), 
                           str_c("unidentified_", genus), species)) %>% 
  mutate(species = str_replace(species, "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{convert OTU column in character from numeric}}
  mutate_if(is.numeric,as.character)

## -------- 5.3. fungi --------

ITS2_tax <- ITS2_OTUs %>% 
  #{{merge taxonomy}}
  left_join(select(read_tsv(data_ITS), OTU, taxonomy), by = "OTU") %>%
  #{{select only OTU number and taxonomy}}
  select(OTU, taxonomy) %>%
  #{{"Fungi" are selected}}
  separate(taxonomy, its_taxonomic_levels, sep = "[|]", extra = "drop") %>%
  #{{remove "k__" bits, etc in taxonomic names and replace missing taxonomic information}}
  mutate_if(is.character, ~sub("^[kpcofgs]__", "", .)) %>%
  mutate_if(is.character, ~sub("\\_phy_.*", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("\\_cls_.*", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("\\_ord_.*", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("\\_fam_.*", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("\\_gen_.*", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("[*]", "unidentified", .)) %>%
  #{{replace unidentified by the higher taxonomic rank}}
  mutate(phylum = if_else(str_detect(phylum, "unidentified"),
                          str_c("unidentified_", kingdom), phylum)) %>%
  mutate(class = if_else(str_detect(class, "unidentified"),
                         str_c("unidentified_",phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "unidentified"),
                         str_c("unidentified_",class), order)) %>%
  mutate(family = if_else(str_detect(family, "unidentified"),
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "unidentified"),
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "unidentified"),
                           str_c("unidentified_",genus), species)) %>%
  #{{homogneize unidentified terms for taxonomic rank}}
  mutate(class = str_replace(class,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{homogenize terms for unidentified or uncharacterized species (sp)}}
  mutate(species = if_else(str_detect(species, "_sp"), 
                           str_c("unidentified_", genus), species)) %>% 
  mutate(species = str_replace(species, "unidentified_unidentified_",
                               "unidentified_")) %>%
  #{{convert OTU column in character from numeric}}
  mutate_if(is.numeric,as.character)

## -------- 5.3. amf --------

amf_tax <- amf_OTUs %>% 
  #{{merge taxonomy}}
  left_join(select(read_tsv(data_18S), OTU, taxonomy), by = "OTU") %>%
  #{{select only OTU number and taxonomy}}
  select(OTU, taxonomy) %>%
  #{{"Fungi" are selected}}
  separate(taxonomy, amf_taxonomic_levels, sep = "[|]", extra = "drop") %>%
  #{{remove "k__" bits, etc in taxonomic names and replace missing taxonomic information}}
  mutate_if(is.character, ~sub("^[kdpcofgs]:", "", .)) %>%
  mutate_if(is.character, ~sub("Glomeromycotina_XX_sp.", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("Glomeromycotina_XX", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("Glomeromycotina_X", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("[*]", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("Diversiporales", "Diversisporales", .)) %>%
  #{{replace unidentified by the higher taxonomic rank}}
  mutate(phylum = if_else(str_detect(phylum, "unidentified"),
                          str_c("unidentified_", kingdom), phylum)) %>%
  mutate(class = if_else(str_detect(class, "unidentified"),
                         str_c("unidentified_",phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "unidentified"),
                         str_c("unidentified_",class), order)) %>%
  mutate(family = if_else(str_detect(family, "unidentified"),
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "unidentified"),
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "unidentified"),
                           str_c("unidentified_",genus), species)) %>%
  #{{homogneize unidentified terms for taxonomic rank}}
  mutate(class = str_replace(class,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{homogenize terms for unidentified or uncharacterized species (sp)}}
  mutate(species = if_else(str_detect(species, "_sp"), 
                           str_c("unidentified_", genus), species)) %>% 
  mutate(species = str_replace(species, "unidentified_unidentified_",
                               "unidentified_")) %>%
  #{{convert OTU column in character from numeric}}
  mutate_if(is.numeric,as.character)

## ---- 6. DECONTAMINATION OF DATA ----

#{{For now, I decide to compute the total number of reads for each OTU
# in control samples. That sum is then substracted from occurrences
# of that OTU in true samples. OTUs present in control samples. The
# rationale is as follows:
# - strong in control samples, weak in true samples (contamination
#   specific to the control samples, will be eliminated by the
#   substraction, i.e, final abundance is 0),
# - present in control samples, present in true samples (systematic
#   contamination, will be mitigated by the substraction)
# - weak in control samples, strong in true samples (cross-talk, will
#   be eliminated/mitigated by the substraction)
# Finally, control samples should be eliminated from the statistical
# analysis now that the substraction is done}}

## -------- 6.1. bacteria --------

#{{Extract control samples}}
d <- replace(bact_OTUs, bact_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  gather("samples", "reads", -OTU) %>%
  filter(!is.na(reads))

#{{Extract samples}}
bact_OTUs_decont <- replace(bact_OTUs, bact_OTUs == 0, NA) %>%
  gather("samples", "reads", -c("OTU", starts_with("CONT"))) %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  spread(samples, reads, fill = 0) %>%
  select(-starts_with("CONT")) 

#{{delete previous data}}
rm(d)
rm(bact_OTUs)

## -------- 6.3. fungi --------

#{{Extract control samples}}
d <- replace(ITS2_OTUs, ITS2_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  gather("samples", "reads", -OTU) %>%
  filter(!is.na(reads))

ITS2_OTUs_decont <- replace(ITS2_OTUs, ITS2_OTUs == 0, NA) %>%
  gather("samples", "reads", -c("OTU", starts_with("CONT"))) %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  spread(samples, reads, fill = 0) %>%
  select(-starts_with("CONT"))

#{{delete previous data}}
rm(d)
rm(ITS2_OTUs)

## -------- 6.3. amf --------

#{{Extract control samples}}
d <- replace(amf_OTUs, amf_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  gather("samples", "reads", -OTU) %>%
  filter(!is.na(reads))

amf_OTUs_decont <- replace(amf_OTUs, amf_OTUs == 0, NA) %>%
  gather("samples", "reads", -c("OTU", starts_with("CONT"))) %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  spread(samples, reads, fill = 0) %>%
  select(-starts_with("CONT"))

#{{delete previous data}}
rm(d)
rm(amf_OTUs)

## ---- 7. FORMAT OTU DATA FOR RAREFACTION ----
## -------- 7.1. bacteria --------

bact_OTUs_t <- bact_OTUs_decont %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

#{{delete previous data}}
rm(bact_OTUs_decont)

## -------- 7.3. fungi --------

ITS2_OTUs_t <- ITS2_OTUs_decont %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

#{{delete previous data}}
rm(ITS2_OTUs_decont)

## -------- 7.4. amf --------

amf_OTUs_t <- amf_OTUs_decont %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

#{{delete previous data}}
rm(amf_OTUs_decont)

## ---- 8. SELECT DATA RAREFACTION THRESHOLD ----
## -------- 8.1. bacteria -------- 

#{{Identification of Low samples and outliers + rarefaction threshold with smallest}}
quantile(rowSums(bact_OTUs_t))
min(rowSums(bact_OTUs_t))
sort(rowSums(bact_OTUs_t))
head(sort(rowSums(bact_OTUs_t)))
smallest_bact <- min(rowSums(bact_OTUs_t))

#{{samples to delete, lower number of sequences}}
#Threshold = 499 séquences = SU_USFLO_CL_02_l
#{{sample to delete, lower number of sequences}
bact_OTUs_t <- bact_OTUs_t %>% 
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{sample to exclude}}
  select(-SU_USFLO_CL_09_l, -SU_USFLO_LS_07_sk, -SU_USFLO_LS_20_l, -SU_USFLO_LS_22_l,
         -SU_USFLO_LS_16_l, -SU_USFLO_LS_07_l, -SU_USFLO_CL_14_sk, -SU_USFLO_CL_12_sk,
         -SU_USFLO_LS_14_sk, -SU_USFLO_CL_20_sk, -SU_USFLO_CL_16_sk, -SU_USFLO_CL_13_l,
         -SU_USFLO_CL_19_sk, -SU_USFLO_CL_14_l, -SU_USFLO_CL_15_sk, -SU_USFLO_CL_24_l,
         -SU_USFLO_LS_19_sk, -SU_USFLO_LS_25_sk, -SU_USFLO_CL_18_sk, -SU_USFLO_CL_22_sk,
         -SU_USFLO_CL_21_sk, -SU_USFLO_LS_12_sk, -SU_USFLO_CL_25_l, -SU_USFLO_CL_23_l,
         -SU_USFLO_CL_22_l, -SU_USFLO_CL_25_sk, -SU_USFLO_CL_20_l, -SU_USFLO_LS_25_l,
         -SU_USFLO_CL_02_r, -SU_USFLO_CL_03_r, -SU_USFLO_CL_06_l, -SU_USFLO_CL_08_sk,
         -SU_USFLO_CL_06_r, -SU_USFLO_CL_21_l, -SU_USFLO_CL_06_sk, -SU_USFLO_LS_24_l,
         -SU_USFLO_CL_07_r, -SU_USFLO_CL_01_sk, -SU_USFLO_LS_23_l, -SU_USFLO_LS_16_sk,
         -SU_USFLO_CL_12_l, -SU_USFLO_LS_17_sk, -SU_USFLO_CL_24_sk, -SU_USFLO_CL_01_r,
         -SU_USFLO_CL_10_r, -SU_USFLO_CL_01_l, -SU_USFLO_LS_20_sk, -SU_USFLO_CL_04_r,
         -SU_USFLO_CL_05_l) %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

quantile(rowSums(bact_OTUs_t))
min(rowSums(bact_OTUs_t))
sort(rowSums(bact_OTUs_t))
head(sort(rowSums(bact_OTUs_t)))
smallest_bact <- min(rowSums(bact_OTUs_t))

## -------- 8.2. fungi --------

#{{Identification of Low samples and outliers + rarefaction threshold with smallest}}
quantile(rowSums(ITS2_OTUs_t))
min(rowSums(ITS2_OTUs_t))
sort(rowSums(ITS2_OTUs_t))
head(sort(rowSums(ITS2_OTUs_t)))
smallest_its2 <- min(rowSums(ITS2_OTUs_t))

#{{sample to delete, lower number of sequences}}
#Threshold = 10835 séquences = ROB_L_s

ITS2_OTUs_t <- ITS2_OTUs_t %>% 
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{sample to exclude}}
  select(-ROB_M4_s, -ROB_M5_s, -ROB_M2_s, -ROB_M3_s, -ROB_B22_s, -ROB_M21_s,
         -ROB_B25_s, -ROB_M1_s, -ROB_M25_s, -ROB_M22_s, -ROB_M24_s, -ROB_M15_s,
         -ROB_M23_s, -ROB_F19_s, -ROB_F18_s, -ROB_N25_s, -ROB_M7_s, -ROB_B11_s,
         -ROB_R11_s, -ROB_N10_s,  -ROB_B9_s, -ROB_R13_s, -ROB_R12_s) %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

quantile(rowSums(ITS2_OTUs_t))
min(rowSums(ITS2_OTUs_t))
sort(rowSums(ITS2_OTUs_t))
head(sort(rowSums(ITS2_OTUs_t)))
smallest_its2 <- min(rowSums(ITS2_OTUs_t))

## -------- 8.3. amf --------

#{{Identification of Low samples and outliers + rarefaction threshold with smallest}}
quantile(rowSums(amf_OTUs_t))
min(rowSums(amf_OTUs_t))
sort(rowSums(amf_OTUs_t))
head(sort(rowSums(amf_OTUs_t)))
smallest_amf <- min(rowSums(amf_OTUs_t))

#{{sample to delete, lower number of sequences}}
#Threshold = 791 séquences = ROB_N7_s

amf_OTUs_t <- amf_OTUs_t %>% 
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{sample to exclude}}
  select(-ROB_M1_s, -ROB_M25_s, -ROB_N14_s,  -ROB_B8_s, -ROB_M21_s, -ROB_M23_s,
         -ROB_M3_s, -ROB_N10_s, -ROB_N11_s, -ROB_N13_s, -ROB_N15_s, -ROB_N16_s,
         -ROB_N25_s, -ROB_R19_s,  -ROB_R2_s, -ROB_M24_s, -ROB_B25_s, -ROB_M20_s,
         -ROB_B2_s, -ROB_M14_s, -ROB_M17_s,  -ROB_R8_s, -ROB_M15_s,  -ROB_B5_s,
         -ROB_H_s, -ROB_F10_s, -ROB_F11_s, -ROB_B24_s, -ROB_M16_s) %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

quantile(rowSums(amf_OTUs_t))
min(rowSums(amf_OTUs_t))
sort(rowSums(amf_OTUs_t))
head(sort(rowSums(amf_OTUs_t)))
smallest_amf <- min(rowSums(amf_OTUs_t))

## ---- 9. DATA RAREFACTION ------
## -------- 9.1. bacteria --------

#{{all data}}
{
set.seed(seed)
#{{rarefaction}}
bact_OTUs_rarefied_t <- bact_OTUs_t %>%
  rrarefy(smallest_bact) %>%
  as_tibble(rownames = NA) %>%
  #{{create a new column with sum of reads for each OTU in all samples}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{delete all sample down to 0 reads in rarefied samples}}
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  filter(total != 0) %>%
  select(-total) %>%
  column_to_rownames(var = "OTU") %>% 
  t() %>%
  as_tibble(rownames = NA)
}

#{{export data}}
{
  output <- "bact_OTUs_all_plant_final.txt"
  
  tmp <- bact_OTUs_rarefied_t %>%
    rownames_to_column(var = "soil_code")
  
  write_delim(tmp, file = output, delim = "\t", col_names = TRUE)
}

#{{delete previous data}}
rm(bact_OTUs_t) 

## -------- 9.3. fungi --------

#{{all data}}
{
  set.seed(seed)
  #{{rarefaction}}
  ITS2_OTUs_rarefied_t <- ITS2_OTUs_t %>%
    rrarefy(smallest_its2) %>%
    as_tibble(rownames = NA) %>%
    #{{create a new column with sum of reads for each OTU in all samples}}
    t() %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "OTU") %>%
    #{{delete all sample down to 0 reads in rarefied samples}}
    mutate(total = rowSums(across(where(is.numeric)))) %>% 
    filter(total != 0) %>%
    select(-total) %>%
    column_to_rownames(var = "OTU") %>% 
    t() %>%
    as_tibble(rownames = NA)
}

#{{export data}}
{
  output <- "fung_OTUs_all_soil_final.txt"
  
  tmp <- ITS2_OTUs_rarefied_t %>%
    rownames_to_column(var = "soil_code")
  
  write_delim(tmp, file = output, delim = "\t", col_names = TRUE)
}

#{{delete previous data}}
rm(ITS2_OTUs_t)  

## -------- 9.3. amf --------

#{{all data}}
{
  set.seed(seed)
  #{{rarefaction}}
  amf_OTUs_rarefied_t <- amf_OTUs_t %>%
    rrarefy(smallest_amf) %>%
    as_tibble(rownames = NA) %>%
    #{{create a new column with sum of reads for each OTU in all samples}}
    t() %>% 
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "OTU") %>%
    #{{delete all sample down to 0 reads in rarefied samples}}
    mutate(total = rowSums(across(where(is.numeric)))) %>% 
    filter(total != 0) %>%
    select(-total) %>%
    column_to_rownames(var = "OTU") %>% 
    t() %>%
    as_tibble(rownames = NA)
}

#{{export data}}
{
  output <- "amf_OTUs_all_soil_final.txt"
  
  tmp <- amf_OTUs_rarefied_t %>%
    rownames_to_column(var = "soil_code")
  
  write_delim(tmp, file = output, delim = "\t", col_names = TRUE)
}

#{{delete previous data}}
rm(amf_OTUs_t)  

## ---- 10. ALPHA-DIVERSITY ----
## -------- 10.1. bacteria --------
## ------------ 10.1.1. estimate usual diversity index -----------

#{{richness}}
R <- bact_OTUs_rarefied_t %>%
  estimateR() %>% 
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>% 
  select(soil_code, S.obs)

#{{Shannon index H; richness + evenness}}
H <- bact_OTUs_rarefied_t %>% 
  diversity(index = "shannon", MARGIN = 1, base =
              exp(1)) %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{Pielou's index of evenness; 0-1, 1 = max. evenness}}
S <- bact_OTUs_rarefied_t %>% 
  specnumber()
J <- H %>% 
  mutate(value = value/log(S))

#{{Simpson's D index; richness + evenness; 0-1; 1 - D rises as evenness increases}}
D <- bact_OTUs_rarefied_t %>% 
  diversity("simpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")
inv_D <- bact_OTUs_rarefied_t %>% 
  diversity("invsimpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 10.1.2. estimate effective number of species - Hill numbers -----------

#This equation has a parameter q that defines its sensitivity to rare species:
#low values of q favor rare species, high values of q favor abundant species. 
#For example, Shannon diversity is of order q = 1, and for Simpson diversity q = 2. 
#When q = 0, diversity = S (richness), because rare species are treated the same as abundant ones.

##{{all data}}
{
##{{richness as hill numbers}}
dR <- bact_OTUs_rarefied_t %>% 
  t() %>% 
  hill_div(qvalue = 0) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{richness + evenness (= shannon diversity) as hill numbers; rare species}}
dREr <- bact_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 1) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

##{{richness + evenness (= inverse Simpson = simpson diversity) as hill numbers; abundant species}}
dREa <- bact_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 2) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{evenness (= shannon evenness) as hill ratio}}
eDRr <- dREr %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)

#{{evenness (= ssimpson evenness) as hill ratio}}
eDRa  <- dREa %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)
}

## ------------ 10.1.3. format diversity table (usual index + Hill numbers) -----------

#{{merge data}}
Indices_bact <- R %>%
  left_join(H, by = "soil_code") %>% 
  left_join(inv_D, by = "soil_code") %>%
  left_join(J, by = "soil_code") %>%
  left_join(dR, by = "soil_code") %>%
  left_join(dREr, by = "soil_code") %>%
  left_join(dREa, by = "soil_code") %>%
  left_join(eDRr, by = "soil_code") %>%
  left_join(eDRa, by = "soil_code")
  
  names(Indices_bact) <- c("soil_code", "Richness", "Shannon", "Inv_Simpson", "Pielou", 
                         "Hill_Richness","Hill_Shannon","Hill_Inv_Simpson", 
                         "Hill_Shannon_evenness", "Hill_Simpson_evenness")

## -------- 10.3. fungi --------
## ------------ 10.3.1. estimate usual diversity index -----------

#{{richness}}
R <- ITS2_OTUs_rarefied_t %>%
  estimateR() %>% 
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>% 
  select(soil_code, S.obs)

#{{Shannon index H; richness + evenness}}
H <- ITS2_OTUs_rarefied_t %>% 
  diversity(index = "shannon", MARGIN = 1, base =
              exp(1)) %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{Pielou's index of evenness; 0-1, 1 = max. evenness}}
S <- ITS2_OTUs_rarefied_t %>% 
  specnumber()
J <- H %>% 
  mutate(value = value/log(S))

#{{Simpson's D index; richness + evenness; 0-1; 1 - D rises as evenness increases}}
D <- ITS2_OTUs_rarefied_t %>% 
  diversity("simpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")
inv_D <- ITS2_OTUs_rarefied_t %>% 
  diversity("invsimpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 10.3.2. estimate effective number of species - Hill numbers -----------

#This equation has a parameter q that defines its sensitivity to rare species:
#low values of q favor rare species, high values of q favor abundant species. 
#For example, Shannon diversity is of order q = 1, and for Simpson diversity q = 2. 
#When q = 0, diversity = S (richness), because rare species are treated the same as abundant ones.

##{all data}
{
##{{richness as hill numbers}}
dR <- ITS2_OTUs_rarefied_t %>% 
  t() %>% 
  hill_div(qvalue = 0) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{richness + evenness (= Shannon entropy) as hill numbers; rare species}}
dREr <- ITS2_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 1) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

##r{{richness + evenness (= inverse Simpson) as hill numbers; abundant species}}
dREa <- ITS2_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 2) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{evenness (= shannon evenness) as hill ratio}}
eDRr <- dREr %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)

#{{evenness (= ssimpson evenness) as hill ratio}}
eDRa  <- dREa %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)
}

## ------------ 10.3.3. format diversity table (usual index + Hill numbers) -----------

#{{merge data}}
Indices_its2 <- R %>%
  left_join(H, by = "soil_code") %>% 
  left_join(inv_D, by = "soil_code") %>%
  left_join(J, by = "soil_code") %>%
  left_join(dR, by = "soil_code") %>%
  left_join(dREr, by = "soil_code") %>%
  left_join(dREa, by = "soil_code") %>%
  left_join(eDRr, by = "soil_code") %>%
  left_join(eDRa, by = "soil_code")

names(Indices_its2) <- c("soil_code", "Richness", "Shannon", "Inv_Simpson", "Pielou", 
                         "Hill_Richness","Hill_Shannon","Hill_Inv_Simpson", 
                         "Hill_Shannon_evenness", "Hill_Simpson_evenness")

## -------- 10.4. amf --------
## ------------ 10.4.1. estimate usual diversity index -----------

#{{richness}}
R <- amf_OTUs_rarefied_t %>%
  estimateR() %>% 
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>% 
  select(soil_code, S.obs)

#{{Shannon index H; richness + evenness}}
H <- amf_OTUs_rarefied_t %>% 
  diversity(index = "shannon", MARGIN = 1, base =
              exp(1)) %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{Pielou's index of evenness; 0-1, 1 = max. evenness}}
S <- amf_OTUs_rarefied_t %>% 
  specnumber()
J <- H %>% 
  mutate(value = value/log(S))

#{{Simpson's D index; richness + evenness; 0-1; 1 - D rises as evenness increases}}
D <- amf_OTUs_rarefied_t %>% 
  diversity("simpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")
inv_D <- amf_OTUs_rarefied_t %>% 
  diversity("invsimpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 10.3.2. estimate effective number of species - Hill numbers -----------

#This equation has a parameter q that defines its sensitivity to rare species:
#low values of q favor rare species, high values of q favor abundant species. 
#For example, Shannon diversity is of order q = 1, and for Simpson diversity q = 2. 
#When q = 0, diversity = S (richness), because rare species are treated the same as abundant ones.

##{all data}
{
  ##{{richness as hill numbers}}
  dR <- amf_OTUs_rarefied_t %>% 
    t() %>% 
    hill_div(qvalue = 0) %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "soil_code")
  
  #{{richness + evenness (= Shannon entropy) as hill numbers; rare species}}
  dREr <- amf_OTUs_rarefied_t %>% 
    t() %>%
    hill_div(qvalue = 1) %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "soil_code")
  
  ##r{{richness + evenness (= inverse Simpson) as hill numbers; abundant species}}
  dREa <- amf_OTUs_rarefied_t %>% 
    t() %>%
    hill_div(qvalue = 2) %>%
    as_tibble(rownames = NA) %>%
    rownames_to_column(var = "soil_code")
  
  #{{evenness (= shannon evenness) as hill ratio}}
  eDRr <- dREr %>%
    left_join(dR, by = "soil_code") %>%
    mutate(value = value.x / value.y) %>%
    select(soil_code,value)
  
  #{{evenness (= ssimpson evenness) as hill ratio}}
  eDRa  <- dREa %>%
    left_join(dR, by = "soil_code") %>%
    mutate(value = value.x / value.y) %>%
    select(soil_code,value)
}

## ------------ 10.3.3. format diversity table (usual index + Hill numbers) -----------

#{{merge data}}
Indices_amf <- R %>%
  left_join(H, by = "soil_code") %>% 
  left_join(inv_D, by = "soil_code") %>%
  left_join(J, by = "soil_code") %>%
  left_join(dR, by = "soil_code") %>%
  left_join(dREr, by = "soil_code") %>%
  left_join(dREa, by = "soil_code") %>%
  left_join(eDRr, by = "soil_code") %>%
  left_join(eDRa, by = "soil_code")

names(Indices_amf) <- c("soil_code", "Richness", "Shannon", "Inv_Simpson", "Pielou", 
                         "Hill_Richness","Hill_Shannon","Hill_Inv_Simpson", 
                         "Hill_Shannon_evenness", "Hill_Simpson_evenness")

## ---- 11. FORMAT FUNCTIONAL MICROBIAL DATA ----
## -------- 11.1. bacteria ----
## ------------ 11.1.5. FAPROTAX ----

#{{format data for faprotax prediction}}
bact_OTUs_tax <- bact_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  left_join(bact_tax, by = "OTU") %>%
  select(-OTU) %>%
  unite(taxonomy, bact_taxonomic_levels, sep = ";")

write_delim(bact_OTUs_tax, file = "tax_table.tsv", delim = "\t", col_names = TRUE)
rm(bact_OTUs_tax)

#{{open a command line}}
#Using FAPROTAX 1.2.7 with classical (tabular) taxon tables
# type d: in a command line to switch in D:\>
# enter the directory path
#To convert a classical taxon table (tab-separated format) into a function table, use:
# type 	collapse_table.py -i tax_table.tsv -o func_table.tsv -g FAPROTAX.txt -d "taxonomy" -c "#" -v 

#{{import faprotax data}}
bact_faprotax <- read_tsv("func_table.tsv")

## -------- 11.2. fungi ----
## ------------ 11.2.1. FUNguild ----

#{{import FUNguild database}}
funguild <- parse_funguild()
attr(funguild, "DownloadDate") #to indicate the database dowloading

#{{format funguild database}}
funguild <- funguild %>%
  mutate_if(is.factor, ~as.character(.)) %>% 
  mutate(taxonomicLevel = case_when(
    taxonomicLevel == "Subspecies" ~ "Species",
    taxonomicLevel == "Variety" ~ "Species",
    taxonomicLevel == "Form" ~ "Species", 
    TRUE ~ taxonomicLevel)) %>% 
  #{{create genus column, including higher taxonomic range}}
  mutate(genus = if_else(taxonomicLevel == "Genus", 
                         taxon, word(taxon, 1, sep = fixed(" ")))) %>%
  mutate(genus = case_when(
    taxonomicLevel == "Family" ~ paste("unidentified_", genus),
    taxonomicLevel == "Order" ~ paste("unidentified_", genus),
    taxonomicLevel == "Phylum" ~ paste("unidentified_", genus),
    TRUE ~ genus)) %>% 
  #{{create species column}}
  mutate(species = if_else(taxonomicLevel == "Species", taxon, 
                           paste("unidentified_", word(taxon, 1, sep = fixed(" ")),
                                 sep = ""))) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>% 
  select(-taxon, -taxonomicLevel) %>% 
  relocate(c(genus,species))

#{{Merge OTU taxonomic data and FUNguild database}}
#{{need to screen all taxonomic range because funguild assignement is based on different taxonomic range }}
#{{some species unassigned or poorly-assigned were re-assigned}}
funcoffee.guilds <- ITS2_tax %>% 
  left_join(select(funguild, genus, trophicMode, guild), by = c("genus" = "genus")) %>%
  left_join(select(funguild, species, trophicMode, guild), by = c("species" = "species")) %>% 
  mutate(trophic_mode = case_when(
    is.na(trophicMode.x) ~ trophicMode.y,
    is.na(trophicMode.y) ~ trophicMode.x,
    species == "Albifimbria_verrucaria" ~ "Pathotroph-Saprotroph-Symbiotroph",
    species == "Akenomyces_costatus" ~ "Saprotroph",
    species == "Alfaria_dandenongensis" ~ "Symbiotroph",
    species == "Allophaeosphaeria cytisi" ~ "Saprotroph",
    species == "Amandinea_punctata" ~ "Symbiotroph",
    species == "Anguillospora_crassa" ~ "Symbiotroph",
    species == "Apseudocercosporella_trigonotidis" ~ "Pathotroph",
    species == "Aquamyces chlorogonii" ~ "Fungal Parasite",
    species == "Aquamyces_chlorogonii" ~ "Pathotroph",
    species == "Backusella_lamprospora" ~ "Symbiotroph",
    species == "Balansia_cyperi" ~ "Pathotroph-Symbiotroph",
    species == "Basidioascus_magus" ~ "Saprotroph",
    species == "Biatriospora_mackinnonii" ~ "Pathotroph",
    species == "Bifiguratus_adelaidae" ~ "Saprotroph-Symbiotroph",
    species == "Bonomyces_sinopicus" ~ "Saprotroph",
    species == "Caloboletus_radicans" ~ "Symbiotroph",
    species == "Chaetospermum_chaetosporum" ~ "Saprotroph",
    species == "Chlorophyllum_brunneum" ~ "Saprotroph",
    species == "Clitocella_popinalis" ~ "Saprotroph",
    species == "Contumyces_rosellus" ~ "Saprotroph",
    species == "Cryphonectria_parasitica" ~ "Pathotroph", #chancre du châtaigner
    species == "Darksidea_alpha" ~ "Symbiotroph",
    species == "Darksidea_beta" ~ "Symbiotroph",
    species == "Darksidea_epsilon" ~ "Symbiotroph",
    species == "Darksidea_zeta" ~ "Symbiotroph",
    species == "Dematiopleospora_fusiformis" ~ "Saprotroph",
    species == "Dermoloma_cuneifolium" ~ "Saprotroph",
    species == "Dermoloma josserandii" ~ "Saprotroph",
    species == "Dermoloma_pseudocuneifolium" ~ "Saprotroph",
    species == "Derxomyces_komagatae" ~ "Saprotroph",
    species == "Desmazierella_acicola" ~ "Symbiotroph",
    species == "Dioszegia_changbaiensis" ~ "Symbiotroph",
    species == "Dioszegia_fristingensis" ~ "Symbiotroph",
    species == "Dioszegia_zsoltii_var._yunnanensis" ~ "Symbiotroph",
    species == "Dioszegia_zsoltii_var._zsoltii" ~ "Symbiotroph",
    species == "Diromma_dirinellum" ~ "Symbiotroph",
    species == "Discinella_terrestris" ~ "Saprotroph",
    species == "Discosia_neofraxinea" ~ "Symbiotroph",
    species == "Dominikia_disticha" ~ "Symbiotroph",
    species == "Dominikia_iranica" ~ "Symbiotroph",
    species == "Drechmeria_balanoides" ~ "Pathotroph",
    species == "Drechmeria_zeospora" ~ "Pathotroph",
    species == "Ectophoma_multirostrata" ~ "Pathotroph",
    species == "Epichloe_amarillans" ~ "Symbiotroph",
    species == "Fitzroyomyces_cyperacearum" ~ "Saprotroph",
    species == "Flavomyces_fulophazii" ~ "Symbiotroph",
    species == "Fusarium_concentricum" ~ "Pathotroph",
    species == "Gastrosporium_simplex" ~ "Pathotroph",
    species == "Geosmithia microcorthyli" ~ "Saprotroph-Symbiotroph",
    species == "Hawksworthiomyces_hibbettii" ~ "Saprotroph",
    species == "Hawksworthiomyces_lignivorus" ~ "Saprotroph",
    species == "Hawksworthiomyces_taylorii" ~ "Saprotroph",
    species == "Heterospora_chenopodii" ~ "Pathotroph-Saprotroph",
    species == "Humicola_fuscoatra" ~ "Pathotroph",
    species == "Humicola_grisea" ~ "Saprotroph",
    species == "Hysterium_angustatum" ~ "Saprotroph",
    species == "Hysterium_pulicare" ~ "Saprotroph",
    species == "Itersonilia_pannonica" ~ "Pathotroph-Symbiotroph",
    species == "Itersonilia_perplexans" ~ "Pathotroph-Symbiotroph",
    species == "Jobellisia_guangdongensis" ~ "Saprotroph",
    species == "Kamienskia_bistrata" ~ "Symbiotroph",
    species == "Kamienskia_perpusilla" ~ "Symbiotroph",
    species == "Leratiomyces_erythrocephalus" ~ "Saprotroph",
    species == "Liberomyces_macrosporus" ~ "Symbiotroph",
    species == "Liberomyces_saliciphilus" ~ "Symbiotroph",
    species == "Lycoperdon_nigrescens" ~ "Saprotroph",
    species == "Magnaporthiopsis_poae" ~ "Pathotroph",
    species == "Mitrula_elegans" ~ "Saprotroph",
    species == "Murilentithecium_clematidis" ~ "Saprotroph",
    species == "Mycenella_trachyspora" ~ "Saprotroph",
    species == "Naganishia_albida" ~ "Pathotroph",
    species == "Naohidea_sebacea" ~ "Pathotroph",
    species == "Neoascochyta_europaea" ~ "Pathotroph", #leaf scorch on wheat
    species == "Neoascochyta_graminicola" ~ "Pathotroph", #leaf scorch on wheat
    species == "Neopestalotiopsis_asiatica" ~ "Pathotroph", #leaf spot on almond
    species == "Neopestalotiopsis_musae" ~ "Pathotroph", #Leaf bright disease on banana
    species == "Oehlia_diaphana" ~ "Symbiotroph",
    species == "Pachyphlodes_citrina" ~ "Symbiotroph",
    species == "Paralepista_flaccida" ~ "Saprotroph",
    species == "Paraboeremia_selaginellae" ~ "Pathotroph",
    species == "Pseudogymnoascus_pannorum" ~ "Pathotroph",
    species == "Pseudopithomyces_rosae" ~ "Pathotroph",
    species == "Pseudoteratosphaeria_flexuosa" ~ "Pathotroph", #leaf spot on eucalypts
    species == "Pseudoteratosphaeria_ohnowa" ~ "Pathotroph", #leaf spot on eucalypts
    species == "Rubroboletus_pulcherrimus" ~ "Symbiotroph",
    species == "Sarcogyne_regularis" ~ "Symbiotroph",
    species == "Sarcoscypha_coccinea" ~ "Saprotroph",
    species == "Septoriella hirta" ~ "Pathotroph",
    species == "Setophoma_terrestris" ~ "Pathotroph", #Corky and Pink Root of Tomato
    species == "Siphula_ceratites" ~ "Symbiotroph",
    species == "Sistotrema_muscicola" ~ "Symbiotroph",
    species == "Spiroplana_centripeta" ~ "Pathotroph",
    species == "Stereocaulon_glabrum" ~ "Symbiotroph",
    species == "Suillellus_queletii" ~ "Symbiotroph",
    species == "Teratosphaeria_viscida" ~ "Pathotroph", #leaf disease on eucalypts
    species == "Teratosphaericola_pseudoafricana" ~ "Pathotroph", #leaf disease on eucalypts
    species == "Thamnidium_elegans" ~ "Saprotroph",
    species == "Thelephora_palmata" ~ "Symbiotroph",
    species == "Tuber_aestivum" ~ "Symbiotroph",
    species == "Tuber_maculatum" ~ "Symbiotroph",
    species == "Tuber_brumale" ~ "Symbiotroph",
    species == "Tuber_excavatum" ~ "Symbiotroph",
    species == "Tuber_gennadii" ~ "Symbiotroph",
    species == "Tuber_huidongense" ~ "Symbiotroph",
    species == "Tuber_melanosporum" ~ "Symbiotroph",
    species == "Tuber_rufum" ~ "Symbiotroph",
    species == "Tuber_vesicoperidium" ~ "Symbiotroph",
    species == "unidentified_Ambisporaceae" ~ "Symbiotroph",
    species == "unidentified_Claroideoglomeraceae" ~ "Symbiotroph",
    species == "unidentified_Glomeromycetes" ~ "Symbiotroph",
    species == "unidentified_GS24" ~ "Symbiotroph",
    species == "unidentified_Tuber" ~ "Symbiotroph",
    species == "Xerocomellus_chrysenteron" ~ "Symbiotroph",
    species == "Xerocomellus_porosporus" ~ "Symbiotroph",
    species == "Xerocomellus_redeuilhii" ~ "Symbiotroph",
    species == "Zygophlyctis_melosirae" ~ "Pathotroph",
    species == "Zymoseptoria_brevis" ~ "Pathotroph",
    TRUE ~ trophicMode.y)) %>% 
  mutate(guild_mode = case_when(
    is.na(guild.x) ~ guild.y,
    is.na(guild.y) ~ guild.x,
    species == "Albifimbria_verrucaria" ~ "Animal Pathogen-Endophyte-Plant Pathogen-Soil Saprotroph",
    species == "Akenomyces_costatus" ~ "Undefined Saprotroph",
    species == "Alfaria_dandenongensis" ~ "Endophyte",
    species == "Allophaeosphaeria cytisi" ~ "Wood Saprotroph",
    species == "Amandinea_punctata" ~ "Lichenized",
    species == "Anguillospora_crassa" ~ "Endophyte",
    species == "Apseudocercosporella_trigonotidis" ~ "Plant Pathogen",
    species == "Aquamyces_chlorogonii" ~ "Fungal Parasite",
    species == "Backusella_lamprospora" ~ "Endophyte",
    species == "Balansia_cyperi" ~ "Epiphyte-Plant Pathogen",
    species == "Basidioascus_magus" ~ "Undefined Saprotroph",
    species == "Biatriospora_mackinnonii" ~ "Animal Pathogen",
    species == "Bifiguratus_adelaidae" ~ "Endophyte-Undefined Saprotroph",
    species == "Bonomyces_sinopicus" ~ "Soil Saprotroph",
    species == "Caloboletus_radicans" ~ "Ectomycorrhizal",
    species == "Chaetospermum_chaetosporum" ~ "Undefined Saprotroph",
    species == "Chlorophyllum_brunneum" ~ "Soil Saprotroph",
    species == "Clitocella_popinalis" ~ "Soil Saprotroph",
    species == "Contumyces_rosellus" ~ "Wood Saprotroph",
    species == "Cryphonectria_parasitica" ~ "Plant pathogen", #chancre du châtaigner
    species == "Darksidea_alpha" ~ "Endophyte",
    species == "Darksidea_beta" ~ "Endophyte",
    species == "Darksidea_epsilon" ~ "Endophyte",
    species == "Darksidea_zeta" ~ "Endophyte",
    species == "Dematiopleospora_fusiformis" ~ "Wood Saprotroph",
    species == "Dermoloma_cuneifolium" ~ "Soil Saprotroph",
    species == "Dermoloma josserandii" ~ "Soil Saprotroph",
    species == "Dermoloma_pseudocuneifolium" ~ "Soil Saprotroph",
    species == "Derxomyces_komagatae" ~ "Leaf Saprotroph",
    species == "Desmazierella_acicola" ~ "Endophyte-Epiphyte",
    species == "Dioszegia_changbaiensis" ~ "Epiphyte",
    species == "Dioszegia_fristingensis" ~ "Epiphyte",
    species == "Dioszegia_zsoltii_var._yunnanensis" ~ "Epiphyte",
    species == "Dioszegia_zsoltii_var._zsoltii" ~ "Epiphyte",
    species == "Diromma_dirinellum" ~ "Lichenized",
    species == "Discinella_terrestris" ~ "Wood Saprotroph",
    species == "Discosia_neofraxinea" ~ "Epiphyte",
    species == "Dominikia_disticha" ~ "Arbuscular Mycorrhizal",
    species == "Dominikia_iranica" ~ "Arbuscular Mycorrhizal",
    species == "Drechmeria_balanoides" ~ "Fungal Parasite",
    species == "Drechmeria_zeospora" ~ "Fungal Parasite",
    species == "Ectophoma_multirostrata" ~ "Plant Pathogen",
    species == "Epichloe_amarillans" ~ "Endophyte",
    species == "Fitzroyomyces_cyperacearum" ~ "Undefined Saprotroph",
    species == "Flavomyces_fulophazii" ~ "Endophyte",
    species == "Fusarium_concentricum" ~ "Plant Pathogen",
    species == "Gastrosporium_simplex" ~ "Fungal parasite",
    species == "Geosmithia microcorthyli" ~ "Animal Endosymbiont-Wood Saprotroph",
    species == "Hawksworthiomyces_hibbettii" ~ "Wood Saprotroph",
    species == "Hawksworthiomyces_lignivorus" ~ "Wood Saprotroph",
    species == "Hawksworthiomyces_taylorii" ~ "Wood Saprotroph",
    species == "Heterospora_chenopodii" ~ "Plant Pathogen-Undefined Saprotroph",
    species == "Humicola_fuscoatra" ~ "Fungal Parasite",
    species == "Humicola_grisea" ~ "Plant Pathogen",
    species == "Hysterium_angustatum" ~ "Wood Saprotroph",
    species == "Hysterium_pulicare" ~ "Wood Saprotroph",
    species == "Itersonilia_pannonica" ~ "Epiphyte-Plant Pathogen",
    species == "Itersonilia_pannonica" ~ "Epiphyte-Plant Pathogen",
    species == "Jobellisia_guangdongensis" ~ "Undefined Saprotroph",
    species == "Kamienskia_bistrata" ~ "Arbuscular Mycorrhizal",
    species == "Kamienskia_perpusilla" ~ "Arbuscular Mycorrhizal",
    species == "Leratiomyces_erythrocephalus" ~ "Litter Saprotroph",
    species == "Liberomyces_macrosporus" ~ "Endophyte",
    species == "Liberomyces_saliciphilus" ~ "Endophyte",
    species == "Lycoperdon_nigrescens" ~ "Soil Saprotroph",
    species == "Magnaporthiopsis_poae" ~ "Plant Pathogen",
    species == "Mitrula_elegans" ~ "Litter Saprotroph",
    species == "Murilentithecium_clematidis" ~ "Undefined Saprotroph",
    species == "Mycenella_trachyspora" ~ "Undefined Saprotroph",
    species == "Naganishia_albida" ~ "Animal Pathogen",
    species == "Naohidea_sebacea" ~ "Fungal Parasite",
    species == "Neoascochyta_europaea" ~ "Plant Pathogen", #leaf scorch on wheat
    species == "Neoascochyta_graminicola" ~ "Plant Pathogen", #leaf scorch on wheat
    species == "Neopestalotiopsis_asiatica" ~ "Plant Pathogen", #leaf spot on almond
    species == "Neopestalotiopsis_musae" ~ "Plant Pathogen", #Leaf bright disease on banana
    species == "Oehlia_diaphana" ~ "Arbuscular Mycorrhizal",
    species == "Pachyphlodes_citrina" ~ "Ectomycorrhizal",
    species == "Paralepista_flaccida" ~ "Litter Saprotroph",
    species == "Paraboeremia_selaginellae" ~ "Plant Pathogen",
    species == "Pseudogymnoascus_pannorum" ~ "Animal pathogen",
    species == "Pseudopithomyces_rosae" ~ "Plant pathogen",
    species == "Pseudoteratosphaeria_flexuosa" ~ "Plant pathogen", #leaf spot on eucalypts
    species == "Pseudoteratosphaeria_ohnowa" ~ "Plant pathogen", #leaf spot on eucalypts
    species == "Rubroboletus_pulcherrimus" ~ "Ectomycorrhizal",
    species == "Sarcogyne_regularis" ~ "Lichenized",
    species == "Sarcoscypha_coccinea" ~ "Litter Saprotroph",
    species == "Septoriella hirta" ~ "Plant pathogen",
    species == "Setophoma_terrestris" ~ "Plant Pathogen", #Corky and Pink Root of Tomato
    species == "Siphula_ceratites" ~ "Lichenized",
    species == "Sistotrema_muscicola" ~ "Ectomycorrhizal",
    species == "Spiroplana_centripeta" ~ "Plant Parasite",
    species == "Stereocaulon_glabrum" ~ "Lichenized",
    species == "Suillellus_queletii" ~ "Ectomycorrhizal",
    species == "Teratosphaeria_viscida" ~ "Plant Pathogen", #leaf disease on eucalypts
    species == "Teratosphaericola_pseudoafricana" ~ "Plant Pathogen", #leaf disease on eucalypts
    species == "Thamnidium_elegans" ~ "Dung Saprotroph",
    species == "Thelephora_palmata" ~ "Ectomycorrhizal",
    species == "Tuber_aestivum" ~ "Ectomycorrhizal",
    species == "Tuber_maculatum" ~ "Ectomycorrhizal",
    species == "Tuber_brumale" ~ "Ectomycorrhizal",
    species == "Tuber_excavatum" ~ "Ectomycorrhizal",
    species == "Tuber_gennadii" ~ "Ectomycorrhizal",
    species == "Tuber_huidongense" ~ "Ectomycorrhizal",
    species == "Tuber_melanosporum" ~ "Ectomycorrhizal",
    species == "Tuber_rufum" ~ "Ectomycorrhizal",
    species == "Tuber_vesicoperidium" ~ "Ectomycorrhizal",
    species == "unidentified_Ambisporaceae" ~ "Arbuscular Mycorrhizal",
    species == "unidentified_Claroideoglomeraceae" ~ "Arbuscular Mycorrhizal",
    species == "unidentified_Glomeromycetes" ~ "Arbuscular Mycorrhizal",
    species == "unidentified_GS24" ~ "Arbuscular Mycorrhizal",
    species == "unidentified_Tuber" ~ "Ectomycorrhizal",
    species == "Xerocomellus_chrysenteron" ~ "Ectomycorrhizal",
    species == "Xerocomellus_porosporus" ~ "Ectomycorrhizal",
    species == "Xerocomellus_redeuilhii" ~ "Ectomycorrhizal",
    species == "Zygophlyctis_melosirae" ~ "Algal Parasite",
    species == "Zymoseptoria_brevis" ~ "Algal Parasite",
    TRUE ~ guild.y)) %>% 
  select(-trophicMode.x, -trophicMode.y, -guild.x, -guild.y)

#{{check for OTU duplicates}}
duplicated_otu <- funcoffee.guilds %>% 
  mutate(dup_otu = duplicated(OTU)) %>% 
  filter(dup_otu)

#{{remove OTU duplicates and homogenize guild and trophic terms}}
funcoffee.guilds <- funcoffee.guilds %>% 
  distinct(OTU, .keep_all = TRUE) %>%
  mutate(trophic_mode = case_when(
    guild_mode == "Endophyte-Plant Pathogen" ~ "Pathotroph",
    guild_mode == "Epiphyte-Plant Pathogen" ~ "Pathotroph",
    guild_mode == "Ericoid Mycorrhizal" ~ "Symbiotroph",
    trophic_mode == "Pathotroph-Saprotroph-Symbiotroph" ~ "unidentified_trophic_mode",
    trophic_mode == "Pathotroph-Saprotroph" & str_detect(species, "unidentified") ~ "uncharacterized_Pathotroph-Saprotroph",
    trophic_mode == "Pathotroph-Symbiotroph" & str_detect(species, "unidentified") ~ "uncharacterized_Pathotroph-Saprotroph",
    trophic_mode == "Saprotroph-Symbiotroph" & str_detect(species, "unidentified") ~ "uncharacterized_Pathotroph-Saprotroph",
    TRUE ~ trophic_mode))

tmp <- filter(funcoffee.guilds, is.na(trophic_mode))
rm(funguild)

## ------------ 11.2.2. FunFun ----

#{{import FunFun database}}
funfun <- fungal_traits()

## ------------ 11.2.3. FungalTraits (Funguild+FunFun+updates) ----

#{{import FunFun database}}
funtraits <- read_tsv("FungalTraits 1.2_ver_16Dec_2020_genera.txt")

#{{format funguild database}}
funtraits <- funtraits %>%
  #{{create genus column}}
  rename_with(tolower)

#{{Merge OTU taxonomic data and FungalTraits database}}
funcoffee.traits <- ITS2_tax %>% 
  left_join(select(funtraits,-phylum, -class, -order, -family), by = "genus")

rm(funtraits)

#{{check for OTU duplicates}}
duplicated_otu <- funcoffee.traits %>% 
  mutate(dup_otu = duplicated(OTU)) %>% 
  filter(dup_otu)

#{{remove OTU duplicates}}
funcoffee.traits <- funcoffee.traits %>% 
  distinct(OTU, .keep_all = TRUE)

#{{some species unassigned or poorly-assigned were re-assigned}}
funcoffee.traits <- funcoffee.traits %>%
  mutate(primary_lifestyle = if_else(str_detect(phylum, "Glomeromycota"), 
                                     str_c("arbuscular_mycorrhizal"), primary_lifestyle))

## ---- 12. TAXONOMIC COUNT ----
## ---- 13. STATICTICS ON ENVIRONMENTAL DATA ----
## ---- 14. A TAXONOMIC OR FUNCTIONAL DISTRIBUTION ----
## -------- 14.1. bacteria ----
## ------------ 14.1.1. taxonomy ----
## ---------------- 14.1.1.1. plant organ & disease ----
## -------------------- 14.2.2.1.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}
tax <- "genus"
tax_target <- "Xanthomonas"
organ <- "root"
fact_o <- "type"
fact_h <- "health_status"

tmp <- bact_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  #{{select the taxonomic rank to analyse}}
  left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
  select(-OTU) %>%
  #{{select the pathogen to analyse}}
  filter(get(tax) == tax_target) %>%
  #{{group by the taxonomic rank}}
  group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
  summarise_all(sum) %>%
  column_to_rownames(var = tax) %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatments to analyse}}
  filter(type == organ) %>%
  select("Xanthomonas") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Xanthomonas)
hist(tmp$Xanthomonas)
ggqqplot(tmp$Xanthomonas, ylab = "abundance")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(log1p(tmp$Xanthomonas))
hist(log1p(tmp$Xanthomonas))
ggqqplot(log1p(tmp$Xanthomonas), ylab = "abundance")

## ------------------------ 14.2.2.1.1.1. anova and tukey tests {stats}---- 
## ------------------------ 14.2.2.1.1.2. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the leaf data}}
{
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "leaf"
  fact_o <- "type"
  fact_h <- "health_status"
  
  #{{Xanthomonas abundance}}
  tmp <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_xan_l <- kruskal.test(log1p(Xanthomonas)~health_status,
                             data = tmp, na.action = na.omit)
} 

#{{select the stalk data}}
{
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "stalk"
  fact_o <- "type"
  fact_h <- "health_status"
  
  #{{Xanthomonas abundance}}
  tmp <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_xan_sk <- kruskal.test(log1p(Xanthomonas)~health_status,
                             data = tmp, na.action = na.omit)
} 

#{{select the root data}}
{
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "root"
  fact_o <- "type"
  fact_h <- "health_status"
  
  #{{Xanthomonas abundance}}
  tmp <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_xan_r <- kruskal.test(log1p(Xanthomonas)~health_status,
                             data = tmp, na.action = na.omit)
} 

## -------------------- 14.1.1.1.1. barchart ----

#{{select file name and corresponding conditions}}
output <- "Barchart_bact_phylum_genus_organ_health.pdf"

#{{jco palette visualisation}}
show_col(pal_simpsons("springfield")(16))

#{{select and re-order the factors}}
flevels_organ <- c("leaf", "stalk", "root")
flevels_health <- c("healthy","disease")
flevels_organ_health <- c("leaf_healthy", "leaf_disease", 
                          "stalk_healthy", "stalk_disease",
                          "root_healthy", "root_disease")

#{{bar1 for phylum}}
{
  tax <- "phylum"
  tax_target <- "phylum"
  fact_o <- "type"
  fact_h <- "health_status"
  
  bar1 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_target) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.002, "< 0.2%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Phylum") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 2)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.2%",
                                   expression(italic("Acidobacteriota")),
                                   expression(italic("Actinobacteriota")),
                                   expression(italic("Bacteroidota")),
                                   expression(italic("Campylobacterota")),
                                   expression(italic("Chloroflexi")),
                                   expression(italic("Cyanobacteria")),
                                   expression(italic("Deinococcota")),
                                   expression(italic("Desulfobacterota")),
                                   expression(italic("Firmicutes")),
                                   expression(italic("Myxococcota")),
                                   expression(italic("Planctomycetota")),
                                   expression(italic("Proteobacteria")),
                                   expression(italic("Verrucomicrobiota")))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar2 for Acidobacteriota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Acidobacteriota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar2 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.005, "< 0.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Acidobacteriota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.5%",
                                   expression(italic("Aridibacter")),
                                   expression(italic("Blastocatella")),
                                   expression(italic("Bryobacter")),
                                   "RB41",
                                   expression(paste("un. ", italic("Blastocatellaceae"))),
                                   expression(paste("un. ", italic("Holophagaceae"))),
                                   expression(paste("un. ", italic("Koribacteraceae"))),
                                   expression(paste("un. ", italic("Solibacteraceae"))),
                                   expression(paste("un. ", italic("Thermoanaerobaculaceae"))),
                                   expression(paste("un. ", italic("Vicinamibacteraceae"))),
                                   expression(italic("Vicinamibacter")))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar3 for Actinobacteriota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Actinobacteriota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar3 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.02, "< 2%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Actinobacteriota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 2%",
                                   expression(italic("Acidothermus")),
                                   expression(italic("Actinoplanes")),
                                   expression(italic("Agromyces")),
                                   expression(italic("Curtobacterium")),
                                   expression(italic("Kineococcus")),
                                   expression(italic("Kineosporia")),
                                   expression(italic("Klenkia")),
                                   expression(italic("Quadrisphaera")),
                                   expression(italic("Saccharopolyspora")),
                                   expression(italic("Streptomyces")),
                                   expression(paste("un. ", italic("Kineosporiaceae"))),
                                   expression(paste("un. ", italic("Microbacteriaceae"))),
                                   expression(paste("un. ", italic("Micromonosporaceae"))),
                                   expression(paste("un. ", italic("Streptomycetaceae"))),
                                   expression(paste("un. ", italic("Thermomonosporaceae"))))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar4 for Bacteroidota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Bacteroidota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar4 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.02, "< 2%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Bacteroidota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 2%",
                                   expression(italic("Chitinophaga")),
                                   expression(italic("Flavobacterium")),
                                   expression(italic("Hymenobacter")),
                                   expression(italic("Niastella")),
                                   expression(italic("Ohtaekwangia")),
                                   expression(italic("Terrimonas")),
                                   expression(paste("un. ", italic("Chinitophagaceae"))),
                                   expression(paste("un. ", italic("Microscillaceae"))),
                                   expression(paste("un. ", italic("Paludibacteraceae"))))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar5 for Campylobacterota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Campylobacterota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar5 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Campylobactorota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   expression(italic("Campylobacter")))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar6 for Chloroflexi}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Chloroflexi"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar6 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.005, "< 0.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Chloroflexi") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.5%",
                                   expression(italic("Chloronema")),
                                   "un. FFCH7168",
                                   expression(italic("Oscillochloris")),
                                   expression(italic("Roseiflexus")),
                                   expression(italic("Thermosporothrix")),
                                   "un. A4b",
                                   expression(paste("un. ", italic("Anaerolineaceae"))),
                                   expression(paste("un. ", italic("Caldilineaceae"))),
                                   expression(paste("un. ", italic("Chloroflexaceae"))),
                                   "un. JG30-KF-AS9",
                                   "un. JG30-KF-CM45",
                                   expression(paste("un. ", italic("Ktedonobacteraceae"))),
                                   expression(paste("un. ", italic("Roseiflexaceae"))))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar7 for Cyanobacteria}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Cyanobacteria"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar7 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Cyanobacteria") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   "un. CENA518",
                                   "un. CENA533",
                                   "un. EcFYyy-200",
                                   "un. JSC-12",
                                   expression(italic("Pantalinema")),
                                   expression(italic("Tolypothrix")),
                                   expression(paste("un. ", italic("Coleofasciculaceae"))),
                                   expression(paste("un. ", italic("Leptolyngbyaceae"))),
                                   expression(paste("un. ", italic("Nostocaceae"))),
                                   expression(paste("un. ", italic("Obscuribacteraceae"))),
                                   expression(paste("un. ", italic("Phormidiaceae"))))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar8 for Deinococcota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Deinococcota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  bar8 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Deinococcota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   expression(italic("Deinococcus")))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar9 for Desulfobacterota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Desulfobacterota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar9 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Desulfobacterota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   expression(italic("Citrifermentans")),
                                   expression(italic("Desulfobulbus")),
                                   expression(italic("Desulfopila")),
                                   expression(italic("DeSulvibrio")),
                                   expression(italic("Geoalkalibacter")),
                                   expression(italic("Geobacter")),
                                   expression(paste("un. ", italic("Desulfocapsaceae"))),
                                   expression(paste("un. ", italic("Desulfovibrionaceae"))),
                                   expression(paste("un. ", italic("Geobacteraceae"))))) +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar10 for Firmicutes}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Firmicutes"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar10 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.03, "< 3%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Firmicutes") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 3%",
                                   expression(italic("Bacillus")),
                                   expression(italic("Exiguobacterium")),
                                   expression(italic("Lactococcus")),
                                   expression(italic("Paenibacillus")),
                                   expression(italic("Pseudogracilibacillus")),
                                   expression(italic("Streptococcus")),
                                   expression(italic("Terribacillus")),
                                   expression(italic("Tumebacillus")),
                                   expression(paste("un. ", italic("Bacillaceae"))),
                                   expression(paste("un. ", italic("Clostridiaceae"))),
                                   expression(paste("un. ", italic("Enterococaceae"))),
                                   expression(paste("un. ", italic("Streptococcaceae"))))) +
    
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar11 for Myxococcota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Myxococcota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar11 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Myxococcota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   expression(italic("Anaeromyxobacter")),
                                   expression(italic("Haliangium")),
                                   "un. P3OB-42",
                                   expression(italic("Pajaroellobacter")),
                                   expression(italic("Phaselicystis")),
                                   "un. Blrii41",
                                   expression(paste("un. ", italic("Polyangiaceae"))),
                                   expression(paste("un. ", italic("Sandaracinaceae"))))) +
    
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar12 for Planctomycetota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Planctomycetota"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar12 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Planctomycetota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1%",
                                   "un. AKYG587",
                                   expression(italic("Aquisphaera")),
                                   "CL500-3",
                                   expression(italic("Fimbriiglobus")),
                                   expression(italic("Gemmata")),
                                   expression(italic("Pirellula")),
                                   "un. SH-PL14",
                                   "un. SM1A02",
                                   expression(paste("un. ", italic("Gemmataceae"))),
                                   expression(paste("un. ", italic("Isosphaeraceae"))),
                                   expression(paste("un. ", italic("Pirelluceae"))),
                                   "un. SG8-4")) +
    
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar13 for Proteobacteria}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Proteobacteria"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  bar13 <- bact_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(bact_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_o), all_of(fact_h)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_o), -all_of(fact_h)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "organ_health", type, health_status, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by(organ_health) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.02, "< 2%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = organ_health, y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Proteobacteria") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_x_discrete(limits = flevels_organ_health,
                     labels = c("Leaf\nhealthy", "Leaf\ndisease",
                                "Stalk\nhealthy", "Stalk\ndisease",
                                "Root\nhealthy", "Root\ndisease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 2%",
                                   expression(italic("Burkholderia and rel.")),
                                   expression(italic("Methylobacterium and rel.")),
                                   expression(italic("Pseudomonas")),
                                   expression(italic("Sphingomonas")),
                                   expression(paste("un. ", italic("Beggiatoaceae"))),
                                   expression(paste("un. ", italic("Comamonadaceae"))),
                                   expression(paste("un. ", italic("Xanthobacteraceae"))),
                                   expression(italic("Xanthomonas")))) +
    
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{assembly of barchart plots}}
{
  (bar1 /
  (bar2 | bar3 | bar4) / (bar5 | bar6 | bar7) /
  (bar8 | bar9 | bar10) / (bar11 | bar12 | bar13))
  
  width <- 25
  height <- 25
  
  ggsave(file = output, width = width , height = height)
}

## -------------------- 14.1.1.1.2. boxplot ----

#{{select file name}}
output <- "Boxplot_Xanthomonas_pathogens_health.pdf"

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
colpalette_health <- c("#46732EFF", "#C80813FF")

#{{select and re-order the factors}}
flevels_organ <- c("leaf", "stalk", "root")
flevels_health <- c("healthy","disease")
flevels_organ_health <- c("leaf_healthy", "leaf_disease", 
                          "stalk_healthy", "stalk_disease",
                          "root_healthy", "root_disease")

#{{boxplot for pathogen in leaf}}
{
  #select data bact_OTUs_l_rarefied_t #or bact_OTUs_rarefied_t
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "leaf"
  fact_o <- "type"
  fact_h <- "health_status"
  
  box_xan_l <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = log1p(get(tax_target)))) + #get() allow to pass a variable as column name
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ggtitle( "Leaf") +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    scale_color_manual(limits = flevels_health,
                      values = colpalette_health) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    xlab("") + 
    ylab("Relative adundance (log[x+1])") +
    annotate("text", x = 1.5, y = 15, size = 5, hjust = 0.5,
             label = paste("P = ", round(krus_xan_l$p.value, digits = 5), sep = ""),
             fontface = "italic", size = 4)
  
}

#{{boxplot for pathogens in stalk}}
{
  #select data bact_OTUs_sk_rarefied_t #or bact_OTUs_rarefied_t
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "stalk"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  box_xan_sk <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = log1p(get(tax_target)))) + #get() allow to pass a variable as column name
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ggtitle( "Stalk") +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    scale_color_manual(limits = flevels_health,
                      values = colpalette_health) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    xlab("") + 
    ylab("Relative adundance (log[x+1])") +
    annotate("text", x = 1.5, y = 15, size = 5, hjust = 0.5,
             label = paste("P = ", round(krus_xan_sk$p.value, digits = 5), sep = ""),
             fontface = "italic", size = 4)
}

#{{boxplot for pathogens in roots}
{
  #select data bact_OTUs_r_rarefied_t #or bact_OTUs_rarefied_t
  data <- bact_OTUs_rarefied_t
  tax <- "genus"
  tax_target <- "Xanthomonas"
  organ <- "root"
  fact_o <- "type"
  fact_h <- "health_status"
  
  
  box_xan_r <- data %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{select the pathogen to analyse}}
    filter(get(tax) == tax_target) %>%
    #{{group by the taxonomic rank}}
    group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = log1p(get(tax_target)))) + #get() allow to pass a variable as column name
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ggtitle( "Root") +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    scale_color_manual(limits = flevels_health,
                      values = colpalette_health) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    xlab("") + 
    ylab("Relative adundance (log[x+1])") +
    annotate("text", x = 1.5, y = 15, size = 5, hjust = 0.5,
             label = paste("P = ", round(krus_xan_r$p.value, digits = 5), sep = ""),
             fontface = "italic", size = 4)
}

#{{assembly of boxplots}}

(box_xan_l | box_xan_sk | box_xan_r) + 
  #plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") +
  #plot_annotation(title = "Xanthomonas abundance in sugarcane organs") & 
  theme(plot.title = element_text(size = 16, hjust = 0.5))

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

## -------- 14.2. fungi ----
## ------------ 14.2.1. taxonomy ----
## ---------------- 14.2.1.1. locality & diversification ----
## -------------------- 14.2.1.1.1. barchart ----

#{{select file name and corresponding conditions}}
output <- "Barchart_fung_phylum_genus_cluster.pdf"


#{{jco palette visualisation}}
show_col(pal_simpsons("springfield")(16))

#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA", "RCF")
flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 

#{{bar1 for phylum}}
{
  tax <- "phylum"
  tax_target <- "phylum"
  fact_d <- "diversification"
  fact_c <- "cluster_code"

bar1 <- ITS2_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  #{{select the taxonomic rank to analyse}}
  left_join(select(ITS2_tax, OTU, all_of(tax)), by = "OTU") %>%
  select(-OTU) %>%
  #{{group by taxonomic rank}}
  group_by_at(tax_target) %>% #group_by_at replace group_by to pass a variable as column name  
  summarise_all(sum) %>%
  column_to_rownames(var = tax) %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  #{{select condition to plot}}
  rownames_to_column(var = "soil_code") %>%
  left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
            by = "soil_code") %>%
  #{{exclude samples}}
  ####filter(plant_type != "root") %>%
  select(-soil_code) %>%
  #{{transform data for plotting}}
  gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
  #{{merge two factors to create a new condition}}
  unite(col = "cluster_diversification", cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{determine relative abundance of each taxonomic rank per conditionl}}
  group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
  mutate(Rel_abund = reads/sum(reads)) %>%
  ungroup() %>%
  #{{all taxa less than x% in abundance for each plant species are renamed}}
  filter(Rel_abund != "NaN") %>%
  mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
  #{{barchart plot}}
  ggplot(aes(x = get(fact_c), y = reads)) +
  #{{proportional stacked visualisation}}
  geom_col(position = "fill", ####color = "black",
           aes(fill = taxonomy)) +
  ggtitle("Phylum") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0, size = 10),
        axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
        axis.title = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3)) +
  scale_x_discrete(limits = flevels_clu,
                   labels = c("N\n(dry)", "M\n(moderatly dry)",
                              "B\n(variable)", "R\n(wet)",
                              "F\n(very wet)")) +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #{{palette automatic}}
  ####scale_fill_nejm() +
  ####scale_fill_npg() +
  ####scale_fill_startrek() +
  ####scale_fill_viridis_d() +
  ####scale_fill_brewer() +
  scale_fill_simpsons(labels = c("< 0.1%",
                                 expression(italic("Ascomycota")),
                                 expression(italic("Basidiomycota")),
                                 expression(italic("Chytridiomycota")),
                                 expression(italic("Mortierellomycota")),
                                 expression(italic("Mucoromycota")),
                                 "un. fungi")) +
  ylab("proportion (%)") +
  xlab("")
}

#{{bar2 for Ascomycota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Ascomycota"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
bar2 <- ITS2_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
  left_join(select(ITS2_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
  filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
  select(-OTU, -all_of(tax)) %>%
  #{{group by taxonomic rank}}
  group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
  summarise_all(sum) %>%
  column_to_rownames(var = tax_lower) %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  #{{select condition to plot}}
  rownames_to_column(var = "soil_code") %>%
  left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
            by = "soil_code") %>%
  #{{exclude samples}}
  ####filter(plant_type != "root") %>%
  select(-soil_code) %>%
  #{{transform data for plotting}}
  gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
  #{{merge two factors to create a new condition}}
  unite(col = "cluster_diversification", cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{determine relative abundance of each taxonomic rank per conditionl}}
  group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
  mutate(Rel_abund = reads/sum(reads)) %>%
  ungroup() %>%
  #{{all taxa less than x% in abundance for each plant species are renamed}}
  filter(Rel_abund != "NaN") %>%
  mutate(taxonomy = if_else(Rel_abund < 0.01, "< 1%", taxonomy)) %>%
  #{{barchart plot}}
  ggplot(aes(x = get(fact_c), y = reads)) +
  #{{proportional stacked visualisation}}
  geom_col(position = "fill", ####color = "black",
           aes(fill = taxonomy)) +
  ggtitle("Ascomycota") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0, size = 10),
        axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
        axis.title = element_text(size = 10)) +
  guides(fill = guide_legend(nrow = 3)) +
  scale_x_discrete(limits = flevels_clu,
                   labels = c("N\n(dry)", "M\n(moderatly dry)",
                              "B\n(variable)", "R\n(wet)",
                              "F\n(very wet)")) +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #{{palette automatic}}
  ####scale_fill_nejm() +
  ####scale_fill_npg() +
  ####scale_fill_startrek() +
  ####scale_fill_viridis_d() +
  ####scale_fill_brewer() +
  scale_fill_simpsons(labels = c("< 1%",
                                 expression(italic("Candida")),
                                 expression(italic("Chordomyces")),
                                 expression(italic("Dimorphiseta")),
                                 expression(italic("Keithomyces")),
                                 expression(italic("Marquandomyces")),
                                 expression(italic("Paraconiothyrium")),
                                 expression(italic("Paramacroventuria")),
                                 expression(italic("Plectosphaerella")),
                                 expression(italic("Stictis")),
                                 expression(italic("Tolypocladium")),
                                 expression(italic("Trichoderma")),
                                 expression(paste("un. ", italic("Ascomycota"))),
                                 expression(paste("un. ", italic("Sordariomycetes"))))) +
  ylab("Proportion (%)") +
  xlab("")
}

#{{bar3 for Basidiomycota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Basidiomycota"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar3 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(ITS2_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.015, "< 1.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Basidiomycota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 3)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1.5%",
                                   expression(italic("Apiotrichum")),
                                   expression(italic("Clitopilus")),
                                   expression(italic("Geastrum")),
                                   expression(italic("Homophron")),
                                   expression(italic("Lindtneria")),
                                   expression(italic("Pterula")),
                                   expression(italic("Saitozyma")),
                                   expression(italic("Termitomyces")),
                                   expression(italic("Trechispora")),
                                   expression(italic("Tremella")),
                                   expression(paste("un. ", italic("Agaricales"))),
                                   expression(paste("un. ", italic("Agaromycetes"))),
                                   expression(paste("un. ", italic("Tremellales"))))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar4 for Chytridiomycota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Chytridiomycota"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar4 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(ITS2_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.015, "< 1.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Chytridiomycota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 3)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1.5%",
                                   expression(italic("Boothiomyces")),
                                   expression(italic("Kochiomyces")),
                                   expression(italic("Rhizophlyctis")),
                                   expression(italic("Sonoraphlyctis")),
                                   expression(italic("Spizellomyces")),
                                   expression(paste("un. ", italic("Chytridiales"))),
                                   expression(paste("un. ", italic("Chytridiomycota"))),
                                   expression(paste("un. ", italic("GS14"))),
                                   expression(paste("un. ", italic("Rhizophlyctidaceae"))),
                                   expression(paste("un. ", italic("Rhizophydiaceae"))),
                                   expression(paste("un. ", italic("Rhizophydiales"))),
                                   expression(paste("un. ", italic("Spizellomycetaceae"))),
                                   expression(paste("un. ", italic("Spizellomycetales"))))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar5 for Mortierellomycota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Mortierellomycota"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar5 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(ITS2_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.015, "< 1.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Mortierellomycota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 3)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1.5%",
                                   expression(italic("Mortierella")),
                                   expression(italic("Podila")),
                                   expression(paste("un. ", italic("Mortierellaceae"))))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar6 for Mucoromycota}}
{
  tax <- "phylum"
  tax_lower <- "genus"
  tax_target <- "Mucoromycota"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar6 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(ITS2_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.015, "< 1.5%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Mucoromycota") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 3)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 1.5%",
                                   expression(italic("Absidia")),
                                   expression(italic("Cunninghamella")),
                                   expression(italic("Gongronella")),
                                   expression(italic("Umbelopsis")),
                                   expression(paste("un. ", italic("Umbelopsidomycetes"))))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{assembly of barchart plots}}
{
  ((bar1) | (bar2 / bar3 / bar4 / bar5 / bar6)) +
    plot_layout(widths = c(1,2)) +
    plot_annotation(tag_levels = 'A') #+
    #plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  width <- 15
  height <- 20
  
  ggsave(file = output, width = width , height = height)
}

## ------------ 14.2.2. guild ----
## ---------------- 14.2.2.1. locality & diversification ----
## -------------------- 14.2.2.1.1. barchart ----

#{{select file name}}
output <- "Barchart_fung_guild_organ_health.pdf"

#{{select and reorder the factor}}
flevels_health <- c("healthy","disease")

#{{bar1 in leaf for health}}
{
  organ <- "leaf"
  fact_h <- "health_status"
  
  bar_guild1 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    #{{select the functional rank to analyse}}
    left_join(select(funsugarcane.guilds, OTU, guild_mode), by = "OTU") %>%
    #{{delete OTU column)
    select(-OTU) %>% 
    #{{group by the trophic rank rank}}
    group_by(guild_mode) %>%
    summarise_all(sum) %>%
    filter(guild_mode != "NA") %>% 
    column_to_rownames(var = "guild_mode") %>%
    #{{transpose data}}
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    #{{#select only trophic data}}
    select(fact_h, where(is.numeric))  %>%
    gather(key = "guild_mode", value = "reads", -fact_h) %>%
    #{{delete taxonomic group with no read}}
    group_by(guild_mode) %>% 
    mutate(sum_group = sum(reads)) %>%
    ####filter(sum_group > 0) %>%
    select(-sum_group) %>% 
    ungroup() %>%
    #{{select specific guild assignment}}
    filter(guild_mode == "Animal Pathogen" |
             guild_mode == "Arbuscular Mycorrhizal" | 
             guild_mode == "Dung Saprotroph" | 
             guild_mode == "Ectomycorrhizal" | 
             guild_mode == "Endophyte" | 
             guild_mode == "Epiphyte" |
             guild_mode == "Ericoid Mycorrhizal" | 
             guild_mode == "Fungal Parasite" |
             guild_mode == "Fungal Parasite-Litter Saprotroph" |
             guild_mode == "Lichen Parasite" |
             guild_mode == "Lichenized" |
             guild_mode == "Insect Pathogen" |
             guild_mode == "Leaf Saprotroph" |
             guild_mode == "Litter Saprotroph" |
             guild_mode == "Orchid Mycorrhizal" |
             guild_mode == "Plant Pathogen" |
             guild_mode == "Endophyte-Plant Pathogen" |
             guild_mode == "Soil Saprotroph" |
             guild_mode == "Wood Saprotroph" |
             guild_mode == "Undefined Saprotroph") %>%
    #{{determine relative abundance of each taxonomic rank per treatment or total}}
    ####group_by(guild_mode) %>%
    ####mutate(Rel_abund = reads/sum(reads)) %>%
    ####ungroup() %>%
    #{{all taxa less than x% in abundance for each land use are renamed}}
    ####mutate(guild_mode = if_else(Rel_abund < 0.1, "< 10%", guild_mode)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = guild_mode)) +
    ggtitle("Leaf") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 12),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 16,
                                     colour = "black"),
          axis.title = element_text(size = 16)) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy","disease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    scale_fill_igv() +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar2 in stalk for health}}
{
  organ <- "stalk"
  fact_h <- "health_status"
  
  bar_guild2 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    #{{select the functional rank to analyse}}
    left_join(select(funsugarcane.guilds, OTU, guild_mode), by = "OTU") %>%
    #{{delete OTU column)
    select(-OTU) %>% 
    #{{group by the trophic rank rank}}
    group_by(guild_mode) %>%
    summarise_all(sum) %>%
    filter(guild_mode != "NA") %>% 
    column_to_rownames(var = "guild_mode") %>%
    #{{transpose data}}
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    #{{#select only trophic data}}
    select(fact_h, where(is.numeric))  %>%
    gather(key = "guild_mode", value = "reads", -fact_h) %>%
    #{{delete taxonomic group with no read}}
    group_by(guild_mode) %>% 
    mutate(sum_group = sum(reads)) %>%
    ####filter(sum_group > 0) %>%
    select(-sum_group) %>% 
    ungroup() %>%
    #{{select specific guild assignment}}
    filter(guild_mode == "Animal Pathogen" |
             guild_mode == "Arbuscular Mycorrhizal" | 
             guild_mode == "Dung Saprotroph" | 
             guild_mode == "Ectomycorrhizal" | 
             guild_mode == "Endophyte" | 
             guild_mode == "Epiphyte" |
             guild_mode == "Ericoid Mycorrhizal" | 
             guild_mode == "Fungal Parasite" |
             guild_mode == "Fungal Parasite-Litter Saprotroph" |
             guild_mode == "Lichen Parasite" |
             guild_mode == "Lichenized" |
             guild_mode == "Insect Pathogen" |
             guild_mode == "Leaf Saprotroph" |
             guild_mode == "Litter Saprotroph" |
             guild_mode == "Orchid Mycorrhizal" |
             guild_mode == "Plant Pathogen" |
             guild_mode == "Endophyte-Plant Pathogen" |
             guild_mode == "Soil Saprotroph" |
             guild_mode == "Wood Saprotroph" |
             guild_mode == "Undefined Saprotroph") %>%
    #{{determine relative abundance of each taxonomic rank per treatment or total}}
    ####group_by(guild_mode) %>%
    ####mutate(Rel_abund = reads/sum(reads)) %>%
    ####ungroup() %>%
    #{{all taxa less than x% in abundance for each land use are renamed}}
    ####mutate(guild_mode = if_else(Rel_abund < 0.1, "< 10%", guild_mode)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = guild_mode)) +
    ggtitle("Stalk") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 12),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 16,
                                     colour = "black"),
          axis.title = element_text(size = 16)) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy","disease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    scale_fill_igv() +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar3 in root for health}}
{
  organ <- "root"
  fact_h <- "health_status"
  
  bar_guild3 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    #{{select the functional rank to analyse}}
    left_join(select(funsugarcane.guilds, OTU, guild_mode), by = "OTU") %>%
    #{{delete OTU column)
    select(-OTU) %>% 
    #{{group by the trophic rank rank}}
    group_by(guild_mode) %>%
    summarise_all(sum) %>%
    filter(guild_mode != "NA") %>% 
    column_to_rownames(var = "guild_mode") %>%
    #{{transpose data}}
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select the treatments to analyse}}
    filter(type == organ) %>%
    #{{#select only trophic data}}
    select(fact_h, where(is.numeric))  %>%
    gather(key = "guild_mode", value = "reads", -fact_h) %>%
    #{{delete taxonomic group with no read}}
    group_by(guild_mode) %>% 
    mutate(sum_group = sum(reads)) %>%
    ####filter(sum_group > 0) %>%
    select(-sum_group) %>% 
    ungroup() %>%
    #{{select specific guild assignment}}
    filter(guild_mode == "Animal Pathogen" |
             guild_mode == "Arbuscular Mycorrhizal" | 
             guild_mode == "Dung Saprotroph" | 
             guild_mode == "Ectomycorrhizal" | 
             guild_mode == "Endophyte" | 
             guild_mode == "Epiphyte" |
             guild_mode == "Ericoid Mycorrhizal" | 
             guild_mode == "Fungal Parasite" |
             guild_mode == "Fungal Parasite-Litter Saprotroph" |
             guild_mode == "Lichen Parasite" |
             guild_mode == "Lichenized" |
             guild_mode == "Insect Pathogen" |
             guild_mode == "Leaf Saprotroph" |
             guild_mode == "Litter Saprotroph" |
             guild_mode == "Orchid Mycorrhizal" |
             guild_mode == "Plant Pathogen" |
             guild_mode == "Endophyte-Plant Pathogen" |
             guild_mode == "Soil Saprotroph" |
             guild_mode == "Wood Saprotroph" |
             guild_mode == "Undefined Saprotroph") %>%
    #{{determine relative abundance of each taxonomic rank per treatment or total}}
    ####group_by(guild_mode) %>%
    ####mutate(Rel_abund = reads/sum(reads)) %>%
    ####ungroup() %>%
    #{{all taxa less than x% in abundance for each land use are renamed}}
    ####mutate(guild_mode = if_else(Rel_abund < 0.1, "< 10%", guild_mode)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = guild_mode)) +
    ggtitle("Root") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 12),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 16,
                                     colour = "black"),
          axis.title = element_text(size = 16)) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy","disease")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    scale_fill_igv() +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{assembly of barcharts}}

(bar_guild1 | bar_guild2 | bar_guild3) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

width <- 15
height <- 9

ggsave(file = output, width = width , height = height)

## ------------ 16.1.4. potential pathogens
## ------------ 14.2.3. traits ----
## ---------------- 14.2.3.1. locality & diversification ----
## -------------------- 14.2.3.1.1. barchart ----

#{{select file name}}
output <- "Barchart_fung_traits_div.pdf"

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))

#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA", "RCF")
flevels_clu <- c("B","F","M","N","R")

#{{select color palette}}
colpalette_div <- c("#D2AF81FF", "#FED439FF", "#91331FFF", "#F05C3BFF", "#D5E4A2FF", "#46732EFF")
colpalette_clu <- c("#71D0F5FF", "#709AE1FF", "#197EC0FF", "#075149FF", "#1A9993FF")

#{{bar1 in soil for diversification}}
{
  
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
    bar_trait1 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    #{{select the functional rank to analyse}}
    left_join(select(funcoffee.traits, OTU, primary_lifestyle), by = "OTU") %>%
    #{{delete OTU column)
    select(-OTU) %>% 
    #{{group by the lifestyle rank}}
    group_by(primary_lifestyle) %>%
    summarise_all(sum) %>%
    filter(primary_lifestyle != "NA") %>% 
    column_to_rownames(var = "primary_lifestyle") %>%
    #{{transpose data}}
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{#select only trophic data}}
    select(fact_d, where(is.numeric))  %>%
    gather(key = "primary_lifestyle", value = "reads", -fact_d) %>%
    #{{delete taxonomic group with no read}}
    group_by(primary_lifestyle) %>% 
    mutate(sum_group = sum(reads)) %>%
    ####filter(sum_group > 0) %>%
    select(-sum_group) %>% 
    ungroup() %>%
    #{{select specific guild assignment}}
    filter(primary_lifestyle == "algal_parasite"  |
             primary_lifestyle == "animal-associated" |
             primary_lifestyle == "animal_endosymbiont" |
             primary_lifestyle == "animal_parasite" |
             primary_lifestyle == "arbuscular_mycorrhizal" |
             primary_lifestyle == "dung_saprotroph" |
             primary_lifestyle == "ectomycorrhizal" |
             primary_lifestyle == "epiphyte" |
             primary_lifestyle == "foliar_endophyte" |
             primary_lifestyle == "lichen_parasite" |
             primary_lifestyle == "lichenized" |
             primary_lifestyle == "litter_saprotroph" |
             primary_lifestyle == "moss_symbiont" |
             primary_lifestyle == "mycoparasite" |
             primary_lifestyle == "nectar/tap_saprotroph" |
             primary_lifestyle == "plant_pathogen" |
             primary_lifestyle == "pollen_saprotroph" |
             primary_lifestyle == "root_endophyte" |
             primary_lifestyle == "soil_saprotroph" |
             primary_lifestyle == "sooty_mold" |
             primary_lifestyle == "unspecified_pathotroph" |
             primary_lifestyle == "unspecified_saprotroph" |
             primary_lifestyle == "unspecified_symbiotroph" |
             primary_lifestyle == "wood_saprotroph") %>%
    #{{determine relative abundance of each taxonomic rank per treatment or total}}
    ####group_by(primary_lifestyle) %>%
    ####mutate(Rel_abund = reads/sum(reads)) %>%
    ####ungroup() %>%
    #{{all taxa less than x% in abundance for each land use are renamed}}
    ####mutate(primary_lifestyle = if_else(Rel_abund < 0.1, "< 10%", primary_lifestyle)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_d), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = primary_lifestyle)) +
    ggtitle("") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 12),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 16,
                                     colour = "black"),
          axis.title = element_text(size = 16)) +
      scale_x_discrete(limits = flevels_div,
                       labels = c("monoculture", "+ banana",
                                  "+ shade tree", "+ banana\n+ shade tree",
                                  "+ banana\n+ shade tree\n+ others", "+ forest")) +
      scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    scale_fill_igv() +
    ylab("Relative abundance (%)") +
    xlab("")
}

#{{bar2 in soil for cluster}}
{
  
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar_trait2 <- ITS2_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>% 
    #{{select the functional rank to analyse}}
    left_join(select(funcoffee.traits, OTU, primary_lifestyle), by = "OTU") %>%
    #{{delete OTU column)
    select(-OTU) %>% 
    #{{group by the lifestyle rank}}
    group_by(primary_lifestyle) %>%
    summarise_all(sum) %>%
    filter(primary_lifestyle != "NA") %>% 
    column_to_rownames(var = "primary_lifestyle") %>%
    #{{transpose data}}
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "soil_code") %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{#select only trophic data}}
    select(fact_c, where(is.numeric))  %>%
    gather(key = "primary_lifestyle", value = "reads", -fact_c) %>%
    #{{delete taxonomic group with no read}}
    group_by(primary_lifestyle) %>% 
    mutate(sum_group = sum(reads)) %>%
    ####filter(sum_group > 0) %>%
    select(-sum_group) %>% 
    ungroup() %>%
    #{{select specific guild assignment}}
    filter(primary_lifestyle == "algal_parasite"  |
             primary_lifestyle == "animal-associated" |
             primary_lifestyle == "animal_endosymbiont" |
             primary_lifestyle == "animal_parasite" |
             primary_lifestyle == "arbuscular_mycorrhizal" |
             primary_lifestyle == "dung_saprotroph" |
             primary_lifestyle == "ectomycorrhizal" |
             primary_lifestyle == "epiphyte" |
             primary_lifestyle == "foliar_endophyte" |
             primary_lifestyle == "lichen_parasite" |
             primary_lifestyle == "lichenized" |
             primary_lifestyle == "litter_saprotroph" |
             primary_lifestyle == "moss_symbiont" |
             primary_lifestyle == "mycoparasite" |
             primary_lifestyle == "nectar/tap_saprotroph" |
             primary_lifestyle == "plant_pathogen" |
             primary_lifestyle == "pollen_saprotroph" |
             primary_lifestyle == "root_endophyte" |
             primary_lifestyle == "soil_saprotroph" |
             primary_lifestyle == "sooty_mold" |
             primary_lifestyle == "unspecified_pathotroph" |
             primary_lifestyle == "unspecified_saprotroph" |
             primary_lifestyle == "unspecified_symbiotroph" |
             primary_lifestyle == "wood_saprotroph") %>%
    #{{determine relative abundance of each taxonomic rank per treatment or total}}
    ####group_by(primary_lifestyle) %>%
    ####mutate(Rel_abund = reads/sum(reads)) %>%
    ####ungroup() %>%
    #{{all taxa less than x% in abundance for each land use are renamed}}
    ####mutate(primary_lifestyle = if_else(Rel_abund < 0.1, "< 10%", primary_lifestyle)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = primary_lifestyle)) +
    ggtitle("") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 12),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 16,
                                     colour = "black"),
          axis.title = element_text(size = 16)) +
    scale_x_discrete(limits = flevels_clu, 
                       labels = c("B","F","M","N","R")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    scale_fill_igv() +
    ylab("Relative abundance (%)") +
    xlab("")
}


#{{assembly of barcharts}}

(bar_trait1 | bar_trait2) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

width <- 15
height <- 9

ggsave(file = output, width = width , height = height)

## -------- 14.3. amf ----
## ------------ 14.2.1. taxonomy ----
## ---------------- 14.2.1.1. locality & diversification ----
## -------------------- 14.2.1.1.1. barchart ----

#{{select file name and corresponding conditions}}
output <- "Barchart_amf_family_genus_cluster.pdf"

#{{jco palette visualisation}}
show_col(pal_simpsons("springfield")(16))

#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA", "RCF")
flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 

#{{bar1 for family}}
{
  tax <- "family"
  tax_target <- "family"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar1 <- amf_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the taxonomic rank to analyse}}
    left_join(select(amf_tax, OTU, all_of(tax)), by = "OTU") %>%
    select(-OTU) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_target) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Family") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.1%",
                                   expression(italic("Archaeosporales")),
                                   expression(italic("Diversisporales")),
                                   expression(italic("Endogonales")),
                                   expression(italic("Glomerales")),
                                   expression(paste("un. ", italic("Glomeromycotina"))))) +
    ylab("proportion (%)") +
    xlab("")
}

#{{bar2 for Archaeosporales}}
{
  tax <- "family"
  tax_lower <- "genus"
  tax_target <- "Archaeosporales"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar2 <- amf_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(amf_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Archaeosporales") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.1%",
                                   expression(italic("Ambispora")))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar3 for Diversisporales}}
{
  tax <- "family"
  tax_lower <- "genus"
  tax_target <- "Diversisporales"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar3 <- amf_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(amf_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Diversisporales") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.1%",
                                   expression(italic("Acaulospora")),
                                   expression(italic("Diversispora")),
                                   expression(italic("Gigaspora")),
                                   expression(italic("Scutellospora")))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar4 for Endogonales}}
{
  tax <- "family"
  tax_lower <- "genus"
  tax_target <- "Endogonales"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar4 <- amf_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(amf_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Endogonales") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.1%",
                                   expression(italic("Endogone")))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{bar5 for Glomerales}}
{
  tax <- "family"
  tax_lower <- "genus"
  tax_target <- "Glomerales"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  bar5 <- amf_OTUs_rarefied_t %>%
    t() %>%
    as_tibble(rownames = NA) %>% 
    rownames_to_column(var = "OTU") %>%
    #{{select the subtaxonomic rank for which the lower taxonomic rank will be analyse}}
    left_join(select(amf_tax, OTU, all_of(tax), all_of(tax_lower)), by = "OTU") %>%
    filter(!!as.name(tax) == tax_target) %>% ##!! are used to unquote the argument making possible to import a column name from tax
    select(-OTU, -all_of(tax)) %>%
    #{{group by taxonomic rank}}
    group_by_at(tax_lower) %>% #group_by_at replace group_by to pass a variable as column name  
    summarise_all(sum) %>%
    column_to_rownames(var = tax_lower) %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    #{{select condition to plot}}
    rownames_to_column(var = "soil_code") %>%
    left_join(select(metadata, soil_code, all_of(fact_d), all_of(fact_c)), 
              by = "soil_code") %>%
    #{{exclude samples}}
    ####filter(plant_type != "root") %>%
    select(-soil_code) %>%
    #{{transform data for plotting}}
    gather(key = "taxonomy", value = "reads", -all_of(fact_d), -all_of(fact_c)) %>%
    #{{merge two factors to create a new condition}}
    unite(col = "cluster_diversification", cluster_code, diversification, 
          sep = "_", remove = FALSE, na.rm = FALSE) %>%
    #{{determine relative abundance of each taxonomic rank per conditionl}}
    group_by_at(fact_c) %>%  #group_by_at replace group_by to import column name from cond
    mutate(Rel_abund = reads/sum(reads)) %>%
    ungroup() %>%
    #{{all taxa less than x% in abundance for each plant species are renamed}}
    filter(Rel_abund != "NaN") %>%
    mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
    #{{barchart plot}}
    ggplot(aes(x = get(fact_c), y = reads)) +
    #{{proportional stacked visualisation}}
    geom_col(position = "fill", ####color = "black",
             aes(fill = taxonomy)) +
    ggtitle("Glomerales") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(hjust = 0, size = 10),
          axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 10, colour = "black"),
          axis.title = element_text(size = 10)) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
    #{{palette automatic}}
    ####scale_fill_nejm() +
    ####scale_fill_npg() +
    ####scale_fill_startrek() +
    ####scale_fill_viridis_d() +
    ####scale_fill_brewer() +
    scale_fill_simpsons(labels = c("< 0.1%",
                                   expression(italic("Glomus")),
                                   expression(italic("Sclerocystis")),
                                   expression(paste("un. ", italic("Glomerales"))))) +
    ylab("Proportion (%)") +
    xlab("")
}

#{{assembly of barchart plots}}
{
  (bar1 | (bar2 / bar3 / bar4 / bar5)) +
    plot_layout(widths = c(1,1)) +
    plot_annotation(tag_levels = 'A') #+
  #plot_layout(guides = "collect") & theme(legend.position = "right")
  
  width <- 12
  height <- 12
  
  ggsave(file = output, width = width , height = height)
}

## ---- 15. ALPHADIVERSITY + VARIANCE ANALYSIS ----
# p-values and fdr (http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx)
# p-values and fdr (https://stats.stackexchange.com/questions/198073/correction-for-multiple-testing-on-a-modest-number-of-tests-10-20-with-fdr)
# p-values and fdr http://www.biostathandbook.com/multiplecomparisons.html
# p-values and fdr http://rcompanion.org/rcompanion/f_01.html

## -------- 15.1. bacteria --------
## ------------ 15.1.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}
tmp <- Indices_bact %>%
  #{{merge with metadata}}
  left_join(metadata, by = "soil_code") %>%
  #{{select subdata}}
  ####filter(plant_type != "No_plant_cover") %>%
  ####filter(salinity == "0_200") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Hill_Richness)
hist(tmp$Hill_Richness)
ggqqplot(tmp$Shannon, ylab = "Richness")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(log1p(tmp$Hill_Richness))
hist(log1p(tmp$Hill_Richness))
ggqqplot(log1p(tmp$Hill_Richness), ylab = "Number of OTUs")

## ------------ 15.1.2. anova and tukey tests {stats}---- 
## ------------ 15.1.3. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the leaf}}
{
  
  #{{health status}}
  tmp <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    filter(type == "leaf") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_health_leaf_dR <- kruskal.test(Hill_Richness~health_status,
                                      data = tmp, na.action = na.omit)
  
  krus_health_leaf_dREr <- kruskal.test(Hill_Shannon~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_leaf_dREa <- kruskal.test(Hill_Inv_Simpson~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_leaf_eDRr <- kruskal.test(Hill_Shannon_evenness~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_leaf_eDRa <- kruskal.test(Hill_Simpson_evenness~health_status,
                                        data = tmp, na.action = na.omit)
} 

#{{select the stalk}}
{
  
  #{{health status}}
  tmp <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    filter(type == "stalk") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_health_stalk_dR <- kruskal.test(Hill_Richness~health_status,
                                       data = tmp, na.action = na.omit)
  
  krus_health_stalk_dREr <- kruskal.test(Hill_Shannon~health_status,
                                         data = tmp, na.action = na.omit)
  
  krus_health_stalk_dREa <- kruskal.test(Hill_Inv_Simpson~health_status,
                                         data = tmp, na.action = na.omit)
  
  krus_health_stalk_eDRr <- kruskal.test(Hill_Shannon_evenness~health_status,
                                         data = tmp, na.action = na.omit)
  
  krus_health_stalk_eDRa <- kruskal.test(Hill_Simpson_evenness~health_status,
                                         data = tmp, na.action = na.omit)
}  

#{{select the root}}
{
  
  #{{health status}}
  tmp <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    filter(type == "root") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_health_root_dR <- kruskal.test(Hill_Richness~health_status,
                                      data = tmp, na.action = na.omit)
  
  krus_health_root_dREr <- kruskal.test(Hill_Shannon~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_root_dREa <- kruskal.test(Hill_Inv_Simpson~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_root_eDRr <- kruskal.test(Hill_Shannon_evenness~health_status,
                                        data = tmp, na.action = na.omit)
  
  krus_health_root_eDRa <- kruskal.test(Hill_Simpson_evenness~health_status,
                                        data = tmp, na.action = na.omit)
}  

## ------------ 15.1.4. boxplot ------------

#{{select file name and corresponding conditions}}
output <- "Boxplot_bact_diversity_organ_health.pdf"

#{{select and re-order the factors}}
flevels_health <- c("healthy", "disease")

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{select color palette}}
colpalette_health <- c("#46732EFF", "#C80813FF")


#{{format kruskall statistic}}
krus_health_leaf_dR_annotation <- paste("Kruskall-Wallis test",
                                        "\nHealth status P = ", round(krus_health_leaf_dR$p.value, digits = 4),
                                        sep = "")
krus_health_leaf_dREr_annotation <- paste("Kruskall-Wallis test",
                                          "\nHealth status P = ", round(krus_health_leaf_dREr$p.value, digits = 4),
                                          sep = "")
krus_health_leaf_dREa_annotation <- paste("Kruskall-Wallis test",
                                          "\nHealth status P = ", round(krus_health_leaf_dREa$p.value, digits = 4),
                                          sep = "")


krus_health_stalk_dR_annotation <- paste("Kruskall-Wallis test",
                                         "\nHealth status P = ", round(krus_health_stalk_dR$p.value, digits = 4),
                                         sep = "")
krus_health_stalk_dREr_annotation <- paste("Kruskall-Wallis test",
                                           "\nHealth Status P = ", round(krus_health_stalk_dREr$p.value, digits = 4),
                                           sep = "")
krus_health_stalk_dREa_annotation <- paste("Kruskall-Wallis test",
                                           "\nHealth Status P = ", round(krus_health_stalk_dREa$p.value, digits = 4),
                                           sep = "")

krus_health_root_dR_annotation <- paste("Kruskall-Wallis test",
                                        "\nHealth status P = ", round(krus_health_root_dR$p.value, digits = 4),
                                        sep = "")
krus_health_root_dREr_annotation <- paste("Kruskall-Wallis test",
                                          "\nHealth status P = ", round(krus_health_root_dREr$p.value, digits = 4),
                                          sep = "")
krus_health_root_dREa_annotation <- paste("Kruskall-Wallis test",
                                          "\nHealth status P = ", round(krus_health_root_dREa$p.value, digits = 4),
                                          sep = "")

#{{boxplot for diversity in leaf per health status}}
{
  
  organ <- "leaf"
  fact_h <- "health_status"
  
  box_health_leaf_dR <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,400)) +
    annotate("text", x = 1.5, y = 350,
             label = krus_health_leaf_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                      values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Leaf") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nspecies richness")
  
  box_health_leaf_dREr <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,300)) +
    annotate("text", x = 1.5, y = 250,
             label = krus_health_leaf_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                      values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Leaf") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nexponential of Shannon's entropy index")
  
  box_health_leaf_dREa <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,200)) +
    annotate("text", x = 1.5, y = 175,
             label = krus_health_leaf_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Leaf") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\ninverse of Simpson's concentration index")
  
}

#{{boxplot for diversity in stalk per health status}}
{
  
  organ <- "stalk"
  fact_h <- "health_status"
  
  box_health_stalk_dR <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,400)) +
    annotate("text", x = 1.5, y = 350,
             label = krus_health_stalk_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Stalk") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nspecies richness")
  
  box_health_stalk_dREr <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,300)) +
    annotate("text", x = 1.5, y = 250,
             label = krus_health_stalk_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Stalk") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nexponential of Shannon's entropy index")
  
  box_health_stalk_dREa <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,200)) +
    annotate("text", x = 1.5, y = 175,
             label = krus_health_stalk_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Stalk") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\ninverse of Simpson's concentration index")
  
}

#{{boxplot for diversity in root per health status}}
{
  
  organ <- "root"
  fact_h <- "health_status"
  
  box_health_root_dR <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,400)) +
    annotate("text", x = 1.5, y = 350,
             label = krus_health_root_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Root") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nspecies richness")
  
  box_health_root_dREr <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,300)) +
    annotate("text", x = 1.5, y = 250,
             label = krus_health_root_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Root") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nexponential of Shannon's entropy index")
  
  box_health_root_dREa <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    filter(type == organ) %>%
    ggplot(aes(x = get(fact_h), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_h)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_h)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,200)) +
    annotate("text", x = 1.5, y = 175,
             label = krus_health_root_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_health,
                     labels = c("healthy", "disease")) +
    scale_color_manual(limits = flevels_health,
                       values = colpalette_health) +
    labs(fill = "Health status") +
    labs(title = "", subtitle = "Root") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\ninverse of Simpson's concentration index")
  
}

#{{assembly of boxplots}}
{ 
  ((box_health_leaf_dR | box_health_stalk_dR | box_health_root_dR) /
     (box_health_leaf_dREr | box_health_stalk_dREr | box_health_root_dREr) /
     (box_health_leaf_dREa | box_health_stalk_dREa | box_health_root_dREa)) + 
    plot_layout(guides = "collect") & theme(legend.position = "none")
  
  width <- 15
  height <- 10
  
  ggsave(file = output, width = width , height = height)
  
}

## -------- 15.2 fungi --------
## ------------ 15.2.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}
tmp <- Indices_its2 %>%
  #{{merge with metadata}}
  left_join(metadata, by = "soil_code") %>%
  #{{select subdata}}
  ####filter(type == "root") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Hill_Richness)
hist(tmp$Hill_Richness)
ggqqplot(tmp$Shannon, ylab = "Richness")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(log1p(tmp$Hill_Richness))
hist(log1p(tmp$Hill_Richness))
ggqqplot(log1p(tmp$Hill_Richness), ylab = "Number of OTUs")

## ------------ 15.2.2. anova and tukey tests {stats}---- 
#{{Results have to be carefully interpreted, because aov() assume nested anova)}}
## ------------ 15.1.3. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the soil + diversification}}
{
  
  #{{diversification status}}
  tmp <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    ####filter(type == "leaf") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_div_soil_dR <- kruskal.test(Hill_Richness~diversification,
                                      data = tmp, na.action = na.omit)
  
  krus_div_soil_dREr <- kruskal.test(Hill_Shannon~diversification,
                                        data = tmp, na.action = na.omit)
  
  krus_div_soil_dREa <- kruskal.test(Hill_Inv_Simpson~diversification,
                                        data = tmp, na.action = na.omit)
  
  krus_div_soil_eDRr <- kruskal.test(Hill_Shannon_evenness~diversification,
                                        data = tmp, na.action = na.omit)
  
  krus_div_soil_eDRa <- kruskal.test(Hill_Simpson_evenness~diversification,
                                        data = tmp, na.action = na.omit)
} 

#{{select the soil + cluster}}
{
  
  #{{health status}}
  tmp <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    ####filter(type == "soil") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_clu_soil_dR <- kruskal.test(Hill_Richness~cluster_code,
                                   data = tmp, na.action = na.omit)
  
  krus_clu_soil_dREr <- kruskal.test(Hill_Shannon~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_dREa <- kruskal.test(Hill_Inv_Simpson~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_eDRr <- kruskal.test(Hill_Shannon_evenness~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_eDRa <- kruskal.test(Hill_Simpson_evenness~cluster_code,
                                     data = tmp, na.action = na.omit)
}  

## ------------ 15.1.4. boxplot ------------

#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA", "RCF")
flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{select color palette}}
colpalette_div <- c("#D2AF81FF", "#FED439FF", "#91331FFF", "#F05C3BFF", "#D5E4A2FF", "#46732EFF")

#{{viridis palette visualisation}}
show_col(viridis(10))
colpalette_clu <- c("#FDE725FF", "#B4DE2CFF", "#35B779FF", "#26828EFF", "#31688EFF")

#{{format kruskall statistic}}
krus_div_soil_dR_annotation <- paste("Diversification P = ", round(krus_div_soil_dR$p.value, digits = 4),
                                        sep = "")
krus_div_soil_dREr_annotation <- paste("Diversification P = ", round(krus_div_soil_dREr$p.value, digits = 4),
                                          sep = "")
krus_div_soil_dREa_annotation <- paste("Diversification P = ",  round(krus_div_soil_dREa$p.value, digits = 4),
                                          sep = "")


krus_clu_soil_dR_annotation <- paste("Cluster P = ", round(krus_clu_soil_dR$p.value, digits = 5),
                                         sep = "")
krus_clu_soil_dREr_annotation <- paste("Cluster P = ", round(krus_clu_soil_dREr$p.value, digits = 5),
                                           sep = "")
krus_clu_soil_dREa_annotation <- paste("Cluster P = ", round(krus_clu_soil_dREa$p.value, digits = 5),
                                           sep = "")

#{{boxplot for diversity in soil per diversification level}}
{
  #{{select file name and corresponding conditions}}
  output <- "Boxplot_fung_soil_diversification.pdf"
  
  organ <- "soil"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  box_div_soil_dR <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,1500)) +
    annotate("text", x = 3.5, y = 1500,
             label = krus_div_soil_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                         "+ shade tree", "+ banana\n+ shade tree",
                                         "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 0")
  
  box_div_soil_dREr <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,350)) +
    annotate("text", x = 3.5, y = 350,
             label = krus_div_soil_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                         "+ shade tree", "+ banana\n+ shade tree",
                                         "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 1")
  
  box_div_soil_dREa <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,350)) +
    annotate("text", x = 3.5, y = 350,
             label = krus_div_soil_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                         "+ shade tree", "+ banana\n+ shade tree",
                                         "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 2")
  
}

#{{assembly of boxplots}}
{ 
  (box_div_soil_dR / box_div_soil_dREr / box_div_soil_dREa) +
    plot_layout(guides = "collect") & theme(legend.position = "none")
  
  width <- 9
  height <- 15
  
  ggsave(file = output, width = width , height = height)
  
}

#{{boxplot for diversity in soil per cluster}}
{
  #{{select file name and corresponding conditions}}
  output <- "Boxplot_fung_soil_cluster.pdf"
  
  organ <- "soil"
  fact_c <- "diversification"
  fact_c <- "cluster_code"
  
  box_clu_soil_dR <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,1500)) +
    annotate("text", x = 3, y = 1500,
             label = krus_clu_soil_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                      labels = c("N\n(dry)", "M\n(moderatly dry)",
                                 "B\n(variable)", "R\n(wet)",
                                 "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 0")
  
  box_clu_soil_dREr <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,350)) +
    annotate("text", x = 3, y = 350,
             label = krus_clu_soil_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 1")
  
  box_clu_soil_dREa <- Indices_its2 %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,350)) +
    annotate("text", x = 3, y = 350,
             label = krus_clu_soil_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 2")
  
}

#{{assembly of boxplots}}
{ 
  (box_clu_soil_dR / box_clu_soil_dREr / box_clu_soil_dREa) +
    plot_layout(guides = "collect") & theme(legend.position = "none")
  
  width <- 9
  height <- 15
  
  ggsave(file = output, width = width , height = height)
  
}

## -------- 15.2 amf --------
## ------------ 15.2.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}
tmp <- Indices_its2 %>%
  #{{merge with metadata}}
  left_join(metadata, by = "soil_code") %>%
  #{{select subdata}}
  ####filter(type == "root") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Hill_Richness)
hist(tmp$Hill_Richness)
ggqqplot(tmp$Shannon, ylab = "Richness")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(log1p(tmp$Hill_Richness))
hist(log1p(tmp$Hill_Richness))
ggqqplot(log1p(tmp$Hill_Richness), ylab = "Number of OTUs")

## ------------ 15.2.2. anova and tukey tests {stats}---- 
#{{Results have to be carefully interpreted, because aov() assume nested anova)}}
## ------------ 15.1.3. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the soil + diversification}}
{
  
  #{{diversification status}}
  tmp <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    ####filter(type == "leaf") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_div_soil_dR <- kruskal.test(Hill_Richness~diversification,
                                   data = tmp, na.action = na.omit)
  
  krus_div_soil_dREr <- kruskal.test(Hill_Shannon~diversification,
                                     data = tmp, na.action = na.omit)
  
  krus_div_soil_dREa <- kruskal.test(Hill_Inv_Simpson~diversification,
                                     data = tmp, na.action = na.omit)
  
  krus_div_soil_eDRr <- kruskal.test(Hill_Shannon_evenness~diversification,
                                     data = tmp, na.action = na.omit)
  
  krus_div_soil_eDRa <- kruskal.test(Hill_Simpson_evenness~diversification,
                                     data = tmp, na.action = na.omit)
} 

#{{select the soil + cluster}}
{
  
  #{{health status}}
  tmp <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{subdata to select}}
    ####filter(type == "soil") %>%
    droplevels()
  
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_clu_soil_dR <- kruskal.test(Hill_Richness~cluster_code,
                                   data = tmp, na.action = na.omit)
  
  krus_clu_soil_dREr <- kruskal.test(Hill_Shannon~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_dREa <- kruskal.test(Hill_Inv_Simpson~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_eDRr <- kruskal.test(Hill_Shannon_evenness~cluster_code,
                                     data = tmp, na.action = na.omit)
  
  krus_clu_soil_eDRa <- kruskal.test(Hill_Simpson_evenness~cluster_code,
                                     data = tmp, na.action = na.omit)
}  

## ------------ 15.1.4. boxplot ------------

#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA", "RCF")
flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{select color palette}}
colpalette_div <- c("#D2AF81FF", "#FED439FF", "#91331FFF", "#F05C3BFF", "#D5E4A2FF", "#46732EFF")

#{{viridis palette visualisation}}
show_col(viridis(10))
colpalette_clu <- c("#FDE725FF", "#B4DE2CFF", "#35B779FF", "#26828EFF", "#31688EFF")

#{{format kruskall statistic}}
krus_div_soil_dR_annotation <- paste("Diversification P = ", round(krus_div_soil_dR$p.value, digits = 4),
                                     sep = "")
krus_div_soil_dREr_annotation <- paste("Diversification P = ", round(krus_div_soil_dREr$p.value, digits = 4),
                                       sep = "")
krus_div_soil_dREa_annotation <- paste("Diversification P = ",  round(krus_div_soil_dREa$p.value, digits = 4),
                                       sep = "")


krus_clu_soil_dR_annotation <- paste("Cluster P = ", round(krus_clu_soil_dR$p.value, digits = 4),
                                     sep = "")
krus_clu_soil_dREr_annotation <- paste("Cluster P = ", round(krus_clu_soil_dREr$p.value, digits = 4),
                                       sep = "")
krus_clu_soil_dREa_annotation <- paste("Cluster P = ", round(krus_clu_soil_dREa$p.value, digits = 4),
                                       sep = "")

#{{boxplot for diversity in soil per diversification level}}
{
  #{{select file name and corresponding conditions}}
  output <- "Boxplot_amf_soil_diversification.pdf"
  
  organ <- "soil"
  fact_d <- "diversification"
  fact_c <- "cluster_code"
  
  box_div_soil_dR <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,200)) +
    annotate("text", x = 3.5, y = 200,
             label = krus_div_soil_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                "+ shade tree", "+ banana\n+ shade tree",
                                "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 0")
  
  box_div_soil_dREr <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,100)) +
    annotate("text", x = 3.5, y = 100,
             label = krus_div_soil_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                "+ shade tree", "+ banana\n+ shade tree",
                                "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 1")
  
  box_div_soil_dREa <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_d), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_d)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_d)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,100)) +
    annotate("text", x = 3.5, y = 100,
             label = krus_div_soil_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_div,
                     labels = c("monoculture", "+ banana",
                                "+ shade tree", "+ banana\n+ shade tree",
                                "+ banana\n+ shade tree\n+ others", "+ forest")) +
    scale_color_manual(limits = flevels_div,
                       values = colpalette_div) +
    labs(fill = "Diversification") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 2")
  
}

#{{assembly of boxplots}}
{ 
  (box_div_soil_dR / box_div_soil_dREr / box_div_soil_dREa) +
    plot_layout(guides = "collect") & theme(legend.position = "none")
  
  width <- 9
  height <- 15
  
  ggsave(file = output, width = width , height = height)
  
}

#{{boxplot for diversity in soil per cluster}}
{
  #{{select file name and corresponding conditions}}
  output <- "Boxplot_amf_soil_cluster.pdf"
  
  organ <- "soil"
  fact_c <- "diversification"
  fact_c <- "cluster_code"
  
  box_clu_soil_dR <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Richness)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,200)) +
    annotate("text", x = 3, y = 200,
             label = krus_clu_soil_dR_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 0")
  
  box_clu_soil_dREr <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Shannon)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,100)) +
    annotate("text", x = 3, y = 100,
             label = krus_clu_soil_dREr_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 1")
  
  box_clu_soil_dREa <- Indices_amf %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{select subdata}}
    ####filter(type == organ) %>%
    ggplot(aes(x = get(fact_c), y = Hill_Inv_Simpson)) +
    geom_boxplot(aes(color = get(fact_c)), outlier.shape = NA) +
    geom_jitter(aes(color = get(fact_c)),
                position = position_jitter(0.1),
                alpha = 1/2, cex = 3.5) +
    theme_bw(base_size = 12) +
    ylim(c(0,100)) +
    annotate("text", x = 3, y = 100,
             label = krus_clu_soil_dREa_annotation,
             fontface = "italic", size = 4) +
    scale_x_discrete(limits = flevels_clu,
                     labels = c("N\n(dry)", "M\n(moderatly dry)",
                                "B\n(variable)", "R\n(wet)",
                                "F\n(very wet)")) +
    scale_color_manual(limits = flevels_clu,
                       values = colpalette_clu) +
    labs(fill = "Climatic cluster") +
    labs(title = "", subtitle = "") +
    theme(plot.title = element_text(hjust = 0, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "right",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    xlab("") + 
    ylab("Effective number of species\nq = 2")
  
}

#{{assembly of boxplots}}
{ 
  (box_clu_soil_dR / box_clu_soil_dREr / box_clu_soil_dREa) +
    plot_layout(guides = "collect") & theme(legend.position = "none")
  
  width <- 9
  height <- 15
  
  ggsave(file = output, width = width , height = height)
  
}

## ---- 16. BETADIVERSITY + MULTIVARIATE ANALYSIS ----
## -------- 16.1. bacteria --------
## ------------ 16.1.1. format data ------------
## ---------------- 16.1.1.1. OTU level ----------------

type <- "organ_health"

#{{select the data}}
com <- bact_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  ####filter(plant_type != "No_plant_cover") %>%
  ####droplevels() %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  droplevels() %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  #{{delete OTUs with no read}}
  gather(key = "OTU", value = "reads", -soil_code) %>%
  group_by(OTU) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(OTU, reads) %>%
  column_to_rownames(var = "soil_code")

## ---------------- 16.1.1.2. Taxonomic level ----------------

type <- "organ_health"

#{{select the data}}
com <- bact_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  ####filter(type == type) %>%
  droplevels() %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  column_to_rownames(var = "soil_code") %>%
  #{{select the taxonomic range}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  ####mutate(OTU = as.numeric(OTU)) %>%
  left_join(select(bact_tax,OTU, genus), by = "OTU") %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  #{{delete taxonomic group with no read}}
  gather(key = "taxonomy", value = "reads", -genus) %>%
  group_by(taxonomy) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(taxonomy, reads) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA)

## ------------ 16.1.2. data transformation & standardization ------------

#{{no data tranformation}}
transcom <- com 

#{{select abundance data tranformation}}
transcom <- log1p(com)  #{{sqrt() or log1p()}}

#{{select abundance data tranformation and standardization}}
transcom <- decostand(com,"hellinger", MARGIN = 1)  #{{standardization method = "total" or "hellinger" or "normalize" or "chi.square")
transcom <- wisconsin(com)  #{{standardiozation "wisconsin"}}
transcom <- decostand(log1p(com),"chi.square") #{{abundance data tranformation and standardization}}

barplot(table(unlist(transcom)),las = 1)

## ------------ 16.1.3. distance matrices ------------
## ---------------- 16.1.3.1. samples ----------------

distcom <- vegdist(transcom, method = "bray")
distcom <- vegdist(transcom, method = "morisita")
distcom <- vegdist(transcom, method = "horn")
distcom <- vegdist(transcom, binary = TRUE, method = "raup")
distcom <- vegdist(transcom, method = "jaccard")
#{{adapted with standardized data}
distcom <- vegdist(transcom, method = "euclidean") 

## ---------------- 16.1.3.2. OTU or Taxa ----------------

distcom_t <- vegdist(t(transcom), method = "bray")
distcom_t <- vegdist(t(transcom), method = "morisita")
distcom_t <- vegdist(t(transcom), binary = TRUE, method = "raup")
distcom_t <- vegdist(t(transcom), method = "jaccard")
#{{adapted with standardized data}
distcom_t <- vegdist(t(transcom), method = "euclidean")

## ------------ 16.1.4. homova ------------

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  select(soil_code, where(is.character)) %>%
  droplevels() %>%
  column_to_rownames(var = "soil_code")

#{{need to transform in vector for betadisper if numeric}}
#### factors <- as.factor(select(factors$Plant_type)) 

#{{select distance}}
#{{select each factor to test}}
set.seed(seed)
Homovaresults <- betadisper(distcom, factors$type)
perm_Homovaresults_organ <- permutest(Homovaresults)
perm_Homovaresults_organ

set.seed(seed)
Homovaresults <- betadisper(distcom, factors$health_status)
perm_Homovaresults_health <- permutest(Homovaresults)
perm_Homovaresults_health


#{{extract homova p-values}}
hom_organ_pvalues <- perm_Homovaresults_organ$tab$`Pr(>F)`
hom_health_status_pvalues <- perm_Homovaresults_health$tab$`Pr(>F)`

## ------------ 18.1.5. permanova ------------

#log1p and Morisita similarity index was selected

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  select(where(is.character)) %>%
  column_to_rownames(var = "soil_code")

#{{select distance}}
#{{select all factors to test}}
set.seed(seed)
Permanovaresults <- adonis2(transcom~type*health_status, data = factors, 
                           method = "morisita", permutation = 9999)
Permanovaresults

#{{extract permanova R2 and p-values}}
perm_R2 <- Permanovaresults$R2
perm_pvalues <- Permanovaresults$`Pr(>F)`

#{{pairwise comparisons}}

set.seed(seed)
PairPermanovaresults <- pairwise.adonis2(transcom~organ_health, data = factors,
                                         method = "morisita", 
                                         sqrt.dist = TRUE, permutation = 9999, 
                                         p.adjust.methods = "fdr",
                                         na.action = na.omit)

PairPermanovaresults

#{{extract permanova R2 and p-values}}
pperm_leaf_healthy_root_healthy_R2 <- PairPermanovaresults$leaf_healthy_vs_root_healthy$R2[1]
pperm_leaf_healthy_root_healthy_pvalues <- PairPermanovaresults$leaf_healthy_vs_root_healthy$`Pr(>F)`[1]
pperm_leaf_healthy_stalk_healthy_R2 <- PairPermanovaresults$leaf_healthy_vs_stalk_healthy$R2[1]
pperm_leaf_healthy_stalk_healthy_pvalues <- PairPermanovaresults$leaf_healthy_vs_stalk_healthy$`Pr(>F)`[1]
pperm_leaf_healthy_leaf_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_leaf_disease$R2[1]
pperm_leaf_healthy_leaf_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_leaf_disease$`Pr(>F)`[1]
pperm_leaf_healthy_stalk_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_stalk_disease$R2[1]
pperm_leaf_healthy_stalk_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_stalk_disease$`Pr(>F)`[1]
pperm_leaf_healthy_root_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_root_disease$R2[1]
pperm_leaf_healthy_root_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_root_disease$`Pr(>F)`[1]
pperm_root_healthy_stalk_healthy_R2 <- PairPermanovaresults$root_healthy_vs_stalk_healthy$R2[1]
pperm_root_healthy_stalk_healthy_pvalues <- PairPermanovaresults$root_healthy_vs_stalk_healthy$`Pr(>F)`[1]
pperm_root_healthy_leaf_disease_R2 <- PairPermanovaresults$root_healthy_vs_leaf_disease$R2[1]
pperm_root_healthy_leaf_disease_pvalues <- PairPermanovaresults$root_healthy_vs_leaf_disease$`Pr(>F)`[1]
pperm_root_healthy_stalk_disease_R2 <- PairPermanovaresults$root_healthy_vs_stalk_disease$R2[1]
pperm_root_healthy_stalk_disease_pvalues <- PairPermanovaresults$root_healthy_vs_stalk_disease$`Pr(>F)`[1]
pperm_root_healthy_root_disease_R2 <- PairPermanovaresults$root_healthy_root_disease$R2[1]
pperm_root_healthy_root_disease_pvalues <- PairPermanovaresults$root_healthy_root_disease$`Pr(>F)`[1]
pperm_stalk_healthy_stalk_disease_R2 <- PairPermanovaresults$stalk_healthy_stalk_disease$R2[1]
pperm_stalk_healthy_stalk_disease_pvalues <- PairPermanovaresults$stalk_healthy_stalk_disease$`Pr(>F)`[1]
pperm_stalk_healthy_root_disease_R2 <- PairPermanovaresults$stalk_healthy_vs_root_disease$R2[1]
pperm_stalk_healthy_root_disease_pvalues <- PairPermanovaresults$stalk_healthy_vs_root_disease$`Pr(>F)`[1]
pperm_leaf_disease_stalk_disease_R2 <- PairPermanovaresults$leaf_disease_vs_stalk_disease$R2[1]
pperm_leaf_disease_stalk_disease_pvalues <- PairPermanovaresults$leaf_disease_vs_stalk_disease$`Pr(>F)`[1]
pperm_leaf_disease_root_disease_R2 <- PairPermanovaresults$leaf_disease_vs_root_disease$R2[1]
pperm_leaf_disease_root_disease_pvalues <- PairPermanovaresults$leaf_disease_vs_root_disease$`Pr(>F)`[1]
pperm_stalk_disease_root_disease_R2 <- PairPermanovaresults$stalk_disease_vs_root_disease$R2[1]
pperm_stalk_disease_root_disease_pvalues <- PairPermanovaresults$stalk_disease_vs_root_disease$`Pr(>F)`[1]

## ------------ 16.1.6. hierarchical clustering ------------

#{{https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html}}
#{{https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2}}
#{{http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning}}

## ---------------- 16.1.6.1. format data ----------------

#{{select distance matrix for hclust or pvclust for bootstrap value}}
set.seed(seed)
hc_sample <- hclust(distcom, method = "ward.D2")
set.seed(seed)
hc_taxa <- hclust(distcom_t, method = "ward.D2")

#{{convert hclust or pvclust object for ggplot2}}
dhc <- dendro_data(hc_sample, type = "rectangle") 

## ---------------- 16.1.6.2. plot hierachical clustering ----------------

#{{select file name}}
output <- "HC_bact_otu_organ_health_hellinger_euclidean_ward.pdf"

#{{format scales}}
y_min <- min(dhc$segments$yend)
y_max <- max(dhc$segments$yend)

hc2 <- ggplot(dhc$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  #{{check the order of label in dhc$labels}}
  geom_text(data = dhc$labels, aes(x, y, fontface = "italic",
                                   label = c("healthy roots", "healthy stalk", 
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "healthy roots", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "disease stalk", "healthy stalk",
                                             "disease stalk", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "healthy leaf", "healthy leaf",
                                             "healthy leaf", "disease leaf",
                                             "healthy leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease roots", "disease roots", 
                                             "healthy stalk", "healthy roots", 
                                             "healthy roots", "healthy roots",
                                             "healthy leaf", "healthy leaf", 
                                             "healthy leaf")),
            hjust = 1.1, vjust = 0, angle = 90, size = 4) +
  ylim(y_min-1.5, y_max+0.5) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        axis.line.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Euclidean distance") +
  xlab("")

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

## ------------ 16.1.5. heatmap {{ComplexHeatmap}} + hierarchical clustering ------------

#{{select file name}}
output <- "Heatmap_bact_phylum_soils_log1p_euclidean_ward.pdf"

#{{generate plot}}
set.seed(seed)
cluster_name <- cutree(hc_taxa, k = 3)

set.seed(seed)
heat1 <- Heatmap(t(transcom), name = "log1p",
                 #{{hierarchical clustering}}
                 #cluster_rows = hc_taxa,
                 cluster_columns = hc_sample,
                 #{{dendogram split}}
                 row_split = cluster_name, column_split = 3,
                 #{{theme format}}
                 border = TRUE,
                 row_names_gp = gpar(fontsize = 4),
                 show_row_names = FALSE,
                 heatmap_legend_param = list(
                   at = c(0, 5, 10),
                   labels = c("0", "5", "10"),
                   title = "log1p(abun)",
                   legend_height = unit(4, "cm"),
                   title_position = "leftcenter-rot"),
                 #{{Viridis palette}}
                 col = inferno(10))

##{{Creates a gTree object compatible with ggplot2}}
grobheat1 <- grid.grabExpr(draw(heat1))
##{{Convert gTree object in ggplot2 object}}
ggheat1 <- as_ggplot(grobheat1)

width <- 9
height <- 15

ggsave(file = output, width = width , height = height)
rm(grobheat1)


#{{extract taxa list for each cluster}}
taxacluster <- cluster_name %>%
  data.frame(OTU = names(cluster_name), cluster = unname(cluster_name)) %>%
  select(OTU,cluster)  %>%
  left_join(select(bactOTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

write_delim(taxacluster, path <- "HC_bact_phylum_soils_log1p_euclidean_ward_list.txt",
            delim="\t", col_names=TRUE)

#{{count abundance per cluster and taxa for a plant_type}}
taxacluster_count <- bactOTUs_rarefied_factor_t %>%
  filter(Type == "Soils", Plant_type == "Ceratonia_siliqua") %>%
  column_to_rownames(var = "Plant_type") %>%
  select(-soil_code,-Station_code, -Country, -Type, -Water_status) %>%
  t() %>%
  log1p() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  filter (Ceratonia_siliqua != "0") %>%
  left_join(select(taxacluster, OTU, cluster, genus), by = "OTU") %>%
  count(cluster, genus, wt = Ceratonia_siliqua, name = "reads_per_genus_per_cluster")

## ------------ 18.4.7. NMDS ------------
## ---------------- 18.4.7.1. Format data ----------------

#{{select distance matrix and generate nmds matrices}}
set.seed(seed)
com.nmds <- metaMDS(distcom, trymax = 100) 

#{{extraction of stress values}}
stress <- com.nmds$stress

#{{extraction of coordinates}}
data.scores <- as_tibble(scores(com.nmds),rownames = NA) %>%
  #{{Merge with sample data}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = "organ_health", type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #new category by merging two factors
  column_to_rownames(var = "soil_code")

#{{format annotation position}}
x_min <- min(data.scores$NMDS1)
x_max <- max(data.scores$NMDS1)
y_min <- min(data.scores$NMDS2)
y_max <- max(data.scores$NMDS2)

#{{format nmds statistic annotation}}
stress_annotation <- paste("Stress: ", round(stress, digits = 4), sep = "")

#{{format homova statictic annotation corresponding for both water status}}
homova_annotation <- paste("Homova\n",
                           "sugarcane organs: P = ", round(hom_organ_pvalues[1], digits = 4), "\n",
                           "sugarcane health: P = ", round(hom_health_status_pvalues[1], digits = 6),
                           sep = "")

#{{format permanova statictic annotation for both water status}}
permanova_annotation <- paste("Permanova\n",
                              "sugarcane organs: ", 
                              "R2 = ", round(perm_R2[1], digits = 4),
                              "; P = ", round(perm_pvalues[1], digits = 4), "\n",
                              "sugarcane health: ",
                              "R2 = ", round(perm_R2[2], digits = 4),
                              "; P = ", round(perm_pvalues[2], digits = 4), "\n",
                              "organs x health: ",
                              "R2 = ", round(perm_R2[3], digits = 4),
                              "; P = ", round(perm_pvalues[3], digits = 4),
                              sep = "")

## ---------------- 18.4.7.2. Plot NMDS ----------------

#{{select file name}}
output <- "nMDS_bact_organ_health_log1p_morisita.pdf"

#{{select and re-order the factors}}
flevels_organ <- c("leaf","stalk","root")
flevels_health <- c("healthy","disease")
flevels_organ_health <- c("leaf_healthy","stalk_healthy","root_healthy",
                          "leaf_disease","stalk_disease","root_disease")


#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{select color palette}}
colpalette_organ_health <- c("#46732EFF", "#46732EFF", "#46732EFF",
                             "#C80813FF", "#C80813FF", "#C80813FF")
colpalette_health <- c("#46732EFF", "#C80813FF")

colpalette_organ <- c("#FED439FF", "#D5E4A2FF", "#91331FFF")

#{{plot}}
{
  nmds <- ggscatter(data = data.scores, x = "NMDS1", y = "NMDS2",
                    color = "organ_health", shape = "type",
                    mean.point = TRUE,
                    mean.point.size = 10,
                    star.plot = TRUE,
                    star.plot.lty = 3,
                    ellipse = TRUE,
                    ellipse.type = "confidence",
                    ellipse.border.remove = TRUE,
                    ellipse.alpha = 0.2) +
    #{{use geom_point and scale shape for global Plant_type and Water use analysis}}
    geom_point(aes(shape = type, color = organ_health), size = 3) +
    scale_shape_discrete(limits = flevels_organ,
                         name = "Sugarcane organs",
                         labels = c("leaf", "stalk", "root"),
                         guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
    theme_bw(base_size = 16) +
    theme(legend.title = element_text(vjust = 0.5, hjust = 0, size = 12),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.text = element_text(hjust = 0, size = 10),
          legend.background = element_blank(),
          aspect.ratio = 1) +
    ggtitle("") +
    ylim(y_min - 0.1,y_max) +
    xlim(x_min,x_max) +
    annotate("text", x = x_min + 0.05, y = y_max, hjust = 0, label = stress_annotation) +
    annotate("text", x = x_min + 0.05, y = y_min + 1.0, hjust = 0, vjust = 0.5, label = homova_annotation) +
    annotate("text", x = x_min + 0.05, y = y_min + 0.2, hjust = 0, vjust = 0.5, label = permanova_annotation) +
    scale_colour_manual(name = "Sugarcane health",
                        limits = flevels_organ_health,
                        labels = c("healthy", "healthy", "healthy",
                                   "disease", "disease", "disease"),
                        values = colpalette_organ_health,
                        guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
    scale_fill_manual(limits = flevels_organ_health,
                      labels = c("healthy", "healthy", "healthy",
                                 "disease", "disease", "disease"),
                      values = colpalette_organ_health) +
    guides(fill = "none") +
    coord_equal()
}

#{add marginal densities}}
{
  #{{along x axis}} 
  xdens <- axis_canvas(nmds, axis = "x") +
    geom_density(data = data.scores, aes(x = NMDS1, fill = type),
                 alpha = 0.7, size = 0.2) +
    scale_fill_manual(limits = flevels_organ,
                      values = colpalette_organ)
  #{{along y axis}}
  #{{Need to set coord_flip = TRUE, if you plan to use coord_flip()}}
  ydens <- axis_canvas(nmds, axis = "y", coord_flip = TRUE) +
    geom_density(data = data.scores, aes(x = NMDS2, fill = health_status),
                 alpha = 0.7, size = 0.2) +
    scale_fill_manual(limits = flevels_health,
                      values = colpalette_health) +
    coord_flip()
}

#{{assembly of density plots to nmds plot}}   
nmds <- insert_xaxis_grob(nmds, xdens, grid::unit(.2, "null"), position = "top")
nmds <- insert_yaxis_grob(nmds, ydens, grid::unit(.2, "null"), position = "right")

ggdraw(nmds)

width <- 12
height <- 12

ggsave(file = output, width = width , height = height)

## -------- 18.3. fungi ------
## ------------ 18.3.1 format data ------------
## ---------------- 18.3.1.1. OTU level ----------------

#{{select the data}}
com <- ITS2_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  droplevels() %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  filter(diversification != "RCTF") %>% #diversification level RCTF is deleted 
  droplevels() %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  #{{delete OTUs with no read}}
  gather(key = "OTU", value = "reads", -soil_code) %>%
  group_by(OTU) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(OTU, reads) %>%
  column_to_rownames(var = "soil_code")

## ---------------- 18.3.1.3. Taxonomic level ----------------

#{{select the data}}
type <- "organ_health"

#{{select the data}}
com <- ITS2_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  ####filter(type == type) %>%
  droplevels() %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  column_to_rownames(var = "soil_code") %>%
  #{{select the taxonomic range}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  ####mutate(OTU = as.numeric(OTU)) %>%
  left_join(select(ITS2_tax,OTU, genus), by = "OTU") %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  #{{delete taxonomic group with no read}}
  gather(key = "taxonomy", value = "reads", -genus) %>%
  group_by(taxonomy) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(taxonomy, reads) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA)

## ------------ 18.3.2. data transformation & standardization ------------

#{{no data tranformation}}
transcom <- com 

#{{select abundance data tranformation}}
transcom <- log1p(com)  #{{sqrt() or log1p()}}

#{{select abundance data tranformation and standardization}}
transcom <- decostand(com,"hellinger", MARGIN = 1)  #{{standardization method = "total" or "hellinger" or "normalize" or "chi.square")
transcom <- wisconsin(com)  #{{standardiozation "wisconsin"}}
transcom <- decostand(log1p(com),"hellinger") #{{abundance data tranformation and standardization}}

barplot(table(unlist(transcom)),las = 1)

## ------------ 18.3.3 distance matrices ------------
## ---------------- 18.3.3.1. samples ----------------

distcom <- vegdist(transcom, method = "bray")
distcom <- vegdist(transcom, method = "horn")
distcom <- vegdist(transcom, method = "morisita")
distcom <- vegdist(transcom, binary = TRUE, method = "raup")
distcom <- vegdist(transcom, method = "jaccard")
#{{adapted with standardized data}
distcom <- vegdist(transcom, method = "euclidean") 

## ---------------- 18.3.3.2. OTU or Taxa ----------------

distcom_t <- vegdist(t(transcom), method = "bray")
distcom_t <- vegdist(t(transcom), method = "horn")
distcom_t <- vegdist(t(transcom), method = "morisita")
distcom_t <- vegdist(t(transcom), binary = TRUE, method = "raup")
distcom_t <- vegdist(t(transcom), method = "jaccard")
#{{adapted with standardized data}
distcom_t <- vegdist(t(transcom), method = "euclidean")

## ------------ 18.3.4. homova ------------

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  droplevels() %>%
  column_to_rownames(var = "soil_code")

#{{need to transform in vector for betadisper if numeric}}
#### factors <- as.factor(select(factors$Plant_type)) 

#{{select distance}}
#{{select each factor to test}}
set.seed(seed)
Homovaresults <- betadisper(distcom, factors$diversification)
perm_Homovaresults_div <- permutest(Homovaresults)
perm_Homovaresults_div

set.seed(seed)
Homovaresults <- betadisper(distcom, factors$cluster_code)
perm_Homovaresults_cluster <- permutest(Homovaresults)
perm_Homovaresults_cluster

#{{extract homova p-values}}
hom_div_pvalues <- perm_Homovaresults_div$tab$`Pr(>F)`
hom_cluster_pvalues <- perm_Homovaresults_cluster$tab$`Pr(>F)`

## ------------ 18.1.5. permanova ------------

#no data transformation, no standardization and Horn similarity index was selected

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  select(where(is.character)) %>%
  column_to_rownames(var = "soil_code")

#{{select distance}}
#{{select all factors to test}}
set.seed(seed)
Permanovaresults <- adonis3(transcom~cluster_code*diversification, data = factors, 
                            method = "horn", permutation = 9999)
Permanovaresults$aov.tab

#{{extract permanova R2 and p-values}}
perm_R2 <- Permanovaresults$aov.tab$R2
perm_pvalues <- Permanovaresults$aov.tab$`Pr(>F)`

#{{pairwise comparisons}}

set.seed(seed)
PairPermanovaresults <- pairwise.adonis2(transcom~cluster_code, data = factors,
                                         method = "horn", 
                                         sqrt.dist = TRUE, permutation = 9999, 
                                         p.adjust.methods = "fdr",
                                         na.action = na.omit)

PairPermanovaresults

#{{extract permanova R2 and p-values}}
pperm_leaf_healthy_root_healthy_R2 <- PairPermanovaresults$leaf_healthy_vs_root_healthy$R2[1]
pperm_leaf_healthy_root_healthy_pvalues <- PairPermanovaresults$leaf_healthy_vs_root_healthy$`Pr(>F)`[1]
pperm_leaf_healthy_stalk_healthy_R2 <- PairPermanovaresults$leaf_healthy_vs_stalk_healthy$R2[1]
pperm_leaf_healthy_stalk_healthy_pvalues <- PairPermanovaresults$leaf_healthy_vs_stalk_healthy$`Pr(>F)`[1]
pperm_leaf_healthy_leaf_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_leaf_disease$R2[1]
pperm_leaf_healthy_leaf_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_leaf_disease$`Pr(>F)`[1]
pperm_leaf_healthy_stalk_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_stalk_disease$R2[1]
pperm_leaf_healthy_stalk_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_stalk_disease$`Pr(>F)`[1]
pperm_leaf_healthy_root_disease_R2 <- PairPermanovaresults$leaf_healthy_vs_root_disease$R2[1]
pperm_leaf_healthy_root_disease_pvalues <- PairPermanovaresults$leaf_healthy_vs_root_disease$`Pr(>F)`[1]
pperm_root_healthy_stalk_healthy_R2 <- PairPermanovaresults$root_healthy_vs_stalk_healthy$R2[1]
pperm_root_healthy_stalk_healthy_pvalues <- PairPermanovaresults$root_healthy_vs_stalk_healthy$`Pr(>F)`[1]
pperm_root_healthy_leaf_disease_R2 <- PairPermanovaresults$root_healthy_vs_leaf_disease$R2[1]
pperm_root_healthy_leaf_disease_pvalues <- PairPermanovaresults$root_healthy_vs_leaf_disease$`Pr(>F)`[1]
pperm_root_healthy_stalk_disease_R2 <- PairPermanovaresults$root_healthy_vs_stalk_disease$R2[1]
pperm_root_healthy_stalk_disease_pvalues <- PairPermanovaresults$root_healthy_vs_stalk_disease$`Pr(>F)`[1]
pperm_root_healthy_root_disease_R2 <- PairPermanovaresults$root_healthy_root_disease$R2[1]
pperm_root_healthy_root_disease_pvalues <- PairPermanovaresults$root_healthy_root_disease$`Pr(>F)`[1]
pperm_stalk_healthy_stalk_disease_R2 <- PairPermanovaresults$stalk_healthy_stalk_disease$R2[1]
pperm_stalk_healthy_stalk_disease_pvalues <- PairPermanovaresults$stalk_healthy_stalk_disease$`Pr(>F)`[1]
pperm_stalk_healthy_root_disease_R2 <- PairPermanovaresults$stalk_healthy_vs_root_disease$R2[1]
pperm_stalk_healthy_root_disease_pvalues <- PairPermanovaresults$stalk_healthy_vs_root_disease$`Pr(>F)`[1]
pperm_leaf_disease_stalk_disease_R2 <- PairPermanovaresults$leaf_disease_vs_stalk_disease$R2[1]
pperm_leaf_disease_stalk_disease_pvalues <- PairPermanovaresults$leaf_disease_vs_stalk_disease$`Pr(>F)`[1]
pperm_leaf_disease_root_disease_R2 <- PairPermanovaresults$leaf_disease_vs_root_disease$R2[1]
pperm_leaf_disease_root_disease_pvalues <- PairPermanovaresults$leaf_disease_vs_root_disease$`Pr(>F)`[1]
pperm_stalk_disease_root_disease_R2 <- PairPermanovaresults$stalk_disease_vs_root_disease$R2[1]
pperm_stalk_disease_root_disease_pvalues <- PairPermanovaresults$stalk_disease_vs_root_disease$`Pr(>F)`[1]


## ------------ 18.3.4. hierarchical clustering ------------

#{{https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html}}
#{{https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2}}
#{{http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning}}

## ---------------- 18.3.4.1. format data ----------------

#{{select distance matrix for hclust or pvclust; bootstrap value}}
set.seed(seed)
hc_sample <- hclust(distcom, method = "ward.D2")
set.seed(seed)
hc_taxa <- hclust(distcom_t, method = "ward.D2")

#{{convert hclust or pvclust object for ggplot2}}
dhc <- dendro_data(hc_sample, type = "rectangle") 

## ---------------- 18.3.4.2. plot hierachical clustering ----------------

#{{select file name}}
output <- "HC_fung_otu_organ_health_hellinger_euclidean_ward.pdf"

#{{format scales}}
y_min <- min(dhc$segments$yend)
y_max <- max(dhc$segments$yend)

hc2 <- ggplot(dhc$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  #{{check the order of label in dhc$labels}}
  geom_text(data = dhc$labels, aes(x, y, fontface = "italic",
                                   label = c("healthy roots", "healthy stalk", 
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "healthy roots", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "disease stalk", "healthy stalk",
                                             "disease stalk", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "healthy leaf", "healthy leaf",
                                             "healthy leaf", "disease leaf",
                                             "healthy leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease roots", "disease roots", 
                                             "healthy stalk", "healthy roots", 
                                             "healthy roots", "healthy roots",
                                             "healthy leaf", "healthy leaf", 
                                             "healthy leaf")),
            hjust = 1.1, vjust = 0, angle = 90, size = 4) +
  ylim(y_min-1.5, y_max+0.5) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        axis.line.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Euclidean distance") +
  xlab("")

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

## ------------ 18.3.5. heatmap {{ComplexHeatmap}} + hierarchical clustering ------------

#{{select file name}}
output <-"Heatmap_fung_phylum_soils_hellinger_euclidean_ward.pdf"

#{{generate plot}}
set.seed(seed)
cluster_name <- cutree(hc_taxa, k = 2)

set.seed(seed)
heat2 <- Heatmap(t(transcom), name = "Hellinger",
                 #{{hierarchical clustering}}
                 #cluster_rows = hc_taxa,
                 cluster_columns = hc_sample,
                 #{{dendogram split}}
                 row_split = cluster_name, column_split = 2,
                 #{{theme format}}
                 border = TRUE,
                 row_names_gp = gpar(fontsize = 4),
                 show_row_names = FALSE,
                 heatmap_legend_param = list(
                   at = c(0, 5, 10),
                   labels = c("0", "5", "10"),
                   title = "log1p(abun)",
                   legend_height = unit(4, "cm"),
                   title_position = "leftcenter-rot"),
                 #{{Viridis palette}}
                 col =inferno(10))

##{{Creates a gTree object compatible with ggplot2}}
grobheat2 <- grid.grabExpr(draw(heat2))
##{{Convert gTree object in ggplot2 object}}
ggheat2 <- as_ggplot(grobheat2)


width <- 9
height <- 15

ggsave(file = output, width = width , height = height)
rm(grobheat2)

#{{extract taxa list for each cluster}}
taxacluster <- cluster_name %>%
  data.frame(OTU = names(cluster_name), cluster = unname(cluster_name)) %>%
  select(OTU,cluster)  %>%
  left_join(select(ITS2OTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

write_delim(taxacluster, path <- "HC_fung_phylum_soils_log1p_euclidean_ward_list.txt",
            delim="\t", col_names=TRUE)


#{{count abundance per cluster and taxa for a plant_type}}
taxacluster_count <- ITS2OTUs_rarefied_factor_t %>%
  filter(Type == "Soils", Plant_type == "Ceratonia_siliqua") %>%
  column_to_rownames(var = "Plant_type") %>%
  select(-soil_code,-Station_code, -Country, -Type, -Water_status) %>%
  t() %>%
  log1p() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  filter (Ceratonia_siliqua != "0") %>%
  left_join(select(taxacluster, OTU, cluster, genus), by = "OTU") %>%
  count(cluster, genus, wt = Ceratonia_siliqua, name = "reads_per_genus_per_cluster")



## ------------ 18.4.7. NMDS ------------
## ---------------- 18.4.7.1. Format data ----------------

#{{select distance matrix and generate nmds matrices}}
set.seed(seed)
com.nmds <- metaMDS(distcom, trymax = 100) 

#{{extraction of stress values}}
stress <- com.nmds$stress

#{{extraction of coordinates}}
data.scores <- as_tibble(scores(com.nmds),rownames = NA) %>%
  #{{Merge with sample data}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = "cluster_div", cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #new category by merging two factors
  column_to_rownames(var = "soil_code")

#{{format annotation position}}
x_min <- min(data.scores$NMDS1)
x_max <- max(data.scores$NMDS1)
y_min <- min(data.scores$NMDS2)
y_max <- max(data.scores$NMDS2)

#{{format nmds statistic annotation}}
stress_annotation <- paste("Stress: ", round(stress, digits = 4), sep = "")

#{{format homova statictic annotation corresponding for both water status}}
homova_annotation <- paste("Homova\n",
                           "cluster: P = ", round(hom_cluster_pvalues[1], digits = 4), "\n",
                           "diversification: P = ", round(hom_div_pvalues[1], digits = 6),
                           sep = "")

#{{format permanova statictic annotation for both water status}}
permanova_annotation <- paste("Permanova\n",
                              "cluster: ", 
                              "R2 = ", round(perm_R2[1], digits = 4),
                              "; P = ", round(perm_pvalues[1], digits = 4), "\n",
                              "diversification: ",
                              "R2 = ", round(perm_R2[2], digits = 4),
                              "; P = ", round(perm_pvalues[2], digits = 4), "\n",
                              "cluster x diversification: ",
                              "R2 = ", round(perm_R2[3], digits = 4),
                              "; P = ", round(perm_pvalues[3], digits = 4),
                              sep = "")

## ---------------- 18.4.7.2. Plot NMDS ----------------

#{{select file name}}
output <- "nMDS_fung_clu_notrans_horn.pdf"
{
  
#{{select and re-order the factors}}
flevels_organ <- c("soil", "root")
flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA")
flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{select color palette}}
colpalette_div <- c("#D2AF81FF", "#FED439FF", "#91331FFF", "#F05C3BFF", "#D5E4A2FF") #diversification level RCTF has been deleted
#{{viridis palette visualisation}}
show_col(viridis(10))
colpalette_clu <- c("#FDE725FF", "#B4DE2CFF", "#35B779FF", "#26828EFF", "#31688EFF")


#{{plot}}
{
  nmds <- ggscatter(data = data.scores, x = "NMDS1", y = "NMDS2",
                    color = "cluster_code",
                    size = 0.5,
                    mean.point = TRUE,
                    mean.point.size = 10,
                    ####star.plot = TRUE,
                    ####star.plot.lty = 3,
                    ellipse = TRUE,
                    ellipse.level = 0.95,
                    ellipse.type = "confidence",
                    ellipse.border.remove = TRUE,
                    ellipse.alpha = 0.2) +
    #{{use geom_point and scale shape for diversification and cluster}}
    geom_point(aes(color = cluster_code), size = 0.5) +
    theme_bw(base_size = 16) +
    theme(legend.title = element_text(vjust = 0.5, hjust = 0, size = 12),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.text = element_text(hjust = 0, size = 10),
          legend.background = element_blank(),
          aspect.ratio = 1) +
    ggtitle("") +
    ylim(y_min,y_max) +
    xlim(x_min,x_max) +
    annotate("text", x = x_min, y = y_max, hjust = 0, label = stress_annotation) +
    annotate("text", x = x_min, y = y_min + 0.5, hjust = 0, vjust = 0.5, label = homova_annotation) +
    annotate("text", x = x_min, y = y_min, hjust = 0, vjust = 0.5, label = permanova_annotation) +
    scale_colour_manual(name = "",
                        limits = flevels_clu,
                        labels = c("N\n(dry)", "M\n(moderatly dry)",
                                   "B\n(variable)", "R\n(wet)",
                                   "F\n(very wet)"),
                        values = colpalette_clu,
                        guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
    scale_fill_manual(limits = flevels_clu,
                      labels = c("N\n(dry)", "M\n(moderatly dry)",
                                 "B\n(variable)", "R\n(wet)",
                                 "F\n(very wet)"),
                      values = colpalette_clu) +
    guides(fill = "none") +
    coord_equal()
}

#{add marginal densities}}
{
  #{{along x axis}} 
  xdens <- axis_canvas(nmds, axis = "x") +
    geom_density(data = data.scores, aes(x = NMDS1, fill = diversification),
                 alpha = 0.7, size = 0.2) +
    scale_fill_manual(limits = flevels_div,
                      values = colpalette_div)
  #{{along y axis}}
  #{{Need to set coord_flip = TRUE, if you plan to use coord_flip()}}
  ydens <- axis_canvas(nmds, axis = "y", coord_flip = TRUE) +
    geom_density(data = data.scores, aes(x = NMDS2, fill = diversification),
                 alpha = 0.7, size = 0.2) +
    scale_fill_manual(limits = flevels_div ,
                      values = colpalette_div) +
    coord_flip()
  nmds <- insert_xaxis_grob(nmds, xdens, grid::unit(.2, "null"), position = "top")
  nmds <- insert_yaxis_grob(nmds, ydens, grid::unit(.2, "null"), position = "right")
}

#{{assembly of density plots to nmds plot}}  
{
  ggdraw(nmds)
  
  width <- 12
  height <- 12
  
  ggsave(file = output, width = width , height = height)
}


}

## -------- 18.2. amf ------
## ------------ 18.2.1 format data ------------
## ---------------- 18.2.1.1. OTU level ----------------

#{{select the data}}
com <- amf_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  droplevels() %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #filter(diversification != "RCTF") %>% #diversification level RCTF is deleted 
  droplevels() %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  #{{delete OTUs with no read}}
  gather(key = "OTU", value = "reads", -soil_code) %>%
  group_by(OTU) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(OTU, reads) %>%
  column_to_rownames(var = "soil_code")

## ---------------- 18.3.1.3. Taxonomic level ----------------

#{{select the data}}
type <- "organ_health"

#{{select the data}}
com <- ITS2_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatment}}
  ####filter(type == type) %>%
  droplevels() %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{select only OTUs}}
  select(soil_code, where(is.numeric)) %>%
  column_to_rownames(var = "soil_code") %>%
  #{{select the taxonomic range}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  ####mutate(OTU = as.numeric(OTU)) %>%
  left_join(select(ITS2_tax,OTU, genus), by = "OTU") %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  #{{delete taxonomic group with no read}}
  gather(key = "taxonomy", value = "reads", -genus) %>%
  group_by(taxonomy) %>% 
  mutate(sum_group = sum(reads)) %>%
  filter(sum_group > 0) %>%
  select(-sum_group) %>% 
  ungroup() %>%
  spread(taxonomy, reads) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA)

## ------------ 18.3.2. data transformation & standardization ------------

#{{no data tranformation}}
transcom <- com 

#{{select abundance data tranformation}}
transcom <- log1p(com)  #{{sqrt() or log1p()}}

#{{select abundance data tranformation and standardization}}
transcom <- decostand(com,"hellinger", MARGIN = 1)  #{{standardization method = "total" or "hellinger" or "normalize" or "chi.square")
transcom <- wisconsin(com)  #{{standardiozation "wisconsin"}}
transcom <- decostand(log1p(com),"hellinger") #{{abundance data tranformation and standardization}}

barplot(table(unlist(transcom)),las = 1)

## ------------ 18.3.3 distance matrices ------------
## ---------------- 18.3.3.1. samples ----------------

distcom <- vegdist(transcom, method = "bray")
distcom <- vegdist(transcom, method = "horn")
distcom <- vegdist(transcom, method = "morisita")
distcom <- vegdist(transcom, binary = TRUE, method = "raup")
distcom <- vegdist(transcom, method = "jaccard")
#{{adapted with standardized data}
distcom <- vegdist(transcom, method = "euclidean") 

## ---------------- 18.3.3.2. OTU or Taxa ----------------

distcom_t <- vegdist(t(transcom), method = "bray")
distcom_t <- vegdist(t(transcom), method = "horn")
distcom_t <- vegdist(t(transcom), method = "morisita")
distcom_t <- vegdist(t(transcom), binary = TRUE, method = "raup")
distcom_t <- vegdist(t(transcom), method = "jaccard")
#{{adapted with standardized data}
distcom_t <- vegdist(t(transcom), method = "euclidean")

## ------------ 18.3.4. homova ------------

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  droplevels() %>%
  column_to_rownames(var = "soil_code")

#{{need to transform in vector for betadisper if numeric}}
#### factors <- as.factor(select(factors$Plant_type)) 

#{{select distance}}
#{{select each factor to test}}
set.seed(seed)
Homovaresults <- betadisper(distcom, factors$diversification)
perm_Homovaresults_div <- permutest(Homovaresults)
perm_Homovaresults_div

set.seed(seed)
Homovaresults <- betadisper(distcom, factors$cluster_code)
perm_Homovaresults_cluster <- permutest(Homovaresults)
perm_Homovaresults_cluster

#{{extract homova p-values}}
hom_div_pvalues <- perm_Homovaresults_div$tab$`Pr(>F)`
hom_cluster_pvalues <- perm_Homovaresults_cluster$tab$`Pr(>F)`

## ------------ 18.1.5. permanova ------------

#no data transformation, no standardization and Horn similarity index was selected

#{{select factors}}
factors <- transcom %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = cluster_div , cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #create a new column
  select(where(is.character)) %>%
  column_to_rownames(var = "soil_code")

#{{select distance}}
#{{select all factors to test}}
set.seed(seed)
Permanovaresults <- adonis3(transcom~cluster_code*diversification, data = factors, 
                            method = "horn", permutation = 9999)
Permanovaresults$aov.tab

#{{extract permanova R2 and p-values}}
perm_R2 <- Permanovaresults$aov.tab$R2
perm_pvalues <- Permanovaresults$aov.tab$`Pr(>F)`

#{{pairwise comparisons}}

set.seed(seed)
PairPermanovaresults <- pairwise.adonis2(transcom~cluster_code, data = factors,
                                         method = "horn", 
                                         sqrt.dist = TRUE, permutation = 9999, 
                                         p.adjust.methods = "fdr",
                                         na.action = na.omit)

PairPermanovaresults

#{{extract permanova R2 and p-values}}
pperm_B_vs_F_R2 <- PairPermanovaresults$B_vs_F$R2[1]
pperm_B_vs_F_pvalues <- PairPermanovaresults$B_vs_F$`Pr(>F)`[1]

pperm_B_vs_M_R2 <- PairPermanovaresults$B_vs_M$R2[1]
pperm_B_vs_M_pvalues <- PairPermanovaresults$B_vs_M$`Pr(>F)`[1]

pperm_B_vs_N_R2 <- PairPermanovaresults$B_vs_N$R2[1]
pperm_B_vs_N_pvalues <- PairPermanovaresults$B_vs_N$`Pr(>F)`[1]

pperm_B_vs_R_R2 <- PairPermanovaresults$B_vs_R$R2[1]
pperm_B_vs_R_pvalues <- PairPermanovaresults$B_vs_R$`Pr(>F)`[1]

pperm_B_vs_R_R2 <- PairPermanovaresults$B_vs_R$R2[1]
pperm_B_vs_R_pvalues <- PairPermanovaresults$B_vs_R$`Pr(>F)`[1]

pperm_F_vs_M_R2 <- PairPermanovaresults$F_vs_M$R2[1]
pperm_F_vs_M_pvalues <- PairPermanovaresults$F_vs_M$`Pr(>F)`[1]

pperm_F_vs_N_R2 <- PairPermanovaresults$F_vs_N$R2[1]
pperm_F_vs_N_pvalues <- PairPermanovaresults$F_vs_N$`Pr(>F)`[1]

pperm_F_vs_R_R2 <- PairPermanovaresults$F_vs_R$R2[1]
pperm_F_vs_R_pvalues <- PairPermanovaresults$F_vs_R$`Pr(>F)`[1]

pperm_M_vs_N_R2 <- PairPermanovaresults$M_vs_N$R2[1]
pperm_M_vs_N_pvalues <- PairPermanovaresults$M_vs_N$`Pr(>F)`[1]

pperm_N_vs_R_R2 <- PairPermanovaresults$N_vs_R$R2[1]
pperm_N_vs_R_pvalues <- PairPermanovaresults$N_vs_R$`Pr(>F)`[1]


## ------------ 18.3.4. hierarchical clustering ------------

#{{https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html}}
#{{https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2}}
#{{http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning}}

## ---------------- 18.3.4.1. format data ----------------

#{{select distance matrix for hclust or pvclust; bootstrap value}}
set.seed(seed)
hc_sample <- hclust(distcom, method = "ward.D2")
set.seed(seed)
hc_taxa <- hclust(distcom_t, method = "ward.D2")

#{{convert hclust or pvclust object for ggplot2}}
dhc <- dendro_data(hc_sample, type = "rectangle") 

## ---------------- 18.3.4.2. plot hierachical clustering ----------------

#{{select file name}}
output <- "HC_fung_otu_organ_health_hellinger_euclidean_ward.pdf"

#{{format scales}}
y_min <- min(dhc$segments$yend)
y_max <- max(dhc$segments$yend)

hc2 <- ggplot(dhc$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  #{{check the order of label in dhc$labels}}
  geom_text(data = dhc$labels, aes(x, y, fontface = "italic",
                                   label = c("healthy roots", "healthy stalk", 
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "healthy stalk", "disease stalk",
                                             "healthy roots", "healthy stalk",
                                             "disease stalk", "healthy roots",
                                             "healthy stalk", "healthy roots",
                                             "disease stalk", "healthy stalk",
                                             "disease stalk", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "disease leaf",
                                             "disease leaf", "healthy leaf",
                                             "healthy leaf", "healthy leaf",
                                             "healthy leaf", "disease leaf",
                                             "healthy leaf", "healthy leaf",
                                             "disease leaf", "disease leaf",
                                             "disease roots", "disease roots", 
                                             "healthy stalk", "healthy roots", 
                                             "healthy roots", "healthy roots",
                                             "healthy leaf", "healthy leaf", 
                                             "healthy leaf")),
            hjust = 1.1, vjust = 0, angle = 90, size = 4) +
  ylim(y_min-1.5, y_max+0.5) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        axis.line.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Euclidean distance") +
  xlab("")

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

## ------------ 18.3.5. heatmap {{ComplexHeatmap}} + hierarchical clustering ------------

#{{select file name}}
output <-"Heatmap_fung_phylum_soils_hellinger_euclidean_ward.pdf"

#{{generate plot}}
set.seed(seed)
cluster_name <- cutree(hc_taxa, k = 2)

set.seed(seed)
heat2 <- Heatmap(t(transcom), name = "Hellinger",
                 #{{hierarchical clustering}}
                 #cluster_rows = hc_taxa,
                 cluster_columns = hc_sample,
                 #{{dendogram split}}
                 row_split = cluster_name, column_split = 2,
                 #{{theme format}}
                 border = TRUE,
                 row_names_gp = gpar(fontsize = 4),
                 show_row_names = FALSE,
                 heatmap_legend_param = list(
                   at = c(0, 5, 10),
                   labels = c("0", "5", "10"),
                   title = "log1p(abun)",
                   legend_height = unit(4, "cm"),
                   title_position = "leftcenter-rot"),
                 #{{Viridis palette}}
                 col =inferno(10))

##{{Creates a gTree object compatible with ggplot2}}
grobheat2 <- grid.grabExpr(draw(heat2))
##{{Convert gTree object in ggplot2 object}}
ggheat2 <- as_ggplot(grobheat2)


width <- 9
height <- 15

ggsave(file = output, width = width , height = height)
rm(grobheat2)

#{{extract taxa list for each cluster}}
taxacluster <- cluster_name %>%
  data.frame(OTU = names(cluster_name), cluster = unname(cluster_name)) %>%
  select(OTU,cluster)  %>%
  left_join(select(ITS2OTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

write_delim(taxacluster, path <- "HC_fung_phylum_soils_log1p_euclidean_ward_list.txt",
            delim="\t", col_names=TRUE)


#{{count abundance per cluster and taxa for a plant_type}}
taxacluster_count <- ITS2OTUs_rarefied_factor_t %>%
  filter(Type == "Soils", Plant_type == "Ceratonia_siliqua") %>%
  column_to_rownames(var = "Plant_type") %>%
  select(-soil_code,-Station_code, -Country, -Type, -Water_status) %>%
  t() %>%
  log1p() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  filter (Ceratonia_siliqua != "0") %>%
  left_join(select(taxacluster, OTU, cluster, genus), by = "OTU") %>%
  count(cluster, genus, wt = Ceratonia_siliqua, name = "reads_per_genus_per_cluster")



## ------------ 18.4.7. NMDS ------------
## ---------------- 18.4.7.1. Format data ----------------

#{{select distance matrix and generate nmds matrices}}
set.seed(seed)
com.nmds <- metaMDS(distcom, trymax = 100) 

#{{extraction of stress values}}
stress <- com.nmds$stress

#{{extraction of coordinates}}
data.scores <- as_tibble(scores(com.nmds),rownames = NA) %>%
  #{{Merge with sample data}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = "cluster_div", cluster_code, diversification, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>% #new category by merging two factors
  column_to_rownames(var = "soil_code")

#{{format annotation position}}
x_min <- min(data.scores$NMDS1)
x_max <- max(data.scores$NMDS1)
y_min <- min(data.scores$NMDS2)
y_max <- max(data.scores$NMDS2)

#{{format nmds statistic annotation}}
stress_annotation <- paste("Stress: ", round(stress, digits = 4), sep = "")

#{{format homova statictic annotation corresponding for both water status}}
homova_annotation <- paste("Homova\n",
                           "cluster: P = ", round(hom_cluster_pvalues[1], digits = 4), "\n",
                           "diversification: P = ", round(hom_div_pvalues[1], digits = 6),
                           sep = "")

#{{format permanova statictic annotation for both water status}}
permanova_annotation <- paste("Permanova\n",
                              "cluster: ", 
                              "R2 = ", round(perm_R2[1], digits = 4),
                              "; P = ", round(perm_pvalues[1], digits = 4), "\n",
                              "diversification: ",
                              "R2 = ", round(perm_R2[2], digits = 4),
                              "; P = ", round(perm_pvalues[2], digits = 4), "\n",
                              "cluster x diversification: ",
                              "R2 = ", round(perm_R2[3], digits = 4),
                              "; P = ", round(perm_pvalues[3], digits = 4),
                              sep = "")

                              

ppermanova_annotation <- paste("Pairwise permanova\n",
                              "B vs F: ", 
                              "R2 = ", round(pperm_B_vs_F_R2[1], digits = 4),
                              "; P = ", round(pperm_B_vs_F_pvalues[1], digits = 4), "\n",
                              "B vs M: ", 
                              "R2 = ", round(pperm_B_vs_M_R2[1], digits = 4),
                              "; P = ", round(pperm_B_vs_M_pvalues[1], digits = 4), "\n",
                              "B vs N: ", 
                              "R2 = ", round(pperm_B_vs_N_R2[1], digits = 4),
                              "; P = ", round(pperm_B_vs_N_pvalues[1], digits = 4), "\n",
                              "B vs R: ", 
                              "R2 = ", round(pperm_B_vs_R_R2[1], digits = 4),
                              "; P = ", round(pperm_B_vs_R_pvalues[1], digits = 4), "\n",
                              "F vs M: ", 
                              "R2 = ", round(pperm_F_vs_M_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_M_pvalues[1], digits = 4), "\n",
                              "F vs N: ", 
                              "R2 = ", round(pperm_F_vs_N_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_N_pvalues[1], digits = 4), "\n",
                              "F vs R: ", 
                              "R2 = ", round(pperm_F_vs_R_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_R_pvalues[1], digits = 4), "\n",
                              "M vs N: ", 
                              "R2 = ", round(pperm_F_vs_N_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_N_pvalues[1], digits = 4), "\n",
                              "M vs R: ", 
                              "R2 = ", round(pperm_F_vs_N_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_N_pvalues[1], digits = 4), "\n",
                              "N vs R: ", 
                              "R2 = ", round(pperm_F_vs_N_R2[1], digits = 4),
                              "; P = ", round(pperm_F_vs_N_pvalues[1], digits = 4), "\n",
                              sep = "")

## ---------------- 18.4.7.2. Plot NMDS ----------------

#{{select file name}}
output <- "nMDS_amf_clu_notrans_horn.pdf"
{
  
  #{{select and re-order the factors}}
  flevels_organ <- c("soil", "root")
  flevels_div <- c("RCM", "RCB", "RCT", "RCTB", "RCTBA","RCF")
  flevels_clu <- c("N","M","B","R","F") #from dry, moderate, moderate-variable, wet, very wet 
  
  #{{simpson palette visualisation}}
  show_col(pal_simpsons("springfield")(16))
  #{{select color palette}}
  colpalette_div <- c("#D2AF81FF", "#FED439FF", "#91331FFF", "#F05C3BFF", "#D5E4A2FF")
  #{{viridis palette visualisation}}
  show_col(viridis(10))
  colpalette_clu <- c("#FDE725FF", "#B4DE2CFF", "#35B779FF", "#26828EFF", "#31688EFF")
  
  
  #{{plot}}
  {
    nmds <- ggscatter(data = data.scores, x = "NMDS1", y = "NMDS2",
                      color = "cluster_code",
                      size = 0.5,
                      mean.point = TRUE,
                      mean.point.size = 10,
                      ####star.plot = TRUE,
                      ####star.plot.lty = 3,
                      ellipse = TRUE,
                      ellipse.level = 0.95,
                      ellipse.type = "confidence",
                      ellipse.border.remove = TRUE,
                      ellipse.alpha = 0.2) +
      #{{use geom_point and scale shape for diversification and cluster}}
      geom_point(aes(color = cluster_code), size = 0.5) +
      theme_bw(base_size = 16) +
      theme(legend.title = element_text(vjust = 0.5, hjust = 0, size = 12),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.text = element_text(hjust = 0, size = 10),
            legend.background = element_blank(),
            aspect.ratio = 1) +
      ggtitle("") +
      ylim(y_min,y_max) +
      xlim(x_min,x_max) +
      annotate("text", x = x_min, y = y_max, hjust = 0, label = stress_annotation) +
      annotate("text", x = x_min, y = y_min + 0.04, hjust = 0, vjust = 0.5, label = homova_annotation) +
      annotate("text", x = x_min, y = y_min + 0.01, hjust = 0, vjust = 0.5, label = permanova_annotation) +
      scale_colour_manual(name = "",
                          limits = flevels_clu,
                          labels = c("N\n(dry)", "M\n(moderatly dry)",
                                     "B\n(variable)", "R\n(wet)",
                                     "F\n(very wet)"),
                          values = colpalette_clu,
                          guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
      scale_fill_manual(limits = flevels_clu,
                        labels = c("N\n(dry)", "M\n(moderatly dry)",
                                   "B\n(variable)", "R\n(wet)",
                                   "F\n(very wet)"),
                        values = colpalette_clu) +
      guides(fill = "none") +
      coord_equal()
  }
  
  #{add marginal densities}}
  {
    #{{along x axis}} 
    xdens <- axis_canvas(nmds, axis = "x") +
      geom_density(data = data.scores, aes(x = NMDS1, fill = diversification),
                   alpha = 0.7, size = 0.2) +
      scale_fill_manual(limits = flevels_div,
                        values = colpalette_div)
    #{{along y axis}}
    #{{Need to set coord_flip = TRUE, if you plan to use coord_flip()}}
    ydens <- axis_canvas(nmds, axis = "y", coord_flip = TRUE) +
      geom_density(data = data.scores, aes(x = NMDS2, fill = diversification),
                   alpha = 0.7, size = 0.2) +
      scale_fill_manual(limits = flevels_div ,
                        values = colpalette_div) +
      coord_flip()
    
    nmds <- insert_xaxis_grob(nmds, xdens, grid::unit(.2, "null"), position = "top")
    nmds <- insert_yaxis_grob(nmds, ydens, grid::unit(.2, "null"), position = "right")
  }
  
  #{{assembly of density plots to nmds plot}}   
  {
    ggdraw(nmds)
    
    width <- 12
    height <- 12
    
    ggsave(file = output, width = width , height = height)
  }
 
}

## ---- 19. VENN DIAGRAM ANALYSIS ----
## -------- 19.1 bacteria --------

#{{create the categories to compare}}
OTUsCS_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Ceratonia_siliqua") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsCS_soil <- OTUsCS_soil[,-(which(colSums(OTUsCS_soil) == 0))]
OTUsCS_soil <- OTUsCS_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsCS_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Ceratonia_siliqua") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsCS_root <- OTUsCS_root[,-(which(colSums(OTUsCS_root) == 0))]
OTUsCS_root <- OTUsCS_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsPL_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Pistacia_lentiscus") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsPL_soil <- OTUsPL_soil[,-(which(colSums(OTUsPL_soil) == 0))]
OTUsPL_soil <- OTUsPL_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsPL_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Pistacia_lentiscus") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsPL_root <- OTUsPL_root[,-(which(colSums(OTUsPL_root) == 0))]
OTUsPL_root <- OTUsPL_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsRM_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Retama_monosperma") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsRM_soil <- OTUsRM_soil[,-(which(colSums(OTUsRM_soil) == 0))]
OTUsRM_soil <- OTUsRM_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsRM_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Retama_monosperma") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsRM_root <- OTUsRM_root[,-(which(colSums(OTUsRM_root) == 0))]
OTUsRM_root <- OTUsRM_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsGA_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Globularia_alypum") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsGA_soil <- OTUsGA_soil[,-(which(colSums(OTUsGA_soil) == 0))]
OTUsGA_soil <- OTUsGA_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsGA_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Globularia_alypum") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsGA_root <- OTUsGA_root[,-(which(colSums(OTUsGA_root) == 0))]
OTUsGA_root <- OTUsGA_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsHV_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Hordeum_vulgare") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsHV_soil <- OTUsHV_soil[,-(which(colSums(OTUsHV_soil) == 0))]
OTUsHV_soil <- OTUsHV_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsHV_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Hordeum_vulgare") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsHV_root <- OTUsHV_root[,-(which(colSums(OTUsHV_root) == 0))]
OTUsHV_root <- OTUsHV_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsVB_soil <- bactOTUs_rarefied_factor_t %>%
  filter (Type=="Soils", Plant_type=="Vicia_faba") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsVB_soil <- OTUsVB_soil[,-(which(colSums(OTUsVB_soil) == 0))]
OTUsVB_soil <- OTUsVB_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsVB_root <- bactOTUs_rarefied_factor_t %>%
  filter (Type=="Roots", Water_status=="no_stress", Plant_type=="Vicia_faba") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsVB_root <- OTUsVB_root[,-(which(colSums(OTUsVB_root) == 0))]
OTUsVB_root <- OTUsVB_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

#{{set up of venn plot with ggvenn}}
#{{Colorblind-friendly palette with grey}}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#{{set up the list of categories}}
cat <- list("C. siliqua roots (CS)" = OTUsCS_root$OTU, "C. siliqua soils" = OTUsCS_soil$OTU,
            "C. siliqua roots (PL)" = OTUsPL_root$OTU, "P. lentiscus soils" = OTUsPL_soil$OTU,
            "C. siliqua roots (RM)" = OTUsRM_root$OTU, "R. monosperma soils" = OTUsRM_soil$OTU,
            "C. siliqua roots (GA)" = OTUsGA_root$OTU, "G. alypum soils" = OTUsGA_soil$OTU,
            "C. siliqua roots (HV)" = OTUsHV_root$OTU, "H. vulgare soils" = OTUsHV_soil$OTU,
            "C. siliqua roots (VB)" = OTUsVB_root$OTU, "V. faba soils" = OTUsVB_soil$OTU) 

#{{select the file name}}
output <- "venn_bact_soil_roots_no_stress.pdf"

venn1 <- ggvenn(cat,c("C. siliqua soils", "C. siliqua roots (CS)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn2 <- ggvenn(cat,c("P. lentiscus soils", "C. siliqua roots (PL)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn3 <- ggvenn(cat,c("R. monosperma soils", "C. siliqua roots (RM)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn4 <- ggvenn(cat,c("G. alypum soils", "C. siliqua roots (GA)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn5 <- ggvenn(cat,c("H. vulgare soils", "C. siliqua roots (HV)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn6 <- ggvenn(cat,c("V. faba soils", "C. siliqua roots (VB)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

#{{assembly of venn plots}}
((venn1 | venn2) / (venn3 | venn4) / (venn5 | venn6)) +
  plot_annotation(tag_levels = 'A')


width <- 10
height <- 15

ggsave(file = output, width = width , height = height)

## -------- 19.2 fungi ----
#{{create the categories to compare}}
OTUsCS_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Ceratonia_siliqua") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsCS_soil <- OTUsCS_soil[,-(which(colSums(OTUsCS_soil) == 0))]
OTUsCS_soil <- OTUsCS_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsCS_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Ceratonia_siliqua") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsCS_root <- OTUsCS_root[,-(which(colSums(OTUsCS_root) == 0))]
OTUsCS_root <- OTUsCS_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsPL_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Pistacia_lentiscus") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsPL_soil <- OTUsPL_soil[,-(which(colSums(OTUsPL_soil) == 0))]
OTUsPL_soil <- OTUsPL_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsPL_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Pistacia_lentiscus") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsPL_root <- OTUsPL_root[,-(which(colSums(OTUsPL_root) == 0))]
OTUsPL_root <- OTUsPL_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsRM_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Retama_monosperma") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsRM_soil <- OTUsRM_soil[,-(which(colSums(OTUsRM_soil) == 0))]
OTUsRM_soil <- OTUsRM_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsRM_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Retama_monosperma") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsRM_root <- OTUsRM_root[,-(which(colSums(OTUsRM_root) == 0))]
OTUsRM_root <- OTUsRM_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsGA_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Globularia_alypum") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsGA_soil <- OTUsGA_soil[,-(which(colSums(OTUsGA_soil) == 0))]
OTUsGA_soil <- OTUsGA_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsGA_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Globularia_alypum") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsGA_root <- OTUsGA_root[,-(which(colSums(OTUsGA_root) == 0))]
OTUsGA_root <- OTUsGA_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsHV_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Soils", Plant_type == "Hordeum_vulgare") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsHV_soil <- OTUsHV_soil[,-(which(colSums(OTUsHV_soil) == 0))]
OTUsHV_soil <- OTUsHV_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsHV_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type == "Roots", Water_status == "no_stress", Plant_type == "Hordeum_vulgare") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsHV_root <- OTUsHV_root[,-(which(colSums(OTUsHV_root) == 0))]
OTUsHV_root <- OTUsHV_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsVB_soil <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type=="Soils", Plant_type=="Vicia_faba") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsVB_soil <- OTUsVB_soil[,-(which(colSums(OTUsVB_soil) == 0))]
OTUsVB_soil <- OTUsVB_soil %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

OTUsVB_root <- ITS2OTUs_rarefied_factor_t %>%
  filter (Type=="Roots", Water_status=="no_stress", Plant_type=="Vicia_faba") %>%
  select(-soil_code,-Station_code,-Country,-Type,-Plant_type,-Water_status) %>%
  droplevels()
OTUsVB_root <- OTUsVB_root[,-(which(colSums(OTUsVB_root) == 0))]
OTUsVB_root <- OTUsVB_root %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>%
  select (OTU)

#{{set up of Venn plot with ggvenn}}
#{{Colorblind-friendly palette with grey}}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#{{set up the list of categories}}
cat <- list("C. siliqua roots (CS)" = OTUsCS_root$OTU, "C. siliqua soils" = OTUsCS_soil$OTU,
            "C. siliqua roots (PL)" = OTUsPL_root$OTU, "P. lentiscus soils" = OTUsPL_soil$OTU,
            "C. siliqua roots (RM)" = OTUsRM_root$OTU, "R. monosperma soils" = OTUsRM_soil$OTU,
            "C. siliqua roots (GA)" = OTUsGA_root$OTU, "G. alypum soils" = OTUsGA_soil$OTU,
            "C. siliqua roots (HV)" = OTUsHV_root$OTU, "H. vulgare soils" = OTUsHV_soil$OTU,
            "C. siliqua roots (VB)" = OTUsVB_root$OTU, "V. faba soils" = OTUsVB_soil$OTU) 

#{{select the file name}}
output <- "venn_fung_soil_roots_no_stress.pdf"

venn1 <- ggvenn(cat,c("C. siliqua soils", "C. siliqua roots (CS)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn2 <- ggvenn(cat,c("P. lentiscus soils", "C. siliqua roots (PL)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn3 <- ggvenn(cat,c("R. monosperma soils", "C. siliqua roots (RM)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn4 <- ggvenn(cat,c("G. alypum soils", "C. siliqua roots (GA)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn5 <- ggvenn(cat,c("H. vulgare soils", "C. siliqua roots (HV)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

venn6 <- ggvenn(cat,c("V. faba soils", "C. siliqua roots (VB)"), fill_color = cbPalette[c(2,4)], set_name_size = 4, text_size = 4)

#{{assembly of venn plots}}
((venn1 | venn2) / (venn3 | venn4) / (venn5 | venn6)) +
  plot_annotation(tag_levels = 'A')


width <- 10
height <- 15

ggsave(file = output, width = width , height = height)

## ---- 20. STATIC UPSET ANALYSIS ----
## -------- 20.1 bacteria ------

#{{select data}}
upcom <- bact_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{select only OTUs}}
  select(organ_health, where(is.numeric)) %>%
  group_by(organ_health) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "organ_health") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  pivot_longer(!OTU, names_to = "organ_health", values_to = "reads") %>%  
  group_by(OTU) %>%
  mutate(sum_group = sum(reads)) %>%
  ungroup() %>%
  pivot_wider(names_from = organ_health, values_from = reads) %>%
  filter(sum_group != 0)  %>%
  select(-sum_group) %>%
  mutate_if(is.numeric, ~as.character(.)) %>%
  column_to_rownames(var = "OTU")

#{{covert numeric to binary}}}}
upcom[upcom != "0"] <- "1"

#{{convert character value in numeric for UpsetR}}}}
upcom <- upcom %>%
  mutate_if(is.character, ~as.numeric(.)) %>%
  rownames_to_column(var = "OTU")
upcom$OTU <- as.factor(upcom$OTU)

#{{create plot with UpSetR package}}}}
upset(upcom, sets = c("leaf_disease", "leaf_healthy", "root_disease", "root_healthy",
                      "stalk_disease", "stalk_healthy"),
      mb.ratio = c(0.55, 0.45), point.size = 3.5, line.size = 2, 
      mainbar.y.label = "number of OTUs\n(intersection)", 
      sets.x.label = "number of OTUs\n(Sugarcane organs and health status)", 
      text.scale = c(1.5, 1.5, 1.5, 1, 1, 1), 
      empty.intersections = "on", order.by = "freq")

## -------- 20.2 fungi ------

#{{select data}}
upcom <- ITS2_OTUs_rarefied_t %>%
  #{{merge with metadata}}
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  unite(col = organ_health , type, health_status, 
        sep = "_", remove = FALSE, na.rm = FALSE) %>%
  #{{select only OTUs}}
  select(organ_health, where(is.numeric)) %>%
  group_by(organ_health) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "organ_health") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  pivot_longer(!OTU, names_to = "organ_health", values_to = "reads") %>%  
  group_by(OTU) %>%
  mutate(sum_group = sum(reads)) %>%
  ungroup() %>%
  pivot_wider(names_from = organ_health, values_from = reads) %>%
  filter(sum_group != 0)  %>%
  select(-sum_group) %>%
  mutate_if(is.numeric, ~as.character(.)) %>%
  column_to_rownames(var = "OTU")

#{{covert numeric to binary}}}}
upcom[upcom != "0"] <- "1"

#{{convert character value in numeric for UpsetR}}}}
upcom <- upcom %>%
  mutate_if(is.character, ~as.numeric(.)) %>%
  rownames_to_column(var = "OTU")
upcom$OTU <- as.factor(upcom$OTU)

#{{create plot with UpSetR package}}}}
upset(upcom, sets = c("leaf_disease", "leaf_healthy", "root_disease", "root_healthy",
                      "stalk_disease", "stalk_healthy"),
      mb.ratio = c(0.55, 0.45), point.size = 3.5, line.size = 2, 
      mainbar.y.label = "number of OTUs\n(intersection)", 
      sets.x.label = "number of OTUs\n(Sugarcane organs and health status)", 
      text.scale = c(1.5, 1.5, 1.5, 1, 1, 1), 
      empty.intersections = "on", order.by = "freq")

## ---- 21. MICROBIAL INDICATORS ----
## ------ 21.1. bacteria ------
## ---------- 21.1.1. health ----------
## -------------- 21.1.1.1 OTU level ----------
## ------------------ 21.1.1.1.1. data --------------

#{{select data}}
tmp <- bact_OTUs_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.1.1.1 multi-group analysis --------------------
set.seed(seed)
multigroup_bact <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_bact$sign) %>%
  add_column(OTU = rownames(multigroup_bact$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_bact$A) %>%
  add_column(OTU = rownames(multigroup_bact$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_bact$B) %>%
  add_column(OTU = rownames(multigroup_bact$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_bact_otu_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(bactOTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_bact_otu_roots_nostress_plant_p.adj_sig  <- filter(indic_bact_otu_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_bact_otu_roots_nostress_plant_p.adj_sig,path <- "indic_bact_otu_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.1.1.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_bact <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_bact)
multispe_bact


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_bact <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_bact)
multispe_bact

## -------------- 21.1.1.2 taxonomic level ----------
## ------------------ 21.1.1.2.3. family ------------------
#{{select data}}
tmp <- bactFamily_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.2.3.1. multi-group analysis --------------------
set.seed(seed)
multigroup_bact <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_bact$sign) %>%
  add_column(OTU = rownames(multigroup_bact$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_bact$A) %>%
  add_column(OTU = rownames(multigroup_bact$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_bact$B) %>%
  add_column(OTU = rownames(multigroup_bact$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_bact_family_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(bactOTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_bact_family_roots_nostress_plant_p.adj_sig  <- filter(indic_bact_family_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_bact_family_roots_nostress_plant_p.adj_sig,path <- "indic_bact_family_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.2.3.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_bact <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_bact)
multispe_bact


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_bact <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_bact)
multispe_bact

## ------------------ 21.1.1.2.4. genus ------------------
#{{select data}}
tmp <- bactGenus_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.2.4.1. multi-group analysis --------------------
set.seed(seed)
multigroup_bact <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_bact$sign) %>%
  add_column(OTU = rownames(multigroup_bact$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_bact$A) %>%
  add_column(OTU = rownames(multigroup_bact$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_bact$B) %>%
  add_column(OTU = rownames(multigroup_bact$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_bact_genus_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(bactOTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_bact_genus_roots_nostress_plant_p.adj_sig  <- filter(indic_bact_genus_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_bact_genus_roots_nostress_plant_p.adj_sig,path <- "indic_bact_genus_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.2.4.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_bact <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_bact)
multispe_bact


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_bact <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_bact, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_bact)
multispe_bact

## ------ 21.2. fungi (roots) ----
## ---------- 21.1.1. water treatment ----------
## -------------- 21.1.1.1 OTU level ----------
## ------------------ 21.1.1.1.1. data --------------
#{{select data}}
tmp <- ITS2OTUs_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.1.1.1 multi-group analysis --------------------
set.seed(seed)
multigroup_fung <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_fung$sign) %>%
  add_column(OTU = rownames(multigroup_fung$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_fung$A) %>%
  add_column(OTU = rownames(multigroup_fung$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_fung$B) %>%
  add_column(OTU = rownames(multigroup_fung$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_fung_otu_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(ITS2OTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_fung_otu_roots_nostress_plant_p.adj_sig  <- filter(indic_fung_otu_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_fung_otu_roots_nostress_plant_p.adj_sig,path <- "indic_fung_otu_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.1.1.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_fung <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_fung)
multispe_fung


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_fung <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_fung)
multispe_fung

## -------------- 21.1.1.2 taxonomic level ----------
## ------------------ 21.1.1.2.3. family ------------------
#{{select data}}
tmp <- ITS2Family_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.2.3.1. multi-group analysis --------------------
set.seed(seed)
multigroup_fung <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_fung$sign) %>%
  add_column(OTU = rownames(multigroup_fung$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_fung$A) %>%
  add_column(OTU = rownames(multigroup_fung$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_fung$B) %>%
  add_column(OTU = rownames(multigroup_fung$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_fung_family_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(ITS2OTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_fung_family_roots_nostress_plant_p.adj_sig  <- filter(indic_fung_family_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_fung_family_roots_nostress_plant_p.adj_sig,path <- "indic_fung_family_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.2.3.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_fung <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_fung)
multispe_fung


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_fung <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_fung)
multispe_fung

## ------------------ 21.1.1.2.4. genus ------------------
#{{select data}}
tmp <- ITS2Genus_rarefied_factor_t %>%
  filter(Type == "Roots", Water_status == "no_stress", Plant_type != "Control") %>%
  droplevels()

#{{select only factor}}
factor <- tmp$Plant_type
#{{select only taxa}}
taxa <- select(tmp, -soil_code, -Station_code, -Country, -Type, -Plant_type, -Water_status)

#{{select abundance data tranformation}}
transtaxa <- log1p(taxa)  #{{sqrt(x) or log1p(x) or decostand(x,"hellinger", MARGIN = 1) or wisconsin(com) or decostand(log1p(com),"normalize")}}

## ------------------------ 21.1.1.2.4.1. multi-group analysis --------------------
set.seed(seed)
multigroup_fung <- multipatt(transtaxa, factor, duleg = TRUE, func = "IndVal.g",control = how(nperm=999))
summary(multigroup_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.5, Bt = 0.5)

#{{extract table of stats + p.value correction}}
indic.sign <-  as_tibble(multigroup_fung$sign) %>%
  add_column(OTU = rownames(multigroup_fung$sign)) %>%
  mutate(p.adjust.bh = p.adjust(p.value, method = "BH")) %>%
  select(OTU,everything())
indic.A <- as_tibble(multigroup_fung$A) %>%
  add_column(OTU = rownames(multigroup_fung$A)) %>%
  rename(A_Ceratonia_siliqua = Ceratonia_siliqua,
         A_Globularia_alypum = Globularia_alypum,
         A_Hordeum_vulgare = Hordeum_vulgare,
         A_Pistacia_lentiscus = Pistacia_lentiscus,
         A_Retama_monosperma = Retama_monosperma,
         A_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())
indic.B <- as.tibble(multigroup_fung$B) %>%
  add_column(OTU = rownames(multigroup_fung$B)) %>%
  rename(B_Ceratonia_siliqua = Ceratonia_siliqua,
         B_Globularia_alypum = Globularia_alypum,
         B_Hordeum_vulgare = Hordeum_vulgare,
         B_Pistacia_lentiscus = Pistacia_lentiscus,
         B_Retama_monosperma = Retama_monosperma,
         B_Vicia_faba = Vicia_faba) %>%
  select(OTU,everything())

#{{format multipatt results + taxonomic assignation}}
indic_fung_genus_roots_nostress_plant_p.adj  <- indic.A %>%
  left_join(indic.B, by = "OTU") %>%
  left_join(indic.sign, by = "OTU") %>%
  left_join(select(ITS2OTUs_tax,OTU,kingdom,phylum,class,order,family,genus,species), by = "OTU")

#{{select indicator only with p.adjust <=0.05}}
indic_fung_genus_roots_nostress_plant_p.adj_sig  <- filter(indic_fung_genus_roots_nostress_plant_p.adj, p.adjust.bh <= 0.05)

write_delim(indic_fung_genus_roots_nostress_plant_p.adj_sig,path <- "indic_fung_genus_roots_nostress_plant_p.adj_sig.txt",
            delim="\t", col_names=TRUE)

## ------------------------ 21.1.1.2.4.2. multi-species analysis --------------------

#{{select "Ceratonia_siliqua", "Globularia_alypum", "Hordeum_vulgare", "Pistacia_lentiscus", "Retama_monosperma", "Vicia_faba")}}
multispe_fung <- indicators(transtaxa, factor, func = "IndVal.g", max.order = 5,
                            group = "Hordeum_vulgare", At = 0.7, Bt = 0.7)
summary(multispe_fung)
multispe_fung


#{{determine group of species + multigroup analysis}}
combspe <- combinespecies(transtaxa, max.order = 3)$XC
multigroupspe_fung <- multipatt(transtaxa, factor, func = "IndVal.g",control = how(nperm=999))
summary(multigroupspe_fung, indvalcomp = TRUE, alpha = 0.05, At = 0.7, Bt = 0.7)

summary(multispe_fung)
multispe_fung

## ---- 22. Differential compositional analysis {ALDEX2}---- 
## ---- 23. Differential abundance analysis {edgeR}----
#{{for more information , see Statistical Analysis of Microbiome Data with R}}
## ------ 23.1. bacteria ------
## ---------- 23.1.1. plant organ & disease ----------
## -------------- 23.1.2.1 OTU level --------------
## ------------------ 23.1.2.1.1. data ------------------

#{{select data}}
data <- bact_OTUs_rarefied_t %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  filter(type == "stalk") %>%
  droplevels()

#{{select only taxa}}
taxa <- data %>%
  column_to_rownames(var = "soil_code") %>%
  select(-station_code, -type, -country, -health_status) %>%
  t() 

#{{select only factor}}
factor <- data$health_status

#{{Generate edgeR object}}
transtaxa <- DGEList(counts = taxa, group = factor)

#{{pre-filtering of data based on cpm filter}}
#{{pre-filtering increase the sensitivity a pre-filtering is needed}}
keep <- rowSums(cpm(transtaxa) > 0.5) >= 3
transtaxa_filter <- transtaxa[keep, ]

#{{normalization of pre-filtered data}}
transtaxa_norm <- calcNormFactors(transtaxa_filter)

#{{Fit a generalized linear model to estimate the dispersion}}
design <- model.matrix(~factor) #{{with intercept}}
design_nointercep <- model.matrix(~0 + factor, data = transtaxa_norm$samples) #{{with no intercept}}
transtaxa_disp <- estimateDisp(transtaxa_norm, design_nointercep, robust = TRUE)
transtaxa_fit <- glmQLFit(transtaxa_disp, design_nointercep, robust = TRUE)

## ----------------------- 23.1.2.1.1.1 perform comparison tests [GLM model with quasi-likelihood F-test] ------------------
#{{likelihood ratio test can be useful glmLRT() datasets with no replicates}} 
transtaxa_diff <- glmQLFTest(transtaxa_fit, contrast = c(-1,1))
#{{visualize top most significant taxa wih FDR}}
topTags(transtaxa_diff)

## ----------------------- 22.1.1.2.1.2. Format tests ------------------

#{{change name}}
output <- "Diff_bact_otu_stalk_health_p.adj_sig.txt"

#{{select taxa only with FDR p.adjust <=0.05 for all tests}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "BH")) %>%
  select(OTU,everything()) %>%
  filter(p.adjust.bh <= 0.05) %>%
  ####type_convert() %>%
  #{{add taxonomy}}
  left_join(select(bact_tax, OTU, domain, phylum, class, order, family, genus, species), by = "OTU")

write_delim(data, path <- output, delim="\t", col_names = TRUE)

## ----------------------- 23.2.2.2.1.3. DI and DSI index -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "fdr")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(OTU,everything()) %>%
  #{{add taxonomy}}
  left_join(select(bact_tax, OTU, domain, phylum, class, order, family, genus, species), by = "OTU")

#{{rename categories}}
data$differential <- data$differential %>%
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{extract values}}
index <- tibble(differential = c("ns", "depleted", "enriched")) %>%
  left_join(count(data, differential), by = "differential") %>%
  mutate(n = replace_na(n, 0))

ns <- index %>% filter(differential =="ns") %>% pull(n)
enriched <- index %>% filter(differential =="enriched") %>% pull(n)
depleted <- index %>% filter(differential =="depleted") %>% pull(n) 

#{{estimate index}}
DI_bact_health <- depleted / enriched
DSI_bact_health <- (enriched + depleted) / (enriched + depleted + ns)

## ----------------------- 23.2.2.2.1.4. Volcano plot -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "BH")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(OTU,everything()) %>%
  left_join(select(bact_tax, OTU, domain, phylum, class, order, family, genus, species), by = "OTU")

#{{rename categories}}
data$differential <- data$differential %>% 
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{rearrange the color according the palette}}
colpalette_health <- c("#46732EFF", "#C80813FF", "#8A9197FF")

flevels_health <- c("enriched","depleted","ns")

#{{plot}}

output <- "Volcano_bact_otu_stalk_health_p.adj_sig.pdf"

volc1 <- data %>%
  ggplot(aes(x = logCPM, y = logFC)) +
  geom_point(aes(color = differential), size = 6) +
  scale_color_manual(name = "", 
                     limits = flevels_health,
                     values = colpalette_health) +
  ggtitle ("Healthy vs Diseases") +
  theme_grey(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom") +
  geom_text_repel(
    data = subset(data, p.adjust.bh <= 0.05 & logCPM > 10 & (logFC > 3 | logFC < -3)),
    aes(label = species),
    size = 2,
    segment.size = 0.1,
    box.padding = unit(3, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 50)) +
  ylim(-15,15) +
  xlim(10,20) +
  xlab ("Log2 (Count per million)") +
  ylab ("Log2 (Fold change)")

width <- 12
height <- 9

ggsave(file = output, width = width , height = height)

## -------------- 23.1.2.1 taxonomic level --------------
## ------------------ 23.1.2.1.1. data ------------------

#{{select taxonomic level and organ}}
data <- bact_OTUs_rarefied_t
tax <- "genus"
organ <- "leaf"
fact_o <- "type"
fact_h <- "health_status"

data <- bact_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  #{{select the taxonomic rank to analyse}}
  left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
  select(-OTU) %>%
  #{{group by the taxonomic rank}}
  group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name
  summarise_all(sum) %>%
  column_to_rownames(var = tax) %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  #{{select the treatments to analyse}}
  filter(type == organ) %>%
  droplevels()

#{{select only taxa}}
taxa <- data %>%
  column_to_rownames(var = "soil_code") %>%
  select(-station_code, -type, -country, -health_status) %>%
  t() 

#{{select only factor}}
factor <- data$health_status

#{{Generate edgeR object}}
transtaxa <- DGEList(counts = taxa, group = factor)

#{{pre-filtering of data based on cpm filter}}
#{{pre-filtering increase the sensitivity a pre-filtering is needed}}
keep <- rowSums(cpm(transtaxa) > 0.5) >= 3
transtaxa_filter <- transtaxa[keep, ]

#{{normalization of pre-filtered data}}
transtaxa_norm <- calcNormFactors(transtaxa_filter)

#{{Fit a generalized linear model to estimate the dispersion}}
design <- model.matrix(~factor) #{{with intercept}}
design_nointercep <- model.matrix(~0 + factor, data = transtaxa_norm$samples) #{{with no intercept}}
transtaxa_disp <- estimateDisp(transtaxa_norm, design_nointercep, robust = TRUE)
transtaxa_fit <- glmQLFit(transtaxa_disp, design_nointercep, robust = TRUE)

## ----------------------- 23.1.2.1.1.1 perform comparison tests [GLM model with quasi-likelihood F-test] ------------------
#{{likelihood ratio test can be useful glmLRT() datasets with no replicates}} 
transtaxa_diff <- glmQLFTest(transtaxa_fit, contrast = c(-1,1))
#{{visualize top most significant taxa wih FDR}}
topTags(transtaxa_diff)

## ----------------------- 22.1.1.2.1.2. Format tests ------------------

#{{change name}}
output <- "Diff_bact_genus_leaf_health_p.adj_sig.txt"

#{{select taxa only with FDR p.adjust <=0.05 for all tests}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(genus = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "fdr")) %>%
  select(genus,everything()) %>%
  filter(p.adjust.bh <= 0.05)
  
write_delim(data, path <- output, delim="\t", col_names = TRUE)

## ----------------------- 23.2.2.2.1.3. DI and DSI index -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(genus = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "fdr")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(genus,everything())
  
#{{rename categories}}
data$differential <- data$differential %>%
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{extract values}}
index <- tibble(differential = c("ns", "depleted", "enriched")) %>%
  left_join(count(data, differential), by = "differential") %>%
  mutate(n = replace_na(n, 0))

ns <- index %>% filter(differential =="ns") %>% pull(n)
enriched <- index %>% filter(differential =="enriched") %>% pull(n)
depleted <- index %>% filter(differential =="depleted") %>% pull(n) 

#{{estimate index}}
DI_bact_health <- depleted / enriched
DSI_bact_health <- (enriched + depleted) / (enriched + depleted + ns)

## ----------------------- 23.2.2.2.1.4. Volcano plot -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(genus = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "BH")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(genus,everything()) 
  
#{{rename categories}}
data$differential <- data$differential %>% 
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{rearrange the color according the palette}}
colpalette_health <- c("#46732EFF", "#C80813FF", "#8A9197FF")

flevels_health <- c("enriched","depleted","ns")

#{{plot}}

output <- "Volcano_bact_genus_leaf_health_p.adj_sig.pdf"

volc1 <- data %>%
  ggplot(aes(x = logCPM, y = logFC)) +
  geom_point(aes(color = differential), size = 6) +
  scale_color_manual(name = "", 
                     limits = flevels_health,
                     values = colpalette_health) +
  ggtitle ("Healthy vs Diseases") +
  theme_grey(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom") +
  geom_text_repel(
    data = subset(data, p.adjust.bh <= 0.05 & logCPM > 10 & (logFC > 3 | logFC < -3)),
    aes(label = genus),
    size = 2,
    segment.size = 0.1,
    box.padding = unit(3, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 50)) +
  ylim(-15,15) +
  xlim(10,20) +
  xlab ("Log2 (Count per million)") +
  ylab ("Log2 (Fold change)")

width <- 12
height <- 9

ggsave(file = output, width = width , height = height)


## ------ 23.1. fungi ------
## ---------- 23.1.1. diversification ----------
## -------------- 23.1.2.1 OTU level --------------
## ------------------ 23.1.2.1.1. data ------------------

#{{select data}}
data <- ITS2_OTUs_rarefied_t %>%
  rownames_to_column(var = "soil_code") %>%
  left_join(metadata, by = "soil_code") %>%
  filter(diversification == "RCM" | diversification == "RCTBA") %>% 
  droplevels()

#{{select only taxa}}
taxa <- data %>%
  column_to_rownames(var = "soil_code") %>%
  select(-replicate, -cluster_code, -farmer_name, -sampling_date,
         -altitude, -GPS_lat, -GPS_lon, -diversification, -details) %>%
  t() 

#{{select only factor}}
factor <- data$diversification

#{{Generate edgeR object}}
transtaxa <- DGEList(counts = taxa, group = factor)

#{{pre-filtering of data based on cpm filter}}
#{{pre-filtering increase the sensitivity a pre-filtering is needed}}
keep <- rowSums(cpm(transtaxa) > 0.5) >= 3
transtaxa_filter <- transtaxa[keep, ]

#{{normalization of pre-filtered data}}
transtaxa_norm <- calcNormFactors(transtaxa_filter)

#{{Fit a generalized linear model to estimate the dispersion}}
design <- model.matrix(~factor) #{{with intercept}}
design_nointercep <- model.matrix(~0 + factor, data = transtaxa_norm$samples) #{{with no intercept}}
transtaxa_disp <- estimateDisp(transtaxa_norm, design_nointercep, robust = TRUE)
transtaxa_fit <- glmQLFit(transtaxa_disp, design_nointercep, robust = TRUE)

## ----------------------- 23.1.2.1.1.1 perform comparison tests [GLM model with quasi-likelihood F-test] ------------------
#{{likelihood ratio test can be useful glmLRT() datasets with no replicates}} 
transtaxa_diff <- glmQLFTest(transtaxa_fit, contrast = c(-1,1))
#{{visualize top most significant taxa wih FDR}}
topTags(transtaxa_diff)

## ----------------------- 22.1.1.2.1.2. Format tests ------------------

#{{change name}}
output <- "Diff_fung_RCM_RCB_p.adj_sig.txt"

#{{select taxa only with FDR p.adjust <=0.05 for all tests}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "BH")) %>%
  select(OTU,everything()) %>%
  filter(p.adjust.bh <= 0.05) %>%
  ####type_convert() %>%
  #{{add taxonomy}}
  left_join(select(ITS2_tax, OTU, kingdom, phylum, class, order, family, genus, species), by = "OTU")

write_delim(data, path <- output, delim="\t", col_names = TRUE)

## ----------------------- 23.2.2.2.1.3. DI and DSI index -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "BH")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(OTU,everything()) %>%
  #{{add taxonomy}}
  left_join(select(ITS2_tax, OTU, kingdom, phylum, class, order, family, genus, species), by = "OTU")

#{{rename categories}}
data$differential <- data$differential %>%
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{extract values}}
index <- tibble(differential = c("ns", "depleted", "enriched")) %>%
  left_join(count(data, differential), by = "differential") %>%
  mutate(n = replace_na(n, 0))

ns <- index %>% filter(differential =="ns") %>% pull(n)
enriched <- index %>% filter(differential =="enriched") %>% pull(n)
depleted <- index %>% filter(differential =="depleted") %>% pull(n) 

#{{estimate index}}
DI_fung_health <- depleted / enriched
DSI_fung_health <- (enriched + depleted) / (enriched + depleted + ns)

## ----------------------- 23.2.2.2.1.4. Volcano plot -----------------------

#{{select all taxa}}
data <- transtaxa_diff$table %>%
  as_tibble() %>%
  add_column(OTU = rownames(transtaxa_diff$table)) %>%
  mutate(p.adjust.bh = p.adjust(PValue, method = "fdr")) %>%
  mutate(direction = ifelse (logFC > 0, "enriched","depleted")) %>%
  mutate(significance = ifelse (p.adjust.bh <= 0.05, "sig","nosig")) %>%
  unite(col = "differential", direction, significance, sep = "_", remove = FALSE, na.rm = FALSE) %>%
  select(OTU,everything()) %>%
  left_join(select(ITS2_tax, OTU, kingdom, phylum, class, order, family, genus, species), by = "OTU")

#{{rename categories}}
data$differential <- data$differential %>% 
  recode(enriched_sig = "enriched", depleted_sig = "depleted", enriched_nosig = "ns", depleted_nosig = "ns")

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))
#{{rearrange the color according the palette}}
colpalette_div <- c("#46732EFF", "#FED439FF", "#8A9197FF")

flevels_div <- c("enriched","depleted","ns")

#{{plot}}

output <- "Volcano_fung_RCM_RCTBA_p.adj_sig.pdf"
{
volc1 <- data %>%
  ggplot(aes(x = logCPM, y = logFC)) +
  geom_point(aes(color = differential), size = 6) +
  scale_color_manual(name = "", 
                     limits = flevels_div,
                     values = colpalette_div) +
  ggtitle ("Monoculture vs Shade tree + Banana + others") +
  theme_grey(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "bottom") +
  geom_text_repel(
    data = subset(data, p.adjust.bh <= 0.05 & logCPM > 5 & (logFC > 2 | logFC < -2)),
    aes(label = species),
    size = 2,
    segment.size = 0.1,
    box.padding = unit(3, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = getOption("ggrepel.max.overlaps", default = 50)) +
  ylim(-15,15) +
  xlim(5,20) +
  xlab ("Log2 (Count per million)") +
  ylab ("Log2 (Fold change)")

width <- 12
height <- 9

ggsave(file = output, width = width , height = height)
}
