library(tidyverse)

data_18S <- "data/18S_AMf_PR2/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
amf_taxonomic_levels <- c("kingdom", "domain", "phylum", "class", "order", "family",
                          "genus", "species")

metadata <- read_tsv("data/metadata.txt")

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
  filter(str_detect(taxonomy,"Glomeromycotina")|str_detect(taxonomy,"Endogonales")) %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB_", ignore.case = FALSE)) %>%
  droplevels()

taxonomy <- data_18S %>%
  #{{import data}}
  read_tsv() %>%
  filter(str_detect(taxonomy,"Glomeromycotina")|str_detect(taxonomy,"Endogonales")) %>%
  select(OTU, taxonomy) %>%
  mutate(taxonomy = str_replace_all(taxonomy, ":","_"),
         taxonomy = str_replace_all(taxonomy, "[|]",";"))


d <- amf_OTUs %>% 
  replace(amf_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  pivot_longer(-OTU, names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads))

amf_OTUs_decont <- amf_OTUs %>%
  replace(amf_OTUs == 0, NA) %>%
  pivot_longer(-c("OTU", starts_with("CONT")), names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n,-starts_with("CONT")) %>%
  pivot_wider(names_from = samples, values_from =  reads, values_fill = 0) %>%
  as_tibble() %>%
  column_to_rownames(var = "OTU") 
colSums(amf_OTUs_decont)

# rarefaction
# rarefaction
quantile(rowSums(amf_OTUs_decont))
sort(colSums(amf_OTUs_decont))
plot(1:length(colSums(amf_OTUs_decont)), log10(sort(colSums(amf_OTUs_decont))+1))
10^3
sort(colSums(amf_OTUs_decont))


amf_OTUs_decont2 <- amf_OTUs_decont %>%
  #{{sample to exclude}}
  select(where(~ sum(.x) >= 10^3))
smallest_amf <- min(colSums(amf_OTUs_decont2))

library(vegan)

rarefaction.func <- function(df,rar_sample, n_rare = 100, rarcur = T){
  if(rarcur){ # {{ if rarefaction curve needed}}
    df %>%
      select(starts_with("ROB")) %>%
      t() %>%
      vegan::rarecurve(step = 20, sample = rar_sample, col = "blue", cex = 0.6)
  }
  df %>%
    select(starts_with("ROB")) %>%
    dim -> D
  # {{rarefaction}}
  df_rarefy = matrix(0, nrow =  D[1], ncol = D[2])
  for(n in 1:n_rare){
    df %>%
      select(starts_with("ROB")) %>%
      t() %>%
      vegan::rrarefy(rar_sample) %>%
      t() + df_rarefy -> df_rarefy
    
  }
  round(df_rarefy/n_rare) -> df_rarefy
  df_rarefy %>%
    bind_cols(
      df %>%
        select(!starts_with("ROB"))
    ) -> df_rarefy
  rownames(df_rarefy) = rownames(df)
  return(df_rarefy)
}
amf_OTUs_decont2 %>% 
  mutate(across(starts_with("ROB"), ~as.integer(.x))) -> amf_OTUs_decont2
set.seed(2121349)
rarefaction.func(df = amf_OTUs_decont2, rar_sample = smallest_amf, rarcur = T) -> amf_OTUs_rare

sort(colSums(amf_OTUs_rare))

compute.hill.diversities<- function(community_matrix){
  div_id = c(richness = 0, shanon = 1, simpson = 2)
  
  community_matrix %>%
    dplyr::select(starts_with("ROB")) -> community_matrix
  
  div_id %>%
    purrr::map(\(x) hilldiv::hill_div(community_matrix, qvalue = x)) %>%
    as.data.frame() %>%
    rownames_to_column("soil_code") %>%
    mutate(eveness_shanon = shanon/richness,
           eveness_simson = simpson/richness) %>%
    pivot_longer(-soil_code, names_to = "hill_diversity_type", values_to = "diversity_values")
}

library(readxl)
library(GGally)
soil_df <- read_xlsx("data/soil_data_dry_season_2023.xlsx") %>%
  rename(coffee_system = "Robusta coffee system", soil_code = "SAMPLE NO", sand = "Sand%", clay = "Clay%", silt = "Silt%", pH ="soil pH", 
         OC = "OC%", MO = "OM%", N = "N%", P = "P/mg/kg", K = "K/Cmol/kg", Na = "Na/Cmol/kg") %>%
  mutate(soil_code = str_c("ROB_", soil_code,"_s"))

soil_df %>%
  select(-c(soil_code, Replicate,Zone, Cluster, coffee_system, Replicate)) %>%
  ggpairs()

amf_OTUs_rare %>%
  
  compute.hill.diversities  %>%
  left_join(metadata %>% select(sampling_date,
                                altitude,
                                soil_code), by = "soil_code") %>%
  left_join(soil_df, by = "soil_code")-> diversities_afm_df

## Richness
diversities_afm_df %>%
  pivot_longer(c(altitude,sand,clay,silt,pH,OC,MO,N,P,K,Na)) %>%
  filter(hill_diversity_type == "richness") %>%
  ggplot() +
  geom_point(aes(value, diversity_values)) +
  geom_smooth(aes(value, diversity_values), method ="glm") +
  facet_wrap(~name, scale = "free") +
  theme_minimal()

diversities_afm_df %>%
  pivot_longer(c(sampling_date, Zone, Cluster, coffee_system)) %>%
  filter(hill_diversity_type == "richness") %>%
  ggplot() +
  geom_boxplot(aes(value, diversity_values)) +
  facet_wrap(~name, scale = "free") +
  theme_minimal()
## Shanon
diversities_afm_df %>%
  pivot_longer(c(altitude,sand,clay,silt,pH,OC,MO,N,P,K,Na)) %>%
  filter(hill_diversity_type == "shanon") %>%
  ggplot() +
  geom_point(aes(value, diversity_values)) +
  geom_smooth(aes(value, diversity_values), method ="glm") +
  facet_wrap(~name, scale = "free") +
  theme_minimal()

diversities_afm_df %>%
  pivot_longer(c(sampling_date, Zone, Cluster, coffee_system)) %>%
  filter(hill_diversity_type == "shanon") %>%
  ggplot() +
  geom_boxplot(aes(value, diversity_values)) +
  facet_wrap(~name, scale = "free") +
  theme_minimal()

#Simpson
diversities_afm_df %>%
  pivot_longer(c(altitude,sand,clay,silt,pH,OC,MO,N,P,K,Na)) %>%
  filter(hill_diversity_type == "simpson") %>%
  ggplot() +
  geom_point(aes(value, diversity_values)) +
  geom_smooth(aes(value, diversity_values), method ="glm") +
  facet_wrap(~name, scale = "free") +
  theme_minimal()

diversities_afm_df %>%
  pivot_longer(c(sampling_date, Zone, Cluster, coffee_system)) %>%
  filter(hill_diversity_type == "simpson") %>%
  ggplot() +
  geom_boxplot(aes(value, diversity_values)) +
  facet_wrap(~name, scale = "free") +
  theme_minimal()

### test phylofactor
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# devtools::install_github('reptalex/phylofactor')
library(phylofactor)

tree_fr
PF <- PhyloFactor(Data = amf_OTUs_rare,X = soil_df$coffee_system, nfactors=3)


### test Hmsc

library(Hmsc)

nChains = 2
test.run = FALSE
if (test.run){
  #with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 5
} else {
  #with this option, the vignette evaluates slow but it reproduces the results of the
  #.pdf version
  thin = 5
  samples = 1000
  transient = 500*thin
  verbose = 500*thin
}

Y = as.matrix(amf_OTUs_rare %>%
                select(-c(ROB_M_s, ROB_L_s)) %>%
                filter(rowSums(across(starts_with("ROB"))) > 100) %>%
                t())
colSums(Y) %>% sort
soil_df %>% 
  filter(soil_code%in% rownames(Y))%>%
  select(coffee_system) %>%
  mutate(coffee_system = as.factor(coffee_system))-> soil_df2
dim(Y)
m1 = Hmsc(Y = Y, XData = soil_df2, XFormula = ~coffee_system)


m_converg = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)


