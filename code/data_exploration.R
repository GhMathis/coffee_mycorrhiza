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
  select(amplicon, taxonomy, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB", ignore.case = FALSE)) %>%
  #{{"Arbuscular mycorrhizal fungi & Mucoromycotina fine root endophytes" are selected}}
  filter(str_detect(taxonomy,"Glomeromycotina")|str_detect(taxonomy,"Endogonales")) %>%
  #{{select only OTUs and controls}}
  select(amplicon, starts_with("CONT", ignore.case = FALSE),
         starts_with("ROB_", ignore.case = FALSE)) %>%
  droplevels()
amf_OTUs
taxonomy <- data_18S %>%
  #{{import data}}
  read_tsv() %>%
  filter(str_detect(taxonomy,"Glomeromycotina")|str_detect(taxonomy,"Endogonales")) %>%
  select(amplicon, taxonomy) %>%
  mutate(taxonomy = str_replace_all(taxonomy, ":","_"),
         taxonomy = str_replace_all(taxonomy, "[|]","; "))


d <- amf_OTUs %>% 
  replace(amf_OTUs == 0, NA) %>%
  select(amplicon, starts_with("CONT")) %>%
  pivot_longer(-amplicon, names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads))

amf_OTUs_decont <- amf_OTUs %>%
  replace(amf_OTUs == 0, NA) %>%
  pivot_longer(-c("amplicon", starts_with("CONT")), names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, amplicon, wt = reads), by = "amplicon") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n,-starts_with("CONT")) %>%
  pivot_wider(names_from = samples, values_from =  reads, values_fill = 0) %>%
  as_tibble() %>%
  column_to_rownames(var = "amplicon") 
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
library(ape)




tree_amf <- read.tree(file = "data/18S_AMf_PR2/Robust_MiSeq_18S_AMf_20241120_125_samples_Glomeromycotina.tree")
dim(amf_OTUs_decont2)
pf.heatmap(tree=tree_amf,Data=log1p(amf_OTUs_decont2),color = NA)
pf.heatmap(tree=tree_amf,Data=log1p(amf_OTUs_rare))
rownames(amf_OTUs_rare) %in% tree_amf$tip.label

plot(tree_amf,show.tip.label = F)

soil_df %>% 
  filter(soil_code%in% colnames(amf_OTUs_rare))%>%
  mutate(coffee_system = as.factor(coffee_system)) %>%
  column_to_rownames("soil_code")-> soil_df2
amf_OTUs_rare%>%
  select(-c(ROB_M_s, ROB_L_s)) -> Data_abund
Data_abund = as.matrix(Data_abund[tree_amf$tip.label,])

t(t(amf_OTUs_decont2 )/colSums(amf_OTUs_decont2)) %>%
  as.data.frame %>%
  select(-c(ROB_M_s, ROB_L_s)) -> amf_OTUs_relative

soil_df2 %>%
  arrange(coffee_system) %>%
  rownames -> arrange_sample
amf_OTUs_relative = as.matrix(amf_OTUs_relative[tree_amf$tip.label,arrange_sample])

soil_df2$coffee_system %>% unique
temp_amf_OTUs_relative <- amf_OTUs_relative

colnames(temp_amf_OTUs_relative) <- soil_df2 %>%
  arrange(coffee_system) %>%
  pull(coffee_system)

PF <- PhyloFactor(temp_amf_OTUs_relative,tree_amf,soil_df2$coffee_system,nfactors=3)
library(phytools)
clr <- function(Matrix) apply(Matrix,MARGIN=2,FUN=function(x) log(x)-mean(log(x)))
phylo.heatmap(tree_amf,clr(PF$Data), ftype = "off")
pf.heatmap(PF,factors=1:3)

bin.projection <- pf.BINprojection(PF,factor=3)   
bin.taxa <- bin.projection$otus %>% 
  lapply(.,FUN=function(otus,tax) tax[otus,],tax=as.data.frame(taxonomy) %>% column_to_rownames("amplicon"))

PF$nfactors
bin.observed.geometricMeans <- pf.BINprojection(PF)$Data
bin.observed.relativeAbundances <- pf.BINprojection(PF,rel.abund=T)$Data
bin.predicted.geometricMeans <- pf.BINprojection(PF,prediction=T)$Data
bin.predicted.relativeAbundances <- pf.BINprojection(PF,prediction=T,rel.abund=T)$Data

library(plotrix)
par(mfrow=c(2,2))
stackpoly(t(bin.observed.geometricMeans),stack=T,main='Observed gMean',xat = c(18,28,46,66,84),xaxlab = c("RCB", "RCM", "RCT", "RCTB", "RCTBA"))
stackpoly(t(bin.predicted.geometricMeans),stack=T,main='Predicted gMean',xat = c(18,28,46,66,84),xaxlab = c("RCB", "RCM", "RCT", "RCTB", "RCTBA"))

stackpoly(t(bin.observed.relativeAbundances),stack=T,main='Observed Rel-Abund.',xat = c(18,28,46,66,84),xaxlab = c("RCB", "RCM", "RCT", "RCTB", "RCTBA"))
stackpoly(t(bin.predicted.relativeAbundances),stack=T,main='Predicted Rel-Abund.',xat = c(18,28,46,66,84),xaxlab = c("RCB", "RCM", "RCT", "RCTB", "RCTBA"))
par(mfrow=c(1,1))
prediction <- pf.predict(PF,factors=3)
phylo.heatmap(tree_amf,clr(prediction), ftype = "off")
# PA_datatable <- matrix.to.phyloframe(Data_abund,soil_df2$coffee_system,data.name = 'Abundance')
# PF <- gpf(Data = Data_abund, 
#           tree =tree_amf, MetaData = soil_df2,
#           frmla.phylo = Data~phylo*coffee_system,
#           family = poisson,
#           nfactors=3, algorithm='mStable')
# PF
phyca <- PhyCA(Data_abund,tree_amf,ncores=8,ncomponents = 5)

pf.heatmap(PF,factors=1:3)
summary(Phylofactorization)
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
colSums(Y) %>% sort
Y = as.matrix(amf_OTUs_rare %>%
                select(-c(ROB_M_s, ROB_L_s)) %>%
                filter(rowSums(across(starts_with("ROB"))) >=1) %>%
                t())
colSums(Y) %>% sort
soil_df %>% 
  filter(soil_code%in% rownames(Y))%>%
  mutate(coffee_system = as.factor(coffee_system)) %>%
  column_to_rownames("soil_code")-> soil_df2
dim(Y)
m1 = Hmsc(Y = Y, XData = soil_df2, XFormula = ~coffee_system)


m_converg = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
####
library(PLNmodels)
str(Y)
str(soil_df2)
afm_data <- prepare_data(Y, soil_df2)
myPLN <- PLN(Abundance ~ coffee_system , afm_data)

plot(myPLN)

library(FactoMineR)
library(factoextra)
PCA(Y) -> pca_amf
fviz_pca_ind(pca_amf, col.ind = soil_df2$coffee_system)

#####
library(picante)

library(vegan)

rowSums(Data_abund)
apply(Data_abund, 1, which.max) -> OTU_optima_index
as.data.frame(OTU_optima_index) %>% 
  mutate(OTU_optima_sample = colnames(Data_abund)[OTU_optima_index],
         OTU_tot_abund = apply(Data_abund, 1, sum))%>%
  filter(OTU_tot_abund !=0) %>%
  rownames_to_column("OTU") %>%
  left_join(soil_df2%>%
              rownames_to_column("soil_code"), by = join_by("OTU_optima_sample" == "soil_code")) -> optima_df
optima_df %>%
  column_to_rownames("OTU") %>%
  select(pH, N, P,K,Na,MO,sand,clay) %>%
  dist -> dist_optima

drop.tip(tree_amf, setdiff(tree_amf$tip.label, optima_df$OTU)) -> tree_amf_subset
cophenetic(tree_amf_subset) -> amf_phylo_dist
dist(Data_abund) -> afm_abund_dist
str(dist_optima)
str(amf_phylo_dist)
vegan::mantel.correlog(dist_optima,amf_phylo_dist,n.class=50, nperm =999) -> mantel_phylo_optima
plot(mantel_phylo_optima)
vegan::mantel.correlog(afm_abund_dist,amf_phylo_dist,n.class=50, nperm =999) -> mantel_phylo
plot(mantel_phylo)

dist(t(amf_OTUs_relative)) -> afm_abund_rela_dist
vegan::mantel.correlog(afm_abund_rela_dist,amf_phylo_dist,n.class=50, nperm =999) -> mantel_phylo_rela
plot(mantel_phylo_rela)

betaMNTD <- picante::comdistnt(amf_OTUs_relative, amf_phylo_dist, abundance.weighted = T)
