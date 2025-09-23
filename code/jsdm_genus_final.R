## ----setup, include=FALSE,purl = T--------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(purl = T)


## ----include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(mapview)
library(leaflet)
library(leafem)
library(rasterVis)
library(readxl)
library(sf)
library(pals)
library(viridis)

library(Hmsc)
library(ape)
library(purrr)

main_theme = theme_minimal() +
 theme(aspect.ratio = 1,
        axis.text.x = element_text(colour = "black", size=13, face="italic",angle =45, vjust =1, hjust =1),
        axis.text.y = element_text(colour = "black", size=13, face="italic"),
        axis.title = element_text(colour = "black", size=12),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.title = element_text(colour = "black", size=12,
                                    hjust =0.5),

        legend.text = element_text(colour = "black", size=13))

main_theme_final = theme_bw() +
 theme(aspect.ratio = 1,
        axis.text.x = element_text(colour = "black", size=7, face="italic", vjust =1, hjust =1),
        axis.text.y = element_text(colour = "black", size=7, face="italic"),
        axis.title = element_text(colour = "black", size=9),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        #panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(colour = "black", size=9,
                                    hjust =0.5),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(colour = "black", size=7,face="italic"),
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(0.05,0,0.05,0, "cm"),size = 7),
        title = element_text(colour = "black", size=9))


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_xlsx(path = "data/Metadata_Doreen_070525.xlsx")
soil_metadata <- read_xlsx(path = "data/soil_data_dry_season_2023.xlsx") %>% 
  select(-c(Zone, Cluster, Replicate,`Robusta coffee system`)) %>%
  rename(Sample_code = "SAMPLE NO",
          clay = "Clay%",
          N = "N%",
          OC = "OC%",
          pH = "soil pH",
          sand= "Sand%",
          silt = "Silt%",
          OM = "OM%",
          P = "P/mg/kg",
          K = "K/Cmol/kg",
          Na = "Na/Cmol/kg")

utm_coords <- metadata %>%
  mutate(
    zone = as.numeric(substr(easting, 1, 2)),
    band = substr(easting, 3, 3),
    easting = as.numeric(sub(".* ", "", easting)),
    northing = as.numeric(northing),

    hemisphere = ifelse(band == "N", "north", "south") # M = South Hemisphere
  )

# Function to convert a single row to lat/lon
convert_to_wgs84 <- function(zone, hemisphere, easting, northing,Sample_code) {
  epsg <- if (hemisphere == "north") {
    32600 + zone  # Northern hemisphere UTM EPSG codes
  } else {
    32700 + zone  # Southern hemisphere UTM EPSG codes
  }

  point <- st_sfc(st_point(c(easting, northing)), crs = epsg)
  point_wgs84 <- st_transform(point, crs = 4326)
  coords <- st_coordinates(point_wgs84)
  tibble( lon = coords[1], lat = coords[2],Sample_code = Sample_code)
}
utm_coords %>%
  mutate(zone = case_when(zone==32~ 36,
                             .default = zone),
          northing = case_when(Sample_code == "N17" ~ 86562,
                             .default = northing)) -> utm_coords
# Apply to each row
wgs84_coords <- pmap_dfr(utm_coords[, c("zone", "hemisphere", "easting", "northing","Sample_code")], convert_to_wgs84) %>%
  left_join(metadata, by = "Sample_code") %>%
  st_as_sf(coords = c("lon", "lat"),crs = 4326) %>%
  mutate(
     lon= st_coordinates(.)[,1],
     lat= st_coordinates(.)[,2]
  )


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ndvi_uganda <- raster("outputs/raster/ndvi_uganda_lower_res.tif")

ndvi_df <- raster::extract(ndvi_uganda, wgs84_coords)%>%
 tibble(ndvi = .,
    Sample_code  = wgs84_coords$Sample_code)



## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
precip_uganda_resampled <- raster("outputs/raster/precip_2023_lower_res.tif")

precip_df <- raster::extract(precip_uganda_resampled, wgs84_coords)%>%
 tibble(precip = .,
    Sample_code  = wgs84_coords$Sample_code)




## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmax_uganda_resampled <- raster("outputs/raster/tmax_2023_lower_res.tif")

tmax_df <- raster::extract(tmax_uganda_resampled, wgs84_coords)%>%
 tibble(tmax = .,
    Sample_code  = wgs84_coords$Sample_code)




## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmin_uganda_resampled <- raster("outputs/raster/tmin_2023_lower_res.tif")

tmin_df <- raster::extract(tmin_uganda_resampled, wgs84_coords)%>%
 tibble(tmin = .,
    Sample_code  = wgs84_coords$Sample_code)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
delta_t_uganda <- tmax_uganda_resampled -tmin_uganda_resampled

tdelta_df <- raster::extract(delta_t_uganda, wgs84_coords)%>%
 tibble(tdelta = .,
    Sample_code  = wgs84_coords$Sample_code)



## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wind_2023 <- raster("outputs/raster/wind_2023.tif")
pal_wind <- colorNumeric("Spectral", domain = c(0, 4))
wind_df <- raster::extract(wind_2023, wgs84_coords)%>%
 tibble(wind = .,
    Sample_code  = wgs84_coords$Sample_code)



## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wgs84_coords %>%
  left_join(soil_metadata, by = "Sample_code") %>%
  left_join(ndvi_df, by = "Sample_code") %>%
  left_join(precip_df, by = "Sample_code") %>%
  left_join(tmax_df, by = "Sample_code") %>%
  left_join(tmin_df, by = "Sample_code") %>%
  left_join(tdelta_df, by ="Sample_code") %>%
  left_join(wind_df, by ="Sample_code") -> wgs84_coords



## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
amf_genus_rare <- read.table("outputs/abundance_tables/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.afm_rarefy_genus.txt",
                             header = TRUE, stringsAsFactors = FALSE)
genus_tree <- read.tree("outputs/phylo_tree/Robust_MiSeq_18S_AMf_20241120_125_samples_Eukaryome_Glomeromycota_genus.tree")
plot(genus_tree)
genus_tree <- keep.tip(phy = genus_tree, tip = amf_genus_rare$genus)
plot(genus_tree)
genus_tree <- root(genus_tree, "g__Acaulospora_88")
plot(genus_tree)
# Clean tip labels (example transformation)
genus_rename <- genus_tree$tip.label %>%
  str_remove("g__") %>%
  str_replace_all("[a-z]__", "Unknown ") %>%
  str_replace("_(\\d+)$", " (\\1 OTU)")
length(genus_rename)
genus_tree$tip.label <- genus_rename
plot(genus_tree)
str(genus_tree)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare Y matrix for families
amf_genus_rare <- amf_genus_rare %>%
  mutate(genus = genus %>%
           str_remove("g__") %>%
           str_replace_all("[a-z]__", "Unknown ") %>%
           str_replace("_(\\d+)$", " (\\1 OTU)")) %>%
  remove_rownames() %>%
  column_to_rownames("genus")

Y_genus <- amf_genus_rare %>%
  filter(rowSums(across(starts_with("ROB")) != 0) > 3) %>%
  t()

# Prepare Y without zeros as a truncated version
Y_genus_trucated <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")

# Filter tree to the families present in Y
genus_tree_filter <- keep.tip(phy = genus_tree, tip = colnames(Y_genus))
Y_genus <- Y_genus[, genus_tree_filter$tip.label]
Y_genus_trucated <- Y_genus_trucated[, genus_tree_filter$tip.label]

# Sample names from Y
genus_matrix_samples <- rownames(Y_genus)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
covariate_df <- wgs84_coords %>% 
  as.data.frame() %>%
  mutate(Sample_code = str_c("ROB_", Sample_code, "_s")) %>%
  filter(Sample_code %in% genus_matrix_samples) %>%
  dplyr::select(Sample_code,ndvi, precip, tmax, tmin, tdelta, wind)

covariate_df %>%
  dplyr::filter(Sample_code %in% rownames(Y_genus)) %>%
  dplyr::slice(match(rownames(Y_genus), Sample_code)) %>%
  column_to_rownames("Sample_code") %>%
  mutate(across(where(is.numeric), ~vegan::decostand(.x, method = "standardize")[,1])) %>%
  rownames_to_column("Sample_code") %>%
  mutate(Sample_code = as.factor(Sample_code)) -> covariate_standa_df

studyDesign <- data.frame(sample = rownames(Y_genus), stringsAsFactors = TRUE)
rL <- HmscRandomLevel(units = studyDesign$sample)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mod_pois_null <- Hmsc(
  Y = log1p(Y_genus),
  XData =  covariate_standa_df,
  XFormula = ~ 1,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)


mod_pois_ndvi <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~ ndvi,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)


mod_pois_precip <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~ precip,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_tmax <- Hmsc(
  Y = log1p(Y_genus),
  XData =  covariate_standa_df,
  XFormula = ~ tmax,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_wind  <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~ wind,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_full <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_interaction <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~ (ndvi + precip + tmax +tdelta + wind)^2,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_poly <- Hmsc(
  Y = log1p(Y_genus),
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)


model_list <- list(
  mod_pois_null = mod_pois_null,
  mod_pois_ndvi      = mod_pois_ndvi,
  mod_pois_precip    = mod_pois_precip,
  mod_pois_tmax      = mod_pois_tmax,
  mod_pois_wind      = mod_pois_wind,
  mod_pois_full = mod_pois_full,
  mod_pois_interaction = mod_pois_interaction,
  mod_pois_poly = mod_pois_poly
)



## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
compute.hmsc <- function(hmsc_mod, nfold = 10, nChains = 12,   test.run = T){

  if (test.run){
    #with this option, the vignette runs fast but results are not reliable
    thin = 5
    samples = 20
    transient = 5*thin
    verbose = 5*thin
  } else {
    #with this option, the vignette evaluates slow but it reproduces the results of the
    #.pdf version
    thin = 10
    samples = 2000
    transient = 1000*thin
    verbose = 1000*thin
  }
  
  m_converg <- sampleMcmc(
    hmsc_mod,
    thin = thin, samples = samples, transient = transient,
    nChains = nChains, verbose = verbose, nParallel = nChains)
  
  set.seed(92828)
  partition = createPartition(m_converg, nfolds = nfold, column = "sample")
  
  expl_power = computePredictedValues(m_converg, thin = thin, nParallel = nChains)
  MF_expl = evaluateModelFit(hM = m_converg, predY = expl_power) 
  
  #Predictive power based on cross-validation  

  pred_power = computePredictedValues(m_converg, partition = partition, thin = thin, nParallel = nChains)
  MF_pred = evaluateModelFit(hM = m_converg, predY = pred_power)
  
  mPredY <- apply(abind::abind(pred_power, along = 3), c(1, 2), matrixStats::mean2)
  
  colnames(mPredY) <- colnames(m_converg$Y)
  rownames(mPredY) <- rownames(m_converg$Y)
  
  list(hmsc_mod = hmsc_mod, m_converg = m_converg , mPredY = mPredY, MF_expl = MF_expl, MF_pred = MF_pred)
  }
# save(m_converg, file = "outputs/hmsc_genus_full.Rdata")
# load(file = "outputs/hmsc_genus_full.Rdata")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1:length(model_list) %>%
  purrr::map(\(x) compute.hmsc(model_list[[x]], nfold = 5, nChains = 6,   test.run = F)) -> list_mods_crossvalidation


save(list_mods_crossvalidation,file =  "outputs/list_mods_crossvalidation_genus_pois_log.RData")
# load(file = "outputs/list_mods_crossvalidation_genus_pois_log.RData")








## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Y_genus_bin <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, 0, 1)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mod_bin_null <- Hmsc(
  Y = Y_genus_bin,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_ndvi <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)


mod_bin_precip <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ precip,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_tmax <- Hmsc(
  Y = Y_genus_bin,
  XData =  covariate_standa_df,
  XFormula = ~ tmax,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_wind  <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ wind,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_full <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_interaction <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ (ndvi + precip + tmax +tdelta + wind)^2,
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_poly <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = genus_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)



model_bin_list <- list(
  mod_bin_null      = mod_bin_null,
  mod_bin_ndvi      = mod_bin_ndvi,
  mod_bin_precip    = mod_bin_precip,
  mod_bin_tmax      = mod_bin_tmax,
  mod_bin_wind      = mod_bin_wind,
  mod_bin_full      = mod_bin_full,
  mod_bin_interaction = mod_bin_interaction,
  mod_bin_poly      = mod_bin_poly
)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1:length(model_bin_list) %>%
  purrr::map(\(x) compute.hmsc(model_bin_list[[x]], nfold = 5, nChains = 6,   test.run = T)) -> list_mods_bin_crossvalidation

save(list_mods_bin_crossvalidation, file = "outputs/list_mods_bin_crossvalidation_genus2.RData")
# load(file = "outputs/list_mods_bin_crossvalidation_genus2.RData")










































## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ROC_df %>%
  left_join(balanced_score_df %>%
              rename(threshold_B = "threshold"), by ="pred_name") %>%
  left_join(max_f1_score_df %>%
              dplyr::select(!F1_score) %>%
              rename(threshold_F1 = "threshold"), by ="pred_name") %>%
  arrange(threshold) %>%
  ggplot() +
  facet_wrap(~pred_name) +
  geom_path(aes(x=FPR,y=TPR), size=1) +
  geom_ribbon(aes(x=FPR,ymin=0,ymax=TPR),
              alpha=0.2) +
  geom_vline(aes(xintercept = FPR_b), col = "blue", linetype = 2, linewidth =  1.5)+
  geom_hline(aes(yintercept = TPR_b), col = "blue", linetype = 2, linewidth =  1.5)+
  #geom_vline(aes(xintercept=threshold_F1), col = "red")+
  coord_cartesian(xlim = c(0,1), ylim=c(0,1.01), expand = FALSE) +
  geom_abline( linetype='dashed') +
  labs(x='False Positive Rate \n (False presence predicted / all true absence)', y='True Positive Rate \n(True presence predicted / all true presence)',
       title = 'ROC curve of the model')+
  main_theme

recap_pred_bin_df %>%
  filter(mod_ID == 2, AUC >0.75| is.na(AUC)) %>% #, genus_names!="Diversisporaceae (36 OTU)"
  pull(genus_names) -> good_occurence_prediction




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggspatial)  # fonctions annotation_map_tile

