## ----setup, include=FALSE,purl = T---------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(purl = T)


## ----include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(mapview)
library(leaflet)
library(leafem)
library(rasterVis)
library(readxl)
library(sf)
library(pals)

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


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
metadata <- read_xlsx(path = "data/Metadata_Doreen_070525.xlsx")

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


leaflet() %>% 
  addTiles() %>%
  addMarkers(data=wgs84_coords,~lon, ~lat, popup = ~paste(round(lat,2), round(lon,2)) , label = ~Sample_code)%>%
  addCircleMarkers(data = wgs84_coords, lng = ~lon,
                  lat = ~lat, radius = 8, color = "blue", fillOpacity = 0.8)


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
land_degradation_raster <- raster("outputs/raster/land_degradation_reclassify.tif")
land_degradation_type <- levels(land_degradation_raster)[[1]]
land_degradation_type <- land_degradation_type %>%
  mutate(land_degra = factor(land_degra, levels = c("Low", "Medium", "High")))
land_degradation_df <- raster::extract(land_degradation_raster,  wgs84_coords)%>%
  tibble(ID = .,
  Sample_code  = wgs84_coords$Sample_code)%>%
  left_join(land_degradation_type,  by = "ID")

myPal1 <- rev(RColorBrewer::brewer.pal('PRGn', n=3))


leaflet() %>%
  addTiles() %>% 
  addRasterImage(land_degradation_raster, opacity = 0.7,
                 color = myPal1) %>%
   addLegend(
    "bottomright",
    pal = colorFactor(palette=myPal1, domain = land_degradation_type$land_degra),
    values = land_degradation_type$land_degra,
    title = "Land Degradation"
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
land_cover_raster <- raster("outputs/raster/land_cover_lower_res.tif")

land_cover_type <- levels(land_cover_raster)[[1]] %>%
  dplyr::select(ID,Land_Cover) %>%
  mutate(Land_Cover = factor(Land_Cover, levels = levels(land_cover_raster)[[1]]$Land_Cover))


land_cover_df <- raster::extract(land_cover_raster, wgs84_coords) %>%
 tibble(ID = .,
    Sample_code  = wgs84_coords$Sample_code) %>%
    left_join(land_cover_type, by = "ID")

myPal2 <- pals::stepped2(n=18)

pal_land_cover <- colorFactor(myPal2[unique(land_cover_df$ID)], 
                   domain = wgs84_coords$Land_Cover)
leaflet() %>%
  addTiles() %>% 
  addRasterImage(land_cover_raster, opacity = 0.7,
                 colors = myPal2) %>%
  addLegend(
    "bottomright",
    pal = colorFactor(palette=myPal2, domain = land_cover_type$Land_Cover),
    values =land_cover_type$Land_Cover,
    title = "Land Cover"
  )


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
cluster_agro_raster <- raster("outputs/raster/cluster_agro_lower_res.tif")
cluster_agro_type <- levels(cluster_agro_raster)[[1]] %>%
  dplyr::select(ID, Cluster) %>%
  mutate(Cluster = factor(Cluster, levels = levels(cluster_agro_raster)[[1]]$Cluster))

cluster_agro_df <- raster::extract(cluster_agro_raster, wgs84_coords) %>%
  tibble(ID = ., Sample_code = wgs84_coords$Sample_code) %>%
  left_join(cluster_agro_type, by = "ID")

myPal3 <- ggsci::pal_simpsons()(8)

leaflet() %>%
  addTiles() %>%
  addRasterImage(cluster_agro_raster, opacity = 0.7, colors = myPal3) %>%
  addLegend(
    "bottomright",
    pal = colorFactor(palette = myPal3, domain = cluster_agro_type$Cluster),
    values = cluster_agro_type$Cluster,
    title = "Land Cover"
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
ndvi_uganda <- raster("outputs/raster/ndvi_uganda_lower_res.tif")

ndvi_df <- raster::extract(ndvi_uganda, wgs84_coords)%>%
 tibble(ndvi = .,
    Sample_code  = wgs84_coords$Sample_code)

pal_ndvi <- colorNumeric("viridis", domain = c(0, 1))
leaflet() %>%
  addTiles() %>% 
  addRasterImage(ndvi_uganda, opacity = 0.7) %>%
  addLegend(
    "bottomright",
    pal = pal_ndvi,
    values = c(0, 1),
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
precip_uganda_resampled <- raster("outputs/raster/precip_2023_lower_res.tif")

precip_df <- raster::extract(precip_uganda_resampled, wgs84_coords)%>%
 tibble(precip = .,
    Sample_code  = wgs84_coords$Sample_code)
pal_precip <- colorNumeric("Spectral", domain = c(0, 1))
leaflet() %>%
  addTiles() %>%
  addRasterImage(precip_uganda_resampled, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_precip,
    values = c(0, 1),
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
tmax_uganda_resampled <- raster("outputs/raster/tmax_2023_lower_res.tif")

tmax_df <- raster::extract(tmax_uganda_resampled, wgs84_coords)%>%
 tibble(tmax = .,
    Sample_code  = wgs84_coords$Sample_code)
pal_t <- colorNumeric("Spectral", domain = c(0, 34))
leaflet() %>%
  addTiles() %>%
  addRasterImage(tmax_uganda_resampled, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_t,
    values = c(0, 34),
  )



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
tmin_uganda_resampled <- raster("outputs/raster/tmin_2023_lower_res.tif")

tmin_df <- raster::extract(tmin_uganda_resampled, wgs84_coords)%>%
 tibble(tmin = .,
    Sample_code  = wgs84_coords$Sample_code)
leaflet() %>%
  addTiles() %>%
  addRasterImage(tmin_uganda_resampled, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_t,
    values = c(0, 34),
  )



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
delta_t_uganda <- tmax_uganda_resampled -tmin_uganda_resampled

tdelta_df <- raster::extract(delta_t_uganda, wgs84_coords)%>%
 tibble(tdelta = .,
    Sample_code  = wgs84_coords$Sample_code)
tdelta_df %>% summary()
leaflet() %>%
  addTiles() %>%
  addRasterImage(delta_t_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_t,
    values = c(9, 12),
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
wind_2023 <- raster("outputs/raster/wind_2023.tif")
pal_wind <- colorNumeric("Spectral", domain = c(0, 4))
wind_df <- raster::extract(wind_2023, wgs84_coords)%>%
 tibble(wind = .,
    Sample_code  = wgs84_coords$Sample_code)
leaflet() %>%
  addTiles() %>%
  addRasterImage(wind_2023, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_wind,
    values = c(0, 4),
  )



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
altitude_uganda <- raster("outputs/raster/altitude_uganda.tif")
pal_alt <- colorNumeric("Spectral", domain = c(0, 4500))
altidute_df <- raster::extract(altitude_uganda, wgs84_coords)%>%
 tibble(altitude = .,
    Sample_code  = wgs84_coords$Sample_code)
leaflet() %>%
  addTiles() %>%
  addRasterImage(altitude_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_alt,
    values = c(0, 4),
  )



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
ph_uganda <- raster("outputs/raster/pH.tif")
pal_ph <- colorNumeric("Spectral", domain = c(0, 75))

ph_df <- raster::extract(ph_uganda, wgs84_coords)%>%
 tibble(pH = .,
    Sample_code  = wgs84_coords$Sample_code)

leaflet() %>%
  addTiles() %>%
  addRasterImage(ph_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_ph,
    values = c(0, 75),
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
sand_uganda <- raster("outputs/raster/sand.tif")
pal_sand <- colorNumeric("Spectral", domain = c(0, 550))

sand_df <- raster::extract(sand_uganda, wgs84_coords)%>%
 tibble(sand = .,
    Sample_code  = wgs84_coords$Sample_code)

leaflet() %>%
  addTiles() %>%
  addRasterImage(sand_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_sand,
    values = c(0, 550),
  )



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
nitrogen_uganda <- raster("outputs/raster/nitrogen.tif")
pal_nitrogen <- colorNumeric("Spectral", domain = c(0, 748))

nitrogen_df <- raster::extract(nitrogen_uganda, wgs84_coords)%>%
 tibble(nitrogen = .,
    Sample_code  = wgs84_coords$Sample_code)

leaflet() %>%
  addTiles() %>%
  addRasterImage(nitrogen_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_nitrogen,
    values = c(0, 748),
  )



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
OC_uganda <- raster("outputs/raster/organic_carbone.tif")
pal_OC <- colorNumeric("Spectral", domain = c(0, 1305))

OC_df <- raster::extract(OC_uganda, wgs84_coords)%>%
 tibble(OC = .,
    Sample_code  = wgs84_coords$Sample_code)

leaflet() %>%
  addTiles() %>%
  addRasterImage(OC_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_OC,
    values = c(0, 1305),
  )



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
clay_uganda <- raster("outputs/raster/clay.tif")

pal_clay <- colorNumeric("Spectral", domain = c(0, 484))

clay_df <- raster::extract(clay_uganda, wgs84_coords)%>%
 tibble(clay = .,
    Sample_code  = wgs84_coords$Sample_code)

leaflet() %>%
  addTiles() %>%
  addRasterImage(clay_uganda, opacity = 0.5)%>%
    addLegend(
    "bottomright",
    pal = pal_clay,
    values = c(0, 484),
  )


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
wgs84_coords %>%
  left_join(cluster_agro_df, by = "Sample_code") %>%
  left_join(land_degradation_df, by = "Sample_code") %>%
  left_join(ndvi_df, by = "Sample_code") %>%
  left_join(precip_df, by = "Sample_code") %>%
  left_join(tmax_df, by = "Sample_code") %>%
  left_join(tmin_df, by = "Sample_code") %>%
  left_join(tdelta_df, by ="Sample_code") %>%
  left_join(wind_df, by ="Sample_code") %>%
  left_join(altidute_df, by ="Sample_code") %>% 
  left_join(ph_df, by ="Sample_code") %>%
  left_join(sand_df, by ="Sample_code") %>%
  left_join(nitrogen_df, by ="Sample_code") %>%
  left_join(OC_df, by ="Sample_code") %>%
  left_join(clay_df, by ="Sample_code") -> wgs84_coords



## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
amf_family_decont <- read.table("outputs/abundance_tables/Robust_MiSeq_18S_AMf_20241120_125_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.rarefy.family.txt", header = TRUE, stringsAsFactors = FALSE)
family_tree <- read.tree("outputs/phylo_tree/Robust_MiSeq_18S_AMf_20241120_125_samples_Eukaryome_Glomeromycota_family.tree")

# Inspect and root tree
family_tree$tip.label
plot(family_tree)
family_tree <- root(family_tree, "c__Archaeosporomycetes_4")

# Clean tip labels (example transformation)
family_rename <- family_tree$tip.label %>%
  str_remove("f__") %>%
  str_replace_all("[a-z]__", "Unknown ") %>%
  str_replace("_(\\d+)$", " (\\1 OTU)")
length(family_rename)
family_tree$tip.label <- family_rename
plot(family_tree)
str(family_tree)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare Y matrix for families
amf_family_decont <- amf_family_decont %>%
  mutate(family = family %>%
           str_remove("f__") %>%
           str_replace_all("[a-z]__", "Unknown ") %>%
           str_replace("_(\\d+)$", " (\\1 OTU)")) %>%
  remove_rownames() %>%
  column_to_rownames("family")

Y_family <- amf_family_decont %>%
  filter(rowSums(across(starts_with("ROB")) != 0) > 3) %>%
  t()

# Prepare Y without zeros as a truncated version
Y_family_trucated <- Y_family %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")

# Filter tree to the families present in Y
family_tree_filter <- keep.tip(phy = family_tree, tip = colnames(Y_family))
Y_family <- Y_family[, family_tree_filter$tip.label]
Y_family_trucated <- Y_family_trucated[, family_tree_filter$tip.label]

# Sample names from Y
family_matrix_samples <- rownames(Y_family)



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(GGally)
library(FactoMineR)
covariate_df <- wgs84_coords %>%
  as.data.frame() %>%
  mutate(Sample_code = str_c("ROB_", Sample_code, "_s")) %>%
  filter(Sample_code %in% family_matrix_samples) %>%
  dplyr::select(Sample_code,ndvi, precip, tmax, tmin,tdelta, wind, altitude, clay,OC,nitrogen,sand,pH)

covariate_df %>%
  dplyr::filter(Sample_code %in% rownames(Y_family)) %>%
  dplyr::slice(match(rownames(Y_family), Sample_code)) %>%
  column_to_rownames("Sample_code") %>%
  mutate(across(where(is.numeric), ~vegan::decostand(.x, method = "standardize")[,1])) %>%
  rownames_to_column("Sample_code") %>%
  mutate(Sample_code = as.factor(Sample_code)) -> covariate_standa_df
# write.table(covariate_df, "outputs/covariate_table.txt")

covariate_standa_df %>%
  select(-Sample_code) %>%
  ggpairs(progress = F)

covariate_standa_df %>%
  select(-Sample_code) %>%
  PCA()

covariate_standa_df %>%
  select(ndvi, precip, tmax, tdelta, wind, OC, nitrogen,clay, pH, sand) %>%
  ggpairs(progress = F)

covariate_standa_df %>%
  select(ndvi, precip, tmax, tdelta, wind, OC, nitrogen,clay, pH, sand) %>%
  PCA()

covariate_standa_df %>%
  filter(sand <= -4)

covariate_standa_df %>%
  select(ndvi, precip, tmax, tdelta, wind, OC, nitrogen,clay, pH, sand) %>%
  ggpairs(progress = F)
covariate_standa_df %>%
  select(ndvi, precip, tmax, tdelta, wind) %>%
  ggpairs(progress = F)


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------

# 
# nChains = 12
# test.run = F
# if (test.run){
#   #with this option, the vignette runs fast but results are not reliable
#   thin = 5
#   samples = 10
#   transient = 5*thin
#   verbose = 5*thin
# } else {
#   #with this option, the vignette evaluates slow but it reproduces the results of the
#   #.pdf version
#   thin = 10
#   samples = 750
#   transient = 1000*thin
#   verbose = 1000*thin
# }
# covariate_standa_df %>%
#   filter(Sample_code != "ROB_M19_s")  -> covariate_standa_filter_df
# Y_family %>%
#   as.data.frame() %>%
#   rownames_to_column("Sample_code") %>%
#   filter(Sample_code != "ROB_M19_s") %>%
#   column_to_rownames("Sample_code") %>%
#   as.matrix() -> Y_family_filter
studyDesign <- data.frame(sample = rownames(Y_family), stringsAsFactors = TRUE)
rL <- HmscRandomLevel(units = studyDesign$sample)


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------

# Y_family_filter %>% 
#   rowSums %>% 
#   as.data.frame() %>%
#   rename(readscount = ".") %>%
#   rownames_to_column("Sample_code") -> readscount_df
# 
# covariate_standa_filter_df %>%
#   left_join(readscount_df, by = "Sample_code") %>%
#   mutate(Sample_code = as.factor(Sample_code)) -> covariate_standa_filter_df


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
Y_family_trucated <- Y_family %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------

mod_pois_trunc_null <- Hmsc(
  Y = Y_family_trucated,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

# mod_pois_ndvi <- Hmsc(
#   Y = Y_family,
#   XData = covariate_standa_df,
#   XFormula = ~ ndvi,
#   #phyloTree = family_tree_filter,
#   studyDesign = studyDesign,
#   ranLevels = list(sample = rL),
#   distr = "poisson"
# )
# 
# 
# mod_pois_precip <- Hmsc(
#   Y = Y_family,
#   XData = covariate_standa_df,
#   XFormula = ~ precip,
#   #phyloTree = family_tree_filter,
#   studyDesign = studyDesign,
#   ranLevels = list(sample = rL),
#   distr = "poisson"
# )
# 
# mod_pois_tmax <- Hmsc(
#   Y = Y_family,
#   XData =  covariate_standa_df,
#   XFormula = ~ tmax,
#   #phyloTree = family_tree_filter,
#   studyDesign = studyDesign,
#   ranLevels = list(sample = rL),
#   distr = "poisson"
# )
# 
# mod_pois_wind  <- Hmsc(
#   Y = Y_family,
#   XData = covariate_standa_df,
#   XFormula = ~ wind,
#   #phyloTree = family_tree_filter,
#   studyDesign = studyDesign,
#   ranLevels = list(sample = rL),
#   distr = "poisson"
# )

mod_pois_trunc_full <- Hmsc(
  Y = Y_family_trucated,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_trunc_poly <- Hmsc(
  Y = Y_family_trucated,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)
mod_pois_null_log <- Hmsc(
  Y = log1p(Y_family),
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_full_log <- Hmsc(
  Y = log1p(Y_family),
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_poly_log <- Hmsc(
  Y = log1p(Y_family),
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_null_trunc_log <- Hmsc(
  Y = log(Y_family_trucated),
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_full_trunc_log <- Hmsc(
  Y = log(Y_family_trucated),
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_pois_poly_trunc_log <- Hmsc(
  Y = log(Y_family_trucated),
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "poisson"
)

mod_lognorm_null <- Hmsc(
  Y = Y_family,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  #phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_full <- Hmsc(
  Y = Y_family,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_poly <- Hmsc(
  Y = Y_family,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_null_log <- Hmsc(
  Y = log1p(Y_family),
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_full_log <- Hmsc(
  Y = log1p(Y_family),
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_poly_log <- Hmsc(
  Y = log1p(Y_family),
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_null_trunc <- Hmsc(
  Y = Y_family_trucated,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_full_trunc <- Hmsc(
  Y = Y_family_trucated,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)

mod_lognorm_poly_trunc <- Hmsc(
  Y = Y_family_trucated,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "lognormal poisson"
)


model_list <- list(
  mod_pois_trunc_null      = mod_pois_trunc_null,
  # mod_pois_ndvi      = mod_pois_ndvi,
  # mod_pois_precip    = mod_pois_precip,
  # mod_pois_tmax      = mod_pois_tmax,
  # mod_pois_wind      = mod_pois_wind,
  mod_pois_trunc_full      = mod_pois_trunc_full,
  mod_pois_trunc_poly      = mod_pois_trunc_poly,
  mod_pois_null_log  = mod_pois_null_log,
  mod_pois_full_log  = mod_pois_full_log,
  mod_pois_poly_log  = mod_pois_poly_log,
  mod_pois_null_trunc_log  = mod_pois_null_trunc_log,
  mod_pois_full_trunc_log  = mod_pois_full_trunc_log,
  mod_pois_poly_trunc_log  = mod_pois_poly_trunc_log,
  mod_lognorm_null   = mod_lognorm_null,
  mod_lognorm_full   = mod_lognorm_full,
  mod_lognorm_poly   = mod_lognorm_poly,
  mod_lognorm_null_log   = mod_lognorm_null_log,
  mod_lognorm_full_log   = mod_lognorm_full_log,
  mod_lognorm_poly_log   = mod_lognorm_poly_log,
  mod_lognorm_null_trunc   = mod_lognorm_null_trunc,
  mod_lognorm_full_trunc   = mod_lognorm_full_trunc,
  mod_lognorm_poly_trunc   = mod_lognorm_poly_trunc
)


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
compute.hmsc <- function(hmsc_mod, nfold = 10, nChains = 12,   test.run = T){

  if (test.run){
    #with this option, the vignette runs fast but results are not reliable
    thin = 5
    samples = 10
    transient = 5*thin
    verbose = 5*thin
  } else {
    #with this option, the vignette evaluates slow but it reproduces the results of the
    #.pdf version
    thin = 10
    samples = 1000
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
# save(m_converg, file = "outputs/hmsc_family_full.Rdata")
# load(file = "outputs/hmsc_family_full.Rdata")



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
1:length(model_list) %>%
  purrr::map(\(x) compute.hmsc(model_list[[x]], nfold = 10, nChains = 16,   test.run = F)) -> list_mods_crossvalidation

save(list_mods_crossvalidation,file =  "list_mods_crossvalidation_family.RData")
#load(file = "list_mods_bin_crossvalidation.RData")










## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
Y_family_bin <- Y_family %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, 0, 1)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
mod_bin_null <- Hmsc(
  Y = Y_family_bin,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_ndvi <- Hmsc(
  Y = Y_family_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)


mod_bin_precip <- Hmsc(
  Y = Y_family_bin,
  XData = covariate_standa_df,
  XFormula = ~ precip,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_tmax <- Hmsc(
  Y = Y_family_bin,
  XData =  covariate_standa_df,
  XFormula = ~ tmax,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_wind  <- Hmsc(
  Y = Y_family_bin,
  XData = covariate_standa_df,
  XFormula = ~ wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_full <- Hmsc(
  Y = Y_family_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_poly <- Hmsc(
  Y = Y_family_bin,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = family_tree_filter,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)



model_bin_list <- list(
  mod_bin_null      = mod_bin_null,
  # mod_bin_ndvi      = mod_bin_ndvi,
  # mod_bin_precip    = mod_bin_precip,
  # mod_bin_tmax      = mod_bin_tmax,
  # mod_bin_wind      = mod_bin_wind,
  mod_bin_full      = mod_bin_full,
  mod_bin_poly      = mod_bin_poly
)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------

1:length(model_bin_list) %>%
  purrr::map(\(x) compute.hmsc(model_bin_list[[x]], nfold = 10, nChains = 14,   test.run = F)) -> list_mods_bin_crossvalidation

save(list_mods_bin_crossvalidation, file = "list_mods_bin_crossvalidation_family.RData")
# load(file = "list_mods_truebin_crossvalidation.RData")

