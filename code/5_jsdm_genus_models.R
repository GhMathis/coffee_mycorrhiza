## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(Hmsc)
library(ape)
library(furrr)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
test_run = T # For test, take few minutes
# test_run = F # True run, need 2 day to compute on 10 physical core with a frequency of approximately 4.6 GHz 
if(test_run){
  nfolds = 2
  nChains = 2
}else{
  nfolds =10
  nChains = 20
}


## ----cars---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
covariate_standa_df <- read.table("outputs/data_for_models/covariate_standa_df.txt", stringsAsFactors = T)
Y_genus <- read.table("outputs/data_for_models/Y_genus.txt")
genus_tree <- ape::read.tree("outputs/data_for_models/genus_tree.tree")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
genus_rename <- genus_tree$tip.label %>%
  str_remove("g__") %>%
  str_replace_all("[a-z]__", "Unknown ") %>%
  str_replace("_(\\d+)$", " (\\1 OTU)")

genus_tree$tip.label <- genus_rename

Y_genus <- Y_genus %>%
  rename_with(. %>%
                str_remove("g__") %>%
                str_replace_all("[a-z]__", "Unknown ") %>%
                str_replace("_(\\d+)$", " (\\1 OTU)"))

genus_tree$tip.label %in% colnames(Y_genus)
Y_genus <- Y_genus[,genus_tree$tip.label ]

Y_genus <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code") %>%
  log()
## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
covariate_standa_df[,c("lon", "lat")] %>% duplicated() -> temp 
covariate_standa_df %>%
  filter(lon == covariate_standa_df[temp,c("lon")],
         lat == covariate_standa_df[temp,c("lat")] )
# ROB_F7_s and ROB_F6_s have the same coordinates due to gps precision.
# Add a very small coordinate variation to ROB_F7_s for hmsc to accept cordiante data (hmcs dont take duplicated coordiante). It is also closer to the reality because ROB_F7_s and ROB_F6_s dont overlap.

covariate_standa_df[temp,c("lon", "lat")] <- covariate_standa_df[temp,c("lon", "lat")] + 0.0001
xy <- covariate_standa_df[,c("lon", "lat", "Sample_code")] %>%
  remove_rownames() %>%
  column_to_rownames("Sample_code")
studyDesign <- data.frame(sample = rownames(Y_genus), stringsAsFactors = TRUE)
rL <- HmscRandomLevel(sData = xy, longlat = T)


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mod_pois_null <- Hmsc(
  Y = Y_genus,
  XData =  covariate_standa_df,
  XFormula = ~ 1,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)


mod_pois_ndvi <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~ ndvi,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)


mod_pois_precip <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~ precip,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)

mod_pois_tmax <- Hmsc(
  Y = Y_genus,
  XData =  covariate_standa_df,
  XFormula = ~ tmax,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)

mod_pois_wind  <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~ wind,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)

mod_pois_full <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)

mod_pois_interaction <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~ (ndvi + precip + tmax +tdelta + wind)^2,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
)

mod_pois_poly <- Hmsc(
  Y = Y_genus,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "normal",
  YScale = TRUE
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


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
compute.hmsc <- function(hmsc_mod, nfold = 10, nChains = 10, samples = 2000,  test.run = T){
  
  if (test.run){
    #with this option, the vignette runs fast but results are not reliable
    thin = 5
    samples = 20
    transient = 5*thin
    verbose = 5*thin
  } else {
    #with this option, the vignette evaluates slow but it reproduces the results of the
    #.pdf version
    thin = 5
    transient = 1000*thin
    verbose = 1000*thin
  }
  
  m_converg <- sampleMcmc(
    hmsc_mod,
    thin = thin, samples = samples, transient = transient,
    nChains = nChains, verbose = verbose, nParallel = nChains)
  
  #  set.seed(92828)
  #  partition = createPartition(m_converg, nfolds = nfold, column = "sample")
  #  
  #  expl_power = computePredictedValues(m_converg, thin = thin, nParallel = nChains)
  #  MF_expl = evaluateModelFit(hM = m_converg, predY = expl_power) 
  #  
  #  #Predictive power based on cross-validation  
  #
  #  pred_power = computePredictedValues(m_converg, partition = partition, thin = thin, nParallel = nChains)
  #  MF_pred = evaluateModelFit(hM = m_converg, predY = pred_power)
  #  
  #  mPredY <- apply(abind::abind(pred_power, along = 3), c(1, 2), matrixStats::mean2)
  #  
  #  colnames(mPredY) <- colnames(m_converg$Y)
  #  rownames(mPredY) <- rownames(m_converg$Y)
  
  list(hmsc_mod = hmsc_mod, m_converg = m_converg) 
  # , mPredY = mPredY, MF_expl = MF_expl, MF_pred = MF_pred)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model_list %>%
  purrr::map(\(x) compute.hmsc(x, nfold = nfolds, nChains = nChains, samples = 2000,  test.run = test_run)
  ) -> list_mods

if(file.exists("outputs/models/list_mods_genus_pois_autocorr.RData")){
  # safe save in case a model list in already present in the folder
  save(list_mods, file = "outputs/models/list_mods_genus_pois_autocorr_safesave.RData")
}else{
  save(list_mods, file = "outputs/models/list_mods_genus_pois_autocorr.RData")
}



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Y_genus_bin <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(is.na(value), 0, 1)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mod_bin_null <- Hmsc(
  Y = Y_genus_bin,
  XData =  covariate_standa_df ,
  XFormula = ~ 1,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_ndvi <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)


mod_bin_precip <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ precip,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_tmax <- Hmsc(
  Y = Y_genus_bin,
  XData =  covariate_standa_df,
  XFormula = ~ tmax,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_wind  <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ wind,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_full <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ ndvi + precip + tmax +tdelta + wind,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_interaction <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~ (ndvi + precip + tmax +tdelta + wind)^2,
  phyloTree = genus_tree,
  studyDesign = studyDesign,
  ranLevels = list(sample = rL),
  distr = "probit"
)

mod_bin_poly <- Hmsc(
  Y = Y_genus_bin,
  XData = covariate_standa_df,
  XFormula = ~poly(ndvi,2) + poly(precip,2) + poly(tmax,2) + poly(tdelta,2) + poly(wind,2),
  phyloTree = genus_tree,
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


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


model_bin_list %>%
  purrr::map(\(x) compute.hmsc(x, nfold = nfolds, nChains = nChains, samples = 1000, test.run = test_run)
             )-> list_mods_bin

if(file.exists("outputs/models/list_mods_bin_autocorr_genus.RData")){
  # safe save in case a model list in already present in the folder
  save(list_mods_bin, file = "outputs/models/list_mods_bin__autocorr_genus_safesave.RData")
}else{
  save(list_mods_bin, file = "outputs/models/list_mods_bin_autocorr_genus.RData")
}


