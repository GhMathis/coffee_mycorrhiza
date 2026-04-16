## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(Hmsc)
library(ape)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# test_run = T # For test, take few minutes
test_run = F # True run, need 2 day to compute on 10 physical core with a frequency of approximately 4.6 GHz
if(test_run){
  nfolds = 2
  nChains = 2
}else{
  nfolds =10
  nChains = 10
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

Y_genus_trunc <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code") 
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
studyDesign <- data.frame(
  sample    = as.factor(rownames(Y_genus)),
  sample_id = as.factor(rownames(Y_genus))
)
rL_spatial <- HmscRandomLevel(sData = xy, longlat = TRUE)
rL_iid     <- HmscRandomLevel(units = studyDesign$sample_id)

library(GGally)
covariate_standa_df %>% 
  select(-c(lon, lat, tmin, Sample_code)) %>%
  GGally::ggpairs()
## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All 2^5 = 32 additive combinations of covariates (null through full), with phylogenetic tree
covariates <- c("ndvi", "precip", "tmax", "tdelta", "wind")

all_cov_sets <- c(
  list(character(0)),
  unlist(lapply(1:5, function(k) combn(covariates, k, simplify = FALSE)), recursive = FALSE)
)

make_gauss_name <- function(covs) {
  if (length(covs) == 0) return("mod_pois_null")
  if (length(covs) == 5) return("mod_pois_full")
  paste0("mod_pois_", paste(covs, collapse = "_"))
}

model_list <- setNames(
  lapply(all_cov_sets, function(covs) {
    f <- if (length(covs) == 0) ~ 1 else {as.formula(paste("~", paste(covs, collapse = " + ")))}
    Hmsc(
      Y           = log(Y_genus_trunc),
      XData       = covariate_standa_df,
      XFormula    = f,
      phyloTree   = genus_tree,
      studyDesign = studyDesign,
      ranLevels   = list(sample = rL_spatial, sample_id = rL_iid),
      distr       = "normal",
      YScale      = TRUE
    )
  }),
  sapply(all_cov_sets, make_gauss_name)
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
    thin = 10
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
  print(hmsc_mod)
  list(hmsc_mod = hmsc_mod, m_converg = m_converg) 
  # , mPredY = mPredY, MF_expl = MF_expl, MF_pred = MF_pred)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model_list %>%
  purrr::imap(\(x, nm) {
    cat(sprintf("\n[Gaussian %d/%d] %s\n", which(names(model_list) == nm), length(model_list), nm))
    compute.hmsc(x, nfold = nfolds, nChains = nChains, samples = 2000, test.run = test_run)
  }) -> list_mods
list_mods$mod_pois_null$hmsc_mod$Y %>% str
if(file.exists("outputs/models/list_mods_genus_gauss_phylo_autocorr.RData")){
  # safe save in case a model list in already present in the folder
  save(list_mods, file = "outputs/models/list_mods_genus_gauss_phylo_autocorr_safesave.RData")
}else{
  save(list_mods, file = "outputs/models/list_mods_genus_gauss_phylo_autocorr.RData")
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Y_genus_bin <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value==0, 0, 1)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")


## ----echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# All 2^5 = 32 additive combinations of covariates (null through full), with phylogenetic tree
make_bin_name <- function(covs) {
  if (length(covs) == 0) return("mod_bin_null")
  if (length(covs) == 5) return("mod_bin_full")
  paste0("mod_bin_", paste(covs, collapse = "_"))
}

model_bin_list <- setNames(
  lapply(all_cov_sets, function(covs) {
    f <- if (length(covs) == 0) ~ 1 else {as.formula(paste("~", paste(covs, collapse = " + ")))}
    Hmsc(
      Y           = Y_genus_bin,
      XData       = covariate_standa_df,
      XFormula    = f,
      phyloTree   = genus_tree,
      studyDesign = studyDesign,
      ranLevels   = list(sample = rL_spatial, sample_id = rL_iid),
      distr       = "probit"
    )
  }),
  sapply(all_cov_sets, make_bin_name)
)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


model_bin_list %>%
  purrr::imap(\(x, nm) {
    cat(sprintf("\n[Probit %d/%d] %s\n", which(names(model_bin_list) == nm), length(model_bin_list), nm))
    compute.hmsc(x, nfold = nfolds, nChains = nChains, samples = 1000, test.run = test_run)
  }) -> list_mods_bin

if(file.exists("outputs/models/list_mods_bin_phylo_autocorr_genus.RData")){
  # safe save in case a model list in already present in the folder
  save(list_mods_bin, file = "outputs/models/list_mods_bin_phylo_autocorr_genus_safesave.RData")
}else{
  save(list_mods_bin, file = "outputs/models/list_mods_bin_phylo_autocorr_genus.RData")
}


