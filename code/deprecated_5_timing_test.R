## -----------------------------------------------------------------------
## 5_timing_test.R
## Timing benchmark for Hmsc Gaussian abundance models
##
## Crosses:
##   - model formula : null (~ 1) vs full (~ ndvi + precip + tmax + tdelta + wind)
##   - phylogenetic tree : with vs without
##
## Varies:
##   - MCMC samples : 5 levels of increasing intensity
##
## Random effects: spatial (sData) + iid sample-level (overdispersion fix)
## Measures elapsed time for: fit, explanatory prediction, CV prediction, map prediction
## Produces a ggplot summary saved to outputs/timing/
## -----------------------------------------------------------------------

library(tidyverse)
library(Hmsc)
library(ape)
library(raster)
library(vegan)
library(sf)

## --- Fixed benchmark parameters -----------------------------------------
nChains      <- 2
thin         <- 1
transient    <- 50
nfolds       <- 5
samples_grid <- c(20, 50, 100, 200, 500)

## --- Load data (identical to script 5) ----------------------------------
covariate_standa_df <- read.table(
  "outputs/data_for_models/covariate_standa_df.txt",
  stringsAsFactors = TRUE
)
Y_genus    <- read.table("outputs/data_for_models/Y_genus.txt")
genus_tree <- ape::read.tree("outputs/data_for_models/genus_tree.tree")

## Rename genera (strip g__ prefix, standardize unknowns)
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

Y_genus <- Y_genus[, genus_tree$tip.label]

## Zeros → NA for log-normal model
Y_genus_trunc <- Y_genus %>%
  as.data.frame() %>%
  rownames_to_column("soil_code") %>%
  pivot_longer(-soil_code) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  column_to_rownames("soil_code")

## Fix duplicate GPS coordinates (ROB_F6_s / ROB_F7_s)
dup <- duplicated(covariate_standa_df[, c("lon", "lat")])
covariate_standa_df[dup, c("lon", "lat")] <-
  covariate_standa_df[dup, c("lon", "lat")] + 0.0001

xy <- covariate_standa_df[, c("lon", "lat", "Sample_code")] %>%
  remove_rownames() %>%
  column_to_rownames("Sample_code")

studyDesign <- data.frame(
  sample    = as.factor(rownames(Y_genus)),
  sample_id = as.factor(rownames(Y_genus))
)
rL_spatial <- HmscRandomLevel(sData = xy, longlat = TRUE)
rL_iid     <- HmscRandomLevel(units = studyDesign$sample_id)

## --- Prepare spatial prediction grid (mirrors section 11 of script 6) ---
ndvi_uganda              <- raster("outputs/raster/ndvi_uganda_lower_res.tif")
precip_uganda_resampled  <- raster("outputs/raster/precip_2023_lower_res.tif")
tmax_uganda_resampled    <- raster("outputs/raster/tmax_2023_lower_res.tif")
tmin_uganda_resampled    <- raster("outputs/raster/tmin_2023_lower_res.tif")
delta_t_uganda           <- tmax_uganda_resampled - tmin_uganda_resampled
wind_2023                <- raster("outputs/raster/wind_2023.tif")

data_spatial_df <- rasterToPoints(ndvi_uganda) %>%
  as.data.frame() %>% rename(ndvi = "ndvi_uganda_lower_res") %>%
  left_join(rasterToPoints(precip_uganda_resampled) %>%
              as.data.frame() %>% rename(precip = "precip_2023_lower_res"),
            by = c("x", "y")) %>%
  left_join(rasterToPoints(tmax_uganda_resampled) %>%
              as.data.frame() %>% rename(tmax = "tmax_2023_lower_res"),
            by = c("x", "y")) %>%
  left_join(rasterToPoints(delta_t_uganda) %>%
              as.data.frame() %>% rename(tdelta = "layer"),
            by = c("x", "y")) %>%
  left_join(rasterToPoints(wind_2023) %>%
              as.data.frame() %>% rename(wind = "wind_2023"),
            by = c("x", "y")) %>%
  na.omit() %>%
  mutate(across(c(ndvi, precip, tmax, tdelta, wind),
                ~ decostand(.x, method = "standardize")[, 1]))

## Filter to limited environmental gradient (within training data range)
data_spatial_limited_df <- data_spatial_df %>%
  filter(
    max(covariate_standa_df$ndvi)   > ndvi   & min(covariate_standa_df$ndvi)   < ndvi,
    max(covariate_standa_df$precip) > precip & min(covariate_standa_df$precip) < precip,
    max(covariate_standa_df$tmax)   > tmax   & min(covariate_standa_df$tmax)   < tmax,
    max(covariate_standa_df$tdelta) > tdelta & min(covariate_standa_df$tdelta) < tdelta,
    max(covariate_standa_df$wind)   > wind   & min(covariate_standa_df$wind)   < wind
  )

xy.grid_limited    <- data_spatial_limited_df %>% dplyr::select(x, y) %>% as.matrix()
XData.grid_limited <- data_spatial_limited_df %>% dplyr::select(-c(x, y))

cat(sprintf("Spatial prediction grid: %d pixels\n", nrow(xy.grid_limited)))

## --- Model configuration grid -------------------------------------------
model_configs <- list(
  # list(label = "null_nophylo", formula = ~ 1,                                    phylo = FALSE),
  # list(label = "null_phylo",   formula = ~ 1,                                    phylo = TRUE),
  # list(label = "full_nophylo", formula = ~ ndvi + precip + tmax + tdelta + wind, phylo = FALSE),
  list(label = "full_phylo",   formula = ~ ndvi + precip + tmax + tdelta + wind, phylo = TRUE)
)

## --- Timing loop --------------------------------------------------------
results <- vector("list", length(model_configs) * length(samples_grid))
idx     <- 1L

for (cfg in model_configs) {
  for (s in samples_grid) {

    cat(sprintf(
      "\n=== %-16s | samples = %d ===\n", cfg$label, s
    ))

    ## Build unfitted Hmsc model
    hmsc_mod <- Hmsc(
      Y           = log(Y_genus_trunc),
      XData       = covariate_standa_df,
      XFormula    = cfg$formula,
      phyloTree   = if (cfg$phylo) genus_tree else NULL,
      studyDesign = studyDesign,
      ranLevels   = list(sample = rL_spatial, sample_id = rL_iid),
      distr       = "normal",
      YScale      = TRUE
    )

    ## 1. Fit ---------------------------------------------------------------
    t_fit <- system.time({
      m_fit <- sampleMcmc(
        hmsc_mod,
        thin      = thin,
        samples   = s,
        transient = transient,
        nChains   = nChains,
        verbose   = 0,
        nParallel = nChains
      )
    })

    ## Partition for cross-validation
    set.seed(42)
    partition <- createPartition(m_fit, nfolds = nfolds, column = "sample")

    ## 2. Explanatory prediction -------------------------------------------
    t_expl <- system.time({
      expl_power <- computePredictedValues(
        m_fit, thin = thin, nParallel = nChains
      )
    })

    ## 3. Map prediction (spatial gradient on new pixels) ---------------------
    ## Uses the limited-gradient raster grid (~2000 new locations)
    model_covariates <- setdiff(colnames(m_fit$X), "(Intercept)")
    XData.grid_sub   <- XData.grid_limited[, model_covariates, drop = FALSE]

    t_map <- system.time({
      gradient <- prepareGradient(
        m_fit,
        XDataNew = XData.grid_sub,
        sDataNew = list(sample = xy.grid_limited)
      )
      predY <- predict(m_fit, Gradient = gradient, predictEtaMean = TRUE)
    })

    cat(sprintf(
      "  fit: %.1f s | expl: %.1f s | map: %.1f s\n",
      t_fit["elapsed"], t_expl["elapsed"], t_map["elapsed"]
    ))

    ## 4. CV prediction -----------------------------------------------------
    # t_cv <- system.time({
    #   pred_power <- computePredictedValues(
    #     m_fit, partition = partition, thin = thin, nParallel = nChains
    #   )
    # })

    results[[idx]] <- data.frame(
      label    = cfg$label,
      model    = if (grepl("null", cfg$label)) "null" else "full",
      phylo    = cfg$phylo,
      samples  = s,
      fit_sec  = t_fit["elapsed"],
      expl_sec = t_expl["elapsed"],
      map_sec  = t_map["elapsed"],
      # cv_sec   = t_cv["elapsed"],
      row.names = NULL
    )
    idx <- idx + 1L
  }
}

timing_df <- bind_rows(results)

## --- Save results -------------------------------------------------------
dir.create("outputs/timing", showWarnings = FALSE, recursive = TRUE)
write.csv(timing_df, "outputs/timing/timing_results.csv", row.names = FALSE)
cat("\nResults saved to outputs/timing/timing_results.csv\n")

## --- ggplot summary -----------------------------------------------------
plot_df <- timing_df %>%
  pivot_longer(
    cols      = c(fit_sec, expl_sec, map_sec),
    names_to  = "step",
    values_to = "seconds"
  ) %>%
  mutate(
    step = factor(
      step,
      levels = c("fit_sec", "expl_sec", "map_sec"),
      labels = c("Fit", "Explanatory", "Map prediction")
    ),
    model       = factor(model, levels = c("null", "full"), labels = c("Null", "Full")),
    phylo_label = if_else(phylo, "With phylo", "Without phylo")
  )

p <- ggplot(plot_df, aes(
  x        = samples,
  y        = seconds,
  color    = phylo_label,
  linetype = model,
  shape    = model
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ step, scales = "free_y", nrow = 1) +
  scale_color_manual(
    values = c("With phylo" = "#2196F3", "Without phylo" = "#E53935")
  ) +
  scale_x_continuous(breaks = samples_grid) +
  labs(
    title    = "Hmsc Gaussian abundance model — timing benchmark",
    subtitle = sprintf(
      "nChains = %d, thin = %d, transient = %d, nfolds = %d",
      nChains, thin, transient, nfolds
    ),
    x        = "MCMC samples",
    y        = "Elapsed time (seconds)",
    color    = "Phylogenetic tree",
    linetype = "Model",
    shape    = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )

print(p)
ggsave("outputs/timing/timing_plot.png", p, width = 11, height = 5, dpi = 150)
cat("Plot saved to outputs/timing/timing_plot.png\n")
