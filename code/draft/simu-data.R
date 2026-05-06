## Simulate Opportunistic Flowering Data that responds to temp and ppt cues

# Load Packages
library(rnaturalearth)
library(sf)
library(ggplot2)
library(terra)
library(dplyr)
library(patchwork)
library(dbscan)
library(mgcv)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/simu-data-fxns.R")

# Load in shape of Nevada
nevada <- ne_states(country = "United States of America", 
                    returnclass = "sf") %>% filter(name == "Nevada")

# Read in preprocessed extracts of Nevada for Daymet Daily Layers 2017-2020
tmin_nev <- terra::rast("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/clim-layers/nevada_tmin_2017_2020_daymet.tif")
tmax_nev <- terra::rast("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/clim-layers/nevada_tmin_2017_2020_daymet.tif")
tmean_nev <- terra::rast("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/clim-layers/nevada_tmean_2017_2020_daymet.tif")
ppt_nev <- terra::rast("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/clim-layers/nevada_prcp_2017_2020_daymet.tif")

# Create Simu Species Range
nevada_bbox <- st_bbox(nevada)

species_range <- st_polygon(list(cbind(
  c(-120, -119, -117, -115, -116, -118, -120), # lons
  c(36, 39, 41, 39.5, 37, 35.5, 36)))) %>% #lats
  st_sfc(crs = st_crs(nevada)) %>% 
  st_sf() %>% 
  st_intersection(nevada) # clip to nev boundary

ggplot() + 
  geom_sf(data = nevada, mapping = aes(), fill = "lightgray") + 
  geom_sf(data = species_range, mapping = aes(), fill = "darkgreen", alpha = 0.5) + 
  theme_bw() + 
  ggtitle("Species S Distribution in Nevada")

# Build True Flowering Response (Data Generating Process)

## extract env thresholds given the simu species range
# reproj species range to match clim data
species_range_clim <- st_transform(species_range, crs = crs(tmean_nev))

# create a point for each 1km cell center within the species range (res of daymet clim layer)
tmean_species <- crop(tmean_nev, vect(species_range_clim), mask = TRUE)
ppt_species <- crop(ppt_nev, vect(species_range_clim), mask = TRUE)

# extract cell centers as pts 
cells <- cells(tmean_species[[1]]) # get cell numbers
cell_coords <- xyFromCell(tmean_species[[1]], cells)
cell_pts <- st_as_sf(data.frame(cell_coords), 
                     coords = c("x", "y"), 
                     crs = crs(tmean_nev))

years <- 2017:2020
days_per_year <- c(365,365,365,366)
spring_doy <- 0:120

# extract using a function that grabs percentiles and slope of temp 
env_thresholds <- extract_env_thresholds(tmean_species, 
                                         ppt_species, 
                                         cell_pts, 
                                         years, 
                                         days_per_year, 
                                         spring_doy, 
                                         ppt_low_pct = 0.10, 
                                         ppt_high_pct = 0.90, 
                                         doy_range = c(60, 120), # absolute bounds of the species (can flowering either beginning of March - End of April)
                                         plot = TRUE)

## parameters
params <- list(
  intercept = env_thresholds$intercept,      # flowering DOY when temp = 0°C
  temp_slope = env_thresholds$temp_slope,    # days earlier per 1°C increase
  sigma_spatial = 2, # variation between sites
  sigma_individual = 5, # variation within populations
  ppt_min = env_thresholds$ppt_min, # min ppt mm needed for flower on temp schedule
  ppt_max = env_thresholds$ppt_max, # ppt above which penalty on flowering begins
  ppt_max_penalty = 60 # maximum DOY penalty for excess ppt
)

## visualize response 
# set up ppt responses
ppt_scenarios <- c() # init a empty vector
ppt_scenarios[1] <- params$ppt_min  # make a drought
ppt_scenarios[2] <- (params$ppt_min + params$ppt_max) / 2 # find optimum
ppt_scenarios[3] <- params$ppt_max # excess rain
ppt_labels <- c(paste0("drought(", ppt_scenarios[1], "mm)"), paste0("optimal(", ppt_scenarios[2], "mm)"), paste0("excess(", ppt_scenarios[3], "mm)"))

# set up temp responses
temp_range <- seq(0, 25, by = 0.5)

# apply through
response_distributions <- expand.grid(temp_range = temp_range, ppt = ppt_scenarios) %>%
  mutate(ppt_label = factor(ppt, levels = ppt_scenarios, labels = ppt_labels)) %>%
  mutate(
    distribution = purrr::pmap(list(temp_range, ppt), function(t, p)
      sim_pop_flowering(t, p, n_individuals = 500, par = params)),
    onset_10  = purrr::map_dbl(distribution, ~quantile(.x, 0.10)),
    median_50 = purrr::map_dbl(distribution, ~quantile(.x, 0.50)),
    offset_90 = purrr::map_dbl(distribution, ~quantile(.x, 0.90))
  )

# make a plot df
response_long <- response_distributions %>%
  tidyr::pivot_longer(cols = c(onset_10, median_50, offset_90),
                      names_to = "phenometric",
                      values_to = "doy") %>%
  mutate(phenometric = factor(phenometric,
                              levels = c("onset_10", "median_50", "offset_90"),
                              labels = c("Onset (10%)", "Median (50%)", "Offset (90%)")))

# plot

ggplot(response_long, aes(x = temp_range, y = doy, color = phenometric, linetype = phenometric)) +
  geom_ribbon(data = response_distributions,
              aes(x = temp_range, ymin = onset_10, ymax = offset_90),
              fill = "yellow", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(size = 1) +
  scale_color_manual(values = c("darkred", "darkorange", "steelblue")) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed")) +
  facet_wrap(~ppt_label) +
  labs(color = "Phenometric", linetype = "Phenometric",
       x = "Mean Spring Temperature (°C)", y = "Flowering DOY",
       title = "Flowering Response to Temperature Across Precipitation Scenarios") +
  theme_bw()

# Simulate Population Phenology Across the Species Distribution
## calc cumulative year positions 
year_starts <- c(1, cumsum(days_per_year[-length(days_per_year)]) + 1)
names(year_starts) <- years
## generate true occurrences for all years
cells_true_flowering_all <- list()
for(i in seq_along(years)){
  # get year and when year starts in raster
  year <- years[i]
  year_start <- year_starts[i]
  # get spring layers for this year
  spring_layers_year <- (year_start - 1) + spring_doy
  # calc spring climate for this year
  tmean_spring <- mean(tmean_species[[spring_layers_year]])
  ppt_spring <- sum(ppt_species[[spring_layers_year]])
  # extract clim vals 
  climate_data <- terra::extract(c(tmean_spring, ppt_spring), vect(cell_pts))
  # generate cell phenometrics per year
  cells_true_flowering__dists_year <- cell_pts %>% 
    st_drop_geometry() %>% 
    bind_cols(st_coordinates(cell_pts)) %>% 
    mutate(
      tmean_spring = climate_data$mean,
      ppt_spring = climate_data$sum, 
      year = year
    ) %>% 
    rowwise() %>% 
    mutate(
      cell_flowering_distribution = list(sim_pop_flowering(tmean_spring, ppt = ppt_spring, n_individuals = 100, par = params)
      ), 
      cell_flowering_onset = quantile(cell_flowering_distribution, 0.10), 
      cell_flowering_median = quantile(cell_flowering_distribution, 0.50), 
      cell_flowering_offset = quantile(cell_flowering_distribution, 0.90)) %>% 
    ungroup() 
  cells_true_flowering_all[[as.character(year)]] <- cells_true_flowering__dists_year
  
}

# combine all years
cells_true_flowering <- bind_rows(cells_true_flowering_all)
cells_true_flowering_sf  <- st_as_sf(cells_true_flowering, 
                                     coords = c("X", "Y"), 
                                     crs = crs(tmean_nev))


true_flowering_plot_all <- ggplot() + 
  geom_sf(data = nevada, fill = "gray90", color = "black") + 
  geom_sf(data = species_range, fill = NA, color = "black", size = 1) + 
  geom_sf(data = cells_true_flowering_sf , mapping = aes(color = cell_flowering_onset), size = 0.5) + 
  scale_color_viridis_c(option = "magma", direction = -1, name = "Flowering\nDOY") + 
  facet_wrap(~year) +
  ggtitle("True Flowering Onset (10%) Across Species Range") +
  theme_bw()

print(true_flowering_plot_all)

# Sampling Layer
## start simple, lets evenly sample cells across space and across flowering distributions
n_cells_per_year <- 500
obs_per_cell_range <- c(50,100)

## sample from true pops 
observed_occs <- cells_true_flowering %>% 
  group_by(year) %>% 
  # randomly sample cells 
  slice_sample(n = n_cells_per_year, replace = FALSE) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate( # randomly determine how many obs within cell
    n_obs_at_cell = sample(obs_per_cell_range[1]:obs_per_cell_range[2], 1), 
    # sample flowering from distribution
    observed_flowering_times = list(sample(cell_flowering_distribution, n_obs_at_cell))
  ) %>% 
  ungroup() %>% 
  tidyr::unnest(observed_flowering_times) %>%  # each row will be an obs
  rename(flowering_doy_obs = observed_flowering_times)

## visualize sampling 
# create a spatial plot: 
observed_occs_sf <- st_as_sf(observed_occs, 
                             coords = c("X", "Y"), 
                             crs = crs(tmean_nev))

observed_occs_sf <- observed_occs_sf %>% 
  distinct(geometry, year, n_obs_at_cell)

obs_flowering_plot <- ggplot() + 
  geom_sf(data = nevada, fill = "gray90", color = "black") + 
  geom_sf(data = species_range, fill = NA, color = "black", size = 1) + 
  geom_sf(data = observed_occs_sf, mapping = aes(color = n_obs_at_cell), size = 0.5) + 
  scale_color_viridis_c(option = "magma", direction = +1, name = "Number of Obs") + 
  ggtitle("Flowering Observation Simu") + 
  facet_wrap(~year)

obs_flowering_plot 

# Group Clusters for Building PhenoEstimates
## use spatial clustering
obs_spatial_clusters <- cluster_spatial(observed_occs, 
                                        distance_km = 25, 
                                        min_obs = 7,     
                                        min_cells = 2)
## use enviromental variable clustering
obs_env_clusters <- cluster_env(observed_occs, 
                                env_range_h = 2, 
                                min_obs = 7, 
                                min_cells = 2)

## visualize clustering 
# Visualize clusters in space
obs_spatial_clusters_sf <- st_as_sf(obs_spatial_clusters, 
                                    coords = c("X", "Y"), 
                                    crs = crs(tmean_nev))
obs_env_clusters_sf <- st_as_sf(obs_env_clusters, 
                                coords = c("X", "Y"), 
                                crs = crs(tmean_nev))

ggplot() +
  geom_sf(data = nevada, fill = "gray90") +
  geom_sf(data = obs_spatial_clusters_sf, aes(color = factor(cluster_id)), size = 1) +
  facet_wrap(~year) +
  labs(color = "Cluster", title = "Spatial Clusters by Year") +
  theme_minimal()

ggplot() +
  geom_sf(data = nevada, fill = "gray90") +
  geom_sf(data = obs_env_clusters_sf, aes(color = factor(cluster_id)), size = 1) +
  facet_wrap(~year) +
  labs(color = "Cluster", title = "Temp Env Clusters by Year") +
  theme_minimal()

# Use Phenesse to Calculate Onset for the observed data
cluster_spatial_onset <- obs_spatial_clusters %>%
  group_by(cluster_full_id, year, cluster_compactness) %>%
  summarise(
    n_obs = n(),
    onset_estimate = phenesse::quantile_ci(observations = flowering_doy_obs, 
                                           percentile = 0.10, 
                                           bootstraps = 250)$estimate,
    # Also get cluster-level environmental predictors for GAM
    tmean_spring_mean = mean(tmean_spring),
    ppt_spring_mean = mean(ppt_spring),
    X_centroid = mean(X),
    Y_centroid = mean(Y),
    .groups = "drop"
  )  %>% 
  mutate(weight = cluster_compactness / mean(cluster_compactness)) # normalized

cluster_env_onset <- obs_env_clusters %>%
  group_by(cluster_full_id, year, cluster_temp_sd) %>%
  summarise(
    n_obs = n(),
    onset_estimate = phenesse::quantile_ci(observations = flowering_doy_obs, 
                                           percentile = 0.10, 
                                           bootstraps = 250)$estimate,
    # Also get cluster-level environmental predictors for GAM
    tmean_spring_mean = mean(tmean_spring),
    ppt_spring_mean = mean(ppt_spring),
    X_centroid = mean(X),
    Y_centroid = mean(Y),
    .groups = "drop"
  )  %>% 
  mutate(weight = mean(cluster_temp_sd)/cluster_temp_sd)


# Model Onset using a GAM 
cluster_spatial_onset$cluster_full_id <- as.factor(cluster_spatial_onset$cluster_full_id)
cluster_env_onset$cluster_full_id <- as.factor(cluster_env_onset$cluster_full_id)
# run GAM, overkill here.
gam_spatial <- gam(onset_estimate ~ s(tmean_spring_mean, k = 5) +
                     s(ppt_spring_mean, k = 10, bs = "ad") +
                     s(cluster_full_id, bs = "re"), 
                   data = cluster_spatial_onset, 
                   weights = weight)

gam_env <- gam(onset_estimate ~ s(tmean_spring_mean, k = 7) + 
                 s(ppt_spring_mean, k = 10, bs = "ad") +
                 s(cluster_full_id, bs = "re"), 
               data = cluster_env_onset, 
               weights = weight)

# predict 
cells_true_flowering$predicted_onset_spatial <- predict(gam_spatial, newdata = data.frame(
  tmean_spring_mean = cells_true_flowering$tmean_spring,
  ppt_spring_mean = cells_true_flowering$ppt_spring), 
  exclude = "s(cluster_full_id)",  # Exclude random effect
  newdata.guaranteed = TRUE)

cells_true_flowering$predicted_onset_env <- predict(gam_env, newdata = data.frame(
  tmean_spring_mean = cells_true_flowering$tmean_spring, 
  ppt_spring_mean = cells_true_flowering$ppt_spring), 
  exclude = "s(cluster_full_id)",  # Exclude random effect
  newdata.guaranteed = TRUE)


# Evaluation of Model Fit
## How did we do with fitting response? 
# make a response pred grid
pred_grid <- expand.grid(
  tmean_spring_mean = seq(min(cells_true_flowering$tmean_spring, na.rm = TRUE),
                          max(cells_true_flowering$tmean_spring, na.rm = TRUE),
                          length.out = 100),
  ppt_spring_mean   = ppt_scenarios  # same ppt scenarios as response curve viz
) %>%
  mutate(ppt_label = factor(ppt_spring_mean,
                            levels = ppt_scenarios,
                            labels = ppt_labels))

# add in fitted predicted data
pred_grid$predicted_spatial <- predict(gam_spatial,
                                       newdata = pred_grid,
                                       exclude = "s(cluster_full_id)",
                                       newdata.guaranteed = TRUE,
                                       se.fit = FALSE)
pred_grid$predicted_env <- predict(gam_env,
                                   newdata = pred_grid,
                                   exclude = "s(cluster_full_id)",
                                   newdata.guaranteed = TRUE,
                                   se.fit = FALSE)
## Build the True Response Curve given these same pred_grid vars
pred_grid <- pred_grid %>%
  rowwise() %>%
  mutate(
    ppt_mid = (params$ppt_min + params$ppt_max) / 2,
    ppt_penalty = if (ppt_spring_mean < params$ppt_min) {
      params$ppt_max_penalty
    } else if (ppt_spring_mean > params$ppt_max) {
      params$ppt_max_penalty
    } else {
      params$ppt_max_penalty * ((ppt_spring_mean - ppt_mid) /
                                  (ppt_mid - params$ppt_min))^2
    },
    mu = params$intercept + params$temp_slope * tmean_spring_mean + ppt_penalty,
    true_onset = qnorm(0.10, mean = mu, sd = params$sigma_individual)
  ) %>%
  ungroup()

## Make plot df
pred_long <- pred_grid %>%
  tidyr::pivot_longer(cols = c(predicted_spatial, predicted_env, true_onset),
                      names_to = "source",
                      values_to = "onset") %>%
  mutate(source = recode(source,
                         "predicted_spatial" = "GAM (Spatial)",
                         "predicted_env"     = "GAM (Environmental)",
                         "true_onset"        = "True Response"))

gam_vs_true <- ggplot(pred_long, aes(x = tmean_spring_mean, y = onset,
                                     color = source, linetype = source)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("GAM (Spatial)"       = "steelblue",
                                "GAM (Environmental)" = "darkorange",
                                "True Response"            = "black")) +
  scale_linetype_manual(values = c("GAM (Spatial)"       = "dashed",
                                   "GAM (Environmental)" = "dashed",
                                   "True Response"            = "solid")) +
  facet_wrap(~ppt_label) +
  labs(x = "Mean Spring Temperature (°C)", y = "Flowering Onset (DOY)",
       title = "GAM Estimated vs True Response Curve",
       color = "Source", linetype = "Source") +
  theme_bw()

print(gam_vs_true)

# Check Interpolation vs Extrapolation Error

## MESS
predictor_vars <- c("tmean_spring", "ppt_spring") 

# set-up training data (empirical data) for the models
training_cells_spatial <- obs_spatial_clusters %>% 
  distinct(X, Y, .keep_all = TRUE)
training_cells_env <- obs_env_clusters %>% 
  distinct(X, Y, .keep_all = TRUE)

# determine bounds and dists
for(v in 1:length(predictor_vars)){
  var_name <- predictor_vars[v]
  
  # Spatial
  low  <- min(training_cells_spatial[[var_name]], na.rm = TRUE)
  high <- max(training_cells_spatial[[var_name]], na.rm = TRUE)
  cells_true_flowering[[paste0("mess_spatial_", var_name)]] <-
    ifelse(cells_true_flowering[[var_name]] < low | cells_true_flowering[[var_name]] > high, 0, 1)
  cells_true_flowering[[paste0("mess_dist_spatial_", var_name)]] <-
    ifelse(cells_true_flowering[[paste0("mess_spatial_", var_name)]] == 1, 0,
           pmin(abs(cells_true_flowering[[var_name]] - low),
                abs(cells_true_flowering[[var_name]] - high)))
  
  # Environmental
  low  <- min(training_cells_env[[var_name]], na.rm = TRUE)
  high <- max(training_cells_env[[var_name]], na.rm = TRUE)
  cells_true_flowering[[paste0("mess_env_", var_name)]] <-
    ifelse(cells_true_flowering[[var_name]] < low | cells_true_flowering[[var_name]] > high, 0, 1)
  cells_true_flowering[[paste0("mess_dist_env_", var_name)]] <-
    ifelse(cells_true_flowering[[paste0("mess_env_", var_name)]] == 1, 0,
           pmin(abs(cells_true_flowering[[var_name]] - low),
                abs(cells_true_flowering[[var_name]] - high)))
}

## Combine MESS interpolations only if in range of both vars
cells_true_flowering <- cells_true_flowering %>%
  mutate(
    # 1 only if interpolating in both temp and ppt
    mess_spatial_combined = as.integer(mess_spatial_tmean_spring == 1 & mess_spatial_ppt_spring == 1),
    mess_env_combined     = as.integer(mess_env_tmean_spring     == 1 & mess_env_ppt_spring     == 1),
    # total extrapolation distance as euclidean distance across both vars
    mess_dist_spatial_combined = sqrt(mess_dist_spatial_tmean_spring^2 + mess_dist_spatial_ppt_spring^2),
    mess_dist_env_combined     = sqrt(mess_dist_env_tmean_spring^2     + mess_dist_env_ppt_spring^2)
  )


## calc error
cells_true_flowering <- cells_true_flowering %>%
  mutate(
    error_spatial = predicted_onset_spatial - cell_flowering_onset,
    error_env     = predicted_onset_env     - cell_flowering_onset,
    abs_error_spatial = abs(error_spatial),
    abs_error_env     = abs(error_env),
    # Use combined status since both vars drive the model
    interp_status_spatial = ifelse(mess_spatial_combined == 1, "Interpolation", "Extrapolation"),
    interp_status_env     = ifelse(mess_env_combined     == 1, "Interpolation", "Extrapolation")
  )

## Absolute error by interpolation status (boxplot)
error_long <- cells_true_flowering %>%
  dplyr::select(year, abs_error_spatial, abs_error_env,
                interp_status_spatial, interp_status_env) %>%
  tidyr::pivot_longer(
    cols = c(abs_error_spatial, abs_error_env),
    names_to = "model",
    values_to = "abs_error"
  ) %>%
  mutate(
    interp_status = ifelse(model == "abs_error_spatial", interp_status_spatial, interp_status_env),
    model = recode(model,
                   "abs_error_spatial" = "Spatial Clustering",
                   "abs_error_env"     = "Environmental Clustering")
  )

error_box_plot <- ggplot(error_long, aes(x = interp_status, y = abs_error, fill = interp_status)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = c("Interpolation" = "darkgreen", "Extrapolation" = "darkred")) +
  facet_grid(model ~ year) +
  labs(x = NULL, y = "Absolute Error (DOY)",
       title = "Prediction Error: Interpolation vs. Extrapolation (Combined MESS)",
       fill = "Status") +
  theme_bw()
print(error_box_plot)

## Make MESS Plot
mess_sf <- st_as_sf(cells_true_flowering, coords = c("X","Y"), crs = crs(tmean_nev))

plot_mess(mess_sf, "mess_dist_spatial_combined",
          "MESS: Combined (Spatial Clustering)",
          "/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/MESS-spatial-combined.png")
