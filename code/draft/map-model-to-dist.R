# Load Libraries 
library(data.table)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/proxy-baseline-fxns.R")

bee_occs <- fread("/blue/soltis/millerjared/beeBDC/Bee-data/OutputData/05_cleaned_database.csv")
# additional data on Anthophora hololeuca provided by Michael Orr
extra_hololeuca_data <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/raw/a-hololeuca-from-Orr.csv")
phenovision_flowering_data <- fread("/blue/guralnick/millerjared/phenovision_data/phenovision_flowers_03_15_2026.csv")
cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")


## Load in North American Taxonomy
na_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")

## Pull target taxa, order with the accepted parent name in lead
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Olneya tesota") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

target_taxa <- "Anthophora hololeuca"

described_phen_months <- fread("/blue/soltis/millerjared/Legume-Specialist-Occupancy/data/raw/hololeuca-and-host-pheno.csv")

taxon_phen_hist(target_taxa, bees_or_plants = "bees", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "bees", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 21, # look at histogram when making this decision
                                      max_doy = 234 )

model_out$plots$hopkins
model_out$plots$enviroment
df <- model_out$data

## Pull in dist model, limit to presence range only, apply modified hopkins linear model to P pixels 
anth_holo <- rast("/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif")
anth_holo[anth_holo == 0] <- NA # remove non-presence pts

# reproj to the native wgs84 coordinate system used to build hopkins expectations
anth_holo_reproj <- project(anth_holo, "EPSG:4326", method = "near")

# make into xy df 
anth_holo_df <- as.data.frame(anth_holo_reproj, xy = TRUE, na.rm = TRUE)
anth_holo_df <- anth_holo_df %>% 
  rename(longitude = x, latitude = y) %>% 
  mutate(rowid = 1:n()) %>% 
  as.data.table()

# convert this into a spatVect to grab elevation vals, reproj to match elev tif
na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
anth_holo_vect <- vect(anth_holo_df, 
                       geom = c("longitude", "latitude"),
                       crs = "EPSG:4326")
anth_holo_vect_reproj <- project(anth_holo_vect, crs(na_elev))
elevation_vals <- terra::extract(na_elev, anth_holo_vect_reproj)

anth_holo_df[, elevation_m := elevation_vals[,2]]

## Using this modified dataframe, recreate hopkins based on the anchor from the empirical data
# bring in min vals 
min_empirical_lat <- min(df$latitude)
min_empirical_elev <- min(df$elevation_m)
# these will be zero, we'll just adjust so we're going into negative space
# create expected # of delay days based on hopkins
anth_holo_df[, hopkins_delay_days :=
            ((latitude - min_lat)* 4) + # latitudnal effect
            ((elevation_m - min_elev))/120 * 4] # elevation effect

# Create bins for phenology data based on Hopkins delay by intervals of 4
anth_holo_df[, hopkins_bin := cut(hopkins_delay_days, 
                               breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
# additionally, we want to know what pixels are outside of the range of empirical data training
min_empirical_hopkins <- min(df$hopkins_delay_days)
max_empirical_hopkins <- max(df$hopkins_delay_days)
# binary
anth_holo_df[, outside_empirical := hopkins_delay_days < min_empirical_hopkins |
               hopkins_delay_days > max_empirical_hopkins]
# numeric dist
anth_holo_df[, dist_from_empirical := fifelse(
  hopkins_delay_days < min_empirical_hopkins, 
  min_empirical_hopkins - hopkins_delay_days, # below range
  fifelse(
    hopkins_delay_days > max_empirical_hopkins, 
    hopkins_delay_days - max_empirical_hopkins, # above range
    0 #if within range
  )
)]

# extract midpoint of each hopkins bin to use as predictor
anth_holo_df <- anth_holo_df %>%
  mutate(
    # Extract numeric values from bin intervals
    bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
    bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
    hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
  )
 # Some pixels lack elevation, which is why we get this warning. Probably try and see if we can fix this at a later time. 
anth_holo_df <- anth_holo_df %>% filter(!is.na(hopkins_bin))
# apply hopkins linear model 
hopkins_linear <- model_out$hopkins_model
anth_holo_predict <- predict(hopkins_linear, newdata = anth_holo_df)
anth_holo_df[, hopkins_predict_onset := anth_holo_predict]

# Now plot!
basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
basemap <- basemap %>% filter(name %in% c("California", "Nevada", "Arizona", "Sonora", "Baja California", "Baja California Sur"))
basemap_reproj <- basemap %>% st_transform(st_crs(anth_holo_reproj))
ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = hopkins_predict_onset)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Hopkins Onset Prediction", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = dist_from_empirical)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Hopkins Empirical Training Distance Eval", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

## Now lets do this with Env space
anth_holo <- rast("/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif")
anth_holo[anth_holo == 0] <- NA # remove non-presence pts

# reproj to the native wgs84 coordinate system used to build hopkins expectations
anth_holo_reproj <- project(anth_holo, "EPSG:4326", method = "near")

# make into xy df 
anth_holo_df <- as.data.frame(anth_holo_reproj, xy = TRUE, na.rm = TRUE)
anth_holo_df <- anth_holo_df %>% 
  rename(longitude = x, latitude = y) %>% 
  mutate(rowid = 1:n()) %>% 
  as.data.table()

## convert this into a spatVect to grab elevation & clim vars, reproj to match elev and clim tifs
# load North America Elev Raster & Chelsa annual clim vars
na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
avg_temp <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio1_1981-2010_V.2.1.tif")
temp_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio4_1981-2010_V.2.1.tif")
tmax <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio5_1981-2010_V.2.1.tif")
tmin <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio6_1981-2010_V.2.1.tif")
annual_precip <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio12_1981-2010_V.2.1.tif")
precip_of_wettest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio13_1981-2010_V.2.1.tif")
precip_of_driest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio14_1981-2010_V.2.1.tif")
precip_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio15_1981-2010_V.2.1.tif")
precip_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio15_1981-2010_V.2.1.tif")
avg_monthly_precip_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio16_1981-2010_V.2.1.tif")
avg_monthly_precip_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio17_1981-2010_V.2.1.tif")
avg_monthly_precip_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio18_1981-2010_V.2.1.tif")
avg_monthly_precip_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio19_1981-2010_V.2.1.tif")
anth_holo_vect <- vect(anth_holo_df, 
                       geom = c("longitude", "latitude"),
                       crs = "EPSG:4326")
anth_holo_vect_reproj1 <- project(anth_holo_vect, crs(na_elev))
anth_holo_vect_reproj2 <- project(anth_holo_vect, crs(avg_temp))

# Stack all climate layers
climate_stack <- c(avg_temp, temp_seas, tmax, tmin, annual_precip, 
                   precip_of_wettest_month, precip_of_driest_month, 
                   precip_seas, avg_monthly_precip_wettest_q, 
                   avg_monthly_precip_driest_q, avg_monthly_precip_warmest_q, 
                   avg_monthly_precip_coldest_q)
# make names clearer
names(climate_stack) <- c("avg_temp", "temp_seas", "tmax", "tmin", 
                          "annual_precip", "precip_wettest_month", 
                          "precip_driest_month", "precip_seas", 
                          "precip_wettest_q", "precip_driest_q", 
                          "precip_warmest_q", "precip_coldest_q")

# Extract all raster values
elevation_vals <- terra::extract(na_elev, anth_holo_vect_reproj1)
climate_vals <- terra::extract(climate_stack, anth_holo_vect_reproj2)

# Add all values to these data
anth_holo_df[, elevation_m := elevation_vals[,2]]
anth_holo_df[, names(climate_vals)[-1] := climate_vals[,-1]]

# Make hopkins bins again as a available predictor 
## Using this modified dataframe, recreate hopkins based on the anchor from the empirical data
df <- model_out$data
# bring in min vals 
min_empirical_lat <- min(df$latitude)
min_empirical_elev <- min(df$elevation_m)
# these will be zero, we'll just adjust so we're going into negative space
# create expected # of delay days based on hopkins
anth_holo_df[, hopkins_delay_days :=
               ((latitude - min_lat)* 4) + # latitudnal effect
               ((elevation_m - min_elev))/120 * 4] # elevation effect

# Create bins for phenology data based on Hopkins delay by intervals of 4
anth_holo_df[, hopkins_bin := cut(hopkins_delay_days, 
                                  breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
# additionally, we want to know what pixels are outside of the range of empirical data training
min_empirical_hopkins <- min(df$hopkins_delay_days)
max_empirical_hopkins <- max(df$hopkins_delay_days)
# binary
anth_holo_df[, outside_empirical := hopkins_delay_days < min_empirical_hopkins |
               hopkins_delay_days > max_empirical_hopkins]
# numeric dist
anth_holo_df[, dist_from_empirical := fifelse(
  hopkins_delay_days < min_empirical_hopkins, 
  min_empirical_hopkins - hopkins_delay_days, # below range
  fifelse(
    hopkins_delay_days > max_empirical_hopkins, 
    hopkins_delay_days - max_empirical_hopkins, # above range
    0 #if within range
  )
)]

# extract midpoint of each hopkins bin to use as predictor
anth_holo_df <- anth_holo_df %>%
  mutate(
    # Extract numeric values from bin intervals
    bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
    bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
    hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
  )
# Some pixels lack elevation, which is why we get this warning. Probably try and see if we can fix this at a later time. 
anth_holo_df <- anth_holo_df %>% filter(!is.na(hopkins_bin))

# Furthermore, calculate env distance from training data 
env_model <- model_out$env_model
anth_holo_predict <- predict(env_model, newdata = anth_holo_df)
anth_holo_df[, env_predict_onset := anth_holo_predict]

# Calculate MESS
anth_holo_climate_proj <- project(anth_holo_reproj, crs(temp_seas), method = "near")
# crop the climate rasters
temp_seas_crop <- crop(temp_seas, anth_holo_climate_proj)
precip_coldest_q_crop <- crop(avg_monthly_precip_coldest_q, anth_holo_climate_proj)

# resample the species raster to match the climate raster exactly
anth_holo_resampled <- resample(anth_holo_climate_proj, temp_seas_crop, method = "near")

# mask with the resampled raster
temp_seas_crop <- mask(temp_seas_crop, anth_holo_resampled)
precip_coldest_q_crop <- mask(precip_coldest_q_crop, anth_holo_resampled)
# convert to old raster format
temp_seas_raster <- raster::raster(temp_seas_crop)
precip_coldest_q_raster <- raster::raster(precip_coldest_q_crop)

# stack and clean up names
pred_stack_raster <- raster::stack(temp_seas_raster, precip_coldest_q_raster)
names(pred_stack_raster) <- c("temp_seas", "precip_coldest_q")

# designate training data based on the empirical data
training_vars <- df %>% 
  ungroup() %>%
  select(temp_seas, precip_coldest_q) %>%
  distinct()  %>%  # remove exact cp cases as these wont be useful 
  as.data.frame()

# use dismo MESS for the full prediction raster compared to training vals
mess_raster_old <- dismo::mess(pred_stack_raster, training_vars, full = FALSE) # use old rast to use this fxn
mess_raster <- rast(mess_raster_old)

# extract and append to df 
# Extract to dataframe
mess_vals <- terra::extract(mess_raster, 
                            anth_holo_df[, .(longitude, latitude)])
anth_holo_df[, mess := mess_vals[[2]]]

# plot
ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = env_predict_onset)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Climatic Onset Prediction", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = mess)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

# Try clamping
# Get training data ranges
training_vars <- df %>% 
  ungroup() %>%
  select(temp_seas, precip_coldest_q)

train_ranges <- data.frame(
  var = c("temp_seas", "precip_coldest_q"),
  min = c(min(training_vars$temp_seas, na.rm = TRUE),
          min(training_vars$precip_coldest_q, na.rm = TRUE)),
  max = c(max(training_vars$temp_seas, na.rm = TRUE),
          max(training_vars$precip_coldest_q, na.rm = TRUE))
)

# Create clamped versions of predictors
anth_holo_df[, temp_seas_clamped := pmax(train_ranges$min[1], 
                                         pmin(train_ranges$max[1], temp_seas))]
anth_holo_df[, precip_coldest_q_clamped := pmax(train_ranges$min[2], 
                                                pmin(train_ranges$max[2], precip_coldest_q))]

# make only predictions w/in clamped data
clamped_data <- anth_holo_df[, .(temp_seas = temp_seas_clamped, 
                                 precip_coldest_q = precip_coldest_q_clamped)]

anth_holo_df[, env_predict_onset_clamped := predict(env_model, newdata = clamped_data)]


ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = env_predict_onset_clamped)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Clamped Climatic Onset Prediction", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = mess)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality
# because I dont want to wait all day. 
ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = mess)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

