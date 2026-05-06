### Title: Approximate a Baseline
### Author: JT Miller

### Purpose: Use Phenometrics binned spatially (both in latitudnal and elevational space) to build baseline expectations 

# Load Libraries 
library(data.table)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/proxy-baseline-fxns.R")

## Load in North American Taxonomy
na_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")

## Pull target taxa, order with the accepted parent name in lead
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Olneya tesota") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

target_taxa <- "Anthophora hololeuca"

described_phen_months <- fread("/blue/soltis/millerjared/Legume-Specialist-Occupancy/data/raw/hololeuca-and-host-pheno.csv")

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 13, # look at histogram when making this decision
                                      max_doy = 224 )

model_out$plots$hopkins
model_out$plots$enviroment




## Pull annotated flowering data from PhenoVision and CCH2 Datasets
# bring in these data
# all bee data from beeBDC
#bee_occs <- fread("/blue/soltis/millerjared/beeBDC/Bee-data/OutputData/05_cleaned_database.csv")
# additional data on Anthophora hololeuca provided by Michael Orr
#extra_hololeuca_data <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/raw/a-hololeuca-from-Orr.csv")
phenovision_flowering_data <- fread("/blue/guralnick/millerjared/phenovision_data/phenovision_flowers_03_15_2026.csv")
cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")
# use fxn built to filter these data and organize a synced up table
#taxa_flowering_data <- grab_flowering_occs(target_taxa, phenovision_flowering_data = phenovision_flowering_data, cch2_annotated_specimens = cch2_annotated_specimens)
taxa_flowering_data <- grab_bee_occs(target_taxa, bee_occs = bee_occs, extra_hololeuca_data = extra_hololeuca_data)

## Attach Elevation to these occs
# add a row ID to track where we're combing elev with occs
taxa_flowering_data[, row_id := .I]
taxa_flowering_coords <- taxa_flowering_data[!is.na(longitude) & !is.na(latitude)]
# load North America Elev Raster
na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
# convert occ data to a spatVect in WGS84
taxa_flowering_vect <- vect(taxa_flowering_coords, 
                           geom = c("longitude", "latitude"),
                           crs = "EPSG:4326")
# reproj to elevation raster projection
taxa_flowering_vect_reproj <- project(taxa_flowering_vect, crs(na_elev))
# extract elev values
elevation_vals <- terra::extract(na_elev, taxa_flowering_vect_reproj)
# add elevation values to occ table using row_id as our key
taxa_flowering_coords[, elevation_m := elevation_vals[,2]]
# merge
taxa_flowering_data <- taxa_flowering_data[taxa_flowering_coords[, .(row_id, elevation_m)],
                                           on = "row_id"]
# remove NA elevation vals
taxa_flowering_data <- taxa_flowering_data[!is.na(elevation_m)]

## Hopkins Expectation: For every degree in latitude or 120 meters of elevation, flowering will be delayed by 4 days. 
# use min vals as a baseline start
min_lat <- min(taxa_flowering_data$latitude)
min_elev <- min(taxa_flowering_data$elevation_m)
## create expected # of delay days based on hopkins
taxa_flowering_data[, hopkins_delay_days :=
                          ((latitude - min_lat)* 4) + # latitudnal effect
                          ((elevation_m - min_elev))/120 * 4] # elevation effect



# Create bins based on Hopkins delay by intervals of 4
taxa_flowering_data[, hopkins_bin := cut(hopkins_delay_days, 
                                             breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]

## Run Phenometric cleaning
# remove outliers (use priors of known monthly flowering?)
ggplot() + geom_bar(taxa_flowering_data, mapping = aes(x = dayOfYear))
taxa_flowering_data <- taxa_flowering_data[dayOfYear > 50 & dayOfYear < 200]
# downsample high intensity sample days (within any particular year)
taxa_flowering_data <- taxa_flowering_data  %>%
  group_by(species,year,hopkins_bin,dayOfYear) %>% dplyr::slice_sample(n = 3)
# remove non contemporary years
# taxa_flowering_data <- taxa_flowering_data %>% 
#   filter(year >= 2017)
# remove dayOfYears that have less than 4 total doys sampling
# taxa_flowering_data  <- taxa_flowering_data %>% group_by(species,year,hopkins_bin) %>%
#   filter(n_distinct(dayOfYear)>4)

## Phenesse An Estimate for Onset (10th percentile) based on Hopkins Bins 
## require at least 7 obs and 3 distinct doys per combination
# Calc OnSet
taxa_flowering_onset <- taxa_flowering_data %>% 
  group_by(species,hopkins_bin) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(dayOfYear) > 3) %>%
  group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$dayOfYear, percentile = 0.10, bootstraps=1000)))


# extract estimates
onset_estimates <- taxa_flowering_onset %>% 
  filter(column == "estimate") %>% 
  select(species, hopkins_bin, mean) %>% 
  rename(onset_doy = mean)

## Remove outlier estimates
onset_estimates <-  onset_estimates %>% 
  group_by(species) %>% 
  mutate(zonset = scale(onset_doy)) %>% 
  filter(between(zonset,-3.25,+3.25)) 


## Use Species Anchoring To Create a Predicted Hopkins Onset
# Use Rob's estimate sp anchors fxn
HOPKINS_LAT_SLOPE <- 4
HOPKINS_ELEV_SLOPE <- 4/120 
estimate_species_anchors <- function(df) {
  df |> group_by(species) |>
    group_modify(~{
      x <- .x
      slat <- median(x$lat_mean, na.rm = TRUE)
      selev <- median(x$mid_elev, na.rm = TRUE)  
      
      if (nrow(x) >= 10 && length(unique(x$lat_mean[!is.na(x$lat_mean)])) > 2) {
        m <- lm(onset ~ lat_mean, data = x)
        sdoy <- as.numeric(predict(m, newdata = data.frame(lat_mean = slat)))
      } else {
        sdoy <- median(x$onset, na.rm = TRUE)
      }
      tibble(start_lat = slat, start_elev = selev, start_doy = sdoy) 
    }) |>
    ungroup()
}

add_hopkins_predictions <- function(df) {
  df |>
    mutate(
      hopkins_lat = start_doy + (lat_mean - start_lat) * HOPKINS_LAT_SLOPE,
      hopkins_lat_elev = start_doy +
        (lat_mean - start_lat) * HOPKINS_LAT_SLOPE +
        (mid_elev - start_elev) * HOPKINS_ELEV_SLOPE 
    )
}

# bin summary (for creating hopkins predictions)
bin_summary <- taxa_flowering_data %>%
  group_by(species, hopkins_bin) %>%
  summarise(
    hopkins_lat = median(latitude, na.rm = TRUE),
    hopkins_elev = median(elevation_m, na.rm = TRUE),
    hopkins_delay = median(hopkins_delay_days, na.rm = TRUE),
    n_obs = n(),
    .groups = 'drop'
  )

# join bins 
onset_estimates <- onset_estimates %>% 
  left_join(bin_summary, by = "hopkins_bin") %>% 
  rename(lat_mean = hopkins_lat, 
         mid_elev = hopkins_elev, 
         onset = onset_doy, 
         species = species.x) # redundant

anchored_sp <- estimate_species_anchors(onset_estimates)

onset_estimate <- onset_estimates %>% 
  left_join(anchored_sp, by = 'species')

onset_estimate <- add_hopkins_predictions(onset_estimate)

ggplot(onset_estimate, aes(x = hopkins_lat_elev, y = onset)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Hopkins Expected Onset (day of year)",
       y = "Observed Onset (day of year)",
       title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
  theme_minimal()



# extract midpoint of each hopkins bin to use as predictor
onset_estimate <- onset_estimate %>%
  mutate(
    # Extract numeric values from bin intervals
    bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
    bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
    hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
  )

# Fit a linear model for bins midpoint and the onset val
hopkins_model <- lm(onset ~ hopkins_delay_midpoint, data = onset_estimate)
# visualize
ggplot(onset_estimate, aes(x = hopkins_delay_midpoint, y = onset)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Hopkins Expected Delay Bin (Middle Value for Bin)",
       y = "Observed Onset (day of year)",
       title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
  theme_minimal()

# Extract slope
# coef_test <- summary(hopkins_model)$coefficients
# slope <- coef_test[2, 1]
# se_slope <- coef_test[2, 2]

# Extract the corrected slope (actual delay rate per Hopkins bin unit)
actual_delay_rate <- coef(hopkins_model)[2]
actual_delay_se <- summary(hopkins_model)$coefficients[2, 2]

cat("Hopkins assumes 4 days per unit, but actual rate is:", 
    round(actual_delay_rate, 3), "±", round(actual_delay_se, 3), "\n")

# Apply the corrected slope with uncertainty bounds
onset_estimate <- onset_estimate %>%
  mutate(
    #cCalculate delay from the anchor point using Hopkins bins
    delay_from_anchor = hopkins_delay_midpoint - 
      ((start_lat - min_lat) * 4 + 
         (start_elev - min_elev) / 120 * 4),
    
    # apply actual observed delay rate
    corrected_onset_prediction = start_doy + (delay_from_anchor * actual_delay_rate),
    
    # Add uncertainty bounds based on SE of slope
    corrected_onset_lower = start_doy + (delay_from_anchor * (actual_delay_rate - actual_delay_se)),
    corrected_onset_upper = start_doy + (delay_from_anchor * (actual_delay_rate + actual_delay_se))
  )

# Visualize with SE ribbon
ggplot(onset_estimate, aes(x = corrected_onset_prediction, y = onset)) +
  geom_ribbon(aes(ymin = corrected_onset_lower, ymax = corrected_onset_upper), 
              alpha = 0.2, fill = "blue") +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Corrected Predicted Onset (anchored + actual slope)",
       y = "Observed Onset (day of year)",
       title = paste("Species-Specific Delay Rate for", target_taxa[1])) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 1.2, color = "black"))

hopkins_corrected_model <- lm(onset ~ corrected_onset_prediction, data = onset_estimate)
summary(hopkins_corrected_model)


ggplot() +
  geom_line(onset_estimate, mapping = aes(x = hopkins_delay_midpoint, y = corrected_onset_prediction)) + 
  geom_point(onset_estimate, mapping = aes(x = hopkins_delay_midpoint, y = onset, color = "darkred")) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = 1.2, color = "black")) + 
  labs(x = "Hopkins MidPoint Bins", y = "Corrected Onset Prediction", title = "Predicted Onset (TrendLine) and Empirical (Dots) against Hopkins Gradient")

# Load in a well sampled taxa (Desert Chicory)
# phenovision_flowering_data <- fread("/blue/guralnick/millerjared/phenovision_data/phenovision_flowers_03_15_2026.csv")
# rafinesquia_neomexicana1 <- phenovision_flowering_data[scientificName == "Psorothamnus spinosus"]
# rafinesquia_neomexicana2 <- phenovision_flowering_data[scientificName == "Psorodendron spinosum"]
# rafinesquia_neomexicana <- rbind(rafinesquia_neomexicana1, rafinesquia_neomexicana2)
# rafinesquia_neomexicana <- rafinesquia_neomexicana %>% 
#   mutate(scientificName = "Psorodendron spinosum")

# cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")
# rafinesquia_neomexicana_cch2 <- cch2_annotated_specimens %>% filter(species == "Psorothamnus spinosus")

# clean
# rafinesquia_neomexicana_cch2 <- rafinesquia_neomexicana_cch2 %>% 
#   filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(eventDate)) %>% 
#   select(-scientificName) %>% 
#   mutate(dayOfYear = lubridate::yday(eventDate), 
#          year = lubridate::year(eventDate)) %>% 
#   rename(longitude = decimalLongitude, 
#          latitude = decimalLatitude, 
#          scientificName = species)  %>% 
#   mutate(scientificName = "Psorodendron spinosum")

# select and combine
#rafinesquia_neomexicana <- rbind(select(rafinesquia_neomexicana, scientificName, longitude, latitude, year, dayOfYear), select(rafinesquia_neomexicana_cch2, scientificName, longitude, latitude, year, dayOfYear))


# 
# # make spatVector
# r_neomexicana_vect <- vect(rafinesquia_neomexicana, 
#                                      geom = c("longitude", "latitude"),
#                                      crs = "EPSG:4326")
# # reproj to elevation raster projection
# r_neomexicana_vect_reproj <- project(r_neomexicana_vect, crs(na_elev))
# 
# # extract elev values
# elevation_vals <- terra::extract(na_elev, r_neomexicana_vect_reproj)
# 
# # add elevation values to table
# rafinesquia_neomexicana[, elevation_m := elevation_vals[,2]]
# rafinesquia_neomexicana <- rafinesquia_neomexicana[!is.na(elevation_m)]
# 
# # Hopkins Expectation: For every degree in latitude or 120 meters of elevation, flowering will be delayed by 4 days. 
# ## use min vals as a baseline start
# min_lat <- min(rafinesquia_neomexicana$latitude)
# min_elev <- min(rafinesquia_neomexicana$elevation_m)
# ## create bins
# rafinesquia_neomexicana[, hopkins_delay_days :=
#                           ((latitude - min_lat)* 4) + # latitudnal effect
#                           ((elevation_m - min_elev))/120 * 4 # elevation effect
#                         ]
# 
# 
# # Bin by hopkins delay 
# rafinesquia_neomexicana[, hopkins_bin := cut(hopkins_delay_days, 
#                                              breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
# 
# # remove spring outliers
# rafinesquia_neomexicana2 <- rafinesquia_neomexicana[dayOfYear > 20 & dayOfYear < 182]
# 
# # constrain to only 2017-2025
# rafinesquia_neomexicana2 <- rafinesquia_neomexicana2 %>% filter(year >= 2017) %>% filter(year <= 2025)
# 
# # downsample high intensity sample days
# rafinesquia_neomexicana2 <- rafinesquia_neomexicana2  %>%
#   group_by(scientificName,year,hopkins_bin,dayOfYear) %>% dplyr::slice_sample(n = 3)
# 
# # remove dayOfYears that have less than 4 total doys sampling
# rafinesquia_neomexicana2  <- rafinesquia_neomexicana2 %>% group_by(scientificName,year,hopkins_bin) %>% 
#   filter(n_distinct(dayOfYear)>4) 
# 
# # bin summary (for creating hopkins predictions)
# bin_summary <- rafinesquia_neomexicana2 %>%
#   group_by(scientificName, hopkins_bin) %>%
#   summarise(
#     hopkins_lat = median(latitude, na.rm = TRUE),
#     hopkins_elev = median(elevation_m, na.rm = TRUE),
#     hopkins_delay = median(hopkins_delay_days, na.rm = TRUE),
#     n_obs = n(),
#     .groups = 'drop'
#   )
# 
# # use phenesse to calc bin phenology estimates for onset 
# ## require at least 7 obs and 3 distinct doys per combination
# # Calc OnSet
# r_neomexicana_onset <- rafinesquia_neomexicana2 %>% 
#   group_by(scientificName,hopkins_bin) %>% 
#   filter(n() >= 7) %>%
#   filter(n_distinct(dayOfYear) > 3) %>%
#   group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$dayOfYear, percentile = 0.05, bootstraps=250)))
# 
# 
# # extract estimates
# onset_estimates <- r_neomexicana_onset %>% 
#   filter(column == "estimate") %>% 
#   select(scientificName, hopkins_bin, mean) %>% 
#   rename(onset_doy = mean)
# 
# # extract midpoint of each hopkins bin to use as predictor
# onset_estimates <- onset_estimates %>%
#   mutate(
#     # Extract numeric values from bin intervals
#     bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
#     bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
#     hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
#   )
# 
# # fit a linear model 
# hopkins_model <- lm(onset_doy ~ hopkins_delay_midpoint, data = onset_estimates)
# ggplot(onset_estimates, aes(x = hopkins_delay_midpoint, y = onset_doy)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(x = "Hopkins Expected Delay (days)",
#        y = "Observed Onset (day of year)",
#        title = "Testing Hopkins' Bioclimatic Law") +
#   theme_minimal()
# 
# # Use Rob's estimate sp anchors fxn
# HOPKINS_LAT_SLOPE <- 4
# HOPKINS_ELEV_SLOPE <- 4/120 
# estimate_species_anchors <- function(df) {
#   df |> group_by(scientificName) |>
#     group_modify(~{
#       x <- .x
#       slat <- median(x$lat_mean, na.rm = TRUE)
#       selev <- median(x$mid_elev, na.rm = TRUE)  
#       
#       if (nrow(x) >= 10 && length(unique(x$lat_mean[!is.na(x$lat_mean)])) > 2) {
#         m <- lm(onset ~ lat_mean, data = x)
#         sdoy <- as.numeric(predict(m, newdata = data.frame(lat_mean = slat)))
#       } else {
#         sdoy <- median(x$onset, na.rm = TRUE)
#       }
#       tibble(start_lat = slat, start_elev = selev, start_doy = sdoy) 
#     }) |>
#     ungroup()
# }
# 
# add_hopkins_predictions <- function(df) {
#   df |>
#     mutate(
#       hopkins_lat = start_doy + (lat_mean - start_lat) * HOPKINS_LAT_SLOPE,
#       hopkins_lat_elev = start_doy +
#         (lat_mean - start_lat) * HOPKINS_LAT_SLOPE +
#         (mid_elev - start_elev) * HOPKINS_ELEV_SLOPE 
#     )
# }
# 
# # join bins 
# onset_estimates2 <- onset_estimates %>% 
#   left_join(bin_summary, by = "hopkins_bin") %>% 
#   rename(lat_mean = hopkins_lat, 
#          mid_elev = hopkins_elev, 
#          onset = onset_doy, 
#          scientificName = scientificName.x)
# 
# anchored_sp <- estimate_species_anchors(onset_estimates2)
# 
# onset_estimate3 <- onset_estimates2 %>% 
#   left_join(anchored_sp, by = 'scientificName')
# 
# onset_estimate3 <- add_hopkins_predictions(onset_estimate3)
# 
# ggplot(onset_estimate3, aes(x = hopkins_lat_elev, y = onset)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(x = "Hopkins Expected Onset (day of year)",
#        y = "Observed Onset (day of year)",
#        title = "Testing Hopkins' Bioclimatic Law") +
#   theme_minimal()
# 
# # Extract coefficients
# coef_test <- summary(hopkins_model)$coefficients
# slope <- coef_test[2, 1]
# se_slope <- coef_test[2, 2]
# 
# # Test H0: slope = 1
# t_stat <- (slope - 1) / se_slope
# p_value <- 2 * pt(abs(t_stat), df = nrow(onset_estimates) - 2, lower.tail = FALSE)

