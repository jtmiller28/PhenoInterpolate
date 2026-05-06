### Title: Approximate a Baseline
### Author: JT Miller

### Purpose: Use Phenometrics binned spatially (both in latitudnal and elevational space) to build baseline expectations 

# Load Libraries 
library(data.table)
library(sf)
library(terra)
library(dplyr)
library(performance)
library(ggplot2)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/proxy-baseline-fxns.R")

## Load in North American Taxonomy
na_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")

## Pull target taxa, order with the accepted parent name in lead
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorothamnus emoryi") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

#target_taxa <- "Anthophora hololeuca"

## Pull annotated flowering data from PhenoVision and CCH2 Datasets
# bring in these data
# all bee data from beeBDC
#bee_occs <- fread("/blue/soltis/millerjared/beeBDC/Bee-data/OutputData/05_cleaned_database.csv")
# additional data on Anthophora hololeuca provided by Michael Orr
#extra_hololeuca_data <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/raw/a-hololeuca-from-Orr.csv")
#phenovision_flowering_data <- fread("/blue/guralnick/millerjared/phenovision_data/phenovision_flowers_03_15_2026.csv")
#cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")
# use fxn built to filter these data and organize a synced up table
taxa_flowering_data <- grab_flowering_occs(target_taxa, phenovision_flowering_data = phenovision_flowering_data, cch2_annotated_specimens = cch2_annotated_specimens)
#taxa_flowering_data <- grab_bee_occs(target_taxa, bee_occs = bee_occs, extra_hololeuca_data = extra_hololeuca_data)

## Attach Climate and Elevation to these occs
# add a row ID to track where we're combing elev with occs
taxa_flowering_data[, row_id := .I]
taxa_flowering_coords <- taxa_flowering_data[!is.na(longitude) & !is.na(latitude)]
# load North America Elev Raster
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


# convert occ data to a spatVect in WGS84
taxa_flowering_vect <- vect(taxa_flowering_coords, 
                            geom = c("longitude", "latitude"),
                            crs = "EPSG:4326")
# reproj to elevation raster projection
taxa_flowering_vect_reproj <- project(taxa_flowering_vect, crs(na_elev))
taxa_flowering_vect_reproj2 <- project(taxa_flowering_vect, crs(avg_temp)) # do the same for chelsa bioclim layers
# Stack all climate layers
climate_stack <- c(avg_temp, temp_seas, tmax, tmin, annual_precip, 
                   precip_of_wettest_month, precip_of_driest_month, 
                   precip_seas, avg_monthly_precip_wettest_q, 
                   avg_monthly_precip_driest_q, avg_monthly_precip_warmest_q, 
                   avg_monthly_precip_coldest_q)

# Give them meaningful names
names(climate_stack) <- c("avg_temp", "temp_seas", "tmax", "tmin", 
                          "annual_precip", "precip_wettest_month", 
                          "precip_driest_month", "precip_seas", 
                          "precip_wettest_q", "precip_driest_q", 
                          "precip_warmest_q", "precip_coldest_q")

# Extract all climate values at once
climate_vals <- terra::extract(climate_stack, taxa_flowering_vect_reproj2)
elevation_vals <- terra::extract(na_elev, taxa_flowering_vect_reproj)
# extract clim values

# Add all values to your data
taxa_flowering_coords[, elevation_m := elevation_vals[,2]]
taxa_flowering_coords[, names(climate_vals)[-1] := climate_vals[,-1]]

# add elevation values to occ table using row_id as our key
taxa_flowering_coords[, elevation_m := elevation_vals[,2]]
# merge
taxa_flowering_data <- taxa_flowering_data[taxa_flowering_coords[, .(row_id, elevation_m, avg_temp, temp_seas, tmax, tmin, annual_precip, 
                                                                     precip_wettest_month, precip_driest_month, precip_seas, precip_wettest_q, 
                                                                     precip_driest_q, precip_warmest_q, precip_coldest_q)],
                                           on = "row_id"]
# remove NA elevation vals
taxa_flowering_data <- taxa_flowering_data[!is.na(elevation_m) & !is.na(avg_temp) & !is.na(temp_seas) & !is.na(tmax) & !is.na(tmin) & !is.na(annual_precip) &
                                             !is.na(precip_wettest_month) & !is.na(precip_driest_month) & !is.na(precip_seas) & !is.na(precip_wettest_q) &
                                             !is.na(precip_driest_q) & !is.na(precip_warmest_q) & !is.na(precip_coldest_q)]

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
taxa_flowering_data <- taxa_flowering_data[dayOfYear > 0 & dayOfYear < 200]
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

# Attach env avg summaries per hopkins bin
hopkins_bin_env_summaries <- taxa_flowering_data %>% 
  group_by(species, hopkins_bin) %>% 
  filter(n() >= 7) %>%
  filter(n_distinct(dayOfYear) > 3) %>% 
  summarize(avg_temp = mean(avg_temp), 
            temp_seas = mean(temp_seas), 
            tmax = mean(tmax), 
            tmin = mean(tmin),
            annual_precip = mean(annual_precip), 
            precip_wettest_month = mean(precip_wettest_month), 
            precip_driest_month = mean(precip_driest_month), 
            precip_seas = mean(precip_seas), 
            precip_wettest_q = mean(precip_wettest_q), 
            precip_driest_q = mean(precip_driest_q), 
            precip_warmest_q = mean(precip_warmest_q), 
            precip_coldest_q = mean(precip_coldest_q)) %>% 
  ungroup() %>% 
  dplyr::select(hopkins_bin, avg_temp, temp_seas, tmax, tmin, annual_precip, 
                precip_wettest_month, precip_driest_month, precip_seas, 
                precip_wettest_q, precip_driest_q, precip_warmest_q, 
                precip_coldest_q)

onset_estimates <- onset_estimates %>% 
  left_join(hopkins_bin_env_summaries, by = "hopkins_bin")


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


# Start with full model
predictors <- c("hopkins_delay_midpoint", "avg_temp", "temp_seas", "tmax", "tmin",
                "annual_precip", "precip_wettest_month", "precip_driest_month",
                "precip_seas", "precip_wettest_q", "precip_driest_q", "precip_warmest_q",
                "precip_coldest_q")

predictors <- c("hopkins_delay_midpoint", "temp_seas", "tmax", "tmin",
                "precip_seas", "precip_wettest_q", "precip_driest_q", "precip_warmest_q",
                "precip_coldest_q")
vif_threshold <- 5
repeat{
  # build formula 
  formula_str <- paste("onset ~", paste(predictors, collapse = " + "))
  current_model <- glm(as.formula(formula_str), data = onset_estimate)
  # check collinearity
  vif_results <- performance::check_collinearity(current_model)
  print(vif_results)
  # IF all VIFs are below the threshold, we can break
  if(all(vif_results$VIF < vif_threshold)){
    cat("All VIFs below threshold\n")
    break
  }
  # Remove variable with highest VIF
  max_vif_var <- vif_results$Term[which.max(vif_results$VIF)]
  cat("Removing:", max_vif_var, "with VIF =", max(vif_results$VIF), "\n")
  predictors <- predictors[predictors != max_vif_var]
}
base_model <- glm(as.formula(formula_str), data = onset_estimate)
best_model <- MASS::stepAIC(base_model, direction = "both", trace = FALSE) # run model selection.

summary(best_model)


# Vis 
onset_estimate$predicted <- predict(best_model, type = "response")
onset_estimate$residuals <- residuals(best_model)

r2_val <- r2(best_model)
rmse_val <- sqrt(mean(onset_estimate$residuals^2))

ggplot(onset_estimate, aes(x = onset, y = predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = min(onset_estimate$onset), y = max(onset_estimate$predicted),
           label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_val$R2, rmse_val),
           hjust = 0, vjust = 1) +
  labs(x = "Observed Onset", y = "Predicted Onset") +
  theme_minimal()
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
