### Functions for Proxy Baseline Script 

### Function to grab taxa from both the phenoVision & cch2 data
grab_flowering_occs <- function(taxa, phenovision_flowering_data, cch2_annotated_specimens){
  # search and organize taxa from phenoVision
  taxa_pv_flowering_occs <- phenovision_flowering_data[scientificName %in% c(taxa)]
  taxa_pv_flowering_occs <- taxa_pv_flowering_occs %>% 
    mutate(species = taxa[1]) %>%  # create a field name that contains the alignedName
    rename(verbatimScientificName = scientificName) %>% # make it clear that this is the original name
    mutate(month = lubridate::year(date), 
           day = lubridate::day(date)) %>% 
    select(species, verbatimScientificName, longitude, latitude, year, month, day, dayOfYear) %>% 
    mutate(source = "phenoVision")
  
  # search and organize taxa from CCH2 annotated specimen data
  taxa_spec_flowering_occs <- cch2_annotated_specimens[species %in% c(taxa)]
  taxa_spec_flowering_occs <- taxa_spec_flowering_occs %>% 
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude), !is.na(eventDate)) %>% # remove
    select(-scientificName) %>% # not useful for us
    mutate(dayOfYear = lubridate::yday(eventDate), 
           year = lubridate::year(eventDate), 
           month = lubridate::month(eventDate), 
           day = lubridate::day(eventDate)) %>% 
    rename(longitude = decimalLongitude, 
           latitude = decimalLatitude, 
           verbatimScientificName = species)  %>% 
    mutate(species = taxa[1]) %>% 
    select(species, verbatimScientificName, longitude, latitude, year, month, day, dayOfYear) %>% 
    mutate(source = "specimen")
  
  # Bring datasets together 
  taxa_flowering_occs <- rbind(taxa_pv_flowering_occs, taxa_spec_flowering_occs)
  
}


## Function to grab bee occs
grab_bee_occs <- function(taxon, bee_occs, extra_hololeuca_data){
  # reformat this data to match one another
  extra_hololeuca_data <- extra_hololeuca_data[, .(species = paste(GenusName, Species), 
                                                   decimalLongitude = DecimalLon, 
                                                   decimalLatitude = DecimalLat, 
                                                   year = Yr0, 
                                                   month = Mon0, 
                                                   day = Day0, 
                                                   associatedTaxa = CONCATENATE)]
  extra_hololeuca_data <- extra_hololeuca_data[, `:=`(coordinateUncertaintyInMeters = NA_real_,
                                                      basisOfRecord = NA_character_,
                                                      hasGeospatialIssues = NA,
                                                      isDuplicateOf = NA_character_)]
  extra_hololeuca_data[, month := match(month, month.abb)]
  extra_hololeuca_data[, eventDate := paste(year, month, day, sep = "-")]
  bee_occs <- bee_occs[, .(species, decimalLongitude, decimalLatitude,coordinateUncertaintyInMeters,
                           eventDate,day, month, year, basisOfRecord, hasGeospatialIssues, 
                           isDuplicateOf, associatedTaxa)]
  
  # deal with A. hololeuca 
  bee_occs_a_holo <- bee_occs %>% filter(species == "Anthophora hololeuca")
  data.table::setkey(extra_hololeuca_data, decimalLongitude, decimalLatitude, day, month, year)
  data.table::setkey(bee_occs_a_holo, decimalLongitude, decimalLatitude, day, month, year)
  
  # find records in bee_occs_a_holo that are not in extra_hololeuca_data
  unique_bee_occs <- bee_occs_a_holo[!extra_hololeuca_data]
  
  # combine datasets
  a_holo_final <- rbindlist(list(extra_hololeuca_data, unique_bee_occs), fill = TRUE)
  
  # remove this species from the rest
  bee_occs <- bee_occs %>% filter(species != "Anthophora hololeuca")
  # also remove blank species 
  bee_occs <- bee_occs %>% filter(species != "")
  
  bee_occs <- rbind(bee_occs, a_holo_final)
  
  taxon_occs <- bee_occs %>% 
    filter(species == taxon)
  
  # reformat fields 
  taxon_occs <- taxon_occs %>%
    filter(!is.na(decimalLongitude) & !is.na(decimalLatitude) & !is.na(eventDate)) %>%
    select(species, decimalLongitude, decimalLatitude, eventDate) %>%
    mutate(
      # Parse eventDate - handles both "YYYY-MM-DD" and "YYYY/MM/DD HH:MM:SS"
      eventDate = lubridate::ymd_hms(eventDate, truncated = 3),  # truncated=3 allows dates without time
      dayOfYear = lubridate::yday(eventDate),
      year = lubridate::year(eventDate), 
      month = lubridate::month(eventDate), 
      day = lubridate::month(eventDate)
    ) %>%
    rename(longitude = decimalLongitude,
           latitude = decimalLatitude) %>%
    filter(!is.na(dayOfYear)) %>% 
    mutate(source = "occurrence")
  return(taxon_occs)
    
}

## Build Baseline Model ########################################################
model_baseline_phenology <- function(target_taxa, 
                                     phenometric = 0.10, # set as onset by default
                                     bees_or_plants = "plants", # target taxon group to pull empirical data for
                                     phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                     cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                     bee_occs, # input bee occurrence data
                                     extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                     min_doy = 0, # look at histogram when making this decision
                                     max_doy = 365 # ditto to above
                                     ){

## Get Phenology Data
if(bees_or_plants == "plants"){
  phen_data <- grab_flowering_occs(target_taxa, phenovision_flowering_data = phenovision_flowering_data, cch2_annotated_specimens = cch2_annotated_specimens)
}
if(bees_or_plants == "bees"){
  phen_data <- grab_bee_occs(target_taxa, bee_occs = bee_occs, extra_hololeuca_data = extra_hololeuca_data)
}

## Attach Elevation & Climatic Variable Data
# add a row ID to track where we're combing vars with occs
  phen_data[, row_id := .I]
  phen_data_coords <- phen_data[!is.na(longitude) & !is.na(latitude)]
  # load North America Elev Raster & Chelsa annual clim vars
  na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
  avg_temp <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio1_1981-2010_V.2.1.tif")
  avg_diurnal_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio02_1981-2010_V.2.1.tif")
  isothermality <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio03_1981-2010_V.2.1.tif")
  temp_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio4_1981-2010_V.2.1.tif")
  tmax <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio5_1981-2010_V.2.1.tif")
  tmin <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio6_1981-2010_V.2.1.tif")
  annual_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio07_1981-2010_V.2.1.tif")
  avg_monthly_temp_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio08_1981-2010_V.2.1.tif")
  avg_monthly_temp_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio09_1981-2010_V.2.1.tif")
  avg_monthly_temp_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio10_1981-2010_V.2.1.tif")
  avg_monthly_temp_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio11_1981-2010_V.2.1.tif")
  annual_precip <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio12_1981-2010_V.2.1.tif")
  precip_wettest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio13_1981-2010_V.2.1.tif")
  precip_driest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio14_1981-2010_V.2.1.tif")
  precip_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio15_1981-2010_V.2.1.tif")
  avg_monthly_precip_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio16_1981-2010_V.2.1.tif")
  avg_monthly_precip_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio17_1981-2010_V.2.1.tif")
  avg_monthly_precip_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio18_1981-2010_V.2.1.tif")
  avg_monthly_precip_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio19_1981-2010_V.2.1.tif")
  
# convert occ data to a spatVect in WGS84
  phen_data_vect <- vect(phen_data_coords, 
                              geom = c("longitude", "latitude"),
                              crs = "EPSG:4326")
  # reproj to elevation raster projection & clim raster projections
  phen_data_vect_reproj1 <- project(phen_data_vect, crs(na_elev))
  phen_data_vect_reproj2 <- project(phen_data_vect, crs(avg_temp)) 
  # Stack all climate layers
  climate_stack <- c(avg_temp, avg_diurnal_temp_range, isothermality,
                     temp_seas, tmax, tmin, annual_temp_range,
                     avg_monthly_temp_wettest_q, avg_monthly_temp_driest_q,
                     avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                     annual_precip, precip_wettest_month, precip_driest_month,
                     precip_seas, avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q,
                     avg_monthly_precip_warmest_q, avg_monthly_precip_coldest_q)
  
  # make names clearer
  names(climate_stack) <- c("avg_temp", "avg_diurnal_temp_range", "isothermality",
                            "temp_seas", "tmax", "tmin", "annual_temp_range",
                            "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                            "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                            "annual_precip", "precip_wettest_month", "precip_driest_month",
                            "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                            "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q")
  
  # Extract all raster values
  elevation_vals <- terra::extract(na_elev, phen_data_vect_reproj1)
  climate_vals <- terra::extract(climate_stack, phen_data_vect_reproj2)

  
  # Add all values to these data
  phen_data_coords[, elevation_m := elevation_vals[,2]]
  phen_data_coords[, names(climate_vals)[-1] := climate_vals[,-1]]
  
  # add elevation values to occ table using row_id as our key
  phen_data_coords[, elevation_m := elevation_vals[,2]]
  # merge
  phen_data <- phen_data[phen_data_coords[, .(row_id, elevation_m, avg_temp, avg_diurnal_temp_range, isothermality,
                                              temp_seas, tmax, tmin, annual_temp_range,
                                              avg_monthly_temp_wettest_q, avg_monthly_temp_driest_q,
                                              avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                                              annual_precip, precip_wettest_month, precip_driest_month,
                                              precip_seas, avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q,
                                              avg_monthly_precip_warmest_q, avg_monthly_precip_coldest_q)],
                                             on = "row_id"]
  # remove NA elevation vals
  phen_data <- phen_data[!is.na(elevation_m) & !is.na(avg_temp) & !is.na(avg_diurnal_temp_range) & !is.na(isothermality) & !is.na(temp_seas) 
                         & !is.na(tmax) & !is.na(tmin) & !is.na(annual_temp_range) & !is.na(avg_monthly_temp_wettest_q) & !is.na(avg_monthly_temp_driest_q) &
                           !is.na(avg_monthly_temp_warmest_q) & !is.na(avg_monthly_temp_coldest_q) & !is.na(annual_precip) &
                                               !is.na(precip_wettest_month) & !is.na(precip_driest_month) & !is.na(precip_seas) & !is.na(avg_monthly_precip_wettest_q) &
                                               !is.na(avg_monthly_precip_driest_q) & !is.na(avg_monthly_precip_warmest_q) & !is.na(avg_monthly_precip_coldest_q)]
  
## Create Hopkins Expectation: For every degree in latitude or 120 meters of elevation, flowering will be delayed by 4 days. 
  # use min vals as a baseline start
  min_lat <- min(phen_data$latitude)
  min_elev <- min(phen_data$elevation_m)
  ## create expected # of delay days based on hopkins
  phen_data[, hopkins_delay_days :=
                        ((latitude - min_lat)* 4) + # latitudnal effect
                        ((elevation_m - min_elev))/120 * 4] # elevation effect
  
  # Create bins for phenology data based on Hopkins delay by intervals of 4
  phen_data[, hopkins_bin := cut(hopkins_delay_days, 
                                           breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]

## Run Phenometric Cleaning 
  # remove outliers based on known flowering 
  phen_data <- phen_data %>%  filter(dayOfYear > min_doy & dayOfYear < max_doy)
  # downsample any dayOfYears that have abnormally high sampling on a dayOfYear within a year (City Nature Challenge and BioBlitzs)
  phen_data <- phen_data %>% 
    group_by(species,year,hopkins_bin,dayOfYear) %>% 
    dplyr::slice_sample(n = 3)
  # if sample is either a specimen (Herbarium data) or an occurrence (unspecified in bee data), then remove first day of the months from these data
  # these are possibly a data quality issue
  phen_data <- phen_data %>% 
    filter(!(source %in% c("occurrence", "specimen") & day == 1))
  

## Create a phenesse estimate for phenometric 
 phen_data_est <- phen_data %>% 
    group_by(species,hopkins_bin) %>% 
    filter(n() >= 7) %>%
    filter(n_distinct(dayOfYear) > 3) %>% 
    group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$dayOfYear, percentile = phenometric, bootstraps=1000)))
 
# extract estimate
 phen_estimates <- phen_data_est %>% 
   filter(column == "estimate") %>% 
   select(species, hopkins_bin, mean) %>% 
   rename(onset_doy = mean) %>% ungroup()
 
 
## Attach env avg summaries per hopkins bin
 hopkins_bin_env_summaries <- phen_data %>% 
   group_by(species, hopkins_bin) %>% 
   filter(n() >= 7) %>%
   filter(n_distinct(dayOfYear) > 3) %>% 
   summarize(avg_temp = mean(avg_temp), 
             avg_diurnal_temp_range = mean(avg_diurnal_temp_range), 
             isothermality = mean(isothermality),
             temp_seas = mean(temp_seas), 
             tmax = mean(tmax), 
             tmin = mean(tmin),
             annual_temp_range = mean(annual_temp_range), 
             avg_monthly_temp_wettest_q = mean(avg_monthly_temp_wettest_q),
             avg_monthly_temp_driest_q = mean(avg_monthly_temp_driest_q), 
             avg_monthly_temp_warmest_q = mean(avg_monthly_temp_warmest_q), 
             avg_monthly_temp_coldest_q = mean(avg_monthly_temp_coldest_q), 
             annual_precip = mean(annual_precip), 
             precip_wettest_month = mean(precip_wettest_month), 
             precip_driest_month = mean(precip_driest_month), 
             precip_seas = mean(precip_seas), 
             avg_monthly_precip_wettest_q = mean(avg_monthly_precip_wettest_q), 
             avg_monthly_precip_driest_q = mean(avg_monthly_precip_driest_q), 
             avg_monthly_precip_warmest_q = mean(avg_monthly_precip_warmest_q), 
             avg_monthly_precip_coldest_q = mean(avg_monthly_precip_coldest_q)) %>% 
   ungroup() %>% 
   dplyr::select(hopkins_bin, avg_temp, avg_diurnal_temp_range, isothermality,
                 temp_seas, tmax, tmin, annual_temp_range, avg_monthly_temp_wettest_q, 
                 avg_monthly_temp_driest_q, avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q, 
                 annual_precip, precip_wettest_month, precip_driest_month, precip_seas, 
                 avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q, avg_monthly_precip_warmest_q, 
                 avg_monthly_precip_coldest_q)
 ## Join 
 phen_estimates <- phen_estimates %>% 
   left_join(hopkins_bin_env_summaries, by = "hopkins_bin")
 
## Remove Outliers 
 phen_estimates <-  phen_estimates %>% 
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

 # bin summary (for creating hopkins predictions)
 bin_summary <- phen_data %>%
   group_by(species, hopkins_bin) %>%
   summarise(
     hopkins_lat = median(latitude, na.rm = TRUE),
     hopkins_elev = median(elevation_m, na.rm = TRUE),
     hopkins_delay = median(hopkins_delay_days, na.rm = TRUE),
     n_obs = n(),
     .groups = 'drop'
   )
 
 # join bins 
phen_estimates <- phen_estimates %>% 
   left_join(bin_summary, by = "hopkins_bin") %>% 
   rename(lat_mean = hopkins_lat, 
          mid_elev = hopkins_elev, 
          onset = onset_doy, 
          species = species.x) # redundant
 
 anchored_sp <- estimate_species_anchors(phen_estimates)
 
 
 phen_estimate <- phen_estimates %>% 
   left_join(anchored_sp, by = 'species')
 
 # extract midpoint of each hopkins bin to use as predictor
 phen_estimate <- phen_estimate %>%
   mutate(
     # Extract numeric values from bin intervals
     bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
     bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
     hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
   )
 
## Run a Hopkins model based on bins
 hopkins_model <- lm(onset ~ hopkins_delay_midpoint, data = phen_estimate)
 r2_hop <- performance::r2(hopkins_model)
 res_hop <- residuals(hopkins_model)
 rmse_val_hop <- sqrt(mean(res_hop^2))
 # Extract the corrected slope (actual delay rate per Hopkins bin unit)
 actual_delay_rate <- coef(hopkins_model)[2]
 actual_delay_se <- summary(hopkins_model)$coefficients[2, 2]
 # Apply the corrected slope with uncertainty bounds
 phen_estimate <- phen_estimate %>%
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

## Run a env model with model selection based on bins 
#  if(nrow(phen_estimate) > 3){ # min req to make env 
#    predictors <- c("hopkins_delay_midpoint", "avg_temp", "temp_seas", 
#                    
#                    "precip_seas",  "precip_warmest_q",
#                    "precip_coldest_q")
#  
#  vif_threshold <- 5
#  repeat{
#    # build formula
#    formula_str <- paste("onset ~", paste(predictors, collapse = " + "))
#    current_model <- glm(as.formula(formula_str), data = phen_estimate)
#    # check collinearity
#    vif_results <- performance::check_collinearity(current_model)
#    print(vif_results)
#    # IF all VIFs are below the threshold, we can break
#    if(all(vif_results$VIF < vif_threshold)){
#      cat("All VIFs below threshold\n")
#      break
#    }
#    # Remove variable with highest VIF
#    max_vif_var <- vif_results$Term[which.max(vif_results$VIF)]
#    cat("Removing:", max_vif_var, "with VIF =", max(vif_results$VIF), "\n")
#    predictors <- predictors[predictors != max_vif_var]
#  }
#  base_model <- glm(as.formula(formula_str), data = phen_estimate)
#  best_model <- MASS::stepAIC(base_model, direction = "both", trace = FALSE) # run model selection.
#  env_model <- best_model # for clarity
# phen_estimate$predicted <- predict(env_model, type = "response")
# phen_estimate$residuals <- residuals(env_model)
# r2_obj <- performance::r2(env_model)
# r2_val <- r2_obj$R2
# res_env <- residuals(env_model)
# rmse_val <- sqrt(mean(res_env^2))
#  } 
 # if(nrow(phen_estimate) > 3){ # min req to make env
 #   predictors <- c("hopkins_delay_midpoint", "avg_temp", "avg_diurnal_temp_range", "isothermality",
 #                   "temp_seas", "tmax", "tmin", "annual_temp_range",
 #                   "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
 #                   "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
 #                   "annual_precip", "precip_wettest_month", "precip_driest_month",
 #                   "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
 #                   "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q")
 #   
 #   # Generate all combinations of 3 predictors
 #   all_combos <- combn(predictors, 3, simplify = FALSE)
 #   
 #   # Fit models with model selection and VIF checks
 #   vif_threshold <- 5
 #   model_results <- list()
 #   
 #   for(pred_combo in all_combos) {
 #     # start with 3-predictor model
 #     formula_str <- paste("onset ~", paste(pred_combo, collapse = " + "))
 #     full_model <- glm(as.formula(formula_str), data = phen_estimate)
 #     
 #     # apply stepwise selection
 #     best_model <- MASS::stepAIC(full_model, direction = "both", trace = FALSE)
 #     
 #     # count number of predictors
 #     n_predictors <- length(coef(best_model)) - 1 # minus one to disclude intercept 
 #     
 #     # Check VIF only if 2+ predictors
 #     passes_vif <- TRUE
 #     max_vif <- NA
 #     
 #     if(n_predictors >= 2){
 #       vif_results <- performance::check_collinearity(best_model)
 #       max_vif <- max(vif_results$VIF)
 #       passes_vif <- all(vif_results$VIF < vif_threshold)
 #       
 #       if(!passes_vif){
 #         cat("Excluding model with predictors:", paste(names(coef(best_model))[-1], collapse = ", "),
 #             "- Max VIF =", round(max_vif, 2), "\n")
 #       }
 #     }
 #     # If 1 predictor, automatically passes (no collinearity possible)
 #     
 #     if(passes_vif){
 #       model_results[[length(model_results) + 1]] <- list(
 #         model = best_model,
 #         aic = AIC(best_model),
 #         predictors = names(coef(best_model))[-1], # exclude intercept
 #         formula = formula(best_model),
 #         max_vif = max_vif,
 #         n_predictors = n_predictors, 
 #         pred_signature = paste(sort(names(coef(best_model))[-1]), collapse = "_") # For deduplication
 #       )
 #     }
 #   }
 #   
 #   # Check if we have any valid models
 #   if(length(model_results) == 0){
 #     cat("No models passed VIF threshold of", vif_threshold, "\n")
 #   } else {
 #     cat("\n", length(model_results), "models passed VIF checks\n")
 #     
 #     # DEDUPLICATE models based on predictor signature (extra introduced post-model selection...)
 #     unique_signatures <- unique(sapply(model_results, function(x) x$pred_signature))
 #     unique_models <- lapply(unique_signatures, function(sig) {
 #       # Find first model with this signature
 #       model_results[[which(sapply(model_results, function(x) x$pred_signature == sig))[1]]]
 #     })
 #     
 #     cat("After deduplication:", length(unique_models), "unique models\n")
 #     
 #     # Extract AICs and find best
 #     aics <- sapply(unique_models, function(x) x$aic)
 #     best_aic <- min(aics)
 #     
 #     # Identify models within 2 AIC units of best
 #     top_models_idx <- which(aics - best_aic < 2)
 #     top_models <- unique_models[top_models_idx]
 #     
 #     # Cap to top 10 models as maximum to ensemble 
 #     max_ensemble_size <- 10
 #     if(length(top_models) > max_ensemble_size){
 #       # Sort by AIC and keep best 10
 #       aic_order <- order(sapply(top_models, function(x) x$aic))
 #       top_models <- top_models[aic_order[1:max_ensemble_size]]
 #       cat("Capping ensemble at", max_ensemble_size, "models\n")
 #     }
 #     
 #     cat("Best AIC:", best_aic, "\n")
 #     cat("Number of models within 2 AIC units:", length(top_models), "\n")
 #     
 #     # Create env_model
 #     if(length(top_models) == 1) {
 #       # Single best model
 #       env_model <- top_models[[1]]$model
 #       cat("Using single best model with", top_models[[1]]$n_predictors, "predictor(s):", 
 #           paste(top_models[[1]]$predictors, collapse = ", "), "\n")
 #     } else {
 #       # Ensemble of top models by AIC weighting (lower is higher weight)
 #       cat("Ensembling", length(top_models), "models\n")
 #       delta_aic <- sapply(top_models, function(x) x$aic - best_aic)
 #       aic_weights <- exp(-0.5 * delta_aic)
 #       aic_weights <- aic_weights / sum(aic_weights)
 #       for(i in seq_along(top_models)) {
 #         vif_str <- if(is.na(top_models[[i]]$max_vif)) "N/A" else round(top_models[[i]]$max_vif, 2)
 #         cat("  Model", i, "(", top_models[[i]]$n_predictors, "pred):", 
 #             paste(top_models[[i]]$predictors, collapse = ", "), 
 #             "(AIC =", round(top_models[[i]]$aic, 2), ", Weight =", round(aic_weights[i], 3), ")\n")
 #       }
 #       
 #       # Create ensemble model object
 #       env_model <- structure(
 #         list(
 #           models = lapply(top_models, function(x) x$model),
 #           weights = aic_weights, # Use AIC weights instead of equal weights
 #           call = call("ensemble"),
 #           coefficients = top_models[[1]]$model$coefficients
 #         ),
 #         class = c("ensemble_glm", "list")
 #       )
 #       
 #       # Custom predict method for ensemble with weights
 #       predict.ensemble_glm <- function(object, newdata = NULL, type = "response", ...) {
 #         predictions <- sapply(object$models, function(model) {
 #           predict(model, newdata = newdata, type = type, ...)
 #         })
 #         if(is.matrix(predictions)) {
 #           # Weighted average across models
 #           predictions %*% object$weights
 #         } else {
 #           sum(predictions * object$weights)
 #         }
 #       }
 #     }
 #     
 #     # Generate predictions and residuals
 #     phen_estimate$predicted <- predict(env_model, type = "response")
 #     phen_estimate$residuals <- phen_estimate$onset - phen_estimate$predicted
 #     
 #     # Calculate metrics
 #     res_env <- phen_estimate$residuals
 #     rmse_val <- sqrt(mean(res_env^2))
 #     
 #     # Calculate R²
 #     ss_res <- sum(res_env^2)
 #     ss_tot <- sum((phen_estimate$onset - mean(phen_estimate$onset))^2)
 #     r2_val <- 1 - (ss_res / ss_tot)
 #     
 #     cat("R² =", round(r2_val, 3), "| RMSE =", round(rmse_val, 3), "\n")
 #   }
 # }
 if(nrow(phen_estimate) > 3){ # min req to make env
   predictors <- c("hopkins_delay_midpoint", "avg_temp", "avg_diurnal_temp_range", "isothermality",
                   "temp_seas", "tmax", "tmin", "annual_temp_range",
                   "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                   "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                   "annual_precip", "precip_wettest_month", "precip_driest_month",
                   "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                   "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q")
   
   # Scale predictors using scale()
   phen_estimate_scaled <- phen_estimate
   scaling_params <- list()
   
   for(pred in predictors) {
     if(pred %in% names(phen_estimate)) {
       scaled_var <- scale(phen_estimate[[pred]])
       phen_estimate_scaled[[pred]] <- as.vector(scaled_var)
       scaling_params[[pred]] <- list(
         center = attr(scaled_var, "scaled:center"),
         scale = attr(scaled_var, "scaled:scale")
       )
     }
   }
   
   # Generate all combinations of 3 predictors
   all_combos <- combn(predictors, 3, simplify = FALSE)
   
   # Fit models with model selection and VIF checks
   vif_threshold <- 5
   model_results <- list()
   
   for(pred_combo in all_combos) {
     # Start with 3-predictor model (using scaled data)
     formula_str <- paste("onset ~", paste(pred_combo, collapse = " + "))
     full_model <- glm(as.formula(formula_str), data = phen_estimate_scaled)
     
     # Apply stepwise selection
     best_model <- MASS::stepAIC(full_model, direction = "both", trace = FALSE)
     
     # Count number of predictors
     n_predictors <- length(coef(best_model)) - 1
     
     # Check VIF only if 2+ predictors
     passes_vif <- TRUE
     max_vif <- NA
     
     if(n_predictors >= 2){
       vif_results <- performance::check_collinearity(best_model)
       max_vif <- max(vif_results$VIF)
       passes_vif <- all(vif_results$VIF < vif_threshold)
       
       if(!passes_vif){
         cat("Excluding model with predictors:", paste(names(coef(best_model))[-1], collapse = ", "),
             "- Max VIF =", round(max_vif, 2), "\n")
       }
     }
     
     if(passes_vif){
       model_results[[length(model_results) + 1]] <- list(
         model = best_model,
         aic = AIC(best_model),
         predictors = names(coef(best_model))[-1],
         formula = formula(best_model),
         max_vif = max_vif,
         n_predictors = n_predictors,
         pred_signature = paste(sort(names(coef(best_model))[-1]), collapse = "_")
       )
     }
   }
   
   # Check if we have any valid models
   if(length(model_results) == 0){
     cat("No models passed VIF threshold of", vif_threshold, "\n")
   } else {
     cat("\n", length(model_results), "models passed VIF checks\n")
     
     # DEDUPLICATE models
     unique_signatures <- unique(sapply(model_results, function(x) x$pred_signature))
     unique_models <- lapply(unique_signatures, function(sig) {
       model_results[[which(sapply(model_results, function(x) x$pred_signature == sig))[1]]]
     })
     
     cat("After deduplication:", length(unique_models), "unique models\n")
     
     # Extract AICs and find best
     aics <- sapply(unique_models, function(x) x$aic)
     best_aic <- min(aics)
     
     # Identify models within 2 AIC units of best
     top_models_idx <- which(aics - best_aic < 2)
     top_models <- unique_models[top_models_idx]
     
     # Cap at top 10 models
     max_ensemble_size <- 10
     if(length(top_models) > max_ensemble_size){
       aic_order <- order(sapply(top_models, function(x) x$aic))
       top_models <- top_models[aic_order[1:max_ensemble_size]]
       cat("Capping ensemble at", max_ensemble_size, "models\n")
     }
     
     cat("Best AIC:", best_aic, "\n")
     cat("Number of unique models within 2 AIC units:", length(top_models), "\n")
     
     # Check for high variability if multiple models
     use_single_model <- FALSE
     if(length(top_models) > 1) {
       cat("\nChecking coefficient stability across ensemble candidates:\n")
       
       # Thresholds for instability
       cv_threshold <- 50  # CV > 50% is concerning
       iqr_threshold <- 0.5  # IQR/mean > 0.5 is concerning
       
       all_predictors <- unique(unlist(lapply(top_models, function(x) x$predictors)))
       high_variability_found <- FALSE
       
       for(pred in all_predictors) {
         coefs <- sapply(top_models, function(x) {
           if(pred %in% names(coef(x$model))) {
             coef(x$model)[pred]
           } else {
             NA
           }
         })
         coefs <- coefs[!is.na(coefs)]
         
         if(length(coefs) > 1) {
           mean_coef <- mean(coefs)
           sd_coef <- sd(coefs)
           iqr_coef <- IQR(coefs)
           q25 <- quantile(coefs, 0.25)
           q75 <- quantile(coefs, 0.75)
           
           # Calculate CV only if mean is not near zero
           cv <- NA
           cv_high <- FALSE
           if(abs(mean_coef) > 0.01) {
             cv <- (sd_coef / abs(mean_coef)) * 100
             cv_high <- cv > cv_threshold
           } else {
             cv_high <- TRUE  # Near-zero mean with variation is unstable
           }
           
           # Calculate relative IQR
           iqr_relative <- ifelse(abs(mean_coef) > 0.01, iqr_coef / abs(mean_coef), NA)
           iqr_high <- !is.na(iqr_relative) && iqr_relative > iqr_threshold
           
           # Report
           cat("  ", pred, ":\n")
           if(!is.na(cv)) {
             cat("    Mean =", round(mean_coef, 3), ", SD =", round(sd_coef, 3), 
                 ", CV =", round(cv, 1), "%")
             if(cv_high) cat(" [HIGH CV!]")
             cat("\n")
           } else {
             cat("    Mean ≈ 0 (", round(mean_coef, 4), "), SD =", round(sd_coef, 3))
             cat(" [UNSTABLE - near-zero mean]\n")
           }
           
           cat("    IQR =", round(iqr_coef, 3), " (Q1 =", round(q25, 3), 
               ", Q3 =", round(q75, 3), ")")
           if(iqr_high) cat(" [HIGH IQR!]")
           cat("\n")
           
           # Flag if either metric shows high variability
           if(cv_high || iqr_high) {
             high_variability_found <- TRUE
           }
         } else if(length(coefs) == 1) {
           cat("  ", pred, ": appears in only 1 model (coef =", round(coefs, 3), ")\n")
         }
       }
       
       # Decision: use single model if high variability detected
       if(high_variability_found) {
         cat("\n⚠ High coefficient variability detected across ensemble.\n")
         cat("Selecting single best model instead of ensembling.\n")
         use_single_model <- TRUE
       } else {
         cat("\n Coefficient stability acceptable across ensemble.\n")
       }
       cat("\n")
     }
     
     # Create env_model
     if(length(top_models) == 1 || use_single_model) {
       # Single best model (either only 1 candidate or forced by instability)
       if(use_single_model) {
         # Select the single best model by AIC
         best_idx <- which.min(sapply(top_models, function(x) x$aic))
         top_models <- list(top_models[[best_idx]])
       }
       
       env_model <- top_models[[1]]$model
       env_model$scaling_params <- scaling_params
       class(env_model) <- c("scaled_glm", class(env_model))
       
       cat("Using single best model with", top_models[[1]]$n_predictors, "predictor(s):",
           paste(top_models[[1]]$predictors, collapse = ", "), "\n")
       cat("AIC =", round(top_models[[1]]$aic, 2), "\n")
       
       # Custom predict for single scaled model
       predict.scaled_glm <- function(object, newdata = NULL, type = "response", ...) {
         if(!is.null(newdata)) {
           newdata_scaled <- newdata
           for(pred in names(object$scaling_params)) {
             if(pred %in% names(newdata)) {
               newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                 object$scaling_params[[pred]]$scale
             }
           }
           predict.glm(object, newdata = newdata_scaled, type = type, ...)
         } else {
           predict.glm(object, type = type, ...)
         }
       }
     } else {
       # Ensemble of top models
       cat("Ensembling", length(top_models), "models\n")
       
       delta_aic <- sapply(top_models, function(x) x$aic - best_aic)
       aic_weights <- exp(-0.5 * delta_aic)
       aic_weights <- aic_weights / sum(aic_weights)
       
       for(i in seq_along(top_models)) {
         vif_str <- if(is.na(top_models[[i]]$max_vif)) "N/A" else round(top_models[[i]]$max_vif, 2)
         cat("  Model", i, "(", top_models[[i]]$n_predictors, "pred):",
             paste(top_models[[i]]$predictors, collapse = ", "),
             "(AIC =", round(top_models[[i]]$aic, 2), ", Weight =", round(aic_weights[i], 3), ")\n")
       }
       
       # Create ensemble model object
       env_model <- structure(
         list(
           models = lapply(top_models, function(x) x$model),
           weights = aic_weights,
           scaling_params = scaling_params,
           call = call("ensemble"),
           coefficients = top_models[[1]]$model$coefficients
         ),
         class = c("ensemble_glm", "list")
       )
       
       # Custom predict method for ensemble
       predict.ensemble_glm <- function(object, newdata = NULL, type = "response", ...) {
         if(!is.null(newdata)) {
           newdata_scaled <- newdata
           for(pred in names(object$scaling_params)) {
             if(pred %in% names(newdata)) {
               newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                 object$scaling_params[[pred]]$scale
             }
           }
         } else {
           newdata_scaled <- NULL
         }
         
         predictions <- sapply(object$models, function(model) {
           predict(model, newdata = newdata_scaled, type = type, ...)
         })
         
         if(is.matrix(predictions)) {
           predictions %*% object$weights
         } else {
           sum(predictions * object$weights)
         }
       }
     }
     
     # Generate predictions and residuals
     phen_estimate$predicted <- predict(env_model, newdata = phen_estimate, type = "response")
     phen_estimate$residuals <- phen_estimate$onset - phen_estimate$predicted
     
     # Calculate metrics
     res_env <- phen_estimate$residuals
     rmse_val <- sqrt(mean(res_env^2))
     
     # Calculate R²
     ss_res <- sum(res_env^2)
     ss_tot <- sum((phen_estimate$onset - mean(phen_estimate$onset))^2)
     r2_val <- 1 - (ss_res / ss_tot)
     
     cat("R² =", round(r2_val, 3), "| RMSE =", round(rmse_val, 3), "\n")
   }
 }
## Visualize Outputs 
# Hopkins Only
hop_plot1 <- ggplot(phen_estimate, aes(x = hopkins_delay_midpoint, y = onset)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = min(phen_estimate$hopkins_delay_midpoint), y = max(phen_estimate$onset),
           label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_hop[[2]], rmse_val_hop),
           hjust = 0, vjust = 1) +
  labs(x = "Hopkins Expected Delay Bin (Middle Value for Bin)",
       y = "Observed Onset (day of year)",
       title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
  theme_minimal()
hop_plot2 <- ggplot(phen_estimate, aes(x = corrected_onset_prediction, y = onset)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  #geom_ribbon(aes(ymin = corrected_onset_lower, ymax = corrected_onset_upper), 
              #alpha = 0.2, fill = "blue") +
  geom_point(size = 3) +
  annotate("text", x = min(phen_estimate$onset), y = max(phen_estimate$corrected_onset_prediction),
           label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_hop[[2]], rmse_val_hop),
           hjust = 0, vjust = 1) +
  labs(x = "Corrected Predicted Onset (anchored + actual slope)",
       y = "Empirical Estimate Onset (day of year)",
       title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
  theme_minimal()
# Env Space
if(nrow(phen_estimate) > 6){
env_plot <- ggplot(phen_estimate, aes(x = onset, y = predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = min(phen_estimate$onset), y = max(phen_estimate$predicted),
           label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_val, rmse_val),
           hjust = 0, vjust = 1) +
  labs(x = "Observed Onset", y = "Predicted Onset", 
       title = paste("Using Climatic Vars to Predict", target_taxa[1])) +
  theme_minimal()

}

if(nrow(phen_estimate) > 3){
# return everything as list
return(list(
  hopkins_model = hopkins_model, 
  env_model = env_model, 
  data = phen_data,
  phen_est_data = phen_estimate, 
  plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = env_plot), 
  metrics = list(
    hopkins_r2 = r2_hop[[2]], 
    hopkins_rsme = rmse_val_hop, 
    env_r2 = r2_val, 
    env_rsme = rmse_val
  )
))
} else{
  # return everything as list
  return(list(
    hopkins_model = hopkins_model, 
    env_model = NULL, 
    data = phen_data,
    phen_est_data = phen_estimate, 
    plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = NULL), 
    metrics = list(
      hopkins_r2 = r2_hop[[2]], 
      hopkins_rsme = rmse_val_hop, 
      env_r2 = NULL, 
      env_rsme = NULL
    )
  ))
  
}



}



################################################################################

### Assess Species Histogram
taxon_phen_hist <- function(target_taxa, bees_or_plants = "bees", described_phen_months = NULL){
  
  ## Get Phenology Data
  if(bees_or_plants == "plants"){
    phen_data <- grab_flowering_occs(target_taxa, 
                                     phenovision_flowering_data = phenovision_flowering_data, 
                                     cch2_annotated_specimens = cch2_annotated_specimens)
  } else {
    phen_data <- grab_bee_occs(target_taxa, 
                               bee_occs = bee_occs, 
                               extra_hololeuca_data = extra_hololeuca_data)
  }
  
  ## KDE Smoothing (circular for day of year)
  days <- phen_data$dayOfYear
  # Extend data circularly for boundary handling
  days_circular <- c(days - 365, days, days + 365)
  
  # Compute density
  kde <- density(days_circular, bw = 15, from = 1, to = 365, n = 365)
  kde_df <- data.frame(day = kde$x, density = kde$y)
  
  # Find local minima
  d <- diff(kde_df$density)
  minima_idx <- which(d[-1] > 0 & d[-length(d)] < 0) + 1
  minima_days <- kde_df$day[minima_idx]
  
  cat("Local minima (day of year):", round(minima_days), "\n")
  
  ## Optional: Known phenology periods
  flowering_rects <- NULL
  if(!is.null(described_phen_months)){
    target_row <- described_phen_months[species == target_taxa[1]]
    if(nrow(target_row) > 0){
      month_starts <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
      month_ends <- c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
      month_cols <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      
      flowering_months <- which(as.numeric(target_row[, ..month_cols]) == 1)
      if(length(flowering_months) > 0){
        flowering_rects <- data.frame(
          xmin = month_starts[flowering_months],
          xmax = month_ends[flowering_months]
        )
      }
    }
  }
  
  ## Scale KDE for overlay
  max_count <- max(table(phen_data$dayOfYear))
  kde_df$density_scaled <- kde_df$density * max_count / max(kde_df$density) * 0.9
  
  # Create minima segments for legend
  minima_df <- data.frame(
    x = minima_days,
    xend = minima_days,
    y = -Inf,
    yend = Inf
  )
  
  ## Plot
  plot <- ggplot(data = phen_data) +
    {if(!is.null(flowering_rects)) 
      geom_rect(data = flowering_rects,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, 
                    fill = "Reported Phenophase"),
                alpha = 0.3)} +
    geom_bar(aes(x = dayOfYear), alpha = 0.6) +
    geom_line(data = kde_df, aes(x = day, y = density_scaled, 
                                 color = "Smoothed Phenophase"), 
              linewidth = 1.2) +
    geom_segment(data = minima_df, 
                 aes(x = x, xend = xend, y = y, yend = yend,
                     linetype = "Local Day of Year Minima"),
                 color = "blue", alpha = 0.7) +
    geom_label(data = data.frame(x = minima_days, y = Inf),
               aes(x = x, y = y, label = round(x)),
               vjust = 1.5, color = "blue", size = 3) +
    scale_color_manual(values = c("Smoothed Phenophase" = "darkred"), name = NULL) +
    scale_fill_manual(values = c("Reported Phenophase" = "goldenrod"), name = NULL) +
    scale_linetype_manual(values = c("Local Day of Year Minima" = "dashed"), name = NULL) +
    labs(y = "Count", x = "Day of Year", title = paste(target_taxa[1], "Cumulative PhenoPhase Diagnostic Plot")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(plot)
}
# taxon_phen_hist <- function(target_taxa, bees_or_plants = "bees", described_phen_months = NULL){
#   
#   ## Get Phenology Data
#   if(bees_or_plants == "plants"){
#     phen_data <- grab_flowering_occs(target_taxa, 
#                                      phenovision_flowering_data = phenovision_flowering_data, 
#                                      cch2_annotated_specimens = cch2_annotated_specimens)
#   } else {
#     phen_data <- grab_bee_occs(target_taxa, 
#                                bee_occs = bee_occs, 
#                                extra_hololeuca_data = extra_hololeuca_data)
#   }
#   
#   ## KDE Smoothing (circular for day of year)
#   days <- phen_data$dayOfYear
#   # Extend data circularly for boundary handling
#   days_circular <- c(days - 365, days, days + 365)
#   
#   # Compute density
#   kde <- density(days_circular, bw = 15, from = 1, to = 365, n = 365)
#   kde_df <- data.frame(day = kde$x, density = kde$y)
#   
#   # Find local minima
#   d <- diff(kde_df$density)
#   minima_idx <- which(d[-1] > 0 & d[-length(d)] < 0) + 1
#   minima_days <- kde_df$day[minima_idx]
#   
#   cat("Local minima (day of year):", round(minima_days), "\n")
#   
#   ## Optional: Known phenology periods
#   flowering_rects <- NULL
#   if(!is.null(described_phen_months)){
#     target_row <- described_phen_months[species == target_taxa[1]]
#     if(nrow(target_row) > 0){
#       month_starts <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
#       month_ends <- c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
#       month_cols <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
#                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#       
#       flowering_months <- which(as.numeric(target_row[, ..month_cols]) == 1)
#       if(length(flowering_months) > 0){
#         flowering_rects <- data.frame(
#           xmin = month_starts[flowering_months],
#           xmax = month_ends[flowering_months]
#         )
#       }
#     }
#   }
#   
#   ## Scale KDE for overlay
#   max_count <- max(table(phen_data$dayOfYear))
#   kde_df$density_scaled <- kde_df$density * max_count / max(kde_df$density) * 0.9
#   
#   ## Plot
#   plot <- ggplot(data = phen_data) +
#     {if(!is.null(flowering_rects)) 
#       geom_rect(data = flowering_rects,
#                 aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, 
#                     fill = "Reported Phenophase"),
#                 alpha = 0.3)} +
#     geom_bar(aes(x = dayOfYear), alpha = 0.6) +
#     geom_line(data = kde_df, aes(x = day, y = density_scaled, 
#                                  color = "Smoothed Phenophase"), 
#               linewidth = 1.2) +
#     geom_vline(xintercept = minima_days, 
#                aes(linetype = "Local Day of Year Minima"),
#                color = "blue", alpha = 0.7) +
#     geom_label(data = data.frame(x = minima_days, y = Inf),
#                aes(x = x, y = y, label = round(x)),
#                vjust = 1.5, color = "blue", size = 3) +
#     scale_color_manual(values = c("Smoothed Phenophase" = "darkred"), name = NULL) +
#     scale_fill_manual(values = c("Reported Phenophase" = "goldenrod"), name = NULL) +
#     scale_linetype_manual(values = c("Local Day of Year Minima" = "dashed"), name = NULL) +
#     labs(y = "Count", x = "Day of Year", title = paste(target_taxa[1], "Cumulative PhenoPhase Diagnostic Plot")) +
#     theme_bw() +
#     theme(legend.position = "bottom")
#   return(plot)
# }


### Create Predictions for full distributions ##################################
predict_pheno_for_dist <- function(model_out, 
                                   path_to_target_taxon_PA_raster = NULL,
                                   hopkins = TRUE,
                                   env = TRUE){
  
  ## Pull in dist model, limit to presence range only, apply modified hopkins linear model to P pixels 
  taxon_rast <- rast(paste0(path_to_target_taxon_PA_raster))
  taxon_rast[taxon_rast  == 0] <- NA # remove non-presence pts
  
  # reproj to the native wgs84 coordinate system used to build hopkins expectations
  taxon_rast_reproj <- project(taxon_rast, "EPSG:4326", method = "near")
  
  # make into xy df 
  taxon_df <- as.data.frame(taxon_rast_reproj, xy = TRUE, na.rm = TRUE)
  taxon_df <- taxon_df %>% 
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
  taxon_vect <- vect(taxon_df, 
                         geom = c("longitude", "latitude"),
                         crs = "EPSG:4326")
  taxon_vect_reproj1 <- project(taxon_vect, crs(na_elev))
  taxon_vect_reproj2 <- project(taxon_vect, crs(avg_temp))
  
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
  elevation_vals <- terra::extract(na_elev, taxon_vect_reproj1)
  climate_vals <- terra::extract(climate_stack, taxon_vect_reproj2)
  
  # Add all values to these data
  taxon_df[, elevation_m := elevation_vals[,2]]
  taxon_df[, names(climate_vals)[-1] := climate_vals[,-1]]
  
  # Make hopkins bins a available predictor 
  ## Using this modified dataframe, recreate hopkins based on the anchor from the empirical data
  df <- model_out$data
  # bring in min vals 
  min_empirical_lat <- min(df$latitude)
  min_empirical_elev <- min(df$elevation_m)
  # these will be zero, we'll just adjust so we're going into negative space
  # create expected # of delay days based on hopkins
  taxon_df[, hopkins_delay_days :=
                 ((latitude - min_lat)* 4) + # latitudnal effect
                 ((elevation_m - min_elev))/120 * 4] # elevation effect
  
  # Create bins for phenology data based on Hopkins delay by intervals of 4
  taxon_df[, hopkins_bin := cut(hopkins_delay_days, 
                                    breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
  # additionally, we want to know what pixels are outside of the range of empirical data training
  min_empirical_hopkins <- min(df$hopkins_delay_days)
  max_empirical_hopkins <- max(df$hopkins_delay_days)
  # binary
  taxon_df[, outside_empirical := hopkins_delay_days < min_empirical_hopkins |
                 hopkins_delay_days > max_empirical_hopkins]
  # numeric dist
  taxon_df[, dist_from_empirical := fifelse(
    hopkins_delay_days < min_empirical_hopkins, 
    min_empirical_hopkins - hopkins_delay_days, # below range
    fifelse(
      hopkins_delay_days > max_empirical_hopkins, 
      hopkins_delay_days - max_empirical_hopkins, # above range
      0 #if within range
    )
  )]
  
  # extract midpoint of each hopkins bin to use as predictor
  taxon_df <- taxon_df %>%
    mutate(
      # Extract numeric values from bin intervals
      bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
      bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
      hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
    )
  # Some pixels lack elevation, which is why we get this warning. Probably try and see if we can fix this at a later time. 
  taxon_df <- taxon_df %>% filter(!is.na(hopkins_bin))
  
  ## Create Predictions 
  # apply hopkins linear model 
  hopkins_linear <- model_out$hopkins_model
  taxon_df_predict <- predict(hopkins_linear, newdata = taxon_df)
  taxon_df[, hopkins_predict_onset := taxon_df_predict]
  # apply env glm 
  env_model <- model_out$env_model
  taxon_df_predict <- predict(env_model, newdata = taxon_df)
  taxon_df[, env_predict_onset := taxon_df_predict]
  # further restrict env predictions by CLAMP
  predictor_vars <- attr(terms(env_model), "term.labels")
  df <- model_out$data
  training_vars <- df %>% 
    ungroup() %>%
    select(all_of(predictor_vars)) %>% 
    distinct() %>% # remove redundancy
    as.data.frame()
  
  # calc ranges (made flexible for different predictor_var combos)
  train_ranges <- data.frame(
    var = predictor_vars, 
    min = sapply(predictor_vars, function(v) min(training_vars[[v]])), 
    max = sapply(predictor_vars, function(v) max(training_vars[[v]]))
  )
  
  # clamp!!
  for(i in 1:nrow(train_ranges)){ # for each var
    var_name <- train_ranges$var[i]
    clamped_name <- paste0(var_name, "_clamped")
    taxon_df[, (clamped_name) := pmax(train_ranges$min[i], 
                                          pmin(train_ranges$max[i], get(var_name)))]
  }
  # create clamped prediction df
  clamped_vars <- paste0(predictor_vars, "_clamped")
  clamped_data <- taxon_df[, ..clamped_vars]
  setnames(clamped_data, clamped_vars, predictor_vars)
  
  # make clamped predictions
  taxon_df[, env_predict_onset_clamped := predict(env_model, newdata = clamped_data)]
  
  # Create plots 
  ## get basemap (will need to edit to be more flexible l8ter)
  basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
  basemap <- basemap %>% filter(name %in% c("California", "Nevada", "Arizona", "Sonora", "Baja California", "Baja California Sur"))
  basemap_reproj <- basemap %>% st_transform(st_crs(taxon_rast_reproj))
  hop_pred_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = hopkins_predict_onset)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Hopkins Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Adjusted Hopkins Onset Predictions for", unique(df$species)))
  
  hop_training_dist_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = dist_from_empirical)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Hopkins Empirical Training Distance Eval", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Training Data Distance on Hopkins Onset Predictions for", unique(df$species)))
  
  env_pred_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_onset)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Climatic Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Climate GLM Onset Prediction Plot for", unique(df$species)))
  
  clamped_env_pred_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_onset_clamped)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Clamped Climatic Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Clamped Climate GLM Onset Prediction Plot for", unique(df$species)))
  
## Calculate MESS 
  if(env == TRUE){
  # crop rasters down to actual P area of taxon
  taxon_climate_proj <- project(taxon_rast_reproj, crs(temp_seas), method = "near") # use temp seas as a placeholder for all clims
  climate_crop <- crop(climate_stack[[predictor_vars]], taxon_climate_proj)
  taxon_climate_resamp <- resample(taxon_climate_proj, climate_crop[[1]], method = "near")
  climate_masked <- mask(climate_crop, taxon_climate_resamp)
  # convert to old raster format for dismo::mess()
  climate_stack_old <- raster::stack(climate_masked)
  names(climate_stack_old) <- predictor_vars
  # calc MESS
  mess_raster_old <- dismo::mess(climate_stack_old, training_vars, full = FALSE)
  # convert to terra spatRaster
  mess_raster <- rast(mess_raster_old)
  # extract and append to df 
  mess_vals <- terra::extract(mess_raster, 
                              taxon_df[, .(longitude, latitude)])
  taxon_df[, mess := mess_vals[[2]]]
  
  # plot
  mess_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = mess)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
    theme_minimal() +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("MESS Map for Climate GLM Onset Prediction Plot for", unique(df$species)))
  } else{
    mess_plot <- NULL
  }
  return(list(
    taxon_df = taxon_df,
    plots = list(
      prediction = list(hop_pred_plot, env_pred_plot, clamped_env_pred_plot), 
      quality = list(hop_training_dist_plot, mess_plot)
    )
  ))
  

}













################################################################################

### Second V. adding offset
model_baseline_phenology_w_offset <- function(target_taxa,
                                     bees_or_plants = "plants", # target taxon group to pull empirical data for
                                     phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                     cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                     bee_occs, # input bee occurrence data
                                     extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                     min_doy = 0, # look at histogram when making this decision
                                     max_doy = 365 # ditto to above
){
  
  ## Get Phenology Data
  if(bees_or_plants == "plants"){
    phen_data <- grab_flowering_occs(target_taxa, phenovision_flowering_data = phenovision_flowering_data, cch2_annotated_specimens = cch2_annotated_specimens)
  }
  if(bees_or_plants == "bees"){
    phen_data <- grab_bee_occs(target_taxa, bee_occs = bee_occs, extra_hololeuca_data = extra_hololeuca_data)
  }
  
  ## Attach Elevation & Climatic Variable Data
  # add a row ID to track where we're combing vars with occs
  phen_data[, row_id := .I]
  phen_data_coords <- phen_data[!is.na(longitude) & !is.na(latitude)]
  # load North America Elev Raster & Chelsa annual clim vars
  na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
  avg_temp <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio1_1981-2010_V.2.1.tif")
  avg_diurnal_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio02_1981-2010_V.2.1.tif")
  isothermality <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio03_1981-2010_V.2.1.tif")
  temp_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio4_1981-2010_V.2.1.tif")
  tmax <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio5_1981-2010_V.2.1.tif")
  tmin <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio6_1981-2010_V.2.1.tif")
  annual_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio07_1981-2010_V.2.1.tif")
  avg_monthly_temp_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio08_1981-2010_V.2.1.tif")
  avg_monthly_temp_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio09_1981-2010_V.2.1.tif")
  avg_monthly_temp_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio10_1981-2010_V.2.1.tif")
  avg_monthly_temp_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio11_1981-2010_V.2.1.tif")
  annual_precip <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio12_1981-2010_V.2.1.tif")
  precip_wettest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio13_1981-2010_V.2.1.tif")
  precip_driest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio14_1981-2010_V.2.1.tif")
  precip_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio15_1981-2010_V.2.1.tif")
  avg_monthly_precip_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio16_1981-2010_V.2.1.tif")
  avg_monthly_precip_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio17_1981-2010_V.2.1.tif")
  avg_monthly_precip_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio18_1981-2010_V.2.1.tif")
  avg_monthly_precip_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio19_1981-2010_V.2.1.tif")
  # Monthly vals for termination cues
  pet_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                               pattern = "CHELSA_pet_.*_1981-2010_V.2.1.tif$",
                               full.names = TRUE))
  tasmax_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                               pattern = "CHELSA_tasmax_.*_1981-2010_V.2.1.tif$",
                               full.names = TRUE))
  tasmin_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_tasmin_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  tasmin_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_tasmin_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  wind_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                pattern = "CHELSA_sfcWind_.*_1981-2010_V.2.1.tif$",
                                full.names = TRUE))
  precip_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                pattern = "CHELSA_pr_.*_1981-2010_V.2.1.tif$",
                                full.names = TRUE))
  
  # convert occ data to a spatVect in WGS84
  phen_data_vect <- vect(phen_data_coords, 
                         geom = c("longitude", "latitude"),
                         crs = "EPSG:4326")
  # reproj to elevation raster projection & clim raster projections
  phen_data_vect_reproj1 <- project(phen_data_vect, crs(na_elev))
  phen_data_vect_reproj2 <- project(phen_data_vect, crs(avg_temp)) 
  # Stack all climate layers
  climate_stack <- c(avg_temp, avg_diurnal_temp_range, isothermality,
                     temp_seas, tmax, tmin, annual_temp_range,
                     avg_monthly_temp_wettest_q, avg_monthly_temp_driest_q,
                     avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                     annual_precip, precip_wettest_month, precip_driest_month,
                     precip_seas, avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q,
                     avg_monthly_precip_warmest_q, avg_monthly_precip_coldest_q)
  
  # make names clearer
  names(climate_stack) <- c("avg_temp", "avg_diurnal_temp_range", "isothermality",
                            "temp_seas", "tmax", "tmin", "annual_temp_range",
                            "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                            "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                            "annual_precip", "precip_wettest_month", "precip_driest_month",
                            "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                            "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q")
  
  names(pet_stack) <- c("pet_01", "pet_02", "pet_03", "pet_04", "pet_05", "pet_06", 
                        "pet_07", "pet_08", "pet_09", "pet_10", "pet_11", "pet_12")
  names(tasmax_stack) <- c("tasmax_01", "tasmax_02", "tasmax_03", "tasmax_04", "tasmax_05", "tasmax_06", 
                            "tasmax_07", "tasmax_08", "tasmax_09", "tasmax_10", "tasmax_11", "tasmax_12")
  names(tasmin_stack) <- c("tasmin_01", "tasmin_02", "tasmin_03", "tasmin_04", "tasmin_05", "tasmin_06", 
                           "tasmin_07", "tasmin_08", "tasmin_09", "tasmin_10", "tasmin_11", "tasmin_12")
  names(wind_stack) <- c("sfcWind_01", "sfcWind_02", "sfcWind_03", "sfcWind_04", "sfcWind_05", "sfcWind_06", 
                          "sfcWind_07", "sfcWind_08", "sfcWind_09", "sfcWind_10", "sfcWind_11", "sfcWind_12")
  names(precip_stack) <- c("pr_01", "pr_02", "pr_03", "pr_04", "pr_05", "pr_06", 
                         "pr_07", "pr_08", "pr_09", "pr_10", "pr_11", "pr_12")
  

  
  # Extract all raster values
  elevation_vals <- terra::extract(na_elev, phen_data_vect_reproj1)
  climate_vals <- terra::extract(climate_stack, phen_data_vect_reproj2)
  pet_vals <- terra::extract(pet_stack, phen_data_vect_reproj2)
  tasmax_vals <- terra::extract(tasmax_stack, phen_data_vect_reproj2)
  tasmin_vals <- terra::extract(tasmin_stack, phen_data_vect_reproj2)
  wind_vals <- terra::extract(wind_stack, phen_data_vect_reproj2)
  precip_vals <- terra::extract(precip_stack, phen_data_vect_reproj2)
  
  # Add all values to these data
  phen_data_coords[, elevation_m := elevation_vals[,2]]
  phen_data_coords[, names(climate_vals)[-1] := climate_vals[,-1]]
  phen_data_coords[, names(pet_vals)[-1] := pet_vals[, -1]]
  phen_data_coords[, names(tasmax_vals)[-1] :=tasmax_vals[, -1]]
  phen_data_coords[, names(tasmin_vals)[-1] := tasmin_vals[, -1]]
  phen_data_coords[, names(wind_vals)[-1] := wind_vals[, -1]]
  phen_data_coords[, names(precip_vals)[-1] := precip_vals[, -1]]

  
  
  # add elevation values to occ table using row_id as our key
  phen_data_coords[, elevation_m := elevation_vals[,2]]
  # merge
  phen_data <- phen_data[phen_data_coords[, .(row_id, elevation_m, avg_temp, avg_diurnal_temp_range, isothermality,
                                              temp_seas, tmax, tmin, annual_temp_range,
                                              avg_monthly_temp_wettest_q, avg_monthly_temp_driest_q,
                                              avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                                              annual_precip, precip_wettest_month, precip_driest_month,
                                              precip_seas, avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q,
                                              avg_monthly_precip_warmest_q, avg_monthly_precip_coldest_q, 
                                              pet_01, pet_02, pet_03, pet_04, pet_05, pet_06, 
                                              pet_07, pet_08, pet_09, pet_10, pet_11, pet_12, 
                                              tasmax_01, tasmax_02, tasmax_03, tasmax_04, tasmax_05, tasmax_06, 
                                              tasmax_07, tasmax_08, tasmax_09, tasmax_10, tasmax_11, tasmax_12,
                                              tasmin_01, tasmin_02, tasmin_03, tasmin_04, tasmin_05, tasmin_06, 
                                              tasmin_07, tasmin_08, tasmin_09, tasmin_10, tasmin_11, tasmin_12,
                                              sfcWind_01, sfcWind_02, sfcWind_03, sfcWind_04, sfcWind_05, sfcWind_06, 
                                              sfcWind_07, sfcWind_08, sfcWind_09, sfcWind_10, sfcWind_11, sfcWind_12, 
                                              pr_01, pr_02, pr_03, pr_04, pr_05, pr_06, 
                                              pr_07, pr_08, pr_09, pr_10, pr_11, pr_12)],
                         on = "row_id"]
  # remove NA elevation vals
  phen_data <- phen_data[!is.na(elevation_m) & !is.na(avg_temp) & !is.na(avg_diurnal_temp_range) & !is.na(isothermality) & !is.na(temp_seas) 
                         & !is.na(tmax) & !is.na(tmin) & !is.na(annual_temp_range) & !is.na(avg_monthly_temp_wettest_q) & !is.na(avg_monthly_temp_driest_q) &
                           !is.na(avg_monthly_temp_warmest_q) & !is.na(avg_monthly_temp_coldest_q) & !is.na(annual_precip) &
                           !is.na(precip_wettest_month) & !is.na(precip_driest_month) & !is.na(precip_seas) & !is.na(avg_monthly_precip_wettest_q) &
                           !is.na(avg_monthly_precip_driest_q) & !is.na(avg_monthly_precip_warmest_q) & !is.na(avg_monthly_precip_coldest_q)& 
                         !is.na(pet_01)& !is.na(pet_02)& !is.na(pet_03)& !is.na(pet_04)& !is.na(pet_05)& !is.na(pet_06)& 
                         !is.na(pet_07)& !is.na(pet_08)& !is.na(pet_09)& !is.na(pet_10)& !is.na(pet_11)& !is.na(pet_12)& 
                         !is.na(tasmax_01)& !is.na(tasmax_02)& !is.na(tasmax_03)& !is.na(tasmax_04)& !is.na(tasmax_05)& !is.na(tasmax_06)& 
                         !is.na(tasmax_07)& !is.na(tasmax_08)& !is.na(tasmax_09)& !is.na(tasmax_10)& !is.na(tasmax_11)& !is.na(tasmax_12)&
                         !is.na(tasmin_01)& !is.na(tasmin_02)& !is.na(tasmin_03)& !is.na(tasmin_04)& !is.na(tasmin_05)& !is.na(tasmin_06)& 
                         !is.na(tasmin_07)& !is.na(tasmin_08)& !is.na(tasmin_09)& !is.na(tasmin_10)& !is.na(tasmin_11)& !is.na(tasmin_12)&
                         !is.na(sfcWind_01)& !is.na(sfcWind_02)& !is.na(sfcWind_03)& !is.na(sfcWind_04)& !is.na(sfcWind_05)& !is.na(sfcWind_06)& 
                         !is.na(sfcWind_07)& !is.na(sfcWind_08)& !is.na(sfcWind_09)& !is.na(sfcWind_10)& !is.na(sfcWind_11)& !is.na(sfcWind_12)&
                         !is.na(pr_01)& !is.na(pr_02)& !is.na(pr_03)& !is.na(pr_04)& !is.na(pr_05)& !is.na(pr_06)& 
                         !is.na(pr_07)& !is.na(pr_08)& !is.na(pr_09)& !is.na(pr_10)& !is.na(pr_11)& !is.na(pr_12),]
  
  ## Create Hopkins Expectation: For every degree in latitude or 120 meters of elevation, flowering will be delayed by 4 days. 
  # use min vals as a baseline start
  min_lat <- min(phen_data$latitude)
  min_elev <- min(phen_data$elevation_m)
  ## create expected # of delay days based on hopkins
  phen_data[, hopkins_delay_days :=
              ((latitude - min_lat)* 4) + # latitudnal effect
              ((elevation_m - min_elev))/120 * 4] # elevation effect
  
  # Create bins for phenology data based on Hopkins delay by intervals of 4
  phen_data[, hopkins_bin := cut(hopkins_delay_days, 
                                 breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
  
  ## Run Phenometric Cleaning 
  # remove outliers based on known flowering 
  phen_data <- phen_data %>%  filter(dayOfYear > min_doy & dayOfYear < max_doy)
  # downsample any dayOfYears that have abnormally high sampling on a dayOfYear within a year (City Nature Challenge and BioBlitzs)
  phen_data <- phen_data %>% 
    group_by(species,year,hopkins_bin,dayOfYear) %>% 
    dplyr::slice_sample(n = 3)
  # if sample is either a specimen (Herbarium data) or an occurrence (unspecified in bee data), then remove first day of the months from these data
  # these are possibly a data quality issue
  phen_data <- phen_data %>% 
    filter(!(source %in% c("occurrence", "specimen") & day == 1))
  
  
  ## Create a phenesse estimate for phenometric 
  phen_data_est_onset <- phen_data %>% 
    group_by(species,hopkins_bin) %>% 
    filter(n() >= 7) %>%
    filter(n_distinct(dayOfYear) > 3) %>% 
    group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$dayOfYear, percentile = 0.10, bootstraps=1000)))
  
  phen_data_est_offset <- phen_data %>% 
    group_by(species,hopkins_bin) %>% 
    filter(n() >= 7) %>%
    filter(n_distinct(dayOfYear) > 3) %>% 
    group_modify(~ broom::tidy(phenesse::quantile_ci(observations = .x$dayOfYear, percentile = 0.90, bootstraps=1000)))
  
  # extract estimate
  phen_onset_estimates <- phen_data_est_onset %>% 
    filter(column == "estimate") %>% 
    select(species, hopkins_bin, mean) %>% 
    rename(onset_doy = mean) %>% ungroup()
  
  phen_offset_estimates <- phen_data_est_offset %>% 
    filter(column == "estimate") %>% 
    select(species, hopkins_bin, mean) %>% 
    rename(offset_doy = mean) %>% ungroup() %>% select(-species)
  
  phen_estimates <- phen_onset_estimates %>% 
    left_join(phen_offset_estimates, by = "hopkins_bin")
  
  ### Based on onset estimates, create a window 1 month before and 1 month after (3 months) of average monthly conditions that could influence termination 
  onset_lookup <- phen_estimates %>%
    mutate(onset_month = lubridate::month(as.Date(onset_doy, origin = "1998-01-01"))) %>%
    dplyr::select(onset_month, hopkins_bin) 

  
  # Function to get relevant months (±1 month, wrapping around year)
  get_month_window <- function(center_month) {
    months <- c((center_month - 2) %% 12 + 1, 
                (center_month - 1) %% 12 + 1,
                center_month,
                center_month %% 12 + 1,
                (center_month + 1) %% 12 + 1)
    unique(months[1:5])  # Get -2, 0, +2
  }
  
  phen_data_with_windows <- phen_data %>%
    inner_join(onset_lookup, by = "hopkins_bin") %>%
    rowwise() %>%
    mutate(
      avg_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      } , 
      total_photoperiod_window = {
        m <- get_month_window(onset_month)
        phi <- latitude * pi / 180
        
        total <- 0
        for (month in m) {
          days_in_month <- as.numeric(format(seq(as.Date(paste0("2023-", month, "-01")), 
                                                 by = "month", length.out = 2)[2] - 1, "%d"))
          for (day in 1:days_in_month) {
            doy <- as.numeric(format(as.Date(paste0("2023-", month, "-", day)), "%j"))
            delta <- 0.409 * sin(2 * pi * (doy - 81) / 365)
            w <- acos(-tan(phi) * tan(delta))
            daylength <- 24 * w / pi
            total <- total + daylength
          }
        }
        total
      }
    ) %>%
    ungroup()
  
  
  
  ## Attach env avg summaries per hopkins bin
  hopkins_bin_env_summaries <- phen_data_with_windows %>%
    group_by(species, hopkins_bin) %>%
    filter(n() >= 7) %>%
    filter(n_distinct(dayOfYear) > 3) %>%
    summarize(avg_temp = mean(avg_temp),
              avg_diurnal_temp_range = mean(avg_diurnal_temp_range),
              isothermality = mean(isothermality),
              temp_seas = mean(temp_seas),
              tmax = mean(tmax),
              tmin = mean(tmin),
              annual_temp_range = mean(annual_temp_range),
              avg_monthly_temp_wettest_q = mean(avg_monthly_temp_wettest_q),
              avg_monthly_temp_driest_q = mean(avg_monthly_temp_driest_q),
              avg_monthly_temp_warmest_q = mean(avg_monthly_temp_warmest_q),
              avg_monthly_temp_coldest_q = mean(avg_monthly_temp_coldest_q),
              annual_precip = mean(annual_precip),
              precip_wettest_month = mean(precip_wettest_month),
              precip_driest_month = mean(precip_driest_month),
              precip_seas = mean(precip_seas),
              avg_monthly_precip_wettest_q = mean(avg_monthly_precip_wettest_q),
              avg_monthly_precip_driest_q = mean(avg_monthly_precip_driest_q),
              avg_monthly_precip_warmest_q = mean(avg_monthly_precip_warmest_q),
              avg_monthly_precip_coldest_q = mean(avg_monthly_precip_coldest_q),
              avg_pet_window = mean(avg_pet_window),
              sd_pet_window = mean(sd_pet_window), 
              slope_pet_window = mean(slope_pet_window),
              avg_tasmax_window = mean(avg_tasmax_window),
              sd_tasmax_window = mean(sd_tasmax_window), 
              slope_tasmax_window = mean(slope_tasmax_window), 
              avg_tasmin_window = mean(avg_tasmin_window),
              sd_tasmin_window = mean(sd_tasmin_window), 
              slope_tasmin_window = mean(slope_tasmin_window), 
              avg_sfcWind_window = mean(avg_sfcWind_window), 
              sd_sfcWind_window = mean(sd_sfcWind_window), 
              slope_sfcWind_window = mean(slope_sfcWind_window), 
              avg_pr_window = mean(avg_pr_window), 
              sd_pr_window = mean(sd_pr_window), 
              slope_pr_window = mean(slope_pr_window), 
              total_photoperiod_window = mean(total_photoperiod_window)
              
    ) %>%
    ungroup() %>%
    dplyr::select(hopkins_bin, avg_temp, avg_diurnal_temp_range, isothermality,
                  temp_seas, tmax, tmin, annual_temp_range, avg_monthly_temp_wettest_q,
                  avg_monthly_temp_driest_q, avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                  annual_precip, precip_wettest_month, precip_driest_month, precip_seas,
                  avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q, avg_monthly_precip_warmest_q,
                  avg_monthly_precip_coldest_q, avg_pet_window, sd_pet_window, slope_pet_window, avg_tasmax_window, sd_tasmax_window, slope_tasmax_window,
                  avg_tasmin_window, sd_tasmin_window, slope_tasmin_window, avg_sfcWind_window, sd_sfcWind_window, slope_sfcWind_window, 
                  avg_pr_window, sd_pr_window, slope_pr_window, total_photoperiod_window)
  ## Join 
  phen_estimates <- phen_estimates %>% 
    inner_join(hopkins_bin_env_summaries, by = "hopkins_bin")
  
  ## Remove Outliers 
  phen_estimates <-  phen_estimates %>% 
    group_by(species) %>% 
    mutate(zonset = scale(onset_doy)) %>% 
    filter(between(zonset,-3.25,+3.25))
  
  phen_estimates <-  phen_estimates %>% 
    group_by(species) %>% 
    mutate(zoffset = scale(offset_doy)) %>% 
    filter(between(zoffset,-3.25,+3.25))
  
  
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
  
  # bin summary (for creating hopkins predictions)
  bin_summary <- phen_data %>%
    group_by(species, hopkins_bin) %>%
    summarise(
      hopkins_lat = median(latitude, na.rm = TRUE),
      hopkins_elev = median(elevation_m, na.rm = TRUE),
      hopkins_delay = median(hopkins_delay_days, na.rm = TRUE),
      n_obs = n(),
      .groups = 'drop'
    )
  
  # join bins 
  phen_estimates <- phen_estimates %>% 
    left_join(bin_summary, by = "hopkins_bin") %>% 
    rename(lat_mean = hopkins_lat, 
           mid_elev = hopkins_elev, 
           onset = onset_doy, 
           offset = offset_doy,
           species = species.x) # redundant
  
  anchored_sp <- estimate_species_anchors(phen_estimates)
  
  
  phen_estimate <- phen_estimates %>% 
    left_join(anchored_sp, by = 'species')
  
  # extract midpoint of each hopkins bin to use as predictor
  phen_estimate <- phen_estimate %>%
    mutate(
      # Extract numeric values from bin intervals
      bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
      bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
      hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
    )
  
  ## Run a Hopkins model based on bins
  hopkins_model <- lm(onset ~ hopkins_delay_midpoint, data = phen_estimate)
  r2_hop <- performance::r2(hopkins_model)
  res_hop <- residuals(hopkins_model)
  rmse_val_hop <- sqrt(mean(res_hop^2))
  # Extract the corrected slope (actual delay rate per Hopkins bin unit)
  actual_delay_rate <- coef(hopkins_model)[2]
  actual_delay_se <- summary(hopkins_model)$coefficients[2, 2]
  # Apply the corrected slope with uncertainty bounds
  phen_estimate <- phen_estimate %>%
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
  
  #### ONSET ENV MODEL ######################################################################################
  if(nrow(phen_estimate) > 3){ # min req to make env
    predictors <- c("avg_temp", "avg_diurnal_temp_range", "isothermality",
                    "temp_seas", "tmax", "tmin", "annual_temp_range",
                    "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                    "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                    "annual_precip", "precip_wettest_month", "precip_driest_month",
                    "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                    "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q") # removed hopkins_delay_midpoint...maybe a bad idea?
    
    # Scale predictors using scale()
    phen_estimate_scaled <- phen_estimate
    scaling_params <- list()
    
    for(pred in predictors) {
      if(pred %in% names(phen_estimate)) {
        scaled_var <- scale(phen_estimate[[pred]])
        phen_estimate_scaled[[pred]] <- as.vector(scaled_var)
        scaling_params[[pred]] <- list(
          center = attr(scaled_var, "scaled:center"),
          scale = attr(scaled_var, "scaled:scale")
        )
      }
    }
    
    # Generate all combinations of 3 predictors
    all_combos <- combn(predictors, 3, simplify = FALSE)
    
    # Fit models with model selection and VIF checks
    vif_threshold <- 5
    model_results <- list()
    
    for(pred_combo in all_combos) {
      # Start with 3-predictor model (using scaled data)
      formula_str <- paste("onset ~", paste(pred_combo, collapse = " + "))
      full_model <- glm(as.formula(formula_str), data = phen_estimate_scaled)
      
      # Apply stepwise selection
      best_model <- MASS::stepAIC(full_model, direction = "both", trace = FALSE)
      
      # Count number of predictors
      n_predictors <- length(coef(best_model)) - 1
      
      # Check VIF only if 2+ predictors
      passes_vif <- TRUE
      max_vif <- NA
      
      if(n_predictors >= 2){
        vif_results <- performance::check_collinearity(best_model)
        max_vif <- max(vif_results$VIF)
        passes_vif <- all(vif_results$VIF < vif_threshold)
        
        if(!passes_vif){ # lots of text, was useful while testing.
          # cat("Excluding model with predictors:", paste(names(coef(best_model))[-1], collapse = ", "), 
          #     "- Max VIF =", round(max_vif, 2), "\n")
        }
      }
      
      if(passes_vif){
        model_results[[length(model_results) + 1]] <- list(
          model = best_model,
          aic = AIC(best_model),
          predictors = names(coef(best_model))[-1],
          formula = formula(best_model),
          max_vif = max_vif,
          n_predictors = n_predictors,
          pred_signature = paste(sort(names(coef(best_model))[-1]), collapse = "_")
        )
      }
    }
    
    # Check if we have any valid models
    if(length(model_results) == 0){
      cat("No models passed VIF threshold of", vif_threshold, "\n")
    } else {
      cat("\n", length(model_results), "models passed VIF checks\n")
      
      # DEDUPLICATE models
      unique_signatures <- unique(sapply(model_results, function(x) x$pred_signature))
      unique_models <- lapply(unique_signatures, function(sig) {
        model_results[[which(sapply(model_results, function(x) x$pred_signature == sig))[1]]]
      })
      
      cat("After deduplication:", length(unique_models), "unique models\n")
      
      # Extract AICs and find best
      aics <- sapply(unique_models, function(x) x$aic)
      best_aic <- min(aics)
      
      # Identify models within 2 AIC units of best
      top_models_idx <- which(aics - best_aic < 2)
      top_models <- unique_models[top_models_idx]
      
      # Cap at top 10 models
      max_ensemble_size <- 10
      if(length(top_models) > max_ensemble_size){
        aic_order <- order(sapply(top_models, function(x) x$aic))
        top_models <- top_models[aic_order[1:max_ensemble_size]]
        cat("Capping ensemble at", max_ensemble_size, "models\n")
      }
      
      cat("Best AIC:", best_aic, "\n")
      cat("Number of unique models within 2 AIC units:", length(top_models), "\n")
      
      # Check for high variability if multiple models
      use_single_model <- FALSE
      if(length(top_models) > 1) {
        cat("\nChecking coefficient stability across ensemble candidates:\n")
        
        # Thresholds for instability
        cv_threshold <- 50  # CV > 50% is concerning
        iqr_threshold <- 0.5  # IQR/mean > 0.5 is concerning
        
        all_predictors <- unique(unlist(lapply(top_models, function(x) x$predictors)))
        high_variability_found <- FALSE
        
        for(pred in all_predictors) {
          coefs <- sapply(top_models, function(x) {
            if(pred %in% names(coef(x$model))) {
              coef(x$model)[pred]
            } else {
              NA
            }
          })
          coefs <- coefs[!is.na(coefs)]
          
          if(length(coefs) > 1) {
            mean_coef <- mean(coefs)
            sd_coef <- sd(coefs)
            iqr_coef <- IQR(coefs)
            q25 <- quantile(coefs, 0.25)
            q75 <- quantile(coefs, 0.75)
            
            # Calculate CV only if mean is not near zero
            cv <- NA
            cv_high <- FALSE
            if(abs(mean_coef) > 0.01) {
              cv <- (sd_coef / abs(mean_coef)) * 100
              cv_high <- cv > cv_threshold
            } else {
              cv_high <- TRUE  # Near-zero mean with variation is unstable
            }
            
            # Calculate relative IQR
            iqr_relative <- ifelse(abs(mean_coef) > 0.01, iqr_coef / abs(mean_coef), NA)
            iqr_high <- !is.na(iqr_relative) && iqr_relative > iqr_threshold
            
            # Report
            cat("  ", pred, ":\n")
            if(!is.na(cv)) {
              cat("    Mean =", round(mean_coef, 3), ", SD =", round(sd_coef, 3), 
                  ", CV =", round(cv, 1), "%")
              if(cv_high) cat(" [HIGH CV!]")
              cat("\n")
            } else {
              cat("    Mean ≈ 0 (", round(mean_coef, 4), "), SD =", round(sd_coef, 3))
              cat(" [UNSTABLE - near-zero mean]\n")
            }
            
            cat("    IQR =", round(iqr_coef, 3), " (Q1 =", round(q25, 3), 
                ", Q3 =", round(q75, 3), ")")
            if(iqr_high) cat(" [HIGH IQR!]")
            cat("\n")
            
            # Flag if either metric shows high variability
            if(cv_high || iqr_high) {
              high_variability_found <- TRUE
            }
          } else if(length(coefs) == 1) {
            cat("  ", pred, ": appears in only 1 model (coef =", round(coefs, 3), ")\n")
          }
        }
        
        # Decision: use single model if high variability detected
        if(high_variability_found) {
          cat("\n⚠ High coefficient variability detected across ensemble.\n")
          cat("Selecting single best model instead of ensembling.\n")
          use_single_model <- TRUE
        } else {
          cat("\n Coefficient stability acceptable across ensemble.\n")
        }
        cat("\n")
      }
      
      # Create env_model
      if(length(top_models) == 1 || use_single_model) {
        # Single best model (either only 1 candidate or forced by instability)
        if(use_single_model) {
          # Select the single best model by AIC
          best_idx <- which.min(sapply(top_models, function(x) x$aic))
          top_models <- list(top_models[[best_idx]])
        }
        
        env_model <- top_models[[1]]$model
        env_model$scaling_params <- scaling_params
        class(env_model) <- c("scaled_glm", class(env_model))
        
        cat("Using single best model with", top_models[[1]]$n_predictors, "predictor(s):",
            paste(top_models[[1]]$predictors, collapse = ", "), "\n")
        cat("AIC =", round(top_models[[1]]$aic, 2), "\n")
        
        # Custom predict for single scaled model
        predict.scaled_glm <- function(object, newdata = NULL, type = "response", ...) {
          if(!is.null(newdata)) {
            newdata_scaled <- newdata
            for(pred in names(object$scaling_params)) {
              if(pred %in% names(newdata)) {
                newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                  object$scaling_params[[pred]]$scale
              }
            }
            predict.glm(object, newdata = newdata_scaled, type = type, ...)
          } else {
            predict.glm(object, type = type, ...)
          }
        }
      } else {
        # Ensemble of top models
        cat("Ensembling", length(top_models), "models\n")
        
        delta_aic <- sapply(top_models, function(x) x$aic - best_aic)
        aic_weights <- exp(-0.5 * delta_aic)
        aic_weights <- aic_weights / sum(aic_weights)
        
        for(i in seq_along(top_models)) {
          vif_str <- if(is.na(top_models[[i]]$max_vif)) "N/A" else round(top_models[[i]]$max_vif, 2)
          cat("  Model", i, "(", top_models[[i]]$n_predictors, "pred):",
              paste(top_models[[i]]$predictors, collapse = ", "),
              "(AIC =", round(top_models[[i]]$aic, 2), ", Weight =", round(aic_weights[i], 3), ")\n")
        }
        
        # Create ensemble model object
        env_model <- structure(
          list(
            models = lapply(top_models, function(x) x$model),
            weights = aic_weights,
            scaling_params = scaling_params,
            call = call("ensemble"),
            coefficients = top_models[[1]]$model$coefficients
          ),
          class = c("ensemble_glm", "list")
        )
        
        # Custom predict method for ensemble
        predict.ensemble_glm <- function(object, newdata = NULL, type = "response", ...) {
          if(!is.null(newdata)) {
            newdata_scaled <- newdata
            for(pred in names(object$scaling_params)) {
              if(pred %in% names(newdata)) {
                newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                  object$scaling_params[[pred]]$scale
              }
            }
          } else {
            newdata_scaled <- NULL
          }
          
          predictions <- sapply(object$models, function(model) {
            predict(model, newdata = newdata_scaled, type = type, ...)
          })
          
          if(is.matrix(predictions)) {
            predictions %*% object$weights
          } else {
            sum(predictions * object$weights)
          }
        }
      }
      
      # Generate predictions and residuals
      phen_estimate$predicted <- predict(env_model, newdata = phen_estimate, type = "response")
      phen_estimate$residuals <- phen_estimate$onset - phen_estimate$predicted
      
      # Calculate metrics
      res_env <- phen_estimate$residuals
      rmse_val <- sqrt(mean(res_env^2))
      
      # Calculate R²
      ss_res <- sum(res_env^2)
      ss_tot <- sum((phen_estimate$onset - mean(phen_estimate$onset))^2)
      r2_val <- 1 - (ss_res / ss_tot)
      
      cat("R² =", round(r2_val, 3), "| RMSE =", round(rmse_val, 3), "\n")
    }
  }
  
  ### OFFSET ENV MODEL ##########################################################################
  if(nrow(phen_estimate) > 3){ # min req to make env
    predictors <- c("onset", "avg_pet_window", "sd_pet_window", "slope_pet_window", 
                    "avg_tasmax_window","sd_tasmax_window", "slope_tasmax_window","avg_tasmin_window",           
                    "sd_tasmin_window", "slope_tasmin_window", "avg_sfcWind_window","sd_sfcWind_window", "slope_sfcWind_window", 
                    "avg_pr_window","sd_pr_window", "slope_pr_window", 
                    "avg_temp", "avg_diurnal_temp_range", "isothermality",
                    "temp_seas", "tmax", "tmin", "annual_temp_range",
                    "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                    "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                    "annual_precip", "precip_wettest_month", "precip_driest_month",
                    "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                    "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q", "total_photoperiod_window") # chopped out "hopkins_delay_midpoint" maybe a bad idea?
    
    # Scale predictors using scale()
    phen_estimate_scaled <- phen_estimate
    scaling_params <- list()
    
    for(pred in predictors) {
      if(pred %in% names(phen_estimate)) {
        scaled_var <- scale(phen_estimate[[pred]])
        phen_estimate_scaled[[pred]] <- as.vector(scaled_var)
        scaling_params[[pred]] <- list(
          center = attr(scaled_var, "scaled:center"),
          scale = attr(scaled_var, "scaled:scale")
        )
      }
    }
    
    # Generate all combinations of 3 predictors
    all_combos <- combn(predictors, 3, simplify = FALSE)
    
    # Fit models with model selection and VIF checks
    vif_threshold <- 5
    model_results <- list()
    
    for(pred_combo in all_combos) {
      # Start with 3-predictor model (using scaled data)
      formula_str <- paste("offset ~", paste(pred_combo, collapse = " + "))
      full_model <- glm(as.formula(formula_str), data = phen_estimate_scaled)
      
      # Apply stepwise selection
      best_model <- MASS::stepAIC(full_model, direction = "both", trace = FALSE)
      
      # Count number of predictors
      n_predictors <- length(coef(best_model)) - 1
      
      # Check VIF only if 2+ predictors
      passes_vif <- TRUE
      max_vif <- NA
      
      if(n_predictors >= 2){
        vif_results <- performance::check_collinearity(best_model)
        max_vif <- max(vif_results$VIF)
        passes_vif <- all(vif_results$VIF < vif_threshold)
        
        if(!passes_vif){ # lots of text like b4
         # cat("Excluding model with predictors:", paste(names(coef(best_model))[-1], collapse = ", "),
              #"- Max VIF =", round(max_vif, 2), "\n")
        }
      }
      
      if(passes_vif){
        model_results[[length(model_results) + 1]] <- list(
          model = best_model,
          aic = AIC(best_model),
          predictors = names(coef(best_model))[-1],
          formula = formula(best_model),
          max_vif = max_vif,
          n_predictors = n_predictors,
          pred_signature = paste(sort(names(coef(best_model))[-1]), collapse = "_")
        )
      }
    }
    
    # Check if we have any valid models
    if(length(model_results) == 0){
      cat("No models passed VIF threshold of", vif_threshold, "\n")
    } else {
      cat("\n", length(model_results), "models passed VIF checks\n")
      
      # DEDUPLICATE models
      unique_signatures <- unique(sapply(model_results, function(x) x$pred_signature))
      unique_models <- lapply(unique_signatures, function(sig) {
        model_results[[which(sapply(model_results, function(x) x$pred_signature == sig))[1]]]
      })
      
      cat("After deduplication:", length(unique_models), "unique models\n")
      
      # Extract AICs and find best
      aics <- sapply(unique_models, function(x) x$aic)
      best_aic <- min(aics)
      
      # Identify models within 2 AIC units of best
      top_models_idx <- which(aics - best_aic < 2)
      top_models <- unique_models[top_models_idx]
      
      # Cap at top 10 models
      max_ensemble_size <- 10
      if(length(top_models) > max_ensemble_size){
        aic_order <- order(sapply(top_models, function(x) x$aic))
        top_models <- top_models[aic_order[1:max_ensemble_size]]
        cat("Capping ensemble at", max_ensemble_size, "models\n")
      }
      
      cat("Best AIC:", best_aic, "\n")
      cat("Number of unique models within 2 AIC units:", length(top_models), "\n")
      
      # Check for high variability if multiple models
      use_single_model <- FALSE
      if(length(top_models) > 1) {
        cat("\nChecking coefficient stability across ensemble candidates:\n")
        
        # Thresholds for instability
        cv_threshold <- 50  # CV > 50% is concerning
        iqr_threshold <- 0.5  # IQR/mean > 0.5 is concerning
        
        all_predictors <- unique(unlist(lapply(top_models, function(x) x$predictors)))
        high_variability_found <- FALSE
        
        for(pred in all_predictors) {
          coefs <- sapply(top_models, function(x) {
            if(pred %in% names(coef(x$model))) {
              coef(x$model)[pred]
            } else {
              NA
            }
          })
          coefs <- coefs[!is.na(coefs)]
          
          if(length(coefs) > 1) {
            mean_coef <- mean(coefs)
            sd_coef <- sd(coefs)
            iqr_coef <- IQR(coefs)
            q25 <- quantile(coefs, 0.25)
            q75 <- quantile(coefs, 0.75)
            
            # Calculate CV only if mean is not near zero
            cv <- NA
            cv_high <- FALSE
            if(abs(mean_coef) > 0.01) {
              cv <- (sd_coef / abs(mean_coef)) * 100
              cv_high <- cv > cv_threshold
            } else {
              cv_high <- TRUE  # Near-zero mean with variation is unstable
            }
            
            # Calculate relative IQR
            iqr_relative <- ifelse(abs(mean_coef) > 0.01, iqr_coef / abs(mean_coef), NA)
            iqr_high <- !is.na(iqr_relative) && iqr_relative > iqr_threshold
            
            # Report
            cat("  ", pred, ":\n")
            if(!is.na(cv)) {
              cat("    Mean =", round(mean_coef, 3), ", SD =", round(sd_coef, 3), 
                  ", CV =", round(cv, 1), "%")
              if(cv_high) cat(" [HIGH CV!]")
              cat("\n")
            } else {
              cat("    Mean ≈ 0 (", round(mean_coef, 4), "), SD =", round(sd_coef, 3))
              cat(" [UNSTABLE - near-zero mean]\n")
            }
            
            cat("    IQR =", round(iqr_coef, 3), " (Q1 =", round(q25, 3), 
                ", Q3 =", round(q75, 3), ")")
            if(iqr_high) cat(" [HIGH IQR!]")
            cat("\n")
            
            # Flag if either metric shows high variability
            if(cv_high || iqr_high) {
              high_variability_found <- TRUE
            }
          } else if(length(coefs) == 1) {
            cat("  ", pred, ": appears in only 1 model (coef =", round(coefs, 3), ")\n")
          }
        }
        
        # Decision: use single model if high variability detected
        if(high_variability_found) {
          cat("\n⚠ High coefficient variability detected across ensemble.\n")
          cat("Selecting single best model instead of ensembling.\n")
          use_single_model <- TRUE
        } else {
          cat("\n Coefficient stability acceptable across ensemble.\n")
        }
        cat("\n")
      }
      
      # Create env_model
      if(length(top_models) == 1 || use_single_model) {
        # Single best model (either only 1 candidate or forced by instability)
        if(use_single_model) {
          # Select the single best model by AIC
          best_idx <- which.min(sapply(top_models, function(x) x$aic))
          top_models <- list(top_models[[best_idx]])
        }
        
        env_model_offset <- top_models[[1]]$model
        env_model_offset$scaling_params <- scaling_params
        class(env_model_offset) <- c("scaled_glm", class(env_model_offset))
        
        cat("Using single best model with", top_models[[1]]$n_predictors, "predictor(s):",
            paste(top_models[[1]]$predictors, collapse = ", "), "\n")
        cat("AIC =", round(top_models[[1]]$aic, 2), "\n")
        
        # Custom predict for single scaled model
        predict.scaled_glm <- function(object, newdata = NULL, type = "response", ...) {
          if(!is.null(newdata)) {
            newdata_scaled <- newdata
            for(pred in names(object$scaling_params)) {
              if(pred %in% names(newdata)) {
                newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                  object$scaling_params[[pred]]$scale
              }
            }
            predict.glm(object, newdata = newdata_scaled, type = type, ...)
          } else {
            predict.glm(object, type = type, ...)
          }
        }
      } else {
        # Ensemble of top models
        cat("Ensembling", length(top_models), "models\n")
        
        delta_aic <- sapply(top_models, function(x) x$aic - best_aic)
        aic_weights <- exp(-0.5 * delta_aic)
        aic_weights <- aic_weights / sum(aic_weights)
        
        for(i in seq_along(top_models)) {
          vif_str <- if(is.na(top_models[[i]]$max_vif)) "N/A" else round(top_models[[i]]$max_vif, 2)
          cat("  Model", i, "(", top_models[[i]]$n_predictors, "pred):",
              paste(top_models[[i]]$predictors, collapse = ", "),
              "(AIC =", round(top_models[[i]]$aic, 2), ", Weight =", round(aic_weights[i], 3), ")\n")
        }
        
        # Create ensemble model object
        env_model_offset <- structure(
          list(
            models = lapply(top_models, function(x) x$model),
            weights = aic_weights,
            scaling_params = scaling_params,
            call = call("ensemble"),
            coefficients = top_models[[1]]$model$coefficients
          ),
          class = c("ensemble_glm", "list")
        )
        
        # Custom predict method for ensemble
        predict.ensemble_glm <- function(object, newdata = NULL, type = "response", ...) {
          if(!is.null(newdata)) {
            newdata_scaled <- newdata
            for(pred in names(object$scaling_params)) {
              if(pred %in% names(newdata)) {
                newdata_scaled[[pred]] <- (newdata[[pred]] - object$scaling_params[[pred]]$center) / 
                  object$scaling_params[[pred]]$scale
              }
            }
          } else {
            newdata_scaled <- NULL
          }
          
          predictions <- sapply(object$models, function(model) {
            predict(model, newdata = newdata_scaled, type = type, ...)
          })
          
          if(is.matrix(predictions)) {
            predictions %*% object$weights
          } else {
            sum(predictions * object$weights)
          }
        }
      }
      
      # Generate predictions and residuals
      phen_estimate$predicted_offset <- predict(env_model_offset, newdata = phen_estimate, type = "response")
      phen_estimate$residuals_offset <- phen_estimate$offset - phen_estimate$predicted_offset
      
      # Calculate metrics
      res_env_offset <- phen_estimate$residuals_offset
      rmse_val_offset <- sqrt(mean(res_env_offset^2))
      
      # Calculate R²
      ss_res_offset <- sum(res_env_offset^2)
      ss_tot_offset <- sum((phen_estimate$offset - mean(phen_estimate$offset))^2)
      r2_val_offset <- 1 - (ss_res_offset / ss_tot_offset)
      
      cat("R² =", round(r2_val_offset, 3), "| RMSE =", round(rmse_val_offset, 3), "\n")
    }
  }
  
  
  # Compute Duration
  phen_estimate <- phen_estimate %>% 
    rename(predicted_onset = predicted) %>% 
    mutate(predicted_duration = predicted_offset - predicted_onset)
  
  # ensure that hopkins_delay_midpoint is available in the original raw data as well
  phen_data <- phen_data %>% 
    mutate(
      # Extract numeric values from bin intervals
      bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
      bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
      hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
    )
  
  
  
  ################################################################################################
  ## Visualize Outputs 
  # Hopkins Only
  hop_plot1 <- ggplot(phen_estimate, aes(x = hopkins_delay_midpoint, y = onset)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE) +
    annotate("text", x = min(phen_estimate$hopkins_delay_midpoint), y = max(phen_estimate$onset),
             label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_hop[[2]], rmse_val_hop),
             hjust = 0, vjust = 1) +
    labs(x = "Hopkins Expected Delay Bin (Middle Value for Bin)",
         y = "Observed Onset (day of year)",
         title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
    theme_minimal()
  hop_plot2 <- ggplot(phen_estimate, aes(x = corrected_onset_prediction, y = onset)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    #geom_ribbon(aes(ymin = corrected_onset_lower, ymax = corrected_onset_upper), 
    #alpha = 0.2, fill = "blue") +
    geom_point(size = 3) +
    annotate("text", x = min(phen_estimate$onset), y = max(phen_estimate$corrected_onset_prediction),
             label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_hop[[2]], rmse_val_hop),
             hjust = 0, vjust = 1) +
    labs(x = "Corrected Predicted Onset (anchored + actual slope)",
         y = "Empirical Estimate Onset (day of year)",
         title = paste("Testing Hopkins Bioclimatic Law on", target_taxa[1])) +
    theme_minimal()
  # Env Space
  if(nrow(phen_estimate) > 3){
    env_plot <- ggplot(phen_estimate, aes(x = onset, y = predicted_onset)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      annotate("text", x = min(phen_estimate$onset), y = max(phen_estimate$predicted_onset),
               label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_val, rmse_val),
               hjust = 0, vjust = 1) +
      labs(x = "Observed Onset", y = "Predicted Onset", 
           title = paste("Using Climatic Vars to Predict", target_taxa[1])) +
      theme_minimal()
    env_plot_offset <- ggplot(phen_estimate, aes(x = offset, y = predicted_offset)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      annotate("text", x = min(phen_estimate$offset), y = max(phen_estimate$predicted_offset),
               label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_val_offset, rmse_val_offset),
               hjust = 0, vjust = 1) +
      labs(x = "Observed Offset", y = "Predicted Offset", 
           title = paste("Using Climatic Vars to Predict", target_taxa[1])) +
      theme_minimal()
    
  }
  
  if(nrow(phen_estimate) > 3){
    # return everything as list
    return(list(
      hopkins_model = hopkins_model, 
      env_model_onset = env_model, 
      env_model_offset = env_model_offset,
      data = phen_data,
      phen_est_data = phen_estimate, 
      data_w_offset_vars = phen_data_with_windows,
      plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = c(env_plot, env_plot_offset)), 
      metrics = list(
        hopkins_r2 = r2_hop[[2]], 
        hopkins_rsme = rmse_val_hop, 
        env_r2_onset = r2_val, 
        env_rsme_onset = rmse_val, 
        env_r2_offset = r2_val_offset, 
        env_rsme_offset = rmse_val_offset
      )
    ))
  } else{
    # return everything as list
    return(list(
      hopkins_model = hopkins_model, 
      env_model_onset = NULL, 
      env_model_offset =NULL,
      data = phen_data,
      phen_est_data = phen_estimate, 
      data_w_offset_vars = phen_data_with_windows,
      plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = NULL), 
      metrics = list(
        hopkins_r2 = r2_hop[[2]], 
        hopkins_rsme = rmse_val_hop, 
        env_r2 = NULL, 
        env_rsme = NULL, 
        env_r2_offset = NULL, 
        env_rsme_offset = NULL
      )
    ))
    
  }
  
  
  
}

### Create Predictions for full distributions V2 (adds in ensemble methods) ##################################
predict_pheno_for_dist <- function(model_out, 
                                   path_to_target_taxon_PA_raster = NULL,
                                   hopkins = TRUE,
                                   env = TRUE){
  
  ## Pull in dist model, limit to presence range only, apply modified hopkins linear model to P pixels 
  taxon_rast <- rast(paste0(path_to_target_taxon_PA_raster))
  taxon_rast[taxon_rast  == 0] <- NA # remove non-presence pts
  
  # reproj to the native wgs84 coordinate system used to build hopkins expectations
  taxon_rast_reproj <- project(taxon_rast, "EPSG:4326", method = "near")
  
  # make into xy df 
  taxon_df <- as.data.frame(taxon_rast_reproj, xy = TRUE, na.rm = TRUE)
  taxon_df <- taxon_df %>% 
    rename(longitude = x, latitude = y) %>% 
    mutate(rowid = 1:n()) %>% 
    as.data.table()
  
  ## convert this into a spatVect to grab elevation & clim vars, reproj to match elev and clim tifs
  # load North America Elev Raster & Chelsa annual clim vars
  na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
  avg_temp <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio1_1981-2010_V.2.1.tif")
  avg_diurnal_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio02_1981-2010_V.2.1.tif")
  isothermality <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio03_1981-2010_V.2.1.tif")
  temp_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio4_1981-2010_V.2.1.tif")
  tmax <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio5_1981-2010_V.2.1.tif")
  tmin <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio6_1981-2010_V.2.1.tif")
  annual_temp_range <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio07_1981-2010_V.2.1.tif")
  avg_monthly_temp_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio08_1981-2010_V.2.1.tif")
  avg_monthly_temp_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio09_1981-2010_V.2.1.tif")
  avg_monthly_temp_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio10_1981-2010_V.2.1.tif")
  avg_monthly_temp_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio11_1981-2010_V.2.1.tif")
  annual_precip <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio12_1981-2010_V.2.1.tif")
  precip_wettest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio13_1981-2010_V.2.1.tif")
  precip_driest_month <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio14_1981-2010_V.2.1.tif")
  precip_seas <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio15_1981-2010_V.2.1.tif")
  avg_monthly_precip_wettest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio16_1981-2010_V.2.1.tif")
  avg_monthly_precip_driest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio17_1981-2010_V.2.1.tif")
  avg_monthly_precip_warmest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio18_1981-2010_V.2.1.tif")
  avg_monthly_precip_coldest_q <- rast("/blue/guralnick/share/Chelsa-climate-data/climatologies/CHELSA_bio19_1981-2010_V.2.1.tif")
  
  # Monthly vals for termination cues
  pet_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                               pattern = "CHELSA_pet_.*_1981-2010_V.2.1.tif$",
                               full.names = TRUE))
  tasmax_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_tasmax_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  tasmin_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_tasmin_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  tasmin_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_tasmin_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  wind_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                pattern = "CHELSA_sfcWind_.*_1981-2010_V.2.1.tif$",
                                full.names = TRUE))
  precip_stack <- rast(list.files("/orange/guralnick/Chelsa_climatologies/", 
                                  pattern = "CHELSA_pr_.*_1981-2010_V.2.1.tif$",
                                  full.names = TRUE))
  
  
  taxon_vect <- vect(taxon_df, 
                     geom = c("longitude", "latitude"),
                     crs = "EPSG:4326")
  taxon_vect_reproj1 <- project(taxon_vect, crs(na_elev))
  taxon_vect_reproj2 <- project(taxon_vect, crs(avg_temp))
  
  # Stack all climate layers
  climate_stack <- c(avg_temp, avg_diurnal_temp_range, isothermality,
                     temp_seas, tmax, tmin, annual_temp_range,
                     avg_monthly_temp_wettest_q, avg_monthly_temp_driest_q,
                     avg_monthly_temp_warmest_q, avg_monthly_temp_coldest_q,
                     annual_precip, precip_wettest_month, precip_driest_month,
                     precip_seas, avg_monthly_precip_wettest_q, avg_monthly_precip_driest_q,
                     avg_monthly_precip_warmest_q, avg_monthly_precip_coldest_q)
  
  # make names clearer
  names(climate_stack) <- c("avg_temp", "avg_diurnal_temp_range", "isothermality",
                            "temp_seas", "tmax", "tmin", "annual_temp_range",
                            "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                            "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                            "annual_precip", "precip_wettest_month", "precip_driest_month",
                            "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                            "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q")
  
  names(pet_stack) <- c("pet_01", "pet_02", "pet_03", "pet_04", "pet_05", "pet_06", 
                        "pet_07", "pet_08", "pet_09", "pet_10", "pet_11", "pet_12")
  names(tasmax_stack) <- c("tasmax_01", "tasmax_02", "tasmax_03", "tasmax_04", "tasmax_05", "tasmax_06", 
                           "tasmax_07", "tasmax_08", "tasmax_09", "tasmax_10", "tasmax_11", "tasmax_12")
  names(tasmin_stack) <- c("tasmin_01", "tasmin_02", "tasmin_03", "tasmin_04", "tasmin_05", "tasmin_06", 
                           "tasmin_07", "tasmin_08", "tasmin_09", "tasmin_10", "tasmin_11", "tasmin_12")
  names(wind_stack) <- c("sfcWind_01", "sfcWind_02", "sfcWind_03", "sfcWind_04", "sfcWind_05", "sfcWind_06", 
                         "sfcWind_07", "sfcWind_08", "sfcWind_09", "sfcWind_10", "sfcWind_11", "sfcWind_12")
  names(precip_stack) <- c("pr_01", "pr_02", "pr_03", "pr_04", "pr_05", "pr_06", 
                           "pr_07", "pr_08", "pr_09", "pr_10", "pr_11", "pr_12")
  # Extract all raster values
  elevation_vals <- terra::extract(na_elev, taxon_vect_reproj1)
  climate_vals <- terra::extract(climate_stack, taxon_vect_reproj2)
  pet_vals <- terra::extract(pet_stack, taxon_vect_reproj2)
  tasmax_vals <- terra::extract(tasmax_stack, taxon_vect_reproj2)
  tasmin_vals <- terra::extract(tasmin_stack, taxon_vect_reproj2)
  wind_vals <- terra::extract(wind_stack, taxon_vect_reproj2)
  precip_vals <- terra::extract(precip_stack, taxon_vect_reproj2)
  
  # Add all values to these data
  taxon_df[, elevation_m := elevation_vals[,2]]
  taxon_df[, names(climate_vals)[-1] := climate_vals[,-1]]
  taxon_df[, names(pet_vals)[-1] := pet_vals[, -1]]
  taxon_df[, names(tasmax_vals)[-1] :=tasmax_vals[, -1]]
  taxon_df[, names(tasmin_vals)[-1] := tasmin_vals[, -1]]
  taxon_df[, names(wind_vals)[-1] := wind_vals[, -1]]
  taxon_df[, names(precip_vals)[-1] := precip_vals[, -1]]
  
  # Make hopkins bins a available predictor 
  ## Using this modified dataframe, recreate hopkins based on the anchor from the empirical data
  df <- model_out$data
  # bring in min vals 
  min_empirical_lat <- min(df$latitude)
  min_empirical_elev <- min(df$elevation_m)
  # these will be zero, we'll just adjust so we're going into negative space
  # create expected # of delay days based on hopkins
  taxon_df[, hopkins_delay_days :=
             ((latitude - min_empirical_lat)* 4) + # latitudnal effect
             ((elevation_m - min_empirical_elev))/120 * 4] # elevation effect
  
  # Create bins for phenology data based on Hopkins delay by intervals of 4
  taxon_df[, hopkins_bin := cut(hopkins_delay_days, 
                                breaks = seq(0, max(hopkins_delay_days, na.rm = TRUE) + 4, by = 4), include.lowest = TRUE)]
  # additionally, we want to know what pixels are outside of the range of empirical data training
  min_empirical_hopkins <- min(df$hopkins_delay_days)
  max_empirical_hopkins <- max(df$hopkins_delay_days)
  # binary
  taxon_df[, outside_empirical := hopkins_delay_days < min_empirical_hopkins |
             hopkins_delay_days > max_empirical_hopkins]
  # numeric dist
  taxon_df[, dist_from_empirical := fifelse(
    hopkins_delay_days < min_empirical_hopkins, 
    min_empirical_hopkins - hopkins_delay_days, # below range
    fifelse(
      hopkins_delay_days > max_empirical_hopkins, 
      hopkins_delay_days - max_empirical_hopkins, # above range
      0 #if within range
    )
  )]
  
  # extract midpoint of each hopkins bin to use as predictor
  taxon_df <- taxon_df %>%
    mutate(
      # Extract numeric values from bin intervals
      bin_lower = as.numeric(sub("\\((.+),.+\\]", "\\1", hopkins_bin)),
      bin_upper = as.numeric(sub("\\(.+,(.+)\\]", "\\1", hopkins_bin)),
      hopkins_delay_midpoint = (bin_lower + bin_upper) / 2
    )
  # Some pixels lack elevation, which is why we get this warning. Probably try and see if we can fix this at a later time. 
  taxon_df <- taxon_df %>% filter(!is.na(hopkins_bin))
  taxon_df <- taxon_df %>% filter(!is.na(hopkins_delay_midpoint)) # not entirely sure whats happening here. 
  
  ## Create Predictions 
  # apply hopkins linear model 
  hopkins_linear <- model_out$hopkins_model
  taxon_df_predict <- predict(hopkins_linear, newdata = taxon_df)
  taxon_df[, hopkins_predict_onset := taxon_df_predict]
  # apply env glm 
  ## If its ensembled, we need to account for this and the weight of each model
  predict_ensemble <- function(ensemble_model, newdata) {
    
    # scale the newdata using the scaling parameters
    scaled_data <- newdata
    for (var in names(ensemble_model$scaling_params)) {
      if (var %in% names(newdata)) {
        center <- ensemble_model$scaling_params[[var]]$center
        scale <- ensemble_model$scaling_params[[var]]$scale
        scaled_data[[var]] <- (newdata[[var]] - center) / scale
      }
    }
    
    # get predictions from each model
    n_models <- length(ensemble_model$models)
    predictions <- matrix(NA, nrow = nrow(newdata), ncol = n_models)
    
    for (i in 1:n_models) {
      predictions[, i] <- predict(ensemble_model$models[[i]], 
                                  newdata = scaled_data, 
                                  type = "response")
    }
    
    # weight the predictions
    weighted_pred <- predictions %*% ensemble_model$weights
    
    return(as.vector(weighted_pred))
  }
  ## Bring in the onset & offset models 
  env_model_onset <- model_out$env_model_onset
  env_model_offset <- model_out$env_model_offset
  
  ## Predict Onset, if ensemble used then use predict_ensemble(), if not just use the single best model...
  if(length(model_out$env_model_onset$models) > 0){
    taxon_df_predict_onset <- predict_ensemble(env_model_onset, taxon_df)
    taxon_df[, env_predict_onset := taxon_df_predict_onset]
    cat("Ensemble Model Methods used to build onset prediction, applying this to full taxon range...")
  } else {
    taxon_df_predict_onset <- predict(env_model_onset, taxon_df)
    taxon_df[, env_predict_onset := taxon_df_predict_onset]
    cat("Single Best Model methods used to build onset prediction, applying this to full taxon range..,")
  }
  
  ## Use Onset 'Anchor' to find window for potential offset termination cues 
  ### Find the month onset was predicted to occur
  taxon_df <- taxon_df %>% 
    mutate(onset_month = lubridate::month(as.Date(env_predict_onset, origin = "1998-01-01")))
  ### Create month winodw fxn 
  get_month_window <- function(center_month) {
    months <- c((center_month - 2) %% 12 + 1, 
                (center_month - 1) %% 12 + 1,
                center_month,
                center_month %% 12 + 1,
                (center_month + 1) %% 12 + 1)
    unique(months[1:5])  # Get -2, 0, +2
  }
  ### Find averages for each window 
  taxon_df <- taxon_df %>% 
    rowwise() %>% 
    mutate(
      avg_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_pet_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pet_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_tasmax_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmax_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_tasmin_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("tasmin_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_sfcWind_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("sfcWind_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      },
      avg_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        mean(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      sd_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        sd(as.numeric(cur_data()[cols]), na.rm = TRUE)
      },
      slope_pr_window = {
        m <- get_month_window(onset_month)
        cols <- paste0("pr_", sprintf("%02d", m))
        vals <- as.numeric(cur_data()[cols])
        if (all(is.na(vals))) NA_real_ else coef(lm(vals ~ seq_along(vals)))[2]
      } , 
      total_photoperiod_window = {
        m <- get_month_window(onset_month)
        phi <- latitude * pi / 180
        
        total <- 0
        for (month in m) {
          days_in_month <- as.numeric(format(seq(as.Date(paste0("2023-", month, "-01")), 
                                                 by = "month", length.out = 2)[2] - 1, "%d"))
          for (day in 1:days_in_month) {
            doy <- as.numeric(format(as.Date(paste0("2023-", month, "-", day)), "%j"))
            delta <- 0.409 * sin(2 * pi * (doy - 81) / 365)
            w <- acos(-tan(phi) * tan(delta))
            daylength <- 24 * w / pi
            total <- total + daylength
          }
        }
        total
      }
    ) %>%
    ungroup() 
  taxon_df <- taxon_df %>% as.data.table()
  
  ## Now Predict Offset, if esembled use predict_ensemble, if not just use single best model
  if(length(model_out$env_model_offset$models) > 0){
  taxon_df_predict_offset <- predict_ensemble(env_model_offset, taxon_df)
  taxon_df[, env_predict_offset := taxon_df_predict_offset]
  cat("Ensemble Model Methods used to build offset prediction, applying this to full taxon range...")
  } else {
    taxon_df_predict_offset <- predict(env_model_offset, taxon_df)
    taxon_df[, env_predict_offset := taxon_df_predict_offset]
    cat("Single Best Model methods used to build offset prediction, applying this to full taxon range..,")
  }
  
  ## Calculate Predicted Duration
  taxon_df <- taxon_df %>% 
    mutate(env_predict_duration = env_predict_offset - env_predict_onset)
  
  # Further use a CLAMP to restrict ourselves to Interpolation 
  ## extract all unique predictor variables from the ensemble
  if(length(model_out$env_model_onset$models) > 0){
  predictor_vars_onset <- unique(unlist(lapply(env_model_onset$models, function(model) {
    attr(terms(model), "term.labels")
  })))
  } else{
    predictor_vars_onset <- attr(terms(env_model_onset), "term.labels")
  }
  if(length(model_out$env_model_offset$models) > 0){
    predictor_vars_offset <- unique(unlist(lapply(env_model_offset$models, function(model) {
      attr(terms(model), "term.labels")
    })))
  } else{
    predictor_vars_offset <- attr(terms(env_model_offset), "term.labels")
  }
  df <- model_out$data
  df_offset <- model_out$data_w_offset_vars
  training_vars_onset <- df %>% 
    ungroup() %>%
    select(all_of(predictor_vars_onset)) %>% 
    distinct() %>% # remove redundancy
    as.data.frame()
  training_vars_offset <- df_offset %>% 
    ungroup() %>%
    select(all_of(predictor_vars_offset)) %>% 
    distinct() %>% # remove redundancy
    as.data.frame()
  
  # calc ranges (made flexible for different predictor_var combos)
  train_ranges_onset <- data.frame(
    var = predictor_vars_onset, 
    min = sapply(predictor_vars_onset, function(v) min(training_vars_onset[[v]])), 
    max = sapply(predictor_vars_onset, function(v) max(training_vars_onset[[v]]))
  )
  train_ranges_offset <- data.frame(
    var = predictor_vars_offset, 
    min = sapply(predictor_vars_offset, function(v) min(training_vars_offset[[v]])), 
    max = sapply(predictor_vars_offset, function(v) max(training_vars_offset[[v]]))
  )
  # clamp!!
  for(i in 1:nrow(train_ranges_onset)){ # for each var
    var_name <- train_ranges_onset$var[i]
    clamped_name <- paste0(var_name, "_onset_clamped")
    taxon_df[, (clamped_name) := pmax(train_ranges_onset$min[i], 
                                      pmin(train_ranges_onset$max[i], get(var_name)))]
  }
  for(i in 1:nrow(train_ranges_offset)){ # for each var
    var_name <- train_ranges_offset$var[i]
    clamped_name <- paste0(var_name, "_offset_clamped")
    taxon_df[, (clamped_name) := pmax(train_ranges_offset$min[i], 
                                      pmin(train_ranges_offset$max[i], get(var_name)))]
  }
  # create clamped prediction df
  clamped_vars_onset <- paste0(predictor_vars_onset, "_onset_clamped")
  clamped_vars_offset <- paste0(predictor_vars_offset, "_offset_clamped")
  
  clamped_data_onset <- taxon_df[, ..clamped_vars_onset]
  setnames(clamped_data_onset, clamped_vars_onset, predictor_vars_onset)
  clamped_data_offset <- taxon_df[, ..clamped_vars_offset]
  setnames(clamped_data_offset, clamped_vars_offset, predictor_vars_offset)
  
  
  # make clamped predictions
  if(length(model_out$env_model_onset$models) > 0){
  taxon_df[, env_predict_onset_clamped := predict_ensemble(env_model_onset, newdata = clamped_data_onset)]
    cat("Ensemble Model Methods used to build onset prediction, applying this with clamp to full taxon range...")
  } else{
    taxon_df[, env_predict_onset_clamped := predict(env_model_onset, newdata = clamped_data_onset)]
    cat("Single Best Model Methods used to build onset prediction, applying this with clamp to full taxon range...")
  }
  
  if(length(model_out$env_model_offset$models) > 0){
    taxon_df[, env_predict_offset_clamped := predict_ensemble(env_model_offset, newdata = clamped_data_offset)]
    cat("Ensemble Model Methods used to build offset prediction, applying this with clamp to full taxon range...")
  } else{
    taxon_df[, env_predict_offset_clamped := predict(env_model_offset, newdata = clamped_data_offset)]
    cat("Single Best Model Methods used to build offset prediction, applying this with clamp to full taxon range...")
  }
  
  # Calc Clamped Duration 
  taxon_df <- taxon_df %>% 
    mutate(env_predict_duration_clamped = env_predict_offset_clamped - env_predict_onset_clamped)
  
  # Create plots 
  ## get basemap (will need to edit to be more flexible l8ter)
  basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
  basemap <- basemap %>% filter(name %in% c("California", "Nevada", "Arizona", "Sonora", "Baja California", "Baja California Sur"))
  basemap_reproj <- basemap %>% st_transform(st_crs(taxon_rast_reproj))
  hop_pred_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = hopkins_predict_onset)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Hopkins Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Adjusted Hopkins Onset Predictions for", unique(df$species)))
  
  hop_training_dist_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = dist_from_empirical)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Hopkins Empirical Training Distance Eval", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Training Data Distance on Hopkins Onset Predictions for", unique(df$species)))
  
  env_pred_onset_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_onset)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Climatic Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Climate GLM Onset Prediction Plot for", unique(df$species)))
  
  env_pred_offset_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_offset)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Climatic Offset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Climate GLM Offset Prediction Plot for", unique(df$species)))
  
  env_pred_duration_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_duration)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Climatic Duration Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Climate GLM Duration Prediction Plot for", unique(df$species)))
  
  clamped_env_pred_onset_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_onset_clamped)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Clamped Climatic Onset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Clamped Climate GLM Onset Prediction Plot for", unique(df$species)))
  
  clamped_env_pred_offset_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_offset_clamped)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Clamped Climatic Offset Prediction", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Clamped Climate GLM Offset Prediction Plot for", unique(df$species)))
  
  min_val <- min(taxon_df$env_predict_duration_clamped, na.rm = TRUE)
  max_val <- max(taxon_df$env_predict_duration_clamped, na.rm = TRUE)
  zero_position <- (0 - min_val) / (max_val - min_val)
  
  clamped_env_pred_duration_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_duration_clamped)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_gradientn(
      colors = c("grey80", "grey80", viridis::plasma(100)),
      values = c(0, zero_position - 0.0001, seq(zero_position, 1, length.out = 100)),
      breaks = c(-100, 0, 50, 100, 150, 200),
      labels = c("Unsuitable Duration \n Modeling", "0", "50", "100", "150", "200"),
      name = "Clamped Climatic\nDuration Prediction"
    ) +
    theme_bw() +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = paste("Clamped Climate GLM Duration Prediction Plot for", unique(df$species)))
  
  ## Calculate MESS 
  if(env == TRUE){
    ## ONSET MESS ##############################################################
    # crop rasters down to actual P area of taxon
    taxon_climate_proj <- project(taxon_rast_reproj, crs(temp_seas), method = "near") # use temp seas as a placeholder for all clims
    climate_onset_crop <- crop(climate_stack[[predictor_vars_onset]], taxon_climate_proj)
    taxon_climate_onset_resamp <- resample(taxon_climate_proj, climate_onset_crop[[1]], method = "near")
    climate_onset_masked <- mask(climate_onset_crop, taxon_climate_onset_resamp)
    # convert to old raster format for dismo::mess()
    climate_onset_stack_old <- raster::stack(climate_onset_masked)
    names(climate_onset_stack_old) <- predictor_vars_onset
    # calc MESS
    onset_mess_raster_old <- dismo::mess(climate_onset_stack_old, training_vars_onset, full = FALSE)
    # convert to terra spatRaster
    onset_mess_raster <- rast(onset_mess_raster_old)
    # extract and append to df 
    onset_mess_vals <- terra::extract(onset_mess_raster, 
                                taxon_df[, .(longitude, latitude)])
    taxon_df[, onset_mess := onset_mess_vals[[2]]]
    ############################################################################
    
    ## OFFSET MESS #############################################################
    taxon_climate_proj <- project(taxon_rast_reproj, crs(temp_seas), method = "near")
    climate_stack_offset_intersection <- intersect(names(climate_stack), predictor_vars_offset)
    offset_monthly_vars <- setdiff(predictor_vars_offset, names(climate_stack))
    
    # Get climate variables as rasters
    climate_offset_crop <- crop(climate_stack[[climate_stack_offset_intersection]], taxon_climate_proj)
    taxon_climate_offset_resamp <- resample(taxon_climate_proj, climate_offset_crop[[1]], method = "near")
    climate_offset_masked <- mask(climate_offset_crop, taxon_climate_offset_resamp)
    
    # Create a data frame with coords and monthly vars, then convert to SpatVector
    monthly_df <- taxon_df[, c("longitude", "latitude", offset_monthly_vars), with = FALSE]
    coords_vect <- vect(monthly_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    
    # Reproject to match climate rasters
    coords_reproj <- project(coords_vect, crs(climate_offset_masked))
    
    # Rasterize monthly variables directly to match climate grid
    monthly_raster_list <- lapply(offset_monthly_vars, function(var){
      rast_var <- rasterize(coords_reproj, climate_offset_masked[[1]], field = var, fun = mean)
      names(rast_var) <- var
      return(rast_var)
    })
    
    monthly_window_rast <- rast(monthly_raster_list)
    
    # Combine all variables in correct order
    climate_offset_full <- c(climate_offset_masked, monthly_window_rast)
    climate_offset_full <- climate_offset_full[[predictor_vars_offset]]  # Reorder to match
    
    # Convert to old raster format
    climate_offset_stack_old <- raster::stack(climate_offset_full)
    
    # Get training data
    training_vars_offset <- model_out$data_w_offset_vars %>%
      select(all_of(predictor_vars_offset)) %>%
      as.data.frame()
    
    # calc MESS
    offset_mess_raster_old <- dismo::mess(climate_offset_stack_old, training_vars_offset, full = FALSE)
    offset_mess_raster <- rast(offset_mess_raster_old)
    offset_mess_vals <- terra::extract(offset_mess_raster, taxon_df[, .(longitude, latitude)])
    taxon_df[, offset_mess := offset_mess_vals[[2]]]
    
    ############################################################################
    # plot
    mess_plot_onset <- ggplot() +
      geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = onset_mess)) +
      geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
      scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
      theme_minimal() +
      theme_bw() +
      labs(
        x = "Longitude", 
        y = "Latitude", 
        title = paste("MESS Map for Climate GLM Onset Prediction Plot for", unique(df$species)))
    
    mess_plot_offset <- ggplot() +
      geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = offset_mess)) +
      geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
      scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
      theme_minimal() +
      theme_bw() +
      labs(
        x = "Longitude", 
        y = "Latitude", 
        title = paste("MESS Map for Climate GLM Offset Prediction Plot for", unique(df$species)))
  } else{
    mess_plot <- NULL
  }
  return(list(
    taxon_df = taxon_df,
    plots = list(
      prediction = list(hop_pred_plot, env_pred_onset_plot, env_pred_offset_plot, clamped_env_pred_onset_plot, clamped_env_pred_offset_plot), 
      quality = list(hop_training_dist_plot, mess_plot_onset, mess_plot_offset)
    )
  ))
  
  
}

