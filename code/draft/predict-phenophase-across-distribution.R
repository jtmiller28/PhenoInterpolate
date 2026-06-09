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

  ## If its not ensembled, we need to predict after normalizing scales
  predict_single <- function(scaled_model, newdata, type = "response") {
    # scale the newdata using the scaling parameters
    scaled_data <- newdata
    for (var in names(scaled_model$scaling_params)) {
      if (var %in% names(newdata)) {
        center <- scaled_model$scaling_params[[var]]$center
        scale <- scaled_model$scaling_params[[var]]$scale
        scaled_data[[var]] <- (newdata[[var]] - center) / scale
      }
    }
    # get predictions from the model
    predictions <- predict.glm(scaled_model, newdata = scaled_data, type = type)
    return(as.vector(predictions))
  }
  ## Bring in the onset & offset models 
  env_model_onset <- model_out$env_model_onset
  env_model_duration <- model_out$env_model_duration
  
  ## Predict Onset, if ensemble used then use predict_ensemble(), if not just use the single best model...
  if(length(model_out$env_model_onset$models) > 0){
    taxon_df_predict_onset <- predict_ensemble(env_model_onset, taxon_df)
    taxon_df[, env_predict_onset := taxon_df_predict_onset]
    cat("Ensemble Model Methods used to build onset prediction, applying this to full taxon range...")
  } else {
    taxon_df_predict_onset <- predict_single(env_model_onset, taxon_df)
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
    unique(months[2:4])  # Get -1, 0, +1
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
  taxon_df <- taxon_df %>% as.data.table() #%>% mutate(onset = env_predict_onset) # we need to add this... # Try removing.
  
  ## Now Predict Duration, if ensembled use predict_ensemble, if not just use single best model
  if(length(model_out$env_model_duration$models) > 0){
    taxon_df_predict_duration <- predict_ensemble(env_model_duration, taxon_df)
    taxon_df[, env_predict_duration := taxon_df_predict_duration]
    cat("Ensemble Model Methods used to build duration prediction, applying this to full taxon range...")
  } else {
    taxon_df_predict_duration <- predict_single(env_model_duration, taxon_df)
    taxon_df[, env_predict_duration := taxon_df_predict_duration]
    cat("Single Best Model methods used to build duration prediction, applying this to full taxon range..,")
  }
  
  ## Calculate Predicted Offset 
  taxon_df <- taxon_df %>% 
    mutate(env_predict_offset = env_predict_onset + env_predict_duration)
  
  # Further use a CLAMP to restrict ourselves to Interpolation 
  ## extract all unique predictor variables from the ensemble
  if(length(model_out$env_model_onset$models) > 0){
    predictor_vars_onset <- unique(unlist(lapply(model_out$env_model_onset$models, function(model) {
      attr(terms(model), "term.labels")
    })))
  } else{
    predictor_vars_onset <- attr(terms(env_model_onset), "term.labels")
  }
  if(length(model_out$env_model_duration$models) > 0){
    predictor_vars_duration <- unique(unlist(lapply(model_out$env_model_duration$models, function(model) {
      attr(terms(model), "term.labels")
    })))
  } else{
    predictor_vars_duration <- attr(terms(env_model_duration), "term.labels")
  }
  df <- model_out$data
  df_duration <- model_out$data_w_duration_vars
  training_vars_onset <- df %>%
    ungroup() %>%
    select(all_of(predictor_vars_onset)) %>%
    distinct() %>% # remove redundancy
    as.data.frame()
  training_vars_duration <- df_duration %>%
    ungroup() %>%
    select(all_of(predictor_vars_duration)) %>%
    distinct() %>% # remove redundancy
    as.data.frame()
  # 
  # # calc ranges (made flexible for different predictor_var combos)
  # train_ranges_onset <- data.frame(
  #   var = predictor_vars_onset, 
  #   min = sapply(predictor_vars_onset, function(v) min(training_vars_onset[[v]])), 
  #   max = sapply(predictor_vars_onset, function(v) max(training_vars_onset[[v]]))
  # )
  # train_ranges_duration <- data.frame(
  #   var = predictor_vars_duration, 
  #   min = sapply(predictor_vars_duration, function(v) min(training_vars_duration[[v]])), 
  #   max = sapply(predictor_vars_duration, function(v) max(training_vars_duration[[v]]))
  # )
  # # clamp!!
  # for(i in 1:nrow(train_ranges_onset)){ # for each var
  #   var_name <- train_ranges_onset$var[i]
  #   clamped_name <- paste0(var_name, "_onset_clamped")
  #   taxon_df[, (clamped_name) := pmax(train_ranges_onset$min[i], 
  #                                     pmin(train_ranges_onset$max[i], get(var_name)))]
  # }
  # for(i in 1:nrow(train_ranges_duration)){ # for each var
  #   var_name <- train_ranges_duration$var[i]
  #   clamped_name <- paste0(var_name, "_duration_clamped")
  #   taxon_df[, (clamped_name) := pmax(train_ranges_duration$min[i], 
  #                                     pmin(train_ranges_duration$max[i], get(var_name)))]
  # }
  # # create clamped prediction df
  # clamped_vars_onset <- paste0(predictor_vars_onset, "_onset_clamped")
  # clamped_vars_duration <- paste0(predictor_vars_duration, "_duration_clamped")
  # 
  # clamped_data_onset <- taxon_df[, ..clamped_vars_onset]
  # setnames(clamped_data_onset, clamped_vars_onset, predictor_vars_onset)
  # clamped_data_duration <- taxon_df[, ..clamped_vars_duration]
  # setnames(clamped_data_duration, clamped_vars_duration, predictor_vars_duration)
  # 
  # 
  # # make clamped predictions
  # if(length(model_out$env_model_onset$models) > 0){
  #   taxon_df[, env_predict_onset_clamped := predict_ensemble(env_model_onset, newdata = clamped_data_onset)]
  #   cat("Ensemble Model Methods used to build onset prediction, applying this with clamp to full taxon range...")
  # } else{
  #   taxon_df[, env_predict_onset_clamped := predict(env_model_onset, newdata = clamped_data_onset)]
  #   cat("Single Best Model Methods used to build onset prediction, applying this with clamp to full taxon range...")
  # }
  # 
  # if(length(model_out$env_model_duration$models) > 0){
  #   taxon_df[, env_predict_duration_clamped := predict_ensemble(env_model_duration, newdata = clamped_data_duration)]
  #   cat("Ensemble Model Methods used to build duration prediction, applying this with clamp to full taxon range...")
  # } else{
  #   taxon_df[, env_predict_duration_clamped := predict(env_model_duration, newdata = clamped_data_duration)]
  #   cat("Single Best Model Methods used to build duration prediction, applying this with clamp to full taxon range...")
  # }
  # 
  # # Calc Clamped Offset 
  # taxon_df <- taxon_df %>% 
  #   mutate(env_predict_offset_clamped = env_predict_onset_clamped + env_predict_duration_clamped)
  # 
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
  
  # clamped_env_pred_onset_plot <- ggplot() +
  #   geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_onset_clamped)) +
  #   geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  #   scale_fill_viridis_c(name = "Clamped Climatic Onset Prediction", option = "plasma") +
  #   theme_bw() +
  #   labs(
  #     x = "Longitude", 
  #     y = "Latitude", 
  #     title = paste("Clamped Climate GLM Onset Prediction Plot for", unique(df$species)))
  # 
  # clamped_env_pred_offset_plot <- ggplot() +
  #   geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_offset_clamped)) +
  #   geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  #   scale_fill_viridis_c(name = "Clamped Climatic Offset Prediction", option = "plasma") +
  #   theme_bw() +
  #   labs(
  #     x = "Longitude", 
  #     y = "Latitude", 
  #     title = paste("Clamped Climate GLM Offset Prediction Plot for", unique(df$species)))
  # 
  # min_val <- min(taxon_df$env_predict_duration_clamped, na.rm = TRUE)
  # max_val <- max(taxon_df$env_predict_duration_clamped, na.rm = TRUE)
  # zero_position <- (0 - min_val) / (max_val - min_val)
  # 
  # clamped_env_pred_duration_plot <- ggplot() +
  #   geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = env_predict_duration_clamped)) +
  #   geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  #   scale_fill_gradientn(
  #     colors = c("grey80", "grey80", viridis::plasma(100)),
  #     values = c(0, zero_position - 0.0001, seq(zero_position, 1, length.out = 100)),
  #     breaks = c(-100, 0, 50, 100, 150, 200),
  #     labels = c("Unsuitable Duration \n Modeling", "0", "50", "100", "150", "200"),
  #     name = "Clamped Climatic\nDuration Prediction"
  #   ) +
  #   theme_bw() +
  #   labs(
  #     x = "Longitude",
  #     y = "Latitude",
  #     title = paste("Clamped Climate GLM Duration Prediction Plot for", unique(df$species)))
  # 
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
    
    ## DURATION MESS #############################################################
    taxon_climate_proj <- project(taxon_rast_reproj, crs(temp_seas), method = "near")
    climate_stack_duration_intersection <- intersect(names(climate_stack), predictor_vars_duration)
    duration_monthly_vars <- setdiff(predictor_vars_duration, names(climate_stack))
    
    # Get climate variables as rasters
    if(length(climate_stack_duration_intersection) > 0){
    climate_duration_crop <- crop(climate_stack[[climate_stack_duration_intersection]], taxon_climate_proj)
    taxon_climate_duration_resamp <- resample(taxon_climate_proj, climate_duration_crop[[1]], method = "near")
    climate_duration_masked <- mask(climate_duration_crop, taxon_climate_duration_resamp)
    }
    
    # Create a data frame with coords and monthly vars, then convert to SpatVector
    monthly_df <- taxon_df[, c("longitude", "latitude", duration_monthly_vars), with = FALSE]
    coords_vect <- vect(monthly_df, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    
    # Reproject to match climate rasters
    if(length(climate_stack_duration_intersection) > 0){
    coords_reproj <- project(coords_vect, crs(climate_duration_masked))
    
    # Rasterize monthly variables directly to match climate grid
    monthly_raster_list <- lapply(duration_monthly_vars, function(var){
      rast_var <- rasterize(coords_reproj, climate_duration_masked[[1]], field = var, fun = mean)
      names(rast_var) <- var
      return(rast_var)
    })
    } else{
      coords_reproj <- project(coords_vect, crs(climate_onset_masked)) # use other if the above isnt available
      # Rasterize monthly variables directly to match climate grid
      monthly_raster_list <- lapply(duration_monthly_vars, function(var){
        rast_var <- rasterize(coords_reproj, climate_onset_masked[[1]], field = var, fun = mean)
        names(rast_var) <- var
        return(rast_var)
      })
    }
    
    monthly_window_rast <- rast(monthly_raster_list)
    
    # Combine all variables in correct order
    if(length(climate_stack_duration_intersection) > 0){
    climate_duration_full <- c(climate_duration_masked, monthly_window_rast)
    } else{
      climate_duration_full <- monthly_window_rast
    }
    climate_duration_full <- climate_duration_full[[predictor_vars_duration]]  # Reorder to match
    
    # Convert to old raster format
    climate_duration_stack_old <- raster::stack(climate_duration_full)
    
    # Get training data
    training_vars_duration <- model_out$data_w_duration_vars %>%
      select(all_of(predictor_vars_duration)) %>%
      as.data.frame()
    
    # calc MESS
    duration_mess_raster_old <- dismo::mess(climate_duration_stack_old, training_vars_duration, full = FALSE)
    duration_mess_raster <- rast(duration_mess_raster_old)
    duration_mess_vals <- terra::extract(duration_mess_raster, taxon_df[, .(longitude, latitude)])
    taxon_df[, duration_mess := duration_mess_vals[[2]]]
    
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
    
    mess_plot_duration <- ggplot() +
      geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = duration_mess)) +
      geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
      scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
      theme_minimal() +
      theme_bw() +
      labs(
        x = "Longitude", 
        y = "Latitude", 
        title = paste("MESS Map for Climate GLM duration Prediction Plot for", unique(df$species)))
  } else{
    mess_plot <- NULL
  }
  return(list(
    taxon_df = taxon_df,
    plots = list(
      prediction = list(hop_pred_plot, env_pred_onset_plot, env_pred_duration_plot ),# , clamped_env_pred_onset_plot, clamped_env_pred_offset_plot), 
      quality = list(hop_training_dist_plot, mess_plot_onset, mess_plot_duration)
    )
  ))
  
  
}

