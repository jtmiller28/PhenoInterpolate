
### Second V. adding offset
model_baseline_phenophase_w_duration <- function(target_taxa,
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
    left_join(phen_offset_estimates, by = "hopkins_bin") %>% 
    mutate(duration_length = offset_doy - onset_doy)
  
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
    unique(months[2:4])  # Get -1, 0, +1
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
  
  phen_estimates <-  phen_estimates %>% 
    group_by(species) %>% 
    mutate(zduration = scale(duration_length)) %>% 
    filter(between(zduration,-3.25,+3.25))
  
  
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
           duration = duration_length,
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
      corrected_onset_prediction = start_doy + (delay_from_anchor * (actual_delay_rate*4)),
      
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
  
  ### Duration ENV MODEL ##########################################################################
  if(nrow(phen_estimate) > 3){ # min req to make env
    predictors <- c("avg_pet_window", "sd_pet_window", "slope_pet_window", 
                    "avg_tasmax_window","sd_tasmax_window", "slope_tasmax_window","avg_tasmin_window",           
                    "sd_tasmin_window", "slope_tasmin_window", "avg_sfcWind_window","sd_sfcWind_window", "slope_sfcWind_window", 
                    "avg_pr_window","sd_pr_window", "slope_pr_window", 
                    "avg_temp", "avg_diurnal_temp_range", "isothermality",
                    "temp_seas", "tmax", "tmin", "annual_temp_range",
                    "avg_monthly_temp_wettest_q", "avg_monthly_temp_driest_q",
                    "avg_monthly_temp_warmest_q", "avg_monthly_temp_coldest_q",
                    "annual_precip", "precip_wettest_month", "precip_driest_month",
                    "precip_seas", "avg_monthly_precip_wettest_q", "avg_monthly_precip_driest_q",
                    "avg_monthly_precip_warmest_q", "avg_monthly_precip_coldest_q", "total_photoperiod_window") # chopped out "hopkins_delay_midpoint" & "onset" maybe a bad idea?
    
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
      formula_str <- paste("duration ~", paste(pred_combo, collapse = " + "))
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
        
        env_model_duration <- top_models[[1]]$model
        env_model_duration$scaling_params <- scaling_params
        class(env_model_duration) <- c("scaled_glm", class(env_model_duration))
        
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
        env_model_duration <- structure(
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
      phen_estimate$predicted_duration <- predict(env_model_duration, newdata = phen_estimate, type = "response")
      phen_estimate$residuals_duration <- phen_estimate$duration - phen_estimate$predicted_duration
      
      # Calculate metrics
      res_env_duration <- phen_estimate$residuals_duration
      rmse_val_duration <- sqrt(mean(res_env_duration^2))
      
      # Calculate R²
      ss_res_duration <- sum(res_env_duration^2)
      ss_tot_duration <- sum((phen_estimate$duration - mean(phen_estimate$duration))^2)
      r2_val_duration <- 1 - (ss_res_duration / ss_tot_duration)
      
      cat("R² =", round(r2_val_duration, 3), "| RMSE =", round(rmse_val_duration, 3), "\n")
    }
  }
  
  
  # Compute Duration
  phen_estimate <- phen_estimate %>% 
    rename(predicted_onset = predicted) %>% 
    mutate(predicted_offset = predicted_onset + predicted_duration)
  
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
    env_plot_duration <- ggplot(phen_estimate, aes(x = duration, y = predicted_duration)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      annotate("text", x = min(phen_estimate$duration), y = max(phen_estimate$predicted_duration),
               label = sprintf("R² = %.3f\nRMSE = %.2f days", r2_val_duration, rmse_val_duration),
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
      env_model_duration = env_model_duration,
      data = phen_data,
      phen_est_data = phen_estimate, 
      data_w_duration_vars = phen_data_with_windows,
      plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = c(env_plot, env_plot_duration)), 
      metrics = list(
        hopkins_r2 = r2_hop[[2]], 
        hopkins_rsme = rmse_val_hop, 
        env_r2_onset = r2_val, 
        env_rsme_onset = rmse_val, 
        env_r2_duration = r2_val_duration, 
        env_rsme_duration = rmse_val_duration
      )
    ))
  } else{
    # return everything as list
    return(list(
      hopkins_model = hopkins_model, 
      env_model_onset = NULL, 
      env_model_duration =NULL,
      data = phen_data,
      phen_est_data = phen_estimate, 
      data_w_duration_vars = phen_data_with_windows,
      plots = list(hopkins = c(hop_plot1, hop_plot2), enviroment = NULL), 
      metrics = list(
        hopkins_r2 = r2_hop[[2]], 
        hopkins_rsme = rmse_val_hop, 
        env_r2 = NULL, 
        env_rsme = NULL, 
        env_r2_duration = NULL, 
        env_rsme_duration = NULL
      )
    ))
    
  }
  
  
  
}
