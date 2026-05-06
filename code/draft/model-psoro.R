### Model Psorothamnus Phenology

## Load packages
library(tidyverse)
library(dbscan)
library(phenesse)
library(mgcv) 
library(sf)
library(terra)
library(data.table)

psoro_data <- fread("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/annotated_psoro_sp.csv")
psoro_data <- psoro_data %>% filter(year >= 2017 & year < 2024) %>% filter(doy <= 200)

psoro_emoryi <- psoro_data %>% filter(species == "Psorothamnus emoryi")
psoro_emoryi <- psoro_emoryi %>% 

extract_daymet_climate <- function(df, daymet_dir = "/orange/guralnick/daymet-daily/") {
  
  df_with_climate <- df %>%
    group_by(year) %>%
    group_split() %>%
    map_dfr(function(yr_df) {
      yr <- unique(yr_df$year)
      
      # Load Daymet rasters for this year
      tmax_file <- file.path(daymet_dir, paste0("daymet_v4_daily_na_tmax_", yr, ".nc"))
      tmin_file <- file.path(daymet_dir, paste0("daymet_v4_daily_na_tmin_", yr, ".nc"))
      prcp_file <- file.path(daymet_dir, paste0("daymet_v4_daily_na_prcp_", yr, ".nc"))
      
      tmax <- rast(tmax_file)
      tmin <- rast(tmin_file)
      prcp <- rast(prcp_file)
      
      # Create spatial points
      coords <- yr_df %>% select(longitude, latitude)
      pts <- vect(coords, geom = c("longitude", "latitude"), crs = "EPSG:4326")
      
      # Extract climate for each observation's window (DOY 1 to observed DOY)
      climate_list <- map2_dfr(1:nrow(yr_df), yr_df$doy, function(i, obs_doy) {
        
        # Get layers from DOY 1 to observed DOY
        window_layers <- 1:obs_doy
        
        tmax_window <- tmax[[window_layers]]
        tmin_window <- tmin[[window_layers]]
        prcp_window <- prcp[[window_layers]]
        
        # Extract for this specific point
        tmax_vals <- extract(tmax_window, pts[i], ID = FALSE)
        tmin_vals <- extract(tmin_window, pts[i], ID = FALSE)
        prcp_vals <- extract(prcp_window, pts[i], ID = FALSE)
        
        tibble(
          mean_temp = mean((as.numeric(tmax_vals) + as.numeric(tmin_vals)) / 2, na.rm = TRUE),
          mean_ppt = sum(as.numeric(prcp_vals), na.rm = TRUE)
        )
      })
      
      bind_cols(yr_df, climate_list)
    })
  
  return(df_with_climate)
}

# Apply climate extraction
psoro_emoryi_clim <- extract_daymet_climate(psoro_emoryi)

## Modified cluster function for real data
cluster_obs <- function(df, eps_km = 75, min_pts = 10) {
  
  # Rename columns to match function expectations
  df <- df %>%
    rename(lon = longitude, lat = latitude, doy_observed = doy)
  
  eps_deg <- eps_km / 111
  cluster_list <- list()
  global_id <- 0
  
  for(yr in unique(df$year)){
    yr_df <- df %>% filter(year == yr)
    
    coords <- yr_df %>% select(lon, lat) %>% as.matrix()
    db <- dbscan::dbscan(coords, eps = eps_deg, minPts = min_pts)
    yr_df$cluster <- db$cluster
    
    for (cl in setdiff(unique(yr_df$cluster), 0)) {
      global_id <- global_id + 1
      cl_df <- yr_df %>% filter(cluster == cl)
      
      area_km2 <- tryCatch({
        pts_sf <- st_as_sf(cl_df, coords = c("lon", "lat"), crs = 4326)
        hull <- st_convex_hull(st_union(pts_sf))
        hull_proj <- st_transform(hull, crs = 5070)
        as.numeric(st_area(hull_proj)) / 1e6
      }, error = function(e) pi * (eps_km/2)^2)
      
      cluster_list[[global_id]] <- tibble(
        cluster_id = global_id,
        year       = yr,
        n_obs      = nrow(cl_df),
        mean_lat   = mean(cl_df$lat),
        mean_lon   = mean(cl_df$lon),
        mean_temp  = mean(cl_df$mean_temp),
        mean_ppt   = mean(cl_df$mean_ppt),
        area_km2   = area_km2,
        doys       = list(cl_df$doy_observed)
        # Removed true_mus - only for simulation
      )
    }
  }
  bind_rows(cluster_list)
}

# Run clustering
clusters <- cluster_obs(psoro_emoryi_clim, eps_km = 10000, min_pts = 7)
