### Function to visualize how the hopkins bins look across the taxa's range

vis_hopkins_delay_bins_across_dist <- function(
                                               path_to_target_taxon_PA_raster = NULL){
  
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
  
  # load North America Elev Raster & Chelsa annual clim vars
  na_elev <- rast("/blue/guralnick/millerjared/PhenoElevation/data/NAelevation4.tif")
  
  
  taxon_vect <- vect(taxon_df, 
                     geom = c("longitude", "latitude"),
                     crs = "EPSG:4326")
  taxon_vect_reproj1 <- project(taxon_vect, crs(na_elev))
  
  # Extract all raster values
  elevation_vals <- terra::extract(na_elev, taxon_vect_reproj1)
  
  # Add all values to these data
  taxon_df[, elevation_m := elevation_vals[,2]]
  
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
  
  # Create map visual showing how hopkins bins actually line up over a map. 
  basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
  basemap <- basemap %>% filter(name %in% c("California", "Nevada", "Arizona", "Sonora", "Baja California", "Baja California Sur"))
  basemap_reproj <- basemap %>% st_transform(st_crs(taxon_rast_reproj))
  
  hop_map_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = hopkins_bin)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    #scale_fill_viridis_c(name = "Hopkins Delay Bins", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Hopkings Delay Binning", unique(df$species)))
  phen_est_data <- model_out$phen_est_data
  taxon_df <- taxon_df %>% 
    mutate(bin_has_empirical_estimate = ifelse(hopkins_bin %in% phen_est_data$hopkins_bin, TRUE, FALSE))
  
  # Create color vector based on your bins
  n_bins <- length(unique(taxon_df$hopkins_bin))
  plasma_colors <- viridis::plasma(n_bins)
  names(plasma_colors) <- sort(unique(taxon_df$hopkins_bin))
  
  # Create a new fill variable that accounts for empirical estimate
  taxon_df <- taxon_df %>%
    mutate(fill_var = ifelse(bin_has_empirical_estimate, 
                             as.character(hopkins_bin), 
                             "No Estimate"))
  
  # Set up color palette
  color_palette <- c(plasma_colors, "No Estimate" = "grey50")
  
  hop_est_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = fill_var)) +
    geom_sf(data = basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_manual(name = "Hopkins Delay Bins", 
                      values = color_palette,
                      breaks = c(as.character(sort(unique(taxon_df$hopkins_bin))), "No Estimate")) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude",
         title = paste("Hopkins Delay Binning", unique(df$species)))
  
  hop_training_dist_plot <- ggplot() +
    geom_tile(data = taxon_df, aes(x = longitude, y = latitude, fill = dist_from_empirical)) +
    geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
    scale_fill_viridis_c(name = "Hopkins Empirical Training Distance Eval", option = "plasma") +
    theme_bw() +
    labs(
      x = "Longitude", 
      y = "Latitude", 
      title = paste("Training Data Distance on Hopkins Onset Predictions for", unique(df$species)))
  
  return(list(
    plots = list(
      hop_map_plot, hop_est_plot, hop_training_dist_plot
    )
  ))
  
}
