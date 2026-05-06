### Extract Env Thresholds #####################################################
#' Extracts environmentally-grounded parameter values from species range climate data
#'
#' @param tmean_species SpatRaster, daily tmean layers cropped to species range
#' @param ppt_species SpatRaster, daily ppt layers cropped to species range
#' @param cell_pts sf object, cell center points within species range
#' @param years integer vector, years to process (e.g. 2017:2020)
#' @param days_per_year integer vector, number of days per year (same length as years)
#' @param spring_doy integer vector, DOY range defining spring (e.g. 30:150)
#' @param ppt_low_pct numeric, lower percentile for ppt_min (default 0.05)
#' @param ppt_high_pct numeric, upper percentile for ppt_max (default 0.95)
#' @param doy_range numeric vector length 2, expected flowering DOY range for slope scaling (default c(120, 180))
#' @return named list compatible with params structure
#' @examples
#' # extract_env_thresholds(tmean_species, ppt_species, cell_pts, years = 2017:2020,
#' #   days_per_year = c(365,365,365,366), spring_doy = 30:150)
extract_env_thresholds <- function(tmean_species,
                                   ppt_species,
                                   cell_pts,
                                   years,
                                   days_per_year,
                                   spring_doy,
                                   ppt_low_pct  = 0.05,
                                   ppt_high_pct = 0.95,
                                   doy_range    = c(120, 180), 
                                   plot = FALSE) {
  
  # cumulative layer start positions per year
  year_starts <- c(1, cumsum(days_per_year[-length(days_per_year)]) + 1)
  
  # storage
  all_tmean <- c()
  all_ppt   <- c()
  
  for (i in seq_along(years)) {
    spring_layers_year <- year_starts[i] + spring_doy - 1
    
    tmean_spring_r <- mean(tmean_species[[spring_layers_year]])
    ppt_spring_r   <- sum(ppt_species[[spring_layers_year]])
    
    all_tmean <- c(all_tmean, terra::extract(tmean_spring_r, vect(cell_pts))[, 2])
    all_ppt   <- c(all_ppt,   terra::extract(ppt_spring_r,   vect(cell_pts))[, 2])
  }
  
  # drop NAs
  all_tmean <- all_tmean[!is.na(all_tmean)]
  all_ppt   <- all_ppt[!is.na(all_ppt)]
  
  # ppt thresholds from user-defined percentiles
  ppt_min <- quantile(all_ppt, ppt_low_pct)
  ppt_max <- quantile(all_ppt, ppt_high_pct)
  
  # temp slope: scale so that the temp range (p05 to p95) spans the doy_range
  temp_p05  <- quantile(all_tmean, 0.05)
  temp_p95  <- quantile(all_tmean, 0.95)
  temp_range <- temp_p95 - temp_p05
  temp_slope <- -diff(doy_range) / temp_range  # negative: warmer = earlier
  
  # intercept: anchored so mean temp produces midpoint of doy_range
  intercept <- mean(doy_range) - temp_slope * mean(all_tmean)
  
  # build output params
  out_params <- list(
    intercept        = round(intercept, 1),
    temp_slope       = round(temp_slope, 2),
    sigma_spatial    = 0,
    sigma_individual = 0,
    ppt_min          = round(unname(ppt_min), 1),
    ppt_max          = round(unname(ppt_max), 1),
    ppt_max_penalty  = 60
  )
  
  # --- optional plots ---
  if (plot) {
    df <- data.frame(tmean_spring = all_tmean, ppt_spring = all_ppt)
    
    stats_tmean <- data.frame(
      mean = mean(all_tmean),
      p05  = unname(temp_p05),
      p95  = unname(temp_p95)
    )
    
    stats_ppt <- data.frame(
      mean = mean(all_ppt),
      p05  = unname(ppt_min),   # already computed at user percentiles
      p95  = unname(ppt_max)
    )
    
    p1 <- ggplot(df, aes(x = tmean_spring)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
      geom_vline(xintercept = stats_tmean$mean, color = "black",     linewidth = 1,   linetype = "solid") +
      geom_vline(xintercept = stats_tmean$p05,  color = "firebrick", linewidth = 0.8, linetype = "dashed") +
      geom_vline(xintercept = stats_tmean$p95,  color = "firebrick", linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = stats_tmean$mean, y = Inf,
               label = paste0("Mean: ", round(stats_tmean$mean, 1)),
               vjust = 2, hjust = -0.1, size = 3.5) +
      annotate("text", x = stats_tmean$p05, y = Inf,
               label = paste0(ppt_low_pct * 100, "th: ", round(stats_tmean$p05, 1)),
               vjust = 2, hjust = -0.1, color = "firebrick", size = 3.5) +
      annotate("text", x = stats_tmean$p95, y = Inf,
               label = paste0(ppt_high_pct * 100, "th: ", round(stats_tmean$p95, 1)),
               vjust = 2, hjust = -0.1, color = "firebrick", size = 3.5) +
      labs(title = "Distribution of Mean Spring Temperature",
           x = "Mean Spring Tmean (°C)", y = "Density") +
      theme_classic()
    
    p2 <- ggplot(df, aes(x = ppt_spring)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 50, fill = "darkgreen", alpha = 0.7, color = "white") +
      geom_vline(xintercept = stats_ppt$mean, color = "black",     linewidth = 1,   linetype = "solid") +
      geom_vline(xintercept = stats_ppt$p05,  color = "firebrick", linewidth = 0.8, linetype = "dashed") +
      geom_vline(xintercept = stats_ppt$p95,  color = "firebrick", linewidth = 0.8, linetype = "dashed") +
      annotate("text", x = stats_ppt$mean, y = Inf,
               label = paste0("Mean: ", round(stats_ppt$mean, 1)),
               vjust = 4, hjust = -0.1, size = 3.5) +
      annotate("text", x = stats_ppt$p05, y = Inf,
               label = paste0(ppt_low_pct * 100, "th (ppt_min): ", round(stats_ppt$p05, 1)),
               vjust = 2, hjust = -0.1, color = "firebrick", size = 3.5) +
      annotate("text", x = stats_ppt$p95, y = Inf,
               label = paste0(ppt_high_pct * 100, "th (ppt_max): ", round(stats_ppt$p95, 1)),
               vjust = 2, hjust = -0.5, color = "firebrick", size = 3.5) +
      labs(title = "Distribution of Total Spring Precipitation",
           x = "Total Spring PPT (mm)", y = "Density") +
      theme_classic()
    
    print(p1 + p2)
  }
  
  return(out_params)
}

################################################################################

### Simulate Flowering Population ############################################
#' A function that builds a normal phenology curve based on env vars temperature and ppt
#' 
#' @param temp numeric, avg temperature value for spring
#' @param ppt numeric, cumulative precipitation value for spring
#' @param n_individuals, number of individuals from which we assume the population consists for building the normal distribution
#' @param par, parameter set to feed into simulation. Requires intercept, temp_slope, sigma_spatial, sigma_individual, ppt_min, ppt_max, ppt_max_penalty
#' @return normal distribution simu data
#' examples 
#' # sim_pop_flowering(t, p, n_individuals = 500, par = params) 
sim_pop_flowering <- function(temp, ppt, n_individuals = 100, par = params){
  # temp effect on slope as a linear y = mx + b
  mu_site <- par$temp_slope * temp + par$intercept 
  # make a 'gate' for ppt, basically a threshold must be met for ppt, but too much results in the same penalty
  ppt_mid <- (par$ppt_min + par$ppt_max) / 2  # peak of parabola (minimum penalty)
  # if else regarding the penalty
  ppt_penalty <- if(ppt < par$ppt_min){
    par$ppt_max_penalty  # below minimum: flat max penalty
  } else if(ppt > par$ppt_max){
    par$ppt_max_penalty  # above maximum: flat max penalty
  } else {
    # parabola between ppt_min and ppt_max, 0 at midpoint, ppt_max_penalty at edges
    par$ppt_max_penalty * ((ppt - ppt_mid) / (ppt_mid - par$ppt_min))^2
  }
  # augment our site phenology based on the penalty value 
  mu_site <- mu_site + ppt_penalty
  # pull a normal distribution based on mu_site for the number of individuals, with some variability for sigma_ind
  flowering_times <- rnorm(n_individuals,
                           mean = mu_site,
                           sd = par$sigma_individual)
  # return these normal pulls
  return(flowering_times)
}

################################################################################

### Cluster Data Spatially #####################################################
#' A function that clusters the obs data based on spatial clustering using dbscan
#' 
#' @param data dataframe, contains simu data with observation layer
#' @param distance_km numeric, distance in km for upper limit on cluster size
#' @param min_obs, required number of observations for a cluster to exist
#' @param min_cells, number of grid cells within distance_km required to cluster to form a cluster.
#' @return normal distribution simu data
#' examples 
#' # cluster_spatial(data, distance_km = 50, min_obs = 7, min_cells = 3)
cluster_spatial <- function(data, distance_km = 50, min_obs = 7, min_cells = 3){
  data_clustered <- data %>% 
    group_by(year) %>% 
    group_modify(~{
      # get unique cell locs for clustering
      unique_cells <- .x %>% 
        distinct(X,Y, .keep_all = TRUE) %>% 
        dplyr::select(X,Y, tmean_spring, ppt_spring)
      coords <- as.matrix(unique_cells[, c("X", "Y")])
      # convert dist to meters
      eps_m <- distance_km * 1000
      # DBSCAN clustering on cells, instead of individual obs
      clusters <- dbscan(coords, eps = eps_m, minPts = min_cells)
    # add cluster IDs to unique cells
    unique_cells$cluster_id <- clusters$cluster
    # join cluster IDs back to all obs
    .x %>% 
      left_join(unique_cells %>% dplyr::select(X,Y,cluster_id), by = c("X","Y")) %>% 
      filter(cluster_id != 0) # removes noise cluster
  }) %>% 
  ungroup() %>% 
  # create unique cluster id across years
  mutate(cluster_full_id = paste0(year, "_", cluster_id)) %>% 
  # filter to clusters with min total obs
  group_by(cluster_full_id) %>% 
  filter(n() >= min_obs) %>% 
  # calc spatial metrics per cluster
  mutate(
    cluster_n_obs = n(), # total obs in cluster
    cluster_n_cells = n_distinct(X,Y), # number of cells in cluster
    # Calc cell level centroids for spatial metrics 
    cluster_mean_dist_to_center = {
      cell_coords <- distinct(pick(X,Y))
      centroid_x <- mean(cell_coords$X)
      centroid_y <- mean(cell_coords$Y)
      sqrt(mean((cell_coords$X - centroid_x)^2 + (cell_coords$Y - centroid_y)^2)) / 1000
    }, 
    # mean nearest neighbor dist between cells (km)
    cluster_mean_nn_dist = {
      cell_coords <- distinct(pick(X, Y))
      if(nrow(cell_coords) > 1) {
        coords_mat <- as.matrix(cell_coords)
        dist_mat <- as.matrix(dist(coords_mat))
        diag(dist_mat) <- NA
        mean(apply(dist_mat, 1, min, na.rm = TRUE)) / 1000
      } else {
        NA_real_
      }
    },
    
    # Inverse for weighting (lower distance = higher weight)
    cluster_compactness = 1 / (cluster_mean_dist_to_center + 0.1)
  ) %>%
  ungroup()

return(data_clustered)
}
### Cluster Data by Enviromental Similarity #####################################################
#' A function that clusters the obs data based on env clustering using dbscan
#' 
#' @param data dataframe, contains simu data with observation layer
#' @param env_range_h numeric, the height at which the dendrogram produced by hierarchical clustering is cut, this height represents the dissimilarity threshold. Lower value = Lower envirometal dist
#' base value of 2 means that cells are grouped together if their within 2 standard deviations of each other in combined env var space
#' @param min_obs, required number of observations for a cluster to exist
#' @param min_cells, number of grid cells within distance_km required to cluster to form a cluster.
#' @return normal distribution simu data
#' examples 
#' # cluster_env(data, env_range_h = 2, min_obs = 7, min_cells = 3)
## Environmental Clustering 
cluster_env <- function(data, env_range_h = 2, min_obs = 7, min_cells = 3){
  data_clustered <- data %>%
    group_by(year) %>%
    group_modify(~{
      # get unique cell locations with their env vars
      unique_cells <- .x %>%
        distinct(X, Y, .keep_all = TRUE) %>%
        dplyr::select(X, Y, tmean_spring, ppt_spring)
      
      # scale both variables so neither dominates by units
      env_scaled <- scale(unique_cells[, c("tmean_spring", "ppt_spring")])
      
      # run hierarchical clustering on scaled env space
      env_dist <- dist(env_scaled)
      hc <- hclust(env_dist, method = "complete")
      clusters <- cutree(hc, h = env_range_h)
      
      unique_cells$cluster_id <- clusters
      
      .x %>%
        left_join(unique_cells %>% dplyr::select(X, Y, cluster_id), by = c("X", "Y"))
    }) %>%
    ungroup() %>%
    mutate(cluster_full_id = paste0(year, "_", cluster_id)) %>%
    group_by(cluster_full_id) %>%
    filter(n_distinct(X, Y) >= min_cells,
           n() >= min_obs) %>%
    mutate(
      cluster_n_obs   = n(),
      cluster_n_cells = n_distinct(X, Y),
      # Temperature metrics
      cluster_temp_mean  = {
        cell_temps <- distinct(pick(X, Y, tmean_spring))$tmean_spring
        mean(cell_temps)
      },
      cluster_temp_sd  = {
        cell_temps <- distinct(pick(X, Y, tmean_spring))$tmean_spring
        sd(cell_temps)
      },
      cluster_temp_range = {
        cell_temps <- distinct(pick(X, Y, tmean_spring))$tmean_spring
        max(cell_temps) - min(cell_temps)
      },
      # Precipitation metrics
      cluster_ppt_mean = {
        cell_ppts <- distinct(pick(X, Y, ppt_spring))$ppt_spring
        mean(cell_ppts)
      },
      cluster_ppt_sd = {
        cell_ppts <- distinct(pick(X, Y, ppt_spring))$ppt_spring
        sd(cell_ppts)
      },
      cluster_ppt_range = {
        cell_ppts <- distinct(pick(X, Y, ppt_spring))$ppt_spring
        max(cell_ppts) - min(cell_ppts)
      },
      # Weight: prefer homogeneous clusters in both dimensions
      # combines temp and ppt sd so clusters tight in both get higher weight
      cluster_env_weight = 1 / (cluster_temp_sd + cluster_ppt_sd + 0.1)
    ) %>%
    ungroup()
  return(data_clustered)
}
################################################################################

### Plot MESS ##################################################################
plot_mess <- function(mess_sf, dist_col, title, file){
  max_val <- max(mess_sf[[dist_col]], na.rm = TRUE)
  p <- ggplot() +
    geom_sf(data = nevada, fill = "gray90", color = "black") +
    geom_sf(data = species_range, fill = NA, color = "black", size = 1) +
    geom_sf(data = mess_sf, mapping = aes(color = .data[[dist_col]]), size = 0.5) +
    scale_color_gradientn(
      colors = c("darkgreen", "yellow", "darkred"),
      values = scales::rescale(c(0, 0.001, max_val)),
      name = "Extrapolation\nDistance"
    ) +
    facet_wrap(~year) +
    ggtitle(title) +
    theme_bw()
  print(p)
  ggsave(plot = p, file = file, width = 15, height = 15)
}
################################################################################
