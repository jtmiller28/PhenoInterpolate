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