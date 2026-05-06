### Extract Climate Data for sample distribution

## load libraries 
library(terra)
library(sf)
library(tidyterra)
library(prism)
library(data.table)

# download dir 
prism::prism_set_dl_dir("/blue/guralnick/millerjared/PhenoInterpolate/data/raw/prism_data/")

# read in cali poppy as our sample annual to extract from
cali_poppy <- fread("/blue/guralnick/millerjared/PhenoInterpolate/data/processed/annotated_cali_poppy.csv")
cali_poppy <- cali_poppy %>% filter(year %in% c(2017:2023))


# Create ALL dates from 2017-2023
all_dates <- seq(as.Date("2017-01-01"), as.Date("2023-12-31"), by = "day")

# Download prism data for all dates
cat("Downloading", length(all_dates), "days of PRISM data...\n")
get_prism_dailys(
  type = "ppt",
  dates = all_dates,
  keepZip = FALSE,
  resolution = "4km"
)
#
# for(chunk in date_chunks) {
#   cat("Downloading", length(chunk), "days...\n")
#   get_prism_dailys(
#     type = "tmean",
#     dates = chunk,
#     keepZip = FALSE,
#     resolution = "800m"
#   )
# }
# download prism data
get_prism_dailys(
  type = "tmean", 
  dates = all_dates, 
  keepZip = FALSE, 
  resolution = "4km"
)

# extract vals
extract_prism_for_occs <- function(dt, prism_type = "tmean") {
  pts <- terra::vect(dt, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  unique_dates <- unique(dt$date)
  result <- rep(NA_real_, nrow(dt))
  
  cat("Processing", length(unique_dates), "unique dates\n")
  
  for(i in seq_along(unique_dates)){
    d <- unique_dates[i]
    date_str <- format(d, "%Y-%m-%d")  # Explicit format
    
    pd <- prism::prism_archive_subset(
      type = prism_type,
      temp_period = "daily",
      dates = date_str,
      resolution = "4km"
    )
    if(length(pd) > 0){
      rast_file <- prism::pd_to_file(pd)
      r <- terra::rast(rast_file)
      idx <- which(dt$date == d)
      values <- terra::extract(r, pts[idx])
      result[idx] <- values[, 2]
      if(i %% 100 == 0) cat("Processed", i, "of", length(unique_dates), "dates\n")
    } else {
      cat("No PRISM data found for", date_str, "\n")
    }
  }
  return(result)
}
# make date col 
cali_poppy[, date := as.Date(paste(year, doy), format = "%Y %j")]
# add prism variables to occurrence data
cali_poppy$tmean <- extract_prism_for_occs(cali_poppy, "tmean")
cali_poppy$ppt <- extract_prism_for_occs(cali_poppy, "ppt")
