### Title: Prep Data
### Author: JT Miller
### Date: 01-01-2026

# Purpose: Build psuedo-flowering occurrence datasets. These are based on real data for clustering, but will have 'artificial' response curves for doy to climate to simulate how well we can cluster and analyze these data

## Libraries 
library(data.table)
library(sf)
library(tidyverse)
library(ggplot2)
library(rgnparser)
library(arrow)

## Set path for rgnparser tool
my_path <- Sys.getenv("PATH") # grab our path
Sys.setenv(PATH = paste0(my_path, "/home/millerjared/gnparser"))

## Bring in plant phenology annotated occurrence datasets, align taxonomy.
# read in
annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/phenovision-data/flower_inference_formatted_edit_9.8.25.csv")
cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")
# extract names
specimen_names <- annotated_specimens %>% select(scientificName) %>% distinct() 
# parse names
specimen_names_parsed <- gn_parse_tidy(specimen_names$scientificName)
specimen_names_parsed <- specimen_names_parsed %>% rename(scientificName = verbatim)
# Merge these back with the specimen dataset, bringing in it in as verbatimSpeciesName 
annotated_specimens2 <- annotated_specimens %>% 
  left_join(select(specimen_names_parsed, scientificName, canonicalfull), by = c("scientificName")) %>% 
  mutate(name = canonicalfull)
# Use my Aligned Taxonomic table (works with Smith trees if we need pglmms later)
bocp_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-2025.csv")
bocp_taxonomy <- bocp_taxonomy %>% # remove ncbi, we dont need this here.
  filter(source != "ncbi")
# join annotated specimens
annotated_specimens3 <- annotated_specimens2 %>% 
  left_join(select(bocp_taxonomy, name, alignedParentName), by = c("name")) %>% 
  distinct(, .keep_all = TRUE) # sometimes theres extra names due to how I built the name relation table, get rid of dupes
# rename
annotated_specimen_flowers <- annotated_specimens3 %>% 
  rename(species = alignedParentName)# these specimens are already prepped see prep-cch2-data.R for details...
annotated_specimens_flowers_cch2 <- cch2_annotated_specimens 
# grab doy from the cch2 data
annotated_specimens_flowers_cch2 <- annotated_specimens_flowers_cch2 %>% 
  mutate(doy = lubridate::yday(ymd(eventDate))) %>% 
  mutate(year = lubridate::year(ymd(eventDate))) %>% 
  mutate(day = lubridate::day(ymd(eventDate)))
# Filter for only flowering data
annotated_specimen_flowers <- annotated_specimen_flowers %>% 
  filter(trait == "flowers present")
# break out day 
annotated_specimen_flowers <- annotated_specimen_flowers %>% 
  mutate(day = lubridate::day(ymd(date)))
# combine annotated herbarium specimen phenovision data & cch2 annotated specimen data
annotated_specimen_flowers2 <- annotated_specimen_flowers %>% 
  select(species, longitude, latitude, year, doy = dayOfYear, day)
annotated_specimens_flowers_cch2 <- annotated_specimens_flowers_cch2 %>% 
  select(species, longitude = decimalLongitude, latitude = decimalLatitude, year, doy, day)
annotated_specimens_flowers <- rbind(annotated_specimen_flowers2, annotated_specimens_flowers_cch2)
# remove first dates 
annotated_specimens_flowers_filtered <- annotated_specimens_flowers %>% 
  filter(!is.na(day) & day != 1)
# Gather up the phenovision inaturalist data...
phenologycv_csv <- open_dataset(
  sources = "/blue/guralnick/millerjared/PhenoElevation/data/phenobase-annotations-headers-added-01-2025.csv", 
  col_types = schema(ISBN = string()),
  format = "csv")
phenologycv_csv |>
  group_by(year) |>
  write_dataset(path = "/blue/guralnick/millerjared/PhenoInterpolate/data/processed/phenologycv_presence_data", format = "parquet")
phenologycv_pd <- open_dataset("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/phenovision-data/phenologycv_presence_data")
# filter for flowering only data
phenologycv_pd2 <- phenologycv_pd  |> 
  filter(trait == "flower")
# collect into memory
phenologycv_pds3 <- phenologycv_pd2 |> collect()
# clean up record-level dupes
annotated_inat_flowers <- phenologycv_pds3 %>% 
  group_by(scientific_name) %>% 
  distinct(observed_metadata_url, .keep_all = TRUE) %>% 
  ungroup() %>% 
  rename(name = scientific_name)
# harmonize with bocp taxonomy
annotated_inat_flowers <- annotated_inat_flowers %>% 
  left_join(select(bocp_taxonomy, name, alignedParentName), by = c("name")) %>% 
  distinct(, .keep_all = TRUE)  # sometimes theres extra names due to how I built the name relation table, get rid of dups
# break out date
annotated_inat_flowers <- annotated_inat_flowers %>% 
  mutate(day = lubridate::day(ymd(date)))
# filter and combine with our previous
annotated_inats_flowers <- annotated_inat_flowers %>%  
  select(species = alignedParentName, longitude, latitude, year, doy = day_of_year, day)
# and bring together
annotated_flowers <- rbind(annotated_specimens_flowers_filtered , annotated_inats_flowers)
# create a model dataset
annotated_cali_poppy <- annotated_flowers %>% 
  filter(species == "Eschscholzia californica")

annotated_psoro_sp <- annotated_flowers %>% 
  filter(species %in% c("Psorodendron arborescens", 
                        "Psorodendron fremontii", 
                        "Psorodendron schottii", 
                        "Psorodendron spinosum", 
                        "Psorothamnus emoryi", 
                        "Psorothamnus polydenius"))



# save these data
fwrite(annotated_flowers, "/blue/guralnick/millerjared/PhenoInterpolate/data/processed/annotated_flowers.csv")
fwrite(annotated_cali_poppy, "/blue/guralnick/millerjared/PhenoInterpolate/data/processed/annotated_cali_poppy.csv")
fwrite(annotated_psoro_sp, "/blue/guralnick/millerjared/PhenoInterpolate/data/processed/annotated_psoro_sp.csv")
