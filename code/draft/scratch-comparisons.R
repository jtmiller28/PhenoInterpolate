# Load Libraries 
library(data.table)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/proxy-baseline-fxns.R")

bee_occs <- fread("/blue/soltis/millerjared/beeBDC/Bee-data/OutputData/05_cleaned_database.csv")
# additional data on Anthophora hololeuca provided by Michael Orr
extra_hololeuca_data <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/raw/a-hololeuca-from-Orr.csv")
phenovision_flowering_data <- fread("/blue/guralnick/millerjared/phenovision_data/phenovision_flowers_03_15_2026.csv")
cch2_annotated_specimens <- fread("/blue/soltis/millerjared/SpatioTemporalTradeoffs/data/processed/cch2-data/cch2-cleaned-and-aligned-data.csv")


## Load in North American Taxonomy
na_taxonomy <- fread("/blue/guralnick/millerjared/PlantSweepeR/data/processed/wcvp-ncbi-alignment-na.csv")


target_taxa <- "Anthophora hololeuca"

described_phen_months <- fread("/blue/soltis/millerjared/Legume-Specialist-Occupancy/data/raw/hololeuca-and-host-pheno.csv")

taxon_phen_hist(target_taxa, bees_or_plants = "bees", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "bees", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 21, # look at histogram when making this decision
                                      max_doy = 234 )

model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      bees_or_plants = "bees", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 21, # look at histogram when making this decision
                                      max_doy = 234 )

model_out$plots$hopkins
model_out$plots$enviroment

model_out$env_model

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality
# because I dont want to wait all day. 
ggplot() +
  geom_tile(data = anth_holo_df, aes(x = longitude, y = latitude, fill = mess)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "MESS eval", option = "plasma") +
  theme_minimal() +
  labs(
    x = "Longitude", 
    y = "Latitude")

## Pull target taxa, order with the accepted parent name in lead
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorothamnus polydenius") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 38, # look at histogram when making this decision
                                      max_doy = 244 )
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_polydenius_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality

target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorothamnus emoryi") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      #phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 12, # look at histogram when making this decision
                                      max_doy = 185 )
model_out$plots$hopkins
model_out$plots$enviroment

# dont predict this taxa currently...
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron schottii") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      #phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 15, # look at histogram when making this decision
                                      max_doy = 215)
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_schottii_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron fremontii") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      #phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 19, # look at histogram when making this decision
                                      max_doy = 238 )
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_fremontii_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron spinosum") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      #phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 86, # look at histogram when making this decision
                                      max_doy = 248 )
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_spinosus_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron arborescens") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 30, # look at histogram when making this decision
                                      max_doy = 228 )

model_out <- model_baseline_phenology_w_offset(target_taxa, 
                                      #phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 30, # look at histogram when making this decision
                                      max_doy = 228 )
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_arborescens_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Rafinesquia neomexicana") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)
#target_taxa <- "Anthophora hololeuca"
model_out <- model_baseline_phenology(target_taxa, 
                                      phenometric = 0.10, 
                                      bees_or_plants = "plants", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 0, # look at histogram when making this decision
                                      max_doy = 217 )
model_out$plots$hopkins
model_out$plots$enviroment

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/guralnick/share/SDM_pipeline/moj_son_SDMs/Models/Rafinesquia_neomexicana/Rafinesquia_neomexicana_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = FALSE)

distribution_pred$plots$prediction
distribution_pred$plots$quality



