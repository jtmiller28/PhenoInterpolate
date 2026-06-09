### Model PhenoPhase Pipeline
# Load Libraries 
library(data.table)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/model-baseline-phenophase.R")
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/predict-phenophase-across-distribution.R")
source("/blue/guralnick/millerjared/PhenoInterpolate/code/draft/taxon-phen-hist.R")

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
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-phen-hists.png", height = 10, width = 10)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                      bees_or_plants = "bees", # target taxon group to pull empirical data for
                                      phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                      cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                      bee_occs, # input bee occurrence data
                                      extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                      min_doy = 21, # look at histogram when making this decision
                                      max_doy = 234 )

# Save Visuals
model_out$plots$hopkins[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hopkins-response.png", height = 10, width = 10)
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/a-holo-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/a-holo-predictions.csv")

vis_hop_samp <- vis_hopkins_delay_bins_across_dist(path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/bees/Anthophora_hololeuca_SDM_PA.tif")
vis_hop_samp[[1]][1]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hopkins-delay-bins-full.png", height = 10, width = 10)
vis_hop_samp[[1]][2]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hopkins-delay-bins-gaps.png", height = 10, width = 10)
vis_hop_samp[[1]][3]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-hopkins-delay-bins-extraps.png", height = 10, width = 10)
### Plants #####################################################################
target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorothamnus polydenius") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 38, # look at histogram when making this decision
                                                  max_doy = 244 )

# Save Visuals
model_out$plots$hopkins[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-hopkins-response.png", height = 10, width = 10)
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-p-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_polydenius_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-p-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-p-predictions.csv")

target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorothamnus emoryi") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 12, # look at histogram when making this decision
                                                  max_doy = 185 )

# Save Visuals
model_out$plots$hopkins[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-hopkins-response.png", height = 10, width = 10)
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-e-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_emoryi_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-e-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-e-predictions.csv")


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron schottii") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 15, # look at histogram when making this decision
                                                  max_doy = 215 )

# Save Visuals
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sc-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_schottii_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sc-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sc-predictions.csv")


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron fremontii") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 19, # look at histogram when making this decision
                                                  max_doy = 238 )

# Save Visuals
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-f-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_fremontii_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-f-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-f-predictions.csv")


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron spinosum") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 86, # look at histogram when making this decision
                                                  max_doy = 248 )

# Save Visuals
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sp-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_spinosus_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-sp-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sp-predictions.csv")


target_taxa <- na_taxonomy %>% 
  filter(alignedParentName == "Psorodendron arborescens") %>% 
  arrange(desc(nameStatus == "Accepted"), desc(alignedParentName == name)) %>% 
  pull(name)

taxon_phen_hist(target_taxa, bees_or_plants = "plants", described_phen_months = described_phen_months)

model_out <- model_baseline_phenophase_w_duration(target_taxa, 
                                                  bees_or_plants = "plants", # target taxon group to pull empirical data for
                                                  phenovision_flowering_data, # input inat flowering annotations for flowering plants
                                                  cch2_annotated_specimens,  # input specimen flowering annotations for flowering plants (CA)
                                                  bee_occs, # input bee occurrence data
                                                  extra_hololeuca_data, # extra A. hololeuca occurrence data provided by M. Orr
                                                  min_doy = 30, # look at histogram when making this decision
                                                  max_doy = 228 )

# Save Visuals
model_out$plots$hopkins
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-hopkins-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-clim-onset-r2-plot.png", height = 10, width = 10)
model_out$plots$enviroment[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-clim-duration-r2-plot.png", height = 10, width = 10)
fwrite(model_out$metrics, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-a-metrics.csv")

distribution_pred <- predict_pheno_for_dist(model_out, 
                                            path_to_target_taxon_PA_raster = "/blue/soltis/millerjared/Legume-Specialist-Occupancy/outputs/sdms/plants/Psorothamnus_arborescens_SDM_PA.tif", 
                                            hopkins = TRUE, 
                                            env = TRUE)
distribution_pred$plots$prediction[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-hop-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-env-onset-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$prediction[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-env-duration-pred-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[1]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-hop-training-dist-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[2]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-onset-mess-plot.png", height = 10, width = 10)
distribution_pred$plots$quality[[3]]
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/psoro-a-offset-mess-plot.png", height = 10, width = 10)
fwrite(distribution_pred$taxon_df, "/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-a-predictions.csv")




### Read in taxa predictions
a_holo_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/a-holo-predictions.csv")
psoro_a_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-a-predictions.csv")
psoro_f_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-f-predictions.csv")
psoro_p_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-p-predictions.csv")
psoro_e_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-e-predictions.csv")
psoro_sc_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sc-predictions.csv")
psoro_sp_pred <- fread("/blue/guralnick/millerjared/PhenoInterpolate/outputs/metrics/psoro-sp-predictions.csv")

# First, make all MESS < 0 NA for duration 
a_holo_pred <- a_holo_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_a_pred <- psoro_a_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_f_pred <- psoro_f_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_p_pred <- psoro_p_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_e_pred <- psoro_e_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_sc_pred <- psoro_sc_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

psoro_sp_pred<- psoro_sp_pred %>% 
  mutate(env_predict_onset = ifelse(onset_mess > 0, env_predict_onset, NA)) %>% 
  mutate(env_predict_duration = ifelse(duration_mess > 0, env_predict_duration, NA)) %>% 
  mutate(env_predict_duration = ifelse(env_predict_duration < 0, NA, env_predict_duration))

basemap <- rnaturalearth::ne_states(c("United States of America", "Mexico"), returnclass = "sf")
basemap <- basemap %>% filter(name %in% c("California", "Nevada", "Arizona", "Sonora", "Baja California", "Baja California Sur"))
basemap_reproj <- basemap %>% st_transform(st_crs(4326))

ggplot() +
  geom_tile(data = a_holo_pred, aes(x = longitude, y = latitude, fill = env_predict_onset)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Climatic Onset Prediction", option = "plasma") +
  ggtitle("Anthophora hololeuca Onset Prediction") +
  theme_bw()

ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-predict-onset-map.png", height = 10, width = 10)

ggplot() +
  geom_tile(data = a_holo_pred, aes(x = longitude, y = latitude, fill = env_predict_duration)) +
  geom_sf(basemap_reproj, mapping = aes(), fill = NA, color = "black") +
  scale_fill_viridis_c(name = "Climatic Duration Prediction", option = "plasma") +
  theme_bw()
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/a-holo-predict-duration-map.png", height = 10, width = 10)


a_holo_eg <- a_holo_pred %>% select(env_predict_onset, env_predict_offset) %>% head(n = 1)
psoro_p_eg <- psoro_p_pred %>% select(env_predict_onset, env_predict_offset) %>% head(n = 1)
psoro_f_eg <- psoro_f_pred %>% select(env_predict_onset, env_predict_offset) %>% head(n = 100)
psoro_f_eg <- psoro_f_eg[96,]
psoro_a_eg <- psoro_a_pred %>% select(env_predict_onset,env_predict_offset) %>% head(n = 1)
# Combine the data
plot_data <- bind_rows(
  a_holo_eg %>% mutate(species = "A. hololeuca"),
  psoro_p_eg %>% mutate(species = "P. polydenius"), 
  psoro_a_eg %>% mutate(species = "P. arborescens"),
  psoro_f_eg %>%  mutate(species = "P. fremontii")
)

# Create the plot
ggplot(plot_data, aes(y = species)) +
  geom_segment(aes(x = env_predict_onset, xend = env_predict_offset, 
                   yend = species, color = species),
               linewidth = 8) +
  scale_color_manual(values = c("A. hololeuca" = "steelblue", 
                                "P. polydenius" = "goldenrod", 
                                "P. arborescens" = "lightpink", 
                                "P. fremontii" = "lightgreen")) +
  scale_x_continuous(limits = c(1, 365), breaks = seq(0, 365, by = 50)) +
  labs(x = "Day of Year", y = NULL) +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Phenology Overlap in a pixel")
ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/pixel-overlap-conceptual-figure.png")

# Generate normal distribution data (Day of Year scale)
x <- seq(0, 365, length.out = 1000)
mean_doy <- 182.5  # mid-year
sd_doy <- 60       # adjust spread as needed
y <- dnorm(x, mean = mean_doy, sd = sd_doy)
df_norm <- data.frame(x = x, y = y)

# Calculate 10th and 90th percentiles
p10 <- qnorm(0.10, mean = mean_doy, sd = sd_doy)
p90 <- qnorm(0.90, mean = mean_doy, sd = sd_doy)

# Plot
ggplot(df_norm, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = p10, linetype = "dashed", color = "goldenrod", linewidth = 0.8) +
  geom_vline(xintercept = p90, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = p10, y = max(y) * 0.9, label = "10th percentile", hjust = -0.1) +
  annotate("text", x = p90, y = max(y) * 0.9, label = "90th percentile", hjust = 1.1) +
  scale_x_continuous(limits = c(0, 365), breaks = seq(0, 365, by = 30)) +
  theme_bw() +
  labs(x = "Day of Year", y = "Number of Flowers or Bees", title = "Example of PhenoMetrics")

ggsave("/blue/guralnick/millerjared/PhenoInterpolate/outputs/figures/phenometrics-conceptual.png")
