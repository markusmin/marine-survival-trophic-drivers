## 10.1_Fig2_SDM_maps

# Description: This script generates Figure 2, which illustrates the SDM model outputs for the main text.

# This script also generates supplemental figures showing the SDM outputs for
# each taxon across years.

# This script uses the Bongo biomass outputs.

#### Load libraries + model outputs ####

## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)
library(TMB)
library(RColorBrewer)

# source the prep scripts for each of the surveys
source("R/functions/sample_var.R")
source("R/functions/rmvnorm_prec.R")
source(here::here("R", "PRS_PWCC_make_mesh.R"))
source(here::here("R", "hake_survey_make_mesh.R"))
source(here::here("R", "JSOES_seabirds_make_mesh.R"))
source(here::here("R", "JSOES_make_mesh.R"))


# load the SDM outputs + estimated uncertainty from the stage 1 models
# 05.1_prey_SDM
load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "jsoes_bongo_biomass_SDM_output.rda"))
jsoes_bongo_biomass_SDM_Obj <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_Obj
jsoes_bongo_biomass_SDM_Opt <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_Opt
jsoes_bongo_biomass_SDM_report <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_report

load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "estimated_SE_prey_model.rda"))
SE_pianka_o_csyif_jsoes_bongo_biomass_SDM_output_t_prey_model <- SE_prey_model$SE_pianka_o_csyif_jsoes_bongo_biomass_SDM_output_t_prey_model
SE_pianka_o_cssif_jsoes_bongo_biomass_SDM_output_t_prey_model <- SE_prey_model$SE_pianka_o_cssif_jsoes_bongo_biomass_SDM_output_t_prey_model

# 05.1.2_rf_SDM
load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "rf_SDM_output.rda"))
rf_SDM_Obj <- rf_SDM_output$rf_SDM_Obj
rf_SDM_Opt <- rf_SDM_output$rf_SDM_Opt
rf_SDM_report <- rf_SDM_output$rf_SDM_report

load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "estimated_SE_rf_model.rda"))
SE_pianka_o_csyif_rf_SDM_output_t_prey_model <- SE_prey_model$SE_pianka_o_csyif_rf_SDM_output_t_prey_model
SE_pianka_o_cssif_rf_SDM_output_t_prey_model <- SE_prey_model$SE_pianka_o_cssif_rf_SDM_output_t_prey_model

# 05.2_seabird_SDM
load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "seabird_SDM_output.rda"))
seabird_SDM_Obj <- seabird_SDM_output$seabird_SDM_Obj
seabird_SDM_Opt <- seabird_SDM_output$seabird_SDM_Opt
seabird_SDM_report <- seabird_SDM_output$seabird_SDM_report

load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "estimated_SE_seabird_model.rda"))
SE_pianka_o_cssif_sosh_t_seabird_model <- SE_seabird_model$SE_pianka_o_cssif_sosh_t_seabird_model
SE_pianka_o_cssif_comu_t_seabird_model <- SE_seabird_model$SE_pianka_o_cssif_comu_t_seabird_model
SE_pianka_o_csyif_sosh_t_seabird_model <- SE_seabird_model$SE_pianka_o_csyif_sosh_t_seabird_model
SE_pianka_o_csyif_comu_t_seabird_model <- SE_seabird_model$SE_pianka_o_csyif_comu_t_seabird_model

# 05.3_hake_SDM
load(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "hake_SDM_output.rda"))
hake_SDM_Obj <- hake_SDM_output$hake_SDM_Obj
hake_SDM_Opt <- hake_SDM_output$hake_SDM_Opt
hake_SDM_report <- hake_SDM_output$hake_SDM_report

load(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "estimated_SE_hake_model.rda"))
SE_pianka_o_csyif_hake_t_hake_model <- SE_hake_model$SE_pianka_o_csyif_hake_t_hake_model
SE_pianka_o_cssif_hake_t_hake_model <- SE_hake_model$SE_pianka_o_cssif_hake_t_hake_model


#### Estimate uncertainty in density for each taxon ####

### Prey taxa + salmon from the JSOES Bongo biomass model
# dyn.load(dynlib(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "05_1_jsoes_bongo_biomass_SDM")))

# If this has already been run before, just load the output instead of running the sample_var code again
# load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "estimated_SE_density_prey_SDMs.rda"))
# SE_ln_d_gt_cancer_crab_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_cancer_crab_larvae
# SE_ln_d_gt_non_cancer_crab_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_non_cancer_crab_larvae
# SE_ln_d_gt_shrimp_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_shrimp_larvae
# SE_ln_d_gt_hyperiid_amphipods <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_hyperiid_amphipods

# SE_ln_d_gt_cancer_crab_larvae = sample_var( obj=jsoes_bongo_biomass_SDM_Obj, var_name="ln_d_gt_cancer_crab_larvae", mu=jsoes_bongo_biomass_SDM_Obj$env$last.par.best, prec=jsoes_bongo_biomass_SDM_Opt$SD$jointPrecision )
# SE_ln_d_gt_non_cancer_crab_larvae = sample_var( obj=jsoes_bongo_biomass_SDM_Obj, var_name="ln_d_gt_non_cancer_crab_larvae", mu=jsoes_bongo_biomass_SDM_Obj$env$last.par.best, prec=jsoes_bongo_biomass_SDM_Opt$SD$jointPrecision )
# SE_ln_d_gt_shrimp_larvae = sample_var( obj=jsoes_bongo_biomass_SDM_Obj, var_name="ln_d_gt_shrimp_larvae", mu=jsoes_bongo_biomass_SDM_Obj$env$last.par.best, prec=jsoes_bongo_biomass_SDM_Opt$SD$jointPrecision )
# SE_ln_d_gt_hyperiid_amphipods = sample_var( obj=jsoes_bongo_biomass_SDM_Obj, var_name="ln_d_gt_hyperiid_amphipods", mu=jsoes_bongo_biomass_SDM_Obj$env$last.par.best, prec=jsoes_bongo_biomass_SDM_Opt$SD$jointPrecision )

# save all of these
# SE_jsoes_bongo_biomass_SDM <- list(SE_ln_d_gt_cancer_crab_larvae = SE_ln_d_gt_cancer_crab_larvae,
#                           SE_ln_d_gt_non_cancer_crab_larvae = SE_ln_d_gt_non_cancer_crab_larvae,
#                           SE_ln_d_gt_shrimp_larvae = SE_ln_d_gt_shrimp_larvae,
#                           SE_ln_d_gt_hyperiid_amphipods = SE_ln_d_gt_hyperiid_amphipods)
# 
# save(SE_jsoes_bongo_biomass_SDM,
#      file = here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "estimated_SE_density_prey_SDMs.rda"))


### rf from the rf model
# dyn.load(dynlib(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "05_1_2_rf_SDM")))

# If this has already been run before, just load the output instead of running the sample_var code again
# load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "estimated_SE_density_rf_SDMs.rda"))
# SE_rf_SDM <- SE_rf_SDM$SE_ln_d_gt_rf


# SE_ln_d_gt_rf = sample_var( obj=rf_SDM_Obj, var_name="ln_d_gt_rf", mu=rf_SDM_Obj$env$last.par.best, prec=rf_SDM_Opt$SD$jointPrecision )
# save all of these
# SE_rf_SDM <- list(SE_ln_d_gt_rf = SE_ln_d_gt_rf)

# save(SE_rf_SDM,
     # file = here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "estimated_SE_density_rf_SDMs.rda"))




### seabirds from the seabird model
# dyn.load(dynlib(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "05_2_seabird_SDM_v2")))
# SE_ln_d_gt_sosh = sample_var( obj=seabird_SDM_Obj, var_name="ln_d_gt_sosh", mu=seabird_SDM_Obj$env$last.par.best, prec=seabird_SDM_Opt$SD$jointPrecision )
# SE_ln_d_gt_comu = sample_var( obj=seabird_SDM_Obj, var_name="ln_d_gt_comu", mu=seabird_SDM_Obj$env$last.par.best, prec=seabird_SDM_Opt$SD$jointPrecision )
# 
# # save all of these
# SE_seabird_SDM <- list(SE_ln_d_gt_csyif = SE_ln_d_gt_csyif,
#                           SE_ln_d_gt_cssif = SE_ln_d_gt_cssif,
#                           SE_ln_d_gt_sosh = SE_ln_d_gt_sosh,
#                           SE_ln_d_gt_comu = SE_ln_d_gt_comu)
# 
# save(SE_seabird_SDM,
     # file = here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "estimated_SE_density_seabird_salmon_SDMs.rda"))

### Hake from the hake model
# dyn.load(dynlib(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "05_3_hake_SDM_v1")))
# SE_ln_d_gt_hake = sample_var( obj=hake_SDM_Obj, var_name="ln_d_gt_hake", mu=hake_SDM_Obj$env$last.par.best, prec=hake_SDM_Opt$SD$jointPrecision )
# 
# # save all of these
# SE_hake_SDM <- list(SE_ln_d_gt_hake = SE_ln_d_gt_hake)
# 
# save(SE_hake_SDM,
#      file = here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "estimated_SE_hake_SDMs.rda"))


#### load uncertainty previously calculated
load(here::here("R", "05_stage1_SDM", "05.2_seabird_SDM", "estimated_SE_density_seabird_salmon_SDMs.rda"))

SE_ln_d_gt_csyif <- SE_seabird_SDM$SE_ln_d_gt_csyif
SE_ln_d_gt_cssif <- SE_seabird_SDM$SE_ln_d_gt_cssif
SE_ln_d_gt_sosh <- SE_seabird_SDM$SE_ln_d_gt_sosh
SE_ln_d_gt_comu <- SE_seabird_SDM$SE_ln_d_gt_comu

load(here::here("R", "05_stage1_SDM", "05.3_hake_SDM", "estimated_SE_hake_SDMs.rda"))
SE_ln_d_gt_hake <- SE_hake_SDM$SE_ln_d_gt_hake

load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "estimated_SE_density_rf_SDMs.rda"))
SE_ln_d_gt_rf <- SE_rf_SDM$SE_ln_d_gt_rf

load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "estimated_SE_density_prey_SDMs.rda"))

SE_ln_d_gt_cancer_crab_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_cancer_crab_larvae
SE_ln_d_gt_non_cancer_crab_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_non_cancer_crab_larvae
SE_ln_d_gt_shrimp_larvae <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_shrimp_larvae
SE_ln_d_gt_hyperiid_amphipods <- SE_jsoes_bongo_biomass_SDM$SE_ln_d_gt_hyperiid_amphipods

#### Load shapefiles ####

### create a basemap that's tailored to bongo data
usa_spdf <- st_read(here::here("Data", "map_files", "USA_adm0.shp"))
# load BC
CAN_spdf <- st_read(here::here("Data", "map_files", "canada", "lpr_000b16a_e.shp"))
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)

## load USA states and territories
us_states <- st_read(here("map_files", "us_states_territories", "s_05mr24.shp"))

# crop them to our desired area
# US_west_coast <- sf::st_crop(usa_spdf,
#                              c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))

US_west_coast <- sf::st_crop(usa_spdf,
                             # c(xmin = -126, ymin = 40.42, xmax = -120, ymax = 48.5))
                             # c(xmin = -126, ymin = 39.42, xmax = -120, ymax = 48.5))
                             c(xmin = -127, ymin = 39.42, xmax = -123, ymax = 48.5))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -127, ymin = 44, xmax = -123, ymax = 48.5))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

us_states_proj <- sf::st_transform(us_states, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
us_states_proj_km <- st_as_sf(us_states_proj$geometry/1000, crs = UTM_zone_10_crs)

cape_mendocino <- st_as_sf(data.frame(lat = 40.42, lon = -124.4), coords = c("lon", "lat"), crs = 4326)

cape_mendocino <- sf::st_transform(cape_mendocino, crs = UTM_zone_10_crs) 
cape_mendocino <- st_as_sf(cape_mendocino$geometry/1000, crs = UTM_zone_10_crs)

jsoes_map_limits <- st_as_sf(data.frame(lat = c(44.3,48.55), lon = c(-125.7,-123.2)), coords = c("lon", "lat"), crs = 4326)
jsoes_map_limits_m <- sf::st_transform(jsoes_map_limits, crs = UTM_zone_10_crs) 
jsoes_map_limits_km <- st_as_sf(jsoes_map_limits_m$geometry/1000, crs = UTM_zone_10_crs)

jsoes_map_limits_km_df <- as.data.frame(st_coordinates(jsoes_map_limits_km))

expanded_map_limits <- st_as_sf(data.frame(lat = c(39.8,48.55), lon = c(-125.7,-123.2)), coords = c("lon", "lat"), crs = 4326)
expanded_map_limits_m <- sf::st_transform(expanded_map_limits, crs = UTM_zone_10_crs) 
expanded_map_limits_km <- st_as_sf(expanded_map_limits_m$geometry/1000, crs = UTM_zone_10_crs)

expanded_map_limits_km_df <- as.data.frame(st_coordinates(expanded_map_limits_km))

#### Generate basemaps ####

# Create two basemaps: One for only the JSOES/projection grid area and one
# for everything north of Cape Mendocino (for hake and PRS data)

fig2_jsoes_basemap <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = us_states_proj_km, fill = "transparent") + 
  geom_sf(data = BC_coast_proj_km) + 
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(ylim = jsoes_map_limits_km_df$Y,  xlim = jsoes_map_limits_km_df$X, expand = FALSE) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(here::here("figures", "paper_figures", "fig2_jsoes_basemap.png"), fig2_jsoes_basemap,  
       height = 8, width = 4)

fig2_expanded_basemap <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = BC_coast_proj_km) + 
  ylab("Latitude")+
  xlab("Longitude")+
  coord_sf(ylim = expanded_map_limits_km_df$Y,  xlim = expanded_map_limits_km_df$X, expand = FALSE) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(here::here("figures", "paper_figures", "fig2_expanded_basemap.png"), fig2_expanded_basemap,  
       height = 10, width = 4)



#### Function to create df of grid-cell level overlap ####
# load taxa data
# extract Yearling Interior Chinook
csyif <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_fa")
# extract subyearling Interior Chinook
cssif <- subset(jsoes_long, species == "chinook_salmon_subyearling_interior_fa")
# extract seabirds
sosh <- subset(birds_long, species == "sooty_shearwater")
comu <- subset(birds_long, species == "common_murre")
# extract bongo data
jsoes_bongo_biomass_cancer_crab_larvae <- read.csv(here::here("model_inputs", "jsoes_bongo_biomass_cancer_crab_larvae.csv"))

common_years_seabird_prey <- intersect(intersect(unique(csyif$year), # jsoes trawl
                                                 unique(jsoes_bongo_biomass_cancer_crab_larvae$year)), # jsoes bongo
                                       intersect(unique(rf$year), # PRS/PWCC
                                                 unique(sosh$year)))

common_years_prey <- intersect(
  unique(jsoes_bongo_biomass_cancer_crab_larvae$year), # jsoes bongo
  unique(rf$year) # PRS/PWCC
)



ln_d_gt_A = seabird_SDM_report$ln_d_gt_csyif
model_years_A = sort(unique(csyif$year))
ln_d_gt_B = seabird_SDM_report$ln_d_gt_sosh
model_years_B = min(sosh$year):max(sosh$year)

create_overlap_df <- function(ln_d_gt_A, ln_d_gt_B,
                              model_years_A, model_years_B){
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest and add the density data
  SDM_predict_grid_A <- subset(survey_predict_grid, year %in% model_years_A)
  SDM_predict_grid_A$ln_d_gt_A <- as.vector(ln_d_gt_A)
  
  SDM_predict_grid_B <- subset(survey_predict_grid, year %in% model_years_B)
  SDM_predict_grid_B$ln_d_gt_B <- as.vector(ln_d_gt_B)
  
  SDM_predict_grid_A %>% 
    left_join(dplyr::select(SDM_predict_grid_B, c(X, Y, year, ln_d_gt_B)), by = c("X", "Y", "year")) -> SDM_predict_grid
  
  # calculate overlap per cell
  
  # split up the dfs by year
  ln_d_gt_A_split_year <- split(SDM_predict_grid$ln_d_gt_A, SDM_predict_grid$year)
  ln_d_gt_B_split_year <- split(SDM_predict_grid$ln_d_gt_B, SDM_predict_grid$year)
  
  ov_gt_A_B_vector <- vector()
  
  for (i in 1:length(ln_d_gt_A_split_year)){
    # calculate denominator
    d_gt_A <- exp(ln_d_gt_A_split_year[[i]])
    d_gt_A_sq <- (d_gt_A/sum(d_gt_A))^2
    
    d_gt_B <- exp(ln_d_gt_B_split_year[[i]])
    d_gt_B_sq <- (d_gt_B/sum(d_gt_B))^2
    
    ov_gt_denom_A_B <- sqrt(sum(d_gt_A_sq) * sum(d_gt_B_sq))
    
    # calculate numerator
    ov_gt_num_A_B <- d_gt_A/sum(d_gt_A) * d_gt_B/sum(d_gt_B)
    
    # calculate index for each cell
    ov_gt_A_B <- ov_gt_num_A_B/ov_gt_denom_A_B
    
    ov_gt_A_B_vector <- c(ov_gt_A_B_vector, ov_gt_A_B)
  }
  
  SDM_predict_grid$ov_A_B <- ov_gt_A_B_vector
  
  return(SDM_predict_grid)
}


#### Functions to generate figures ####
generate_prediction_maps_jsoes_domain <- function(ln_d_gt, SE_ln_d_gt, max_dens, model_years,
                                                  plotting_years, density_units, SE_density_units){
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid, year %in% model_years)
  
  # add the density and SE of density
  SDM_predict_grid$ln_d_gt <- as.vector(ln_d_gt)
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  # create a new column to facilitate plotting
  
  SDM_predict_grid %>% 
    mutate(d_gt = exp(ln_d_gt)*max_dens/1000) -> SDM_predict_grid
  
  quantile(SDM_predict_grid$d_gt, 0.5) -> quantile_50
  
  SDM_predict_grid %>% 
    mutate(d_gt_plotting = ifelse(d_gt <= quantile_50, NA, d_gt)) -> SDM_predict_grid
  
  # keep only the years you want to plot
  SDM_predict_grid <- subset(SDM_predict_grid, year %in% plotting_years)
  
  ### Visualize predicted density
  
  
  density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
                          name = density_units) +
    facet_wrap(~year, nrow = 3) +
    theme(axis.text = element_blank(),
          legend.position = c(0.96, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'))

  
  ### Visualize SE in predicted density
  
  # Need to figure out how to appropriately transform the standard error of density
  
  SE_density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = SE_density_units) +
    facet_wrap(~year, nrow = 3) +
    theme(axis.text = element_blank(),
          legend.position = c(0.93, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'))
  
  return(list(density_predict_map, SE_density_predict_map))
}


# a new version of the plotting code for the paper figure - designed for comparing two years
generate_prediction_maps_jsoes_domain_for_fig <- function(ln_d_gt, SE_ln_d_gt, max_dens, model_years,
                                                  plotting_years, density_units, SE_density_units, taxon_name){
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid, year %in% model_years)
  
  # add the density and SE of density
  SDM_predict_grid$ln_d_gt <- as.vector(ln_d_gt)
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  # create a new column to facilitate plotting
  
  SDM_predict_grid %>% 
    mutate(d_gt = exp(ln_d_gt)*max_dens/1000) -> SDM_predict_grid
  
  quantile(SDM_predict_grid$d_gt, 0.5) -> quantile_50
  
  SDM_predict_grid %>% 
    mutate(d_gt_plotting = ifelse(d_gt <= quantile_50, NA, d_gt)) -> SDM_predict_grid
  
  # keep only the years you want to plot
  SDM_predict_grid <- subset(SDM_predict_grid, year %in% plotting_years)
  
  ### Visualize predicted density
  # plotting_breaks <- quantile(SDM_predict_grid$d_gt, c(0, 0.5, 0.75, 0.999))
  # names(plotting_breaks) <- NULL
  # plotting_breaks <- round(plotting_breaks, 1)
  
  plotting_breaks <- c(0, round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.25,1), 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.998,1))
  names(plotting_breaks) <- NULL
  
  density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = d_gt),
              width = 7, height = 7) +
    # scale_fill_viridis_c( trans = "sqrt",
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
    #                       name = density_units,
    #                       breaks = plotting_breaks) +
    
    scale_fill_distiller(palette = "BuPu",
                         trans = "sqrt",
                         direction = 1,
                         # trim extreme high values to make spatial variation more visible
                         na.value = "#54278f", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
                         name = density_units,
                         breaks = plotting_breaks) +
    
    facet_wrap(~year, nrow = 2) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8, margin = margin(r = 0.5, unit = "cm")),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    ggtitle(taxon_name)
  
  
  ### Visualize SE in predicted density
  
  # Need to figure out how to appropriately transform the standard error of density
  
  SE_density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = SE_density_units) +
    facet_wrap(~year, nrow = 3) +
    theme(legend.position = c(0.93, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    ggtitle(taxon_name)
  
  return(list(density_predict_map, SE_density_predict_map))
}


generate_overlap_map_jsoes_domain <- function(overlap_df, plotting_years){
  
  # keep only the years you want to plot
  overlap_df_for_plot <- subset(overlap_df, year %in% plotting_years)
  
  ### Visualize predicted density
  
  # plotting_breaks <- quantile(overlap_df_for_plot$ov_A_B, c(0, 0.5, 1))
  # names(plotting_breaks) <- NULL
  # plotting_breaks <- round(plotting_breaks, 3)
  
  plotting_breaks <- c(0, round(max(overlap_df_for_plot$ov_A_B)*0.5,3), 
                       round(max(overlap_df_for_plot$ov_A_B)*0.99,3))
  names(plotting_breaks) <- NULL
  
  overlap_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = overlap_df_for_plot, aes(x = X, y = Y, fill = ov_A_B),
              width = 7, height = 7) +
    # scale_fill_viridis_c(option = "magma",
    #                      breaks = plotting_breaks,
    #                      lim = c(0, round(max(overlap_df_for_plot$ov_A_B),3)),
    #                      name = "Pianka's O") +
    
    scale_fill_distiller(palette = "OrRd",
                         breaks = plotting_breaks,
                         lim = c(0, round(max(overlap_df_for_plot$ov_A_B),3)),
                         direction = 1,
                         name = "Pianka's O") +
    
    # scale_fill_viridis_c( trans = "sqrt",
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
    #                       name = density_units) +
    facet_wrap(~year, nrow = 2) +
    # guides(fill = guide_legend(title = "Pianka's O")) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8, margin = margin(r = 0.5, unit = "cm")),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    ggtitle("Overlap")
  
  
  return(overlap_predict_map)
}


#### Create df of overlap in a year ####

csyif_sosh_overlap_df <- create_overlap_df(ln_d_gt_A = seabird_SDM_report$ln_d_gt_csyif, 
                                           ln_d_gt_B = seabird_SDM_report$ln_d_gt_sosh,
                                           model_years_A = sort(unique(csyif$year)),
                                           model_years_B = min(sosh$year):max(sosh$year))



#### Paper figure ####

# Here, create a figure that shows distributions and estimated overlap in two years: high overlap and low overlap

# check which years would be good to compare
sosh_overlap_df <- data.frame(years = sort(unique(sosh$year)), csyif_sosh_overlap = seabird_SDM_report$pianka_o_csyif_sosh_t)
arrange(sosh_overlap_df, csyif_sosh_overlap)

# lowest overlap year is 2005, highest overlap is 2011

csyif_prediction_maps_fig2 <- generate_prediction_maps_jsoes_domain_for_fig(ln_d_gt = seabird_SDM_report$ln_d_gt_csyif,
                                                                         SE_ln_d_gt = SE_ln_d_gt_csyif,
                                                                         model_years = sort(unique(csyif$year)),
                                                                         plotting_years = c(2005, 2011),
                                                                         max_dens = max(csyif$n_per_km),
                                                                         density_units = "Predicted\ndensity\n(N per km)",
                                                                         SE_density_units = "SE [log-\ndensity]\n(N per km)",
                                                                         taxon_name = "Yearlings")

sosh_prediction_maps_fig2 <- generate_prediction_maps_jsoes_domain_for_fig(ln_d_gt = seabird_SDM_report$ln_d_gt_sosh,
                                                                   SE_ln_d_gt = SE_ln_d_gt_sosh,
                                                                   model_years = min(sosh$year):max(sosh$year),
                                                                   plotting_years = c(2005, 2011),
                                                                   max_dens = max(sosh$n_per_km2),
                                                                   density_units = "Predicted\ndensity\n(N per km^2)",
                                                                   SE_density_units = "SE [log-\ndensity]\n(N per km^2)",
                                                                   taxon_name = "Sooty Shearwaters")

csyif_sosh_overlap_map_fig2 <- generate_overlap_map_jsoes_domain(overlap_df = csyif_sosh_overlap_df, 
                                                                  plotting_years = c(2005, 2011))



fig2_SDMs_overlap <- ggarrange(csyif_prediction_maps_fig2[[1]], 
                               sosh_prediction_maps_fig2[[1]], 
                               csyif_sosh_overlap_map_fig2, 
                               ncol = 3, nrow = 1)
                               # labels = c("(A)", "(B)", "(C)"),
                               # label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                               # hjust = 0, vjust = 0)

ggsave(here::here("figures", "paper_figures", "fig2_SDMs_overlap.png"), fig2_SDMs_overlap,
       height = 8, width = 6)


#### Make a version of the overlap figure for presentation ####

generate_prediction_maps_jsoes_domain_for_pres <- function(ln_d_gt, SE_ln_d_gt, max_dens, model_years,
                                                          plotting_years, density_units, SE_density_units, taxon_name){
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid, year %in% model_years)
  
  # add the density and SE of density
  SDM_predict_grid$ln_d_gt <- as.vector(ln_d_gt)
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  # create a new column to facilitate plotting
  
  SDM_predict_grid %>% 
    mutate(d_gt = exp(ln_d_gt)*max_dens/1000) -> SDM_predict_grid
  
  quantile(SDM_predict_grid$d_gt, 0.5) -> quantile_50
  
  SDM_predict_grid %>% 
    mutate(d_gt_plotting = ifelse(d_gt <= quantile_50, NA, d_gt)) -> SDM_predict_grid
  
  # keep only the years you want to plot
  SDM_predict_grid <- subset(SDM_predict_grid, year %in% plotting_years)
  
  ### Visualize predicted density
  # plotting_breaks <- quantile(SDM_predict_grid$d_gt, c(0, 0.5, 0.75, 0.999))
  # names(plotting_breaks) <- NULL
  # plotting_breaks <- round(plotting_breaks, 1)
  
  plotting_breaks <- c(0, round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.25,1), 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.998,1))
  names(plotting_breaks) <- NULL
  
  density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = d_gt),
              width = 7, height = 7) +
    # scale_fill_viridis_c( trans = "sqrt",
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
    #                       name = density_units,
    #                       breaks = plotting_breaks) +
    
    scale_fill_distiller(palette = "BuPu",
                         trans = "sqrt",
                         direction = 1,
                         # trim extreme high values to make spatial variation more visible
                         na.value = "#54278f", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
                         name = density_units,
                         breaks = plotting_breaks) +
    
    facet_wrap(~year, nrow = 1) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 10, margin = margin(r = 0.5, b = 0.5, unit = "cm")),
          legend.text = element_text(size = 8),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5))
  
  
  ### Visualize SE in predicted density
  
  # Need to figure out how to appropriately transform the standard error of density
  
  SE_density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = SE_density_units) +
    facet_wrap(~year, nrow = 3) +
    theme(legend.position = "right",
      # legend.position = c(0.93, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          plot.title = element_text(size = 15, hjust = 0.5))
  
  return(list(density_predict_map, SE_density_predict_map))
}


generate_overlap_map_jsoes_domain_for_pres <- function(overlap_df, plotting_years){
  
  # keep only the years you want to plot
  overlap_df_for_plot <- subset(overlap_df, year %in% plotting_years)
  
  ### Visualize predicted density
  
  # plotting_breaks <- quantile(overlap_df_for_plot$ov_A_B, c(0, 0.5, 1))
  # names(plotting_breaks) <- NULL
  # plotting_breaks <- round(plotting_breaks, 3)
  
  plotting_breaks <- c(0, round(max(overlap_df_for_plot$ov_A_B)*0.5,3), 
                       round(max(overlap_df_for_plot$ov_A_B)*0.99,3))
  names(plotting_breaks) <- NULL
  
  overlap_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = overlap_df_for_plot, aes(x = X, y = Y, fill = ov_A_B),
              width = 7, height = 7) +
    # scale_fill_viridis_c(option = "magma",
    #                      breaks = plotting_breaks,
    #                      lim = c(0, round(max(overlap_df_for_plot$ov_A_B),3)),
    #                      name = "Pianka's O") +
    
    scale_fill_distiller(palette = "OrRd",
                         breaks = plotting_breaks,
                         lim = c(0, round(max(overlap_df_for_plot$ov_A_B),3)),
                         direction = 1,
                         name = "Pianka's O") +
    
    # scale_fill_viridis_c( trans = "sqrt",
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
    #                       name = density_units) +
    facet_wrap(~year, nrow = 1) +
    # guides(fill = guide_legend(title = "Pianka's O")) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 10, margin = margin(r = 0.5, b = 0.5, unit = "cm")),
          legend.text = element_text(size = 8),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5))
  
  
  return(overlap_predict_map)
}


csyif_prediction_maps_for_pres <- generate_prediction_maps_jsoes_domain_for_pres(ln_d_gt = seabird_SDM_report$ln_d_gt_csyif,
                                                                            SE_ln_d_gt = SE_ln_d_gt_csyif,
                                                                            model_years = sort(unique(csyif$year)),
                                                                            plotting_years = 2005:2011,
                                                                            max_dens = max(csyif$n_per_km),
                                                                            density_units = "Predicted\ndensity\n(N per km)",
                                                                            SE_density_units = "SE [log-\ndensity]\n(N per km)",
                                                                            taxon_name = "Yearlings")

sosh_prediction_maps_for_pres <- generate_prediction_maps_jsoes_domain_for_pres(ln_d_gt = seabird_SDM_report$ln_d_gt_sosh,
                                                                           SE_ln_d_gt = SE_ln_d_gt_sosh,
                                                                           model_years = min(sosh$year):max(sosh$year),
                                                                           plotting_years = 2005:2011,
                                                                           max_dens = max(sosh$n_per_km2),
                                                                           density_units = "Predicted\ndensity\n(N per\nkm^2)",
                                                                           SE_density_units = "SE [log-\ndensity]\n(N per\nkm^2)",
                                                                           taxon_name = "Sooty Shearwaters")

csyif_sosh_overlap_map_for_pres <- generate_overlap_map_jsoes_domain_for_pres(overlap_df = csyif_sosh_overlap_df, 
                                                                 plotting_years = 2005:2011)



for_pres_SDMs_overlap <- ggarrange(csyif_prediction_maps_for_pres[[1]], 
                               sosh_prediction_maps_for_pres[[1]], 
                               csyif_sosh_overlap_map_for_pres, 
                               ncol = 1, nrow = 3)
# labels = c("(A)", "(B)", "(C)"),
# label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
# hjust = 0, vjust = 0)


ggsave(here::here("figures", "presentation_figures", "SDMs_overlap_for_pres.png"), for_pres_SDMs_overlap,
       height = 10, width = 14)








#### Supplemental figures ####

# Define function to generate figures for supplements
generate_prediction_maps_jsoes_domain_for_supp <- function(ln_d_gt, SE_ln_d_gt, max_dens, model_years,
                                                          density_units, SE_density_units, taxon_name){
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid, year %in% model_years)
  
  # add the density and SE of density
  SDM_predict_grid$ln_d_gt <- as.vector(ln_d_gt)
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  # create a new column to facilitate plotting
  
  SDM_predict_grid %>% 
    mutate(d_gt = exp(ln_d_gt)*max_dens/1000) -> SDM_predict_grid
  
  quantile(SDM_predict_grid$d_gt, 0.5) -> quantile_50
  
  SDM_predict_grid %>% 
    mutate(d_gt_plotting = ifelse(d_gt <= quantile_50, NA, d_gt)) -> SDM_predict_grid
  
  
  ### Visualize predicted density
  # plotting_breaks <- quantile(SDM_predict_grid$d_gt, c(0, 0.5, 0.75, 0.999))
  # names(plotting_breaks) <- NULL
  # plotting_breaks <- round(plotting_breaks, 1)
  
  plotting_breaks <- c(0, 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.1,1), 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.25,1), 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.6,1), 
                       round(quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)*0.998,1))
  names(plotting_breaks) <- NULL
  
  density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = d_gt),
              width = 7, height = 7) +
    # scale_fill_viridis_c( trans = "sqrt",
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
    #                       name = density_units,
    #                       breaks = plotting_breaks) +
    
    scale_fill_distiller(palette = "BuPu",
                         trans = "sqrt",
                         direction = 1,
                         # trim extreme high values to make spatial variation more visible
                         na.value = "#54278f", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt)*max_dens/1000, 0.999)),
                         name = density_units,
                         breaks = plotting_breaks) +
    
    facet_wrap(~year, nrow = 3) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(size = 8, margin = margin(r = 0.5, unit = "cm")),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    ggtitle(taxon_name)
  
  
  ### Visualize SE in predicted density
  
  # Need to figure out how to appropriately transform the standard error of density
  
  SE_density_predict_map <- fig2_jsoes_basemap +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = SE_density_units) +
    facet_wrap(~year, nrow = 3) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(size = 8, margin = margin(r = 0.5, unit = "cm")),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.15, 'cm'),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5)) +
    ggtitle(taxon_name)
  
  return(list(density_predict_map, SE_density_predict_map))
}



### CSYIF
csyif_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = seabird_SDM_report$ln_d_gt_csyif,
                                                               SE_ln_d_gt = SE_ln_d_gt_csyif,
                                                               model_years = sort(unique(csyif$year)),
                                                               max_dens = max(csyif$n_per_km),
                                                               density_units = "Predicted\ndensity\n(N per km)",
                                                               SE_density_units = "SE [log-\ndensity]\n(N per km)",
                                                               taxon_name = "Yearling Interior Fall Chinook")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "csyif_predict_map.png"), csyif_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "csyif_predict_SE_map.png"), csyif_prediction_maps[[2]],  
       height = 8, width = 10)

### CSSIF
cssif_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = seabird_SDM_report$ln_d_gt_cssif,
                                                                        SE_ln_d_gt = SE_ln_d_gt_cssif,
                                                                        model_years = sort(unique(cssif$year)),
                                                                        max_dens = max(cssif$n_per_km),
                                                                        density_units = "Predicted\ndensity\n(N per km)",
                                                                        SE_density_units = "SE [log-\ndensity]\n(N per km)",
                                                                        taxon_name = "Subyearling Interior Fall Chinook")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "cssif_predict_map.png"), cssif_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "cssif_predict_SE_map.png"), cssif_prediction_maps[[2]],  
       height = 8, width = 10)


### SOSH
sosh_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = seabird_SDM_report$ln_d_gt_sosh,
                                                              SE_ln_d_gt = SE_ln_d_gt_sosh,
                                                              taxon_name = "Sooty Shearwaters",
                                                              model_years = min(sosh$year):max(sosh$year),
                                                              max_dens = max(sosh$n_per_km2),
                                                              density_units = "Predicted\ndensity\n(N per km^2)",
                                                              SE_density_units = "SE [log-\ndensity]\n(N per km^2)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "sosh_predict_map.png"), sosh_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "sosh_predict_SE_map.png"), sosh_prediction_maps[[2]],  
       height = 8, width = 10)

### COMU
comu_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = seabird_SDM_report$ln_d_gt_comu,
                                                                       SE_ln_d_gt = SE_ln_d_gt_comu,
                                                                       taxon_name = "Common Murres",
                                                                       model_years = unique(comu$year),
                                                                       max_dens = max(comu$n_per_km2),
                                                                       density_units = "Predicted\ndensity\n(N per km^2)",
                                                                       SE_density_units = "SE [log-\ndensity]\n(N per km^2)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "comu_predict_map.png"), comu_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "comu_predict_SE_map.png"), comu_prediction_maps[[2]],  
       height = 8, width = 10)

### Hake
hake_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = hake_SDM_report$ln_d_gt_hake,
                                                                       SE_ln_d_gt = SE_ln_d_gt_hake,
                                                                       taxon_name = "Pacific Hake",
                                                                       model_years = min(hake$year):max(hake$year),
                                                                       max_dens = max(hake$NASC),
                                                                       density_units = "Predicted\ndensity\n(NASC)",
                                                                       SE_density_units = "SE [log-\ndensity]\n(NASC)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "hake_predict_map.png"), hake_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "hake_predict_SE_map.png"), hake_prediction_maps[[2]],  
       height = 8, width = 10)

### Non-Cancer Crab Larvae
non_cancer_crab_larvae_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = jsoes_bongo_biomass_SDM_report$ln_d_gt_non_cancer_crab_larvae,
                                                                                         SE_ln_d_gt = SE_ln_d_gt_non_cancer_crab_larvae,
                                                                                         taxon_name = "Non-Cancer Crab Larvae",
                                                                                         model_years = min(jsoes_bongo_biomass_non_cancer_crab_larvae$year):max(jsoes_bongo_biomass_non_cancer_crab_larvae$year),
                                                                                         max_dens = max(jsoes_bongo_biomass_non_cancer_crab_larvae$total_sum_of_final_carbon_mg_m3),
                                                                                         density_units = "Predicted\ncarbon\n(mg/m^3)",
                                                                                         SE_density_units = "SE [log-\ncarbon]\n(mg/m^3)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "non_cancer_crab_larvae_predict_map.png"), non_cancer_crab_larvae_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "non_cancer_crab_larvae_predict_SE_map.png"), non_cancer_crab_larvae_prediction_maps[[2]],  
       height = 8, width = 10)

### Cancer Crab Larvae
cancer_crab_larvae_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = jsoes_bongo_biomass_SDM_report$ln_d_gt_cancer_crab_larvae,
                                                                                     SE_ln_d_gt = SE_ln_d_gt_cancer_crab_larvae,
                                                                                     taxon_name = "Cancer Crab Larvae",
                                                                                     model_years = min(jsoes_bongo_biomass_cancer_crab_larvae$year):max(jsoes_bongo_biomass_cancer_crab_larvae$year),
                                                                                     max_dens = max(jsoes_bongo_biomass_cancer_crab_larvae$total_sum_of_final_carbon_mg_m3),
                                                                                     density_units = "Predicted\ncarbon\n(mg/m^3)",
                                                                                     SE_density_units = "SE [log-\ncarbon]\n(mg/m^3)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "cancer_crab_larvae_predict_map.png"), cancer_crab_larvae_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "cancer_crab_larvae_predict_SE_map.png"), cancer_crab_larvae_prediction_maps[[2]],  
       height = 8, width = 10)

### Shrimp Larvae
shrimp_larvae_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = jsoes_bongo_biomass_SDM_report$ln_d_gt_shrimp_larvae,
                                                                                SE_ln_d_gt = SE_ln_d_gt_shrimp_larvae,
                                                                                taxon_name = "Shrimp Larvae",
                                                                                model_years = min(jsoes_bongo_biomass_shrimp_larvae$year):max(jsoes_bongo_biomass_shrimp_larvae$year),
                                                                                max_dens = max(jsoes_bongo_biomass_shrimp_larvae$total_sum_of_final_carbon_mg_m3),
                                                                                density_units = "Predicted\ncarbon\n(mg/m^3)",
                                                                                SE_density_units = "SE [log-\ncarbon]\n(mg/m^3)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "shrimp_larvae_predict_map.png"), shrimp_larvae_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "shrimp_larvae_predict_SE_map.png"), shrimp_larvae_prediction_maps[[2]],  
       height = 8, width = 10)

### Hyperiid Amphipods
hyperiid_amphipod_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = jsoes_bongo_biomass_SDM_report$ln_d_gt_hyperiid_amphipod,
                                                                                    SE_ln_d_gt = SE_ln_d_gt_hyperiid_amphipods,
                                                                                    taxon_name = "Hyperiid Amphipods",
                                                                                    model_years = min(jsoes_bongo_biomass_hyperiid_amphipods$year):max(jsoes_bongo_biomass_hyperiid_amphipods$year),
                                                                                    max_dens = max(jsoes_bongo_biomass_hyperiid_amphipods$total_sum_of_final_carbon_mg_m3),
                                                                                    density_units = "Predicted\ncarbon\n(mg/m^3)",
                                                                                    SE_density_units = "SE [log-\ncarbon]\n(mg/m^3)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "hyperiid_amphipod_predict_map.png"), hyperiid_amphipod_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "hyperiid_amphipod_predict_SE_map.png"), hyperiid_amphipod_prediction_maps[[2]],  
       height = 8, width = 10)

### YOY Rockfishes
rf_prediction_maps <- generate_prediction_maps_jsoes_domain_for_supp(ln_d_gt = rf_SDM_report$ln_d_gt_rf,
                                                                     SE_ln_d_gt = SE_ln_d_gt_rf,
                                                                     taxon_name = "YOY Rockfishes",
                                                                     model_years = unique(rf$year),
                                                                     max_dens = max(rf$total),
                                                                     density_units = "Predicted\ndensity\n(N/km)",
                                                                     SE_density_units = "SE [log-\ndensity]\n(N/km)")


ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "rf_predict_map.png"), rf_prediction_maps[[1]],  
       height = 8, width = 10)

ggsave(here::here("figures", "paper_figures", "biomass_models", "SDM_predict_maps", "rf_predict_SE_map.png"), rf_prediction_maps[[2]],  
       height = 8, width = 10)