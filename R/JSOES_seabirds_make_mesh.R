### JSOES: Create meshes for use with sdmTMB or TMB SDMs for seabird data

# Author: Markus Min
# Last updated: 2025-09-22

# Description: This script will create the mesh necessary for SPDE (used in packages such as sdmTMB)
# for the June survey for the JSOES seabird data

# We will need the following data to construct the mesh:
# 1. The JSOES Seabird Data (for sampling locations)
# 2. Shapefiles representing the West Coast (for defining the survey domain)
# 3. Files for covariates for SST (this ensures that we can align the mesh with 
# the resolution of the covariates used for projecting density)

# We will then use the above files to construct a mesh

#### Load libraries ####

library(readxl)
library(janitor)
library(tidyverse)
library(here)
library(sf)
library(sdmTMB)
# library(TMB)
# library(rnaturalearth)
library(fmesher)
# library(lubridate)
# library(geosphere)
# library(geos)


#### Data 1: load and reformat trawl data ####
birds <- clean_names(read_excel(here::here("Data", "Bird_Density_qy1_Crosstab_Markus_Min_COMU_SOSH_June_8.7.25.xlsx")))

# temporal cut off - 2021 is the last year included
birds <- subset(birds, year <= 2021)


# change nmi to km to keep units consistent
birds %>% 
  mutate(km_from_shore = distance_nm * 1.852) -> birds

# pivot the data longer - this separates out common mures and sooty shearwaters
birds %>% 
  pivot_longer(., cols = c("common_murre", "sooty_shearwater"), values_to = "n_per_km2", names_to = "species") %>% 
  mutate(n = n_per_km2*area_sq_km) -> birds_long


# Add UTM coordinates to get an equal distance projection and get lat/longs
utm_crs <- get_crs(birds_long, c("dec_long", "dec_lat")) # check that

birds_long <- add_utm_columns(
  birds_long,
  ll_names = c("dec_long", "dec_lat"),
  ll_crs = 4326,
  utm_names = c("Lon.km", "Lat.km"),
  utm_crs = utm_crs,
  units = c("km")
)


# Change Lon.km and Lat.km columns to X and Y to avoid problems later
birds_long %>% 
  dplyr::rename(X = Lon.km, Y = Lat.km) -> birds_long

# some data filtering for Markus's purposes: 
# drop 2022:2025 (adult returns aren't complete for those years yet)
birds_long <- filter(birds_long, !(year %in% c(2022:2025)))

#### Data 3: load and reformat spatial data ####
usa_spdf <- st_read(here::here("Data", "map_files", "USA_adm0.shp"))
# load BC
CAN_spdf <- st_read(here::here("Data", "map_files", "canada", "lpr_000b16a_e.shp"))
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)


#### create base map for visualizing data
survey_area_basemap <- ggplot(US_west_coast) +
  geom_sf() +
  geom_sf(data = BC_coast) +
  coord_sf(ylim = c(44,48.5),  xlim = c(-126, -123)) +
  scale_x_continuous(breaks = c(124,125,126)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

survey_area_basemap_km <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = BC_coast_proj_km) + 
  ylab("Latitude")+
  xlab("Longitude")+
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















#### Data 4: load and reformat covariate data ####

#### Covariate 1: SST
# Please note that currently the SST data only goes through 2023
# This data is the Sea Surface Temperature, NOAA Coral Reef Watch Daily Global 5km 
# Satellite SST (CoralTemp) dataset was downloaded from the [NOAA ERDDAP server](https://coastwatch.noaa.gov/erddap/griddap/noaacrwsstDaily.html)
SST1 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_9bdb_17ce_0aa0.csv"), skip = 1)
SST2 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_36e3_7cb3_1109.csv"), skip = 1)
SST3 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_519c_f955_0fea.csv"), skip = 1)
SST4 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_727c_5c10_e6a0.csv"), skip = 1)
SST5 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_1176_fc5f_53d2.csv"), skip = 1)
SST6 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_ae0c_9197_083f.csv"), skip = 1)
SST7 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_fe4f_9c9c_442e.csv"), skip = 1)

SST1 %>% 
  bind_rows(., SST2, SST3, SST4, SST5, SST6, SST7) %>% 
  dplyr::rename(time = UTC, latitude = degrees_north, longitude = degrees_east, SST = degree_C)  %>% 
  mutate(time = ymd_hms(time)) %>% 
  mutate(day_of_month = mday(time)) %>% 
  mutate(year = year(time)) %>% 
  mutate(julian = yday(time)) -> SST


# Convert SST to sf and then to UTM zone 10
st_as_sf(SST, coords = c("longitude", "latitude"), crs = 4326) -> SST_sf
sf::st_transform(SST_sf, crs = UTM_zone_10_crs) -> SST_sf_proj

# convert to km
SST_sf_proj_km <- SST_sf_proj
SST_sf_proj_km$geometry <- SST_sf_proj$geometry/1000

SST_sf_proj_km <- st_set_crs(SST_sf_proj_km, UTM_zone_10_crs)

# Dealing with covariate effects:
# If we want to include SST as a covariate and then predict based on SST, we'll need to decide
# how we want to do that, because there in situ temperature data from JSOES and satellite
# SST data, and they don't exactly align.
# The simpler way to do this is just to take the remotely sensed SST data and align
# it to the JSOES seabird samples

# For this, birds needs to be in SF format as well
st_as_sf(birds_long, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> birds_long_sf

birds_long_sf_split <- split(birds_long_sf, birds_long_sf$year)
SST_sf_proj_km_split <- split(SST_sf_proj_km, SST_sf_proj_km$year) 

birds_years <- as.character(unique(birds_long_sf$year))

# Find the nearest remotely sensed SST measurement for each JSOES sample
birds_long_remote_SST_list <- map(birds_years, function(year) {
  birds_one_year <- birds_long_sf_split[[year]]
  SST_one_year <- SST_sf_proj_km_split[[year]]
  
  # Find nearest features
  SST_index_for_birds <- st_nearest_feature(birds_one_year, SST_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  birds_one_year %>%
    mutate(matched_id = SST_index_for_birds) %>%
    bind_cols(
      st_drop_geometry(SST_one_year[SST_index_for_birds, "SST"]) %>%
        dplyr::rename(remote_SST = SST)
    ) %>% 
    dplyr::select(-matched_id) -> output
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
birds_long <- list_rbind(birds_long_remote_SST_list)

#### Calculate distance to shore

dplyr::select(birds_long, station, km_from_shore, X, Y) -> birds_long_spatial
st_as_sf(birds_long_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> birds_long_sf
st_distance(birds_long_sf, US_west_coast_proj_km) -> birds_long_dist_shore
birds_long$sf_dist_shore <- as.numeric(birds_long_dist_shore)

# confirm that these are reasonable, based on points from survey
# they're similar - but not identical
# hist(as.numeric(birds_long$sf_dist_shore)-birds_long_spatial$km_from_shore)
# and they're pretty biased - the distances calculated using sf tend to be smaller
# Like for SST, let's just use the sf distance from shore calculated for both 
# prediction grid and the samples themselves

# Rescale data for fitting to covariates
birds_long %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> birds_long



#### Create a survey prediction grid ####
# Create a survey grid from scratch, using SST data

# only keep one year, since that's all we need to make a grid
SST_sf_proj %>% 
  filter(year == 2010) -> SST_sf_proj_2010

st_as_sf(SST_sf_proj_2010$geometry/1000, crs = UTM_zone_10_crs) -> SST_sf_proj_km_2010

SST_sf_proj_km_2010$dist_shore <- as.numeric(st_distance(SST_sf_proj_km_2010, US_west_coast_proj_km))


# KEY DECISION: To what area do you want to project the SDM to?
# Some observations of the data:
# The vast majority of stations are within 60 km of shore, but there are just a couple 
# that are further out in the Columbia River Plume up to 83 km from shore
# a total of 15 tows (out of over 2000) are >= 35 nmi (65 km) offshore; many more are 30 or 31 nmi (55.5 or 57.4 km) offshore.
# We also want to exclude values that are within 0.5 km of shore, because these
# don't have measureable SST values (they're basically on shore)
# this can also be seen in the following histogram:
# hist(jsoes_long$km_from_shore)

# For the N/S dimension, the most northerly transect is typically the FS transect and 
# the most southerly transect is the NH transect.
# this corresponds to 44.67 to 48.23 N
# Let's take these as our N/S boundaries and add 10 km as a buffer.
min_Y <- mean(subset(birds_long, grepl("NH", station))$Y)
max_Y <- mean(subset(birds_long, grepl("FS", station))$Y)

# Based on this, I'm going to propose that we generate a prediction that extends:
# from 0.5 to 65 km offshore
# from 44.67 to 48.23 N (+ 10 km on either side)

# trim longitudinally
SST_sf_proj_km_2010 %>% 
  filter(dist_shore <= 65 & dist_shore >= 0.5) -> grid_within_65km_shore

# trim latitudinally
survey_domain <- st_crop(grid_within_65km_shore, xmin = 0, xmax = 100000, ymin=min_Y-10, ymax=max_Y+10)

# let's inspect the grid that we created
# ggplot(survey_domain) +
#   geom_sf()

# We will need to manually trim out some of these areas

# manually remove sections for strait, hood canal (and other inland waters), Grays Harbor, Willapa Bay, Columbia River estuary
strait <- st_as_sfc(st_bbox(c(xmin=-124.65, xmax=-122, ymin=47.9375, ymax=49), crs = "WGS84"))
strait_proj <- sf::st_transform(strait, crs = UTM_zone_10_crs)
strait_proj_km <- st_as_sf(strait_proj/1000, crs = UTM_zone_10_crs)

inland_waters <- st_as_sfc(st_bbox(c(xmin=-123.5, xmax=-120, ymin=44, ymax=49.5), crs = "WGS84"))
inland_waters_proj <- sf::st_transform(inland_waters, crs = UTM_zone_10_crs)
inland_waters_proj_km <- st_as_sf(inland_waters_proj/1000, crs = UTM_zone_10_crs)

grays_harbor <- st_as_sfc(st_bbox(c(xmin=-124.15, xmax=-123.6, ymin=46.83, ymax=47.09), crs = "WGS84"))
grays_harbor_proj <- sf::st_transform(grays_harbor, crs = UTM_zone_10_crs)
grays_harbor_proj_km <- st_as_sf(grays_harbor_proj/1000, crs = UTM_zone_10_crs)

willapa_bay <- st_as_sfc(st_bbox(c(xmin=-124.06, xmax=-123.5, ymin=46.34, ymax=46.8), crs = "WGS84"))
willapa_bay_proj <- sf::st_transform(willapa_bay, crs = UTM_zone_10_crs)
willapa_bay_proj_km <- st_as_sf(willapa_bay_proj/1000, crs = UTM_zone_10_crs)

estuary <- st_as_sfc(st_bbox(c(xmin=-124.02, xmax=-123.5, ymin=46.12, ymax=46.35), crs = "WGS84"))
estuary_proj <- sf::st_transform(estuary, crs = UTM_zone_10_crs)
estuary_proj_km <- st_as_sf(estuary_proj/1000, crs = UTM_zone_10_crs)

survey_domain %>% 
  st_difference(., strait_proj_km) %>% 
  st_difference(., inland_waters_proj_km) %>% 
  st_difference(., grays_harbor_proj_km) %>% 
  st_difference(., willapa_bay_proj_km) %>% 
  st_difference(., estuary_proj_km) -> survey_domain

# Create a concave hull around points
st_concave_hull(st_union(survey_domain), ratio = 0.1) -> survey_domain_polygon

# ggplot(survey_domain_polygon) +
#   geom_sf()

# ggplot(survey_domain) +
#   geom_sf()

# let's inspect our revised survey domain
survey_area_basemap_km +
  geom_sf(data = survey_domain_polygon, fill = "blue", alpha = 0.2)

# Take SST grid and intersect with survey domain
SST_sf_proj %>% 
  mutate(geometry = geometry/1000) %>% 
  st_set_crs(UTM_zone_10_crs) %>% 
  st_intersection(survey_domain_polygon) %>% 
  dplyr::select(SST, year, geometry)-> SST_survey_domain

survey_area_basemap_km +
  geom_sf(data = SST_survey_domain, fill = "blue", alpha = 0.2)

# add in distance from shore - this is now your survey domain, with covariates
SST_survey_domain %>% 
  mutate(sf_dist_shore = as.numeric(st_distance(geometry, US_west_coast_proj_km))) -> survey_domain_cov

# split by year
# take the output from the SST code above as input here
# survey_domain_cov_split <- split(survey_domain_cov, survey_domain_cov$year)


# # add in salinity, by grabbing the closest measurement from the MOM6 hindcast
# # Find the nearest salinity measurement for each cell in the prediction grid
# survey_domain_cov_remote_sos_list <- map(jsoes_years, function(year) {
#   survey_domain_cov_one_year <- survey_domain_cov_split[[year]]
#   sos_one_year <- jsoes_sos_sf_proj_km_split[[year]]
#   
#   # Find nearest features
#   sos_index_for_survey_domain_cov <- st_nearest_feature(survey_domain_cov_one_year, sos_one_year)
#   
#   # Combine with matched points (optional: include matched geometry or attributes)
#   survey_domain_cov_one_year %>%
#     mutate(matched_id = sos_index_for_survey_domain_cov) %>%
#     bind_cols(
#       st_drop_geometry(sos_one_year[sos_index_for_survey_domain_cov, "sos"]) %>%
#         dplyr::rename(remote_sos = sos)
#     ) %>% 
#     dplyr::select(-matched_id) -> output
#   
#   # transform this to a df
#   st_drop_geometry(output) %>% 
#     bind_cols(st_coordinates(output)) -> output
#   
#   return(output)
# })
# 
# # take the list and turn it back into an sf object
# survey_domain_cov <- list_rbind(survey_domain_cov_remote_sos_list)
# 
# survey_domain_cov <- st_as_sf(survey_domain_cov, coords = c("X", "Y"), crs = UTM_zone_10_crs)



# because the survey_domain_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_cov %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_cov

# %>% 
#   # new choice: truncate it to only 65 km from shore. There are only a handful of 
#   # tows that are greater than 57 km offshore (15 out of over 2000)
#   filter(sf_dist_shore <= 65 & sf_dist_shore >= 0.5) -> survey_domain_cov

# plot it - looks good!
# ggplot(survey_domain_cov) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_cov)) -> survey_domain_cov_coords
survey_predict_grid <- data.frame(X = survey_domain_cov_coords$X,
                                  Y = survey_domain_cov_coords$Y,
                                  SST = survey_domain_cov$SST,
                                  # sal = survey_domain_cov$remote_sos,
                                  year = survey_domain_cov$year,
                                  dist_shore = as.numeric(survey_domain_cov$sf_dist_shore))

# rescale prediction grid
survey_predict_grid %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         # sal_scaled = as.numeric(scale(sal)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid

#### Trim to the years between max and min birds years ####
survey_domain_cov %>% 
  filter(year >= min(birds$year) & year <= max(birds$year)) -> survey_domain_cov_birds_allyears

# because the survey_domain_jsoes_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_cov_birds_allyears %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_cov_birds_allyears

# plot it - looks good!
ggplot(survey_domain_cov_birds_allyears) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_cov_birds_allyears)) -> survey_domain_cov_birds_allyears_coords
survey_predict_grid_birds_allyears <- data.frame(X = survey_domain_cov_birds_allyears_coords$X,
                                                    Y = survey_domain_cov_birds_allyears_coords$Y,
                                                    SST = survey_domain_cov_birds_allyears$SST,
                                                    year = survey_domain_cov_birds_allyears$year,
                                                    dist_shore = as.numeric(survey_domain_cov_birds_allyears$sf_dist_shore))


# rescale prediction grid
survey_predict_grid_birds_allyears %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid_birds_allyears

# Ok - now we are finally ready to create our mesh!

#### Create the meshes for SDM ####

#### Create the trawl mesh

# get the sampling coordinates
birds_long %>% 
  distinct(year, station, .keep_all = TRUE) %>% 
  dplyr::select(-species) -> birds_samples

# In order to get our model to predict throughout the survey domain, we'll need to rbind
# the projection grid and sampling coordinates

# first extract the grid for the projection area
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_cov %>% 
  filter(year == 2010) %>% 
  dplyr::select(geometry) -> survey_domain_cov_grid

# extract only coordinates
projection_coords <- as.data.frame(st_coordinates(survey_domain_cov_grid))
colnames(projection_coords) <- c("X", "Y")

# combine the samples and the projection grid
dplyr::select(birds_samples, X, Y) %>% 
  bind_rows(projection_coords) -> birds_all_points

# From this object, we then can create a mesh that will be used by the SPDE method.

# use a boundary to ensure that our mesh matches our survey grid
# in the boundary, remove the rarely sampled far offshore samples that cause bulges in the mesh
bnd <- INLA::inla.nonconvex.hull(cbind(subset(birds_samples, distance_nm < 45)$X, subset(birds_samples, distance_nm < 45)$Y), convex = -0.1)

# make versions of the mesh using the boundary object at varying resolution

birds_inla_mesh_cutoff10 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(birds_all_points$X, birds_all_points$Y),
  cutoff = 10,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "birds_boundary_mesh_cutoff10.png"), width=4, height=6, res=200, units="in")
plot(birds_inla_mesh_cutoff10)
points(x = birds_samples$X, y = birds_samples$Y, cex = 0.3)
dev.off()

birds_inla_mesh_cutoff15 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(birds_all_points$X, birds_all_points$Y),
  cutoff = 15,
  boundary = bnd
)

# make the same mesh using sdmTMB for compatibility with sdmTMB() function
birds_inla_mesh_cutoff15_sdmTMB <- make_mesh(
  data = birds_all_points,
  xy_cols = c("X", "Y"),
  cutoff = 15,
  boundary = bnd
)
# confirm that they are the same
# plot(inla_mesh_cutoff15_sdmTMB)
# plot(inla_mesh_cutoff15)

png(here::here("two_stage_models", "figures", "birds_boundary_mesh_cutoff15.png"), width=4, height=6, res=200, units="in")
plot(birds_inla_mesh_cutoff15)
points(x = birds_samples$X, y = birds_samples$Y, cex = 0.3)
dev.off()


#### Visualize alignment between mesh and projection grid ####

# transform mesh to sf
fm_as_sfc(birds_inla_mesh_cutoff15) %>% 
   st_set_crs(st_crs(survey_domain_cov_grid)) -> birds_inla_mesh_cutoff15_sf

birds_mesh_plus_points_plot <- survey_area_basemap_km +
  geom_sf(data = birds_inla_mesh_cutoff15_sf, fill = NA) +
  geom_sf(data = survey_domain_cov_grid) +
  geom_point(data = birds_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("two_stage_models", "figures", "birds_mesh_plus_points_plot.png"), birds_mesh_plus_points_plot,  height = 10, width = 6)

