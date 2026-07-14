### JSOES: Create meshes for use with sdmTMB or TMB SDMs

# Author: Markus Min
# Last updated: 2025-07-23

# Description: This script will create the mesh necessary for SPDE (used in packages such as sdmTMB)
# for the June survey for both the JSOES Bongo and JSOES Midwater Trawl Data

# We will need the following data to construct the mesh:
# 1. The JSOES Trawl Data (for sampling locations)
# 2. The JSOES Bongo Data (for sampling locations)
# 3. Shapefiles representing the West Coast (for defining the survey domain)
# 4. Files for covariates for SST (this ensures that we can align the mesh with 
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
jsoes <- clean_names(read_excel(here::here("Data", "Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.15.25_FINAL.xlsx")))

# compare May and June temporal coverage
table(subset(jsoes, month == "May")$year)
table(subset(jsoes, month == "June")$year)
# only 2013, 2014, and 2020 are missing for May

# Keep only June data
jsoes <- filter(jsoes, month == "June")

# extract species column names
species_column_names <- colnames(jsoes)[32:ncol(jsoes)]

# change nmi to km to keep units consistent
jsoes %>% 
  mutate(km_from_shore = nmi_from_shore * 1.852) -> jsoes

# there are two tows that are missing lat/lon info, but we have the station code.
# use the lat/lon info from another two conducted at that station in another year
# to population the lat/lon field
jsoes %>% 
  filter(station == "CR35") %>% 
  summarize(CR35_mean_lat = mean(mid_lat, na.rm = T),
            CR35_mean_long = mean(mid_long, na.rm = T)) -> CR35_coords

jsoes %>% 
  filter(station == "QR14") %>% 
  summarize(QR14_mean_lat = mean(mid_lat, na.rm = T),
            QR14_mean_long = mean(mid_long, na.rm = T)) -> QR14_coords

jsoes %>% 
  mutate(mid_lat = ifelse(is.na(mid_lat) & station == "QR14", QR14_coords$QR14_mean_lat,
                          ifelse(is.na(mid_lat) & station == "CR35", CR35_coords$CR35_mean_lat, mid_lat))) %>% 
  mutate(mid_long = ifelse(is.na(mid_long) & station == "QR14", QR14_coords$QR14_mean_long,
                           ifelse(is.na(mid_long) & station == "CR35", CR35_coords$CR35_mean_long, mid_long))) -> jsoes

# pivot the data longer
jsoes %>% 
  pivot_longer(., cols = all_of(species_column_names), values_to = "n_per_km", names_to = "species") %>% 
  mutate(n = n_per_km*trawl_dist_km) -> jsoes_long


# Add UTM coordinates to get an equal distance projection and get lat/longs
utm_crs <- get_crs(jsoes_long, c("mid_long", "mid_lat")) # check that

jsoes_long <- add_utm_columns(
  jsoes_long,
  ll_names = c("mid_long", "mid_lat"),
  ll_crs = 4326,
  utm_names = c("Lon.km", "Lat.km"),
  utm_crs = utm_crs,
  units = c("km")
)


# Change Lon.km and Lat.km columns to X and Y to avoid problems later
jsoes_long %>% 
  dplyr::rename(X = Lon.km, Y = Lat.km) -> jsoes_long

# some data filtering for Markus's purposes: 
# drop 1998 (bongo doesn't have this year) and 2022:2025 (adult returns aren't complete for those years yet)
jsoes_long <- filter(jsoes_long, !(year %in% c(1998, 2022:2025)))

#### Data 2: load and reformat bongo data ####
bongo <- clean_names(read_excel(here::here("Data", "Markus_Min_Plankton Density for Trophic Summary Query_10.11.24.xlsx")))

# there's one sample that's missing spatial data: 20 June   2004 062904NH30
# filter(bongo, (is.na(dec_long)))

# use the spatial information for the same tow from the trawl data to fix this
bongo %>% 
  mutate(dec_lat = ifelse(station_code == "062904NH30", subset(jsoes, station_code == "062904NH30")$mid_lat, dec_lat)) %>% 
  mutate(dec_long = ifelse(station_code == "062904NH30", subset(jsoes, station_code == "062904NH30")$mid_long, dec_long)) -> bongo

# turn this into an sf object
st_as_sf(bongo, coords = c("dec_long", "dec_lat"), crs = 4326) -> bongo_sf

# change CRS to UTM zone 10 (to work in meters)
UTM_zone_10_crs <- 32610
bongo_sf_proj <- sf::st_transform(bongo_sf, crs = UTM_zone_10_crs)

# make this projection into kilometers to help with interpretability
bongo_sf_proj_km <- st_as_sf(bongo_sf_proj$geometry/1000, crs = UTM_zone_10_crs)

# extract geometry
as.data.frame(st_coordinates(bongo_sf_proj_km)) -> bongo_km
# add this back to bongo (X and Y now represent eastings and northings in km)
bind_cols(bongo, bongo_km) -> bongo


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
# it to the JSOES trawl samples

# For this, jsoes needs to be in SF format as well
st_as_sf(jsoes_long, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_long_sf

jsoes_long_sf_split <- split(jsoes_long_sf, jsoes_long_sf$year)
SST_sf_proj_km_split <- split(SST_sf_proj_km, SST_sf_proj_km$year) 

jsoes_years <- as.character(unique(jsoes_long_sf$year))

# Find the nearest remotely sensed SST measurement for each JSOES sample
jsoes_long_remote_SST_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_long_sf_split[[year]]
  SST_one_year <- SST_sf_proj_km_split[[year]]
  
  # Find nearest features
  SST_index_for_jsoes <- st_nearest_feature(jsoes_one_year, SST_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  jsoes_one_year %>%
    mutate(matched_id = SST_index_for_jsoes) %>%
    bind_cols(
      st_drop_geometry(SST_one_year[SST_index_for_jsoes, "SST"]) %>%
        dplyr::rename(remote_SST = SST)
    ) %>% 
    dplyr::select(-matched_id) -> output
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
jsoes_long <- list_rbind(jsoes_long_remote_SST_list)

#### Calculate distance to shore

dplyr::select(jsoes_long, station, km_from_shore, X, Y) -> jsoes_long_spatial
st_as_sf(jsoes_long_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_long_sf
st_distance(jsoes_long_sf, US_west_coast_proj_km) -> jsoes_long_dist_shore
jsoes_long$sf_dist_shore <- as.numeric(jsoes_long_dist_shore)

# confirm that these are reasonable, based on points from survey
# they're similar - but not identical
# hist(as.numeric(jsoes_long$sf_dist_shore)-jsoes_long_spatial$km_from_shore)
# and they're pretty biased - the distances calculated using sf tend to be smaller
# Like for SST, let's just use the sf distance from shore calculated for both 
# prediction grid and the samples themselves

# Rescale data for fitting to covariates
jsoes_long %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> jsoes_long

#### Load salinity data and match to the samples and the prediction grid

# take the output from the SST code above as input here
st_as_sf(jsoes_long, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_long_sf
jsoes_long_sf_split <- split(jsoes_long_sf, jsoes_long_sf$year)

jsoes_sos <- read.csv(here::here("model_inputs", "MOM6_sos_june_survey_domain.csv"))
jsoes_sos <- subset(jsoes_sos, !(is.na(sos)))
# get rid of NAs (these are on shore)
# visualize the missing data
# subset(jsoes_sos, is.na(sos)) -> missing_sos
# 
# ggplot(missing_sos, aes(x = lon, y = lat)) +
#   geom_point() +
#   facet_wrap(~year)



# reformat longitude data
jsoes_sos$lon <- jsoes_sos$lon - 360

# Convert salinity to sf and then to UTM zone 10
st_as_sf(jsoes_sos, coords = c("lon", "lat"), crs = 4326) -> jsoes_sos_sf
sf::st_transform(jsoes_sos_sf, crs = UTM_zone_10_crs) -> jsoes_sos_sf_proj

# convert to km
jsoes_sos_sf_proj_km <- jsoes_sos_sf_proj
jsoes_sos_sf_proj_km$geometry <- jsoes_sos_sf_proj$geometry/1000

jsoes_sos_sf_proj_km <- st_set_crs(jsoes_sos_sf_proj_km, UTM_zone_10_crs)

# visualize the data
ggplot(jsoes_sos_sf_proj_km, aes(color = sos)) + 
  geom_sf() +
  facet_wrap(~year)

jsoes_sos_sf_proj_km_split <- split(jsoes_sos_sf_proj_km, jsoes_sos_sf_proj_km$year) 

# Find the nearest salinity measurement for each JSOES sample
jsoes_long_remote_sos_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_long_sf_split[[year]]
  sos_one_year <- jsoes_sos_sf_proj_km_split[[year]]
  
  # Find nearest features
  sos_index_for_jsoes <- st_nearest_feature(jsoes_one_year, sos_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  jsoes_one_year %>%
    mutate(matched_id = sos_index_for_jsoes) %>%
    bind_cols(
      st_drop_geometry(sos_one_year[sos_index_for_jsoes, "sos"]) %>%
        dplyr::rename(remote_sos = sos)
    ) %>% 
    dplyr::select(-matched_id) -> output
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
jsoes_long <- list_rbind(jsoes_long_remote_sos_list)

# compare remote salinity with in situ salinity

ggplot(data = jsoes_long, aes(x = x3m_sal, y = remote_sos)) +
         geom_point()
cor(subset(jsoes_long, !(is.na(x3m_sal)))$x3m_sal, subset(jsoes_long, !(is.na(x3m_sal)))$remote_sos)
# it's okay but not great


#### Create a survey prediction grid ####
# Create a survey grid from scratch, using SST data

# only keep one year, since that's all we need to make a grid
SST_sf_proj %>% 
  filter(year == 2000) -> SST_sf_proj_2000

st_as_sf(SST_sf_proj_2000$geometry/1000, crs = UTM_zone_10_crs) -> SST_sf_proj_km_2000

SST_sf_proj_km_2000$dist_shore <- as.numeric(st_distance(SST_sf_proj_km_2000, US_west_coast_proj_km))


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
min_Y <- mean(subset(jsoes_long, grepl("NH", station))$Y)
max_Y <- mean(subset(jsoes_long, grepl("FS", station))$Y)

# Based on this, I'm going to propose that we generate a prediction that extends:
# from 0.5 to 65 km offshore
# from 44.67 to 48.23 N (+ 10 km on either side)

# trim longitudinally
SST_sf_proj_km_2000 %>% 
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
survey_domain_cov_split <- split(survey_domain_cov, survey_domain_cov$year)


# add in salinity, by grabbing the closest measurement from the MOM6 hindcast
# Find the nearest salinity measurement for each cell in the prediction grid
survey_domain_cov_remote_sos_list <- map(jsoes_years, function(year) {
  survey_domain_cov_one_year <- survey_domain_cov_split[[year]]
  sos_one_year <- jsoes_sos_sf_proj_km_split[[year]]
  
  # Find nearest features
  sos_index_for_survey_domain_cov <- st_nearest_feature(survey_domain_cov_one_year, sos_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  survey_domain_cov_one_year %>%
    mutate(matched_id = sos_index_for_survey_domain_cov) %>%
    bind_cols(
      st_drop_geometry(sos_one_year[sos_index_for_survey_domain_cov, "sos"]) %>%
        dplyr::rename(remote_sos = sos)
    ) %>% 
    dplyr::select(-matched_id) -> output
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
survey_domain_cov <- list_rbind(survey_domain_cov_remote_sos_list)

survey_domain_cov <- st_as_sf(survey_domain_cov, coords = c("X", "Y"), crs = UTM_zone_10_crs)



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
                                  sal = survey_domain_cov$remote_sos,
                                  year = survey_domain_cov$year,
                                  dist_shore = as.numeric(survey_domain_cov$sf_dist_shore))

# rescale prediction grid
survey_predict_grid %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         sal_scaled = as.numeric(scale(sal)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid

# Ok - now we are finally ready to create our mesh!

#### Create the meshes for SDM ####

#### Create the trawl mesh

# get the sampling coordinates
jsoes_long %>% 
  distinct(year, station_code, .keep_all = TRUE) %>% 
  dplyr::select(-species) -> jsoes_samples

# In order to get our model to predict throughout the survey domain, we'll need to rbind
# the projection grid and sampling coordinates

# first extract the grid for the projection area
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_cov %>% 
  filter(year == 2000) %>% 
  dplyr::select(geometry) -> survey_domain_cov_grid

# extract only coordinates
projection_coords <- as.data.frame(st_coordinates(survey_domain_cov_grid))
colnames(projection_coords) <- c("X", "Y")

# combine the samples and the projection grid
dplyr::select(jsoes_samples, X, Y) %>% 
  bind_rows(projection_coords) -> jsoes_all_points

# From this object, we then can create a mesh that will be used by the SPDE method.

# use a boundary to ensure that our mesh matches our survey grid
# in the boundary, remove the rarely sampled far offshore samples that cause bulges in the mesh
bnd <- INLA::inla.nonconvex.hull(cbind(subset(jsoes_samples, nmi_from_shore < 45)$X, subset(jsoes_samples, nmi_from_shore < 45)$Y), convex = -0.1)

# make versions of the mesh using the boundary object at varying resolution

inla_mesh_cutoff5 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_all_points$X, jsoes_all_points$Y),
  cutoff = 5,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "trawl_boundary_mesh_cutoff5.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff5)
points(x = jsoes_samples$X, y = jsoes_samples$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff10 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_all_points$X, jsoes_all_points$Y),
  cutoff = 10,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "trawl_boundary_mesh_cutoff10.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff10)
points(x = jsoes_samples$X, y = jsoes_samples$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff15 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_all_points$X, jsoes_all_points$Y),
  cutoff = 15,
  boundary = bnd
)

# make the same mesh using sdmTMB for compatibility with sdmTMB() function
inla_mesh_cutoff15_sdmTMB <- make_mesh(
  data = jsoes_all_points,
  xy_cols = c("X", "Y"),
  cutoff = 15,
  boundary = bnd
)
# confirm that they are the same
# plot(inla_mesh_cutoff15_sdmTMB)
# plot(inla_mesh_cutoff15)

png(here::here("two_stage_models", "figures", "trawl_boundary_mesh_cutoff15.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff15)
points(x = jsoes_samples$X, y = jsoes_samples$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff20 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_all_points$X, jsoes_all_points$Y),
  cutoff = 20,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "trawl_boundary_mesh_cutoff20.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff20)
points(x = jsoes_samples$X, y = jsoes_samples$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff25 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_all_points$X, jsoes_all_points$Y),
  cutoff = 25,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "trawl_boundary_mesh_cutoff25.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff25)
points(x = jsoes_samples$X, y = jsoes_samples$Y, cex = 0.3)
dev.off()



#### Create the bongo mesh

# get the sampling coordinates
bongo %>% 
  distinct(year, station, .keep_all = TRUE) %>% 
  dplyr::select(year, station, X, Y) -> bongo_samples

# In order to get our model to predict throughout the survey domain, we'll need to rbind
# the projection grid and sampling coordinates

# first extract the grid for the projection area
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_cov %>% 
  filter(year == 2000) %>% 
  dplyr::select(geometry) -> survey_domain_cov_grid

# extract only coordinates
projection_coords <- as.data.frame(st_coordinates(survey_domain_cov_grid))
colnames(projection_coords) <- c("X", "Y")

# combine the samples and the projection grid
dplyr::select(bongo_samples, X, Y) %>% 
  bind_rows(projection_coords) -> bongo_all_points

# From this object, we then can create a mesh that will be used by the SPDE method.

# mesh with a boundary 
# bnd <- INLA::inla.nonconvex.hull(cbind(bongo_samples$X, bongo_samples$Y), convex = -0.1)
# instead of using the bongo samples to generate the boundary (as above), use 
# the trawl samples so that we have the same projection area. Note 
# that the mesh will be slightly different because the location of bongo
# samples differs from the midwater trawl samples
# remove the rarely sampled far offshore samples that cause bulges in the mesh
bnd <- INLA::inla.nonconvex.hull(cbind(subset(jsoes_samples, nmi_from_shore < 45)$X, subset(jsoes_samples, nmi_from_shore < 45)$Y), convex = -0.1)

# make versions of this mesh at varying resolution

bongo_inla_mesh_cutoff5 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(bongo_all_points$X, bongo_all_points$Y),
  cutoff = 5,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "bongo_boundary_mesh_cutoff5.png"), width=4, height=6, res=200, units="in")
plot(bongo_inla_mesh_cutoff5)
points(x = bongo_samples$X, y = bongo_samples$Y, cex = 0.3)
dev.off()

bongo_inla_mesh_cutoff10 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(bongo_all_points$X, bongo_all_points$Y),
  cutoff = 10,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "bongo_boundary_mesh_cutoff10.png"), width=4, height=6, res=200, units="in")
plot(bongo_inla_mesh_cutoff10)
points(x = bongo_samples$X, y = bongo_samples$Y, cex = 0.3)
dev.off()

bongo_inla_mesh_cutoff15 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(bongo_all_points$X, bongo_all_points$Y),
  cutoff = 15,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "bongo_boundary_mesh_cutoff15.png"), width=4, height=6, res=200, units="in")
plot(bongo_inla_mesh_cutoff15)
points(x = bongo_samples$X, y = bongo_samples$Y, cex = 0.3)
dev.off()

bongo_inla_mesh_cutoff20 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(bongo_all_points$X, bongo_all_points$Y),
  cutoff = 20,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "bongo_boundary_mesh_cutoff20.png"), width=4, height=6, res=200, units="in")
plot(bongo_inla_mesh_cutoff20)
points(x = bongo_samples$X, y = bongo_samples$Y, cex = 0.3)
dev.off()

bongo_inla_mesh_cutoff25 <- fmesher::fm_mesh_2d_inla(
  loc = cbind(bongo_all_points$X, bongo_all_points$Y),
  cutoff = 25,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "bongo_boundary_mesh_cutoff25.png"), width=4, height=6, res=200, units="in")
plot(bongo_inla_mesh_cutoff25)
points(x = bongo_samples$X, y = bongo_samples$Y, cex = 0.3)
dev.off()



#### Visualize alignment between mesh and projection grid ####

# transform mesh to sf
fm_as_sfc(inla_mesh_cutoff15) %>% 
   st_set_crs(st_crs(survey_domain_cov_grid)) -> inla_mesh_cutoff15_sf

mesh_plus_points_plot <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA) +
  geom_sf(data = survey_domain_cov_grid) +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("two_stage_models", "figures", "mesh_plus_points_plot.png"), mesh_plus_points_plot,  height = 10, width = 6)

# build this up sequentially for presentation

# just the survey basemap
survey_map_pres  <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA, color = NA)

ggsave(here::here("figures", "presentation_figures", "survey_area_basemap_km.png"), survey_map_pres,  height = 10, width = 6)

# survey area and samples
survey_samples_map <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA, color = NA) +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "survey_samples_map.png"), survey_samples_map,  height = 10, width = 6)

# survey area and mesh and samples
survey_samples_mesh_map <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA) +
  # geom_sf(data = survey_domain_cov_grid) +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "survey_samples_mesh_map.png"), survey_samples_mesh_map,  height = 10, width = 6)


# projection grid
projection_grid_map <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA) +
  geom_sf(data = survey_domain_cov_grid, size = 1) +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "projection_grid_map.png"), projection_grid_map,  height = 10, width = 6)









