### PRS/PWCC surveys: Create meshes for use with sdmTMB or TMB SDMs

# Author: Markus Min
# Last updated: 2025-08-28

# Description: This script will create the mesh necessary for SPDE (used in packages such as sdmTMB)
# for the Pre-Recruit and Pacific Whiting Conservation Cooperative surveys

# We will need the following data to construct the mesh:
# 1. The PRS and PWCC Data (for sampling locations)
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

# Make sure that the data processing script has been run before you run this script
# source(here::here("R", "PRS_PWCC_data_processing.R"))

#### Data 1: load the processed rockfish data ####
prs_rf <- read.csv(here::here("model_inputs", "PRS_rockfishes.csv"))
pwcc_rf <- read.csv(here::here("model_inputs", "PWCC_rockfishes.csv"))

# spatial cut off - keep only samples north of Cape Mendocino (biogeographic break
# that is considered the southern edge of the NCC) - 40.42 latitude

prs_rf <- subset(prs_rf, start_latitude >= 40.42)
pwcc_rf <- subset(pwcc_rf, net_fishing_lat_dd >= 40.42)


# remove the egg mass from PRS data
prs_rf %>% 
  mutate(to_remove = ifelse(comments == "Egg mass", 1, 0)) %>% 
  filter(is.na(to_remove) | to_remove == 0) %>% 
  dplyr:: select(-to_remove) -> prs_rf


# reformat the PWCC data to match PRS data
# note that when a trawl starts before midnight and ends the next day, the date-times lack a time.
# so, we have to account for that
pwcc_rf %>% 
  mutate(net_fishing_time = ifelse(nchar(net_fishing_time) == 10, paste0(net_fishing_time, " 00:00:00"), net_fishing_time)) %>% 
  mutate(net_fishing_time = ymd_hms(net_fishing_time)) %>% 
  mutate(year = year(net_fishing_time),
         month = month(net_fishing_time),
         day = day(net_fishing_time)) -> pwcc_rf


# combine the prs and pwcc data
mutate(dplyr::select(pwcc_rf, year, month, day, X, Y, total), survey = "PWCC") %>% 
  bind_rows(mutate(dplyr::select(prs_rf, year, month, day, X, Y, total), survey = "PRS")) -> rf

# temporal cut off - 2021 is the last year included
rf <- subset(rf, year <= 2021)






#### Data 2: load and reformat spatial data ####
usa_spdf <- st_read(here::here("Data", "map_files", "USA_adm0.shp"))
# load BC
CAN_spdf <- st_read(here::here("Data", "map_files", "canada", "lpr_000b16a_e.shp"))
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             # c(xmin = -126, ymin = 40.42, xmax = -120, ymax = 48.5))
# c(xmin = -126, ymin = 39.42, xmax = -120, ymax = 48.5))
c(xmin = -126, ymin = 39.42, xmax = -123, ymax = 48.5))

# US_west_coast <- sf::st_crop(usa_spdf,
#                              c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))

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

# second version for the full domain
survey_area_basemap_km_PRS_PWCC <- ggplot(US_west_coast_proj_km) +
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


#### Data 3: load and reformat covariate data ####

#### Covariate 1: SST
# Please note that currently the SST data only goes through 2023
# This data is the Sea Surface Temperature, NOAA Coral Reef Watch Daily Global 5km 
# Satellite SST (CoralTemp) dataset was downloaded from the [NOAA ERDDAP server](https://coastwatch.noaa.gov/erddap/griddap/noaacrwsstDaily.html)
SST_1998 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-1998.csv"), skip = 1)
SST_1999 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-1999.csv"), skip = 1)
SST_2000 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2000.csv"), skip = 1)
SST_2001 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2001.csv"), skip = 1)
SST_2002 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2002.csv"), skip = 1)
SST_2003 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2003.csv"), skip = 1)
SST_2004 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2004.csv"), skip = 1)
SST_2005 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2005.csv"), skip = 1)
SST_2006 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2006.csv"), skip = 1)
SST_2007 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2007.csv"), skip = 1)
SST_2008 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2008.csv"), skip = 1)
SST_2009 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2009.csv"), skip = 1)
SST_2010 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2010.csv"), skip = 1)
SST_2011 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2011.csv"), skip = 1)
SST_2012 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2012.csv"), skip = 1)
SST_2013 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2013.csv"), skip = 1)
SST_2014 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2014.csv"), skip = 1)
SST_2015 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2015.csv"), skip = 1)
SST_2016 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2016.csv"), skip = 1)
SST_2017 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2017.csv"), skip = 1)
SST_2018 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2018.csv"), skip = 1)
SST_2019 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2019.csv"), skip = 1)
SST_2020 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2020.csv"), skip = 1)
SST_2021 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2021.csv"), skip = 1)
SST_2022 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2022.csv"), skip = 1)
SST_2023 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2023.csv"), skip = 1)
SST_2024 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsst-2024.csv"), skip = 1)

SST_1998 %>% 
  bind_rows(., SST_1999, SST_2000, SST_2001, SST_2002, SST_2003, SST_2004, SST_2005,
            SST_2006, SST_2007, SST_2008, SST_2009, SST_2010, SST_2011, SST_2012,
            SST_2013, SST_2014, SST_2015, SST_2016, SST_2017, SST_2018, SST_2019,
            SST_2020, SST_2021, SST_2022, SST_2023, SST_2024) %>% 
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
# match satellite data to the trawl data

# For this, rf needs to be in SF format as well
st_as_sf(rf, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> rf_sf

rf_sf_split <- split(rf_sf, rf_sf$year)
SST_sf_proj_km_split <- split(SST_sf_proj_km, SST_sf_proj_km$year) 

rf_years <- as.character(unique(rf_sf$year))

# Find the nearest remotely sensed SST measurement for each rf sample
rf_remote_SST_list <- map(rf_years, function(year) {
  rf_one_year <- rf_sf_split[[year]]
  SST_one_year <- SST_sf_proj_km_split[[year]]
  
  # Find nearest features
  SST_index_for_rf <- st_nearest_feature(rf_one_year, SST_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  rf_one_year %>%
    mutate(matched_id = SST_index_for_rf) %>%
    bind_cols(
      st_drop_geometry(SST_one_year[SST_index_for_rf, "SST"]) %>%
        dplyr::rename(remote_SST = SST)
    ) %>% 
    dplyr::select(-matched_id) -> output
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
rf <- list_rbind(rf_remote_SST_list)

#### Calculate distance to shore

dplyr::select(rf, X, Y) -> rf_spatial
st_as_sf(rf_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> rf_sf
st_distance(rf_sf, US_west_coast_proj_km) -> rf_dist_shore
rf$sf_dist_shore <- as.numeric(rf_dist_shore)


# Rescale data for fitting to covariates
rf %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> rf


#### Create a survey prediction grid ####
# Create a survey grid from scratch, using SST data

# only keep one year, since that's all we need to make a grid
SST_sf_proj %>% 
  filter(year == 2000) -> SST_sf_proj_2000

st_as_sf(SST_sf_proj_2000$geometry/1000, crs = UTM_zone_10_crs) -> SST_sf_proj_km_2000

SST_sf_proj_km_2000$dist_shore <- as.numeric(st_distance(SST_sf_proj_km_2000, US_west_coast_proj_km))


# KEY DECISION: To what area do you want to project the SDM to?
# the southern edge will be the same: The bottom edge of PRS.
# the northern edge will vary
# Here, we will do this two ways:
# 1) Northern edge to the full extent of the JSOES domain
# 2) Northern edge limited to the extent of the PRS survey

# visualize the overlap
survey_area_basemap_km_PRS_PWCC +
  geom_point(data = rf, aes(x = X, y = Y, color = survey))

ggplot(subset(rf, survey == "PRS"), aes(x = X, y = Y)) +
  geom_point() + 
  facet_wrap(~year)

# 2011 extends much further south than the other years

survey_area_basemap_km_PRS_PWCC +
  geom_point(data = rf, aes(x = X, y = Y, color = as.factor(year)))

# Missing 2010 and 2012 data

# 1) use the JSOES projection grid
jsoes <- clean_names(read_excel(here::here("Data", "Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.15.25_FINAL.xlsx")))

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

# For the N/S dimension, the most northerly transect is typically the FS transect and 
# the most southerly transect is the NH transect.
# this corresponds to 44.67 to 48.23 N
# Let's take these as our N/S boundaries and add 10 km as a buffer.
min_Y_jsoes <- mean(subset(jsoes_long, grepl("NH", station))$Y)
max_Y_jsoes <- mean(subset(jsoes_long, grepl("FS", station))$Y)

# min_Y_PRS <- min(subset(rf, year == 2016)$Y)
# max_Y_PRS <- max(subset(rf, year == 2016)$Y)
min_Y_PRS <- min(rf$Y)
max_Y_PRS <- max(rf$Y)


# Based on this, I'm going to propose that we generate a prediction that extends:
# from 0.5 to 65 km offshore
# from 44.67 to 48.23 N (+ 10 km on either side)

# trim longitudinally
SST_sf_proj_km_2000 %>% 
  filter(dist_shore <= 65 & dist_shore >= 0.5) -> grid_within_65km_shore

# trim latitudinally
survey_domain_jsoes <- st_crop(grid_within_65km_shore, xmin = 0, xmax = 100000, ymin=min_Y_jsoes-10, ymax=max_Y_jsoes+10)

# let's inspect the grid that we created
ggplot(survey_domain_jsoes) +
  geom_sf()

# 2) Limited to the extent of the PRS survey

# ok, let's try cutting off by the 2016 PRS survey. Looks like the typical survey coverage for PRS.
survey_area_basemap_km_PRS_PWCC +
  geom_point(data = subset(rf, year == 2016), aes(x = X, y = Y))


# trim longitudinally
SST_sf_proj_km_2000 %>% 
  filter(dist_shore <= 65 & dist_shore >= 0.5) -> grid_within_65km_shore

# trim latitudinally
survey_domain_PRS <- st_crop(grid_within_65km_shore, xmin = 0, xmax = 100000, ymin=min_Y_PRS-10, ymax=max_Y_jsoes+10)

# let's inspect the grid that we created
ggplot(survey_domain_PRS) +
  geom_sf()

### Trim 1 and 2 to manually trim out some of these areas

# manually remove sections for strait, hood canal (and other inland waters), Grays Harbor, Willapa Bay, Columbia River estuary, Arcata Bay
arcata_bay <- st_as_sfc(st_bbox(c(xmin=-124.265, xmax=-124.06525, ymin=40.6765, ymax=40.87462), crs = "WGS84"))
arcata_bay_proj <- sf::st_transform(arcata_bay, crs = UTM_zone_10_crs)
arcata_bay_proj_km <- st_as_sf(arcata_bay_proj/1000, crs = UTM_zone_10_crs)

strait <- st_as_sfc(st_bbox(c(xmin=-124.65, xmax=-122, ymin=47.9375, ymax=49), crs = "WGS84"))
strait_proj <- sf::st_transform(strait, crs = UTM_zone_10_crs)
strait_proj_km <- st_as_sf(strait_proj/1000, crs = UTM_zone_10_crs)

inland_waters <- st_as_sfc(st_bbox(c(xmin=-123.5, xmax=-120, ymin=40, ymax=49.5), crs = "WGS84"))
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

survey_domain_jsoes %>% 
  st_difference(., arcata_bay_proj_km) %>% 
  st_difference(., strait_proj_km) %>% 
  st_difference(., inland_waters_proj_km) %>% 
  st_difference(., grays_harbor_proj_km) %>% 
  st_difference(., willapa_bay_proj_km) %>% 
  st_difference(., estuary_proj_km) -> survey_domain_jsoes

# Create a concave hull around points
st_concave_hull(st_union(survey_domain_jsoes), ratio = 0.1) -> survey_domain_jsoes_polygon

survey_domain_PRS %>% 
  st_difference(., arcata_bay_proj_km) %>% 
  st_difference(., strait_proj_km) %>% 
  st_difference(., inland_waters_proj_km) %>% 
  st_difference(., grays_harbor_proj_km) %>% 
  st_difference(., willapa_bay_proj_km) %>% 
  st_difference(., estuary_proj_km) -> survey_domain_PRS

# Create a concave hull around points
st_concave_hull(st_union(survey_domain_PRS), ratio = 0.1) -> survey_domain_PRS_polygon

# let's compare our survey domains
survey_area_basemap_km_PRS_PWCC +
  geom_sf(data = survey_domain_jsoes_polygon, fill = "blue", alpha = 0.2)

survey_area_basemap_km_PRS_PWCC +
  geom_sf(data = survey_domain_PRS_polygon, fill = "blue", alpha = 0.2)

### add covariates to both survey domains ###
# Take SST grid and intersect with survey domain
SST_sf_proj %>% 
  mutate(geometry = geometry/1000) %>% 
  st_set_crs(UTM_zone_10_crs) %>% 
  st_intersection(survey_domain_jsoes_polygon) %>% 
  dplyr::select(SST, year, geometry)-> SST_survey_domain_jsoes

# add in distance from shore - this is now your survey domain, with covariates
SST_survey_domain_jsoes %>% 
  mutate(sf_dist_shore = as.numeric(st_distance(geometry, US_west_coast_proj_km))) -> survey_domain_jsoes_cov

# Take SST grid and intersect with survey domain
SST_sf_proj %>% 
  mutate(geometry = geometry/1000) %>% 
  st_set_crs(UTM_zone_10_crs) %>% 
  st_intersection(survey_domain_PRS_polygon) %>% 
  dplyr::select(SST, year, geometry)-> SST_survey_domain_PRS

# add in distance from shore - this is now your survey domain, with covariates
SST_survey_domain_PRS %>% 
  mutate(sf_dist_shore = as.numeric(st_distance(geometry, US_west_coast_proj_km))) -> survey_domain_PRS_cov

#### Trim to the years between max and min rf years ####
survey_domain_jsoes_cov %>% 
  filter(year >= min(rf$year) & year <= max(rf$year)) -> survey_domain_jsoes_cov_rf_allyears

# because the survey_domain_jsoes_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_jsoes_cov_rf_allyears %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_jsoes_cov_rf_allyears

# plot it - looks good!
ggplot(survey_domain_jsoes_cov_rf_allyears) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_jsoes_cov_rf_allyears)) -> survey_domain_jsoes_cov_rf_allyears_coords
survey_predict_grid_jsoes_rf_allyears <- data.frame(X = survey_domain_jsoes_cov_rf_allyears_coords$X,
                                                    Y = survey_domain_jsoes_cov_rf_allyears_coords$Y,
                                                    SST = survey_domain_jsoes_cov_rf_allyears$SST,
                                                    year = survey_domain_jsoes_cov_rf_allyears$year,
                                                    dist_shore = as.numeric(survey_domain_jsoes_cov_rf_allyears$sf_dist_shore))


# rescale prediction grid
survey_predict_grid_jsoes_rf_allyears %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid_jsoes_rf_allyears

# Trim to only the years of data
survey_domain_PRS_cov %>% 
  filter(year >= min(rf$year) & year <= max(rf$year)) -> survey_domain_PRS_cov_rf_allyears

# because the survey_domain_PRS_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_PRS_cov_rf_allyears %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_PRS_cov_rf_allyears

# plot it - looks good!
ggplot(survey_domain_PRS_cov_rf_allyears) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_PRS_cov_rf_allyears)) -> survey_domain_PRS_cov_rf_allyears_coords
survey_predict_grid_PRS <- data.frame(X = survey_domain_PRS_cov_rf_allyears_coords$X,
                                      Y = survey_domain_PRS_cov_rf_allyears_coords$Y,
                                      SST = survey_domain_PRS_cov_rf_allyears$SST,
                                      year = survey_domain_PRS_cov_rf_allyears$year,
                                      dist_shore = as.numeric(survey_domain_PRS_cov_rf_allyears$sf_dist_shore))

# rescale prediction grid
survey_predict_grid_PRS %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid_PRS

# check that we don't have any points with missing data
survey_area_basemap_km_PRS_PWCC +
  geom_point(data = subset(survey_predict_grid_PRS, is.na(SST)), aes(x = X, y = Y))


#### Trim to only the years of data in rf ####
survey_domain_jsoes_cov %>% 
  filter(year %in% unique(rf$year)) -> survey_domain_jsoes_cov_rf_years

# because the survey_domain_jsoes_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_jsoes_cov_rf_years %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_jsoes_cov_rf_years

# plot it - looks good!
ggplot(survey_domain_jsoes_cov_rf_years) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_jsoes_cov_rf_years)) -> survey_domain_jsoes_cov_rf_years_coords
survey_predict_grid_jsoes_rf_years <- data.frame(X = survey_domain_jsoes_cov_rf_years_coords$X,
                                        Y = survey_domain_jsoes_cov_rf_years_coords$Y,
                                        SST = survey_domain_jsoes_cov_rf_years$SST,
                                        year = survey_domain_jsoes_cov_rf_years$year,
                                        dist_shore = as.numeric(survey_domain_jsoes_cov_rf_years$sf_dist_shore))


# rescale prediction grid
survey_predict_grid_jsoes_rf_years %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid_jsoes_rf_years

# Trim to only the years of data
survey_domain_PRS_cov %>% 
  filter(year %in% unique(rf$year)) -> survey_domain_PRS_cov_rf_years

# because the survey_domain_PRS_polygon is constructed from a concave hull, some of the borders
# are a little coarse - and as a result some of the super nearshore data snuck back it.
# filter it out
survey_domain_PRS_cov_rf_years %>% 
  filter(sf_dist_shore >= 0.5) -> survey_domain_PRS_cov_rf_years

# plot it - looks good!
ggplot(survey_domain_PRS_cov_rf_years) + geom_sf()

# Now reformat this object to be a data frame that we can use to predict, rather than an sf object
as.data.frame(st_coordinates(survey_domain_PRS_cov_rf_years)) -> survey_domain_PRS_cov_rf_years_coords
survey_predict_grid_PRS <- data.frame(X = survey_domain_PRS_cov_rf_years_coords$X,
                                      Y = survey_domain_PRS_cov_rf_years_coords$Y,
                                      SST = survey_domain_PRS_cov_rf_years$SST,
                                      year = survey_domain_PRS_cov_rf_years$year,
                                      dist_shore = as.numeric(survey_domain_PRS_cov_rf_years$sf_dist_shore))

# rescale prediction grid
survey_predict_grid_PRS %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid_PRS

# check that we don't have any points with missing data
survey_area_basemap_km_PRS_PWCC +
  geom_point(data = subset(survey_predict_grid_PRS, is.na(SST)), aes(x = X, y = Y))

# Ok - now we are finally ready to create our mesh!

#### Create the meshes for SDM ####

#### Create the trawl mesh

# In order to get our model to predict throughout the survey domain, we'll need to rbind
# the projection grid and sampling coordinates

# To construct our mesh, we want to restrict our analysis to the samples north of Cape Mendocino (40.42 N)

cape_mendocino <- st_as_sf(data.frame(lat = 40.42, lon = -124.4), coords = c("lon", "lat"), crs = 4326)

cape_mendocino <- sf::st_transform(cape_mendocino, crs = UTM_zone_10_crs) 
cape_mendocino <- st_as_sf(cape_mendocino$geometry/1000, crs = UTM_zone_10_crs)


rf_samples_study_domain <- subset(rf, Y >= as.data.frame(st_coordinates(cape_mendocino))$Y)



## For JSOES survey domain:
# first extract the grid for the projection area
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_jsoes_cov_rf_years %>% 
  filter(year == 2016) %>% 
  dplyr::select(geometry) -> survey_domain_jsoes_cov_rf_years_grid

# extract only coordinates
jsoes_projection_coords <- as.data.frame(st_coordinates(survey_domain_jsoes_cov_rf_years_grid))
colnames(jsoes_projection_coords) <- c("X", "Y")

# combine the samples and the projection grid
dplyr::select(rf_samples_study_domain, X, Y) %>% 
  bind_rows(jsoes_projection_coords) -> rf_jsoes_all_points

## For PRS survey domain:
# first extract the grid for the projection area
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_PRS_cov_rf_years %>% 
  filter(year == 2016) %>% 
  dplyr::select(geometry) -> survey_domain_PRS_cov_rf_years_grid

# extract only coordinates
PRS_projection_coords <- as.data.frame(st_coordinates(survey_domain_PRS_cov_rf_years_grid))
colnames(PRS_projection_coords) <- c("X", "Y")

# combine the samples and the projection grid
dplyr::select(rf_samples_study_domain, X, Y) %>% 
  bind_rows(PRS_projection_coords) -> rf_PRS_all_points

# From these objects, we then can create a mesh that will be used by the SPDE method.

# use a boundary to ensure that our mesh matches our survey grid
# remove some samples to ensure no bulges in the mesh
# bnd <- INLA::inla.nonconvex.hull(cbind(rf_samples_study_domain$X, rf_samples_study_domain$Y), convex = -0.1)
bnd <- INLA::inla.nonconvex.hull(cbind(subset(rf_samples_study_domain, sf_dist_shore < 45*1.852)$X, subset(rf_samples_study_domain, sf_dist_shore < 45*1.852)$Y), convex = -0.1)

# make versions of the mesh using the boundary object at varying resolution

inla_mesh_cutoff5_PRS_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_PRS_all_points$X, rf_PRS_all_points$Y),
  cutoff = 5,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff5_PRS_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff5_PRS_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff5_jsoes_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_jsoes_all_points$X, rf_jsoes_all_points$Y),
  cutoff = 5,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff5_jsoes_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff5_jsoes_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()

# I think it makes more sense to proceed with the full JSOES domain. You can always
# trim the domain to calculate overlap post-hoc


inla_mesh_cutoff10_PRS_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_PRS_all_points$X, rf_PRS_all_points$Y),
  cutoff = 10,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff10_PRS_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff10_PRS_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()


inla_mesh_cutoff10_jsoes_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_jsoes_all_points$X, rf_jsoes_all_points$Y),
  cutoff = 10,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff10_jsoes_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff10_jsoes_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff15_PRS_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_PRS_all_points$X, rf_PRS_all_points$Y),
  cutoff = 15,
  boundary = bnd
)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff15_PRS_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff15_PRS_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()

inla_mesh_cutoff15_jsoes_domain <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_jsoes_all_points$X, rf_jsoes_all_points$Y),
  cutoff = 15,
  boundary = bnd
)

# make the same mesh using sdmTMB for compatibility with sdmTMB() function
inla_mesh_cutoff15_jsoes_domain_sdmTMB <- make_mesh(
  data = rf_jsoes_all_points,
  xy_cols = c("X", "Y"),
  cutoff = 15,
  boundary = bnd
)
# confirm that they are the same
plot(inla_mesh_cutoff15_jsoes_domain_sdmTMB)
plot(inla_mesh_cutoff15_jsoes_domain)

png(here::here("two_stage_models", "figures", "inla_mesh_cutoff15_jsoes_domain.png"), width=4, height=6, res=200, units="in")
plot(inla_mesh_cutoff15_jsoes_domain)
points(x = rf$X, y = rf$Y, cex = 0.3)
dev.off()


#### Visualize alignment between mesh and projection grid ####

# transform mesh to sf
fm_as_sfc(inla_mesh_cutoff15_jsoes_domain) %>% 
   st_set_crs(st_crs(survey_domain_jsoes_cov_rf_years_grid)) -> inla_mesh_cutoff15_jsoes_domain_sf

rockfish_mesh_plus_points_plot <- survey_area_basemap_km_PRS_PWCC +
  geom_sf(data = inla_mesh_cutoff15_jsoes_domain_sf, fill = NA) +
  geom_sf(data = survey_domain_jsoes_cov_rf_years_grid) +
  geom_point(data = rf, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("two_stage_models", "figures", "rockfish_mesh_plus_points_plot.png"), rockfish_mesh_plus_points_plot,  height = 10, width = 6)


# build this up sequentially for presentation

# survey area and mesh and samples
rockfish_samples_mesh_map <- survey_area_basemap_km_PRS_PWCC +
  geom_sf(data = inla_mesh_cutoff15_jsoes_domain_sf, fill = NA) +
  # geom_sf(data = survey_domain_cov_grid) +
  # geom_point(data = subset(rf, Y >= 4800), aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)
  geom_point(data = rf, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "rockfish_samples_mesh_map.png"), rockfish_samples_mesh_map,  height = 10, width = 6)
