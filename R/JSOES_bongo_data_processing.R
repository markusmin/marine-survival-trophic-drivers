## JSOES bongo data processing

# source this script to load and reformat the JSOES bongo data

library(readxl)
library(here)
library(janitor)
library(tidyverse)
library(kableExtra)
library(sf)
library(gstat)
library(sp)
library(ggthemes)
library(ape)


bongo <- clean_names(read_excel(here::here("Data", "Markus_Min_Plankton Density for Trophic Summary Query_10.11.24.xlsx")))

# there's one sample that's missing spatial data: 20 June   2004 062904NH30
# filter(bongo, (is.na(dec_long)))

jsoes <- clean_names(read_excel(here::here("Data", "Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.15.25_FINAL.xlsx")))
# Keep only June data
jsoes <- filter(jsoes, month == "June")

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

trophic_groupings <- clean_names(read_excel(here::here("Data", "Markus_Min_Trophic Groupings_10.2.24.xlsx")))

# change genus species to be consistent across these two files, and usual 
# latin name and common name capitalization
bongo$genus_species <- str_to_sentence(bongo$genus_species)
trophic_groupings$genus_species <- str_to_sentence(trophic_groupings$genus_species)
trophic_groupings$common_name <- str_to_title(trophic_groupings$common_name)

trophic_groupings %>% 
  dplyr::select(genus_species, common_name) %>% 
  filter(!(duplicated(genus_species))) -> common_name_df

# get all unique tows
bongo %>% 
  mutate(sampleID = paste0(station, "_", sample_date)) %>% 
  relocate(sampleID) %>% 
  dplyr::select(-c(genus_species, life_history_stage, id_code, count, sum_of_density_number_m3)) %>% 
  filter(!(duplicated(sampleID))) -> bongo_tows

# complete bongo data, so that we have zeros for each taxon in tows where they weren't caught
bongo %>% 
  mutate(sampleID = paste0(station, "_", sample_date)) %>% 
  relocate(sampleID) %>% 
  # drop tow information; we will add this back later
  dplyr::select(-c(id_code, cruise_number, month, year, station_code, net_type, sample_date,
                   sample_time, dec_lat, dec_long, station_depth_m, transect_name, station_number, 
                   station, n_cr_s, day_night, volume_filtered_m3)) %>% 
  complete(sampleID, 
           nesting(genus_species, life_history_stage), 
           fill = list(count = 0, sum_of_density_number_m3 = 0)) -> bongo_catch

bongo_catch %>% 
  left_join(bongo_tows, by = "sampleID") -> bongo_catch


# add common name
bongo_catch %>% left_join(common_name_df, by = "genus_species") -> bongo_catch

# look into different life history stages


# create a version of this where we collapse all life history stages together
bongo_catch %>% 
  dplyr::select(sampleID, genus_species, count, sum_of_density_number_m3) %>% 
  group_by(sampleID, genus_species) %>% 
  summarise_if(is.numeric, sum) -> bongo_catch_LH_collapse

bongo_catch_LH_collapse %>% 
  left_join(bongo_tows, by = "sampleID") %>% 
  ungroup() -> bongo_catch_LH_collapse


#### Processing bongo data into functional groups ####

# The JSOES bongo data contains data on both the direct prey of Pacific salmon as well as the prey items of the prey of Pacific salmon. The direct salmon prey that is caught in this survey includes Euphausiids, Cancer crab larvae, hyperiid amphipods, and non-cancer crab larvae. The abundance of copepods has also been correlated with salmon marine survival, and as such we will also examine these taxa.

# Let's ID everything that's seen in at least 0.5% of samples. This accounts for 118 of 215 total taxa identified in this survey.

# get count of each genus_species
bongo_catch_LH_collapse %>% 
  filter(sum_of_density_number_m3 > 0) %>% 
  count(genus_species) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_samples = n/nrow(bongo_tows)) %>% 
  left_join(common_name_df, by = "genus_species") %>% 
  relocate(common_name, .after = "genus_species") -> bongo_freq_occ


bongo_taxon_groups <- data.frame(genus_species = bongo_freq_occ$genus_species,
                                 group = c(
                                   "krill", "krill", "krill",
                                   "cold-water copepod", "Cancer crab larvae", "Chaetognaths",
                                   "Hyperiid amphipod", "Barnacle larvae", "sea butterflies",
                                   "Crangon shrimp", "Hyperiid amphipod", "fish eggs",
                                   "Caridean shrimp", "Cnidarians", "copepod",
                                   "warm-water copepod", "rockfishes", "northern anchovy",
                                   "ghost shrimp", "Non-Cancer crab larvae", "Cancer crab larvae",
                                   "copepod", "Non-Cancer crab larvae", "Non-Cancer crab larvae",
                                   "Gammarid amphipod", "Non-Cancer crab larvae", "Cnidarians",
                                   "Cnidarians", "Cancer crab larvae", "Cnidarians",
                                   "Barnacle larvae", "copepod", "Non-Cancer crab larvae",
                                   "Cnidarians", "Hyperiid amphipod", "flatfishes",
                                   "Polychaetes", "Caridean shrimp", "squat lobster",
                                   "insects", "myctophids", "flatfishes",
                                   "flatfishes", "sculpins", "flatfishes",
                                   "Non-Cancer crab larvae", "fishes", "squids",
                                   "Gammarid amphipod", "Osmerids", "Non-Cancer crab larvae",
                                   "Non-Cancer crab larvae", "mysids", "Pacific Tomcod",
                                   "Hyperiid amphipod", "sea butterflies", "sculpins",
                                   "Non-Cancer crab larvae", "copepod", "flatfishes", 
                                   "Non-Cancer crab larvae", "krill", "shrimp larvae",
                                   "snailfishes", "insects", "ronquil/prickleback",
                                   "sea butterflies", "Non-Cancer crab larvae", "snailfishes",
                                   "snailfishes", "sculpins", "salps",
                                   "Hyperiid amphipod", "flatfishes", "salps",
                                   
                                   "snailfishes", "isopods", "goby",
                                   "Hyperiid amphipod", "snailfishes", "Cnidarians",
                                   
                                   "rockfishes", "sculpins", "bivalves",
                                   "Non-Cancer crab larvae", "insects", "myctophids",
                                   "Gammarid amphipod", "insects", "mysids",
                                   "flatfishes", "copepods", "insects",
                                   "goby", "sculpins", "shrimp",
                                   "Cnidarians", "sculpins", "mysids",
                                   "mysids", "Hyperiid amphipod", "Non-Cancer crab larvae",
                                   "Pacific sardine", "myctophid", "mysids",
                                   "sculpins", "Cnidarians", "Caprellid amphipod",
                                   "sculpins", "unknown decapod", "fish",
                                   "Non-Cancer crab larvae", "mysids", "Cnidarians",
                                   "siphonophores", "insects", "Hyperiid amphipod",
                                   "Cnidarians", rep(NA, 97)
                                 ))

# combine the shrimps together
bongo_taxon_groups %>% 
  mutate(group = ifelse(grepl("shrimp|mysids", group), "shrimp larvae", group)) -> bongo_taxon_groups

bongo_catch_LH_collapse %>% 
  left_join(bongo_taxon_groups, by = "genus_species") -> bongo_catch_LH_collapse

#### Reformat and prepare data for modeling ####

# In the last section of this script, we will export the following functional groups 
# for inclusion in spatiotemporal models. The data will be formatted as the total 
# density of each functional group at each station.

# - Cancer crab larvae
# - Non-cancer crab larvae
# - shrimp larvae
# - Hyperiid amphipods


reformat_group_data <- function(data, group_name){
  # keep only this taxon
  group_data <- subset(data, group == group_name)
  
  # collapse by this group - also collapse all maturity stages.
  # we may not want to do this in the future but for now it's okay
  group_data %>%
    group_by(across(c(-genus_species, -count, -sum_of_density_number_m3))) %>%
    summarise(total_sum_of_density_number_m3 = sum(sum_of_density_number_m3)) %>% 
    ungroup() -> group_data
  
  return(group_data)
}

jsoes_bongo_cancer_crab_larvae <- reformat_group_data(data = bongo_catch_LH_collapse, group_name = "Cancer crab larvae")

jsoes_bongo_non_cancer_crab_larvae <- reformat_group_data(data = bongo_catch_LH_collapse, group_name = "Non-Cancer crab larvae")

jsoes_bongo_shrimp_larvae <- reformat_group_data(data = bongo_catch_LH_collapse, group_name = "shrimp larvae")

jsoes_bongo_hyperiid_amphipods <- reformat_group_data(data = bongo_catch_LH_collapse, group_name = "Hyperiid amphipod")

# drop 1998 (bongo doesn't have this year) and 2022:2025 (adult returns aren't complete for those years yet)
jsoes_bongo_cancer_crab_larvae <- subset(jsoes_bongo_cancer_crab_larvae, !(year %in% c(1998, 2022:2025)))
jsoes_bongo_hyperiid_amphipods <- subset(jsoes_bongo_hyperiid_amphipods, !(year %in% c(1998, 2022:2025)))
jsoes_bongo_non_cancer_crab_larvae <- subset(jsoes_bongo_non_cancer_crab_larvae, !(year %in% c(1998, 2022:2025)))
jsoes_bongo_shrimp_larvae <- subset(jsoes_bongo_shrimp_larvae, !(year %in% c(1998, 2022:2025)))



# add a nmi and km from shore from shore column based on the station codes
add_nmi_km <- function(data){
  data %>% 
    mutate(nmi_from_shore = parse_number(station)) %>% 
    # change nmi to km to keep units consistent
    mutate(km_from_shore = nmi_from_shore * 1.852) -> output
  
  return(output)
}

jsoes_bongo_cancer_crab_larvae <- add_nmi_km(jsoes_bongo_cancer_crab_larvae)
jsoes_bongo_hyperiid_amphipods <- add_nmi_km(jsoes_bongo_hyperiid_amphipods)
jsoes_bongo_non_cancer_crab_larvae <- add_nmi_km(jsoes_bongo_non_cancer_crab_larvae)
jsoes_bongo_shrimp_larvae <- add_nmi_km(jsoes_bongo_shrimp_larvae)

#### Add covariates to the bongo data ####
#### load and reformat SST data ####
SST1 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_9bdb_17ce_0aa0.csv"), skip = 1)
SST2 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_36e3_7cb3_1109.csv"), skip = 1)
SST3 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_519c_f955_0fea.csv"), skip = 1)
SST4 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_727c_5c10_e6a0.csv"), skip = 1)
SST5 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_1176_fc5f_53d2.csv"), skip = 1)
SST6 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_ae0c_9197_083f.csv"), skip = 1)
SST7 <- read.csv(here::here("Data", "NOAA_SST", "noaacrwsstDaily_fe4f_9c9c_442e.csv"), skip = 1)

SST1 %>% 
  bind_rows(., SST2, SST3, SST4, SST5, SST6, SST7) -> SST

SST %>% 
  dplyr::rename(time = UTC, latitude = degrees_north, longitude = degrees_east, SST = degree_C) -> SST
SST %>% 
  mutate(time = ymd_hms(time)) %>% 
  mutate(day_of_month = mday(time)) %>% 
  mutate(year = year(time)) %>% 
  mutate(julian = yday(time))-> SST

# Looks good to me! Now just need to resolve CRS conflicts

# Convert SST to UTM zone 10
# st_as_sf(SST, coords = c("latitude", "longitude"), crs = "4326") -> SST_sf
st_as_sf(SST, coords = c("longitude", "latitude"), crs = 4326) -> SST_sf
US_west_coast

sf::st_transform(SST_sf, crs = UTM_zone_10_crs) -> SST_sf_proj

# convert to km
SST_sf_proj_km <- SST_sf_proj
SST_sf_proj_km$geometry <- SST_sf_proj$geometry/1000

SST_sf_proj_km <- st_set_crs(SST_sf_proj_km, UTM_zone_10_crs)

# Match this SST data to the JSOES bongo data
# jsoes needs to be in SF format as well

st_as_sf(jsoes_bongo_cancer_crab_larvae, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_cancer_crab_larvae_sf
st_as_sf(jsoes_bongo_hyperiid_amphipods, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_hyperiid_amphipods_sf
st_as_sf(jsoes_bongo_non_cancer_crab_larvae, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_non_cancer_crab_larvae_sf
st_as_sf(jsoes_bongo_shrimp_larvae, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_shrimp_larvae_sf

jsoes_bongo_cancer_crab_larvae_sf_split <- split(jsoes_bongo_cancer_crab_larvae_sf, jsoes_bongo_cancer_crab_larvae_sf$year)
jsoes_bongo_hyperiid_amphipods_sf_split <- split(jsoes_bongo_hyperiid_amphipods_sf, jsoes_bongo_hyperiid_amphipods_sf$year)
jsoes_bongo_non_cancer_crab_larvae_sf_split <- split(jsoes_bongo_non_cancer_crab_larvae_sf, jsoes_bongo_non_cancer_crab_larvae_sf$year)
jsoes_bongo_shrimp_larvae_sf_split <- split(jsoes_bongo_shrimp_larvae_sf, jsoes_bongo_shrimp_larvae_sf$year)
SST_sf_proj_km_split <- split(SST_sf_proj_km, SST_sf_proj_km$year) 

jsoes_years <- as.character(unique(jsoes_bongo_cancer_crab_larvae_sf$year))

jsoes_bongo_cancer_crab_larvae_remote_SST_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_bongo_cancer_crab_larvae_sf_split[[year]]
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
jsoes_bongo_cancer_crab_larvae <- list_rbind(jsoes_bongo_cancer_crab_larvae_remote_SST_list)

jsoes_bongo_hyperiid_amphipods_sf_split_remote_SST_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_bongo_hyperiid_amphipods_sf_split[[year]]
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
jsoes_bongo_hyperiid_amphipods <- list_rbind(jsoes_bongo_hyperiid_amphipods_sf_split_remote_SST_list)


jsoes_bongo_non_cancer_crab_larvae_remote_SST_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_bongo_non_cancer_crab_larvae_sf_split[[year]]
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
jsoes_bongo_non_cancer_crab_larvae <- list_rbind(jsoes_bongo_non_cancer_crab_larvae_remote_SST_list)


jsoes_bongo_shrimp_larvae_remote_SST_list <- map(jsoes_years, function(year) {
  jsoes_one_year <- jsoes_bongo_shrimp_larvae_sf_split[[year]]
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
jsoes_bongo_shrimp_larvae <- list_rbind(jsoes_bongo_shrimp_larvae_remote_SST_list)

#### Calculate distance to shore ####
# We will use the distance from shore calculated by sf rather than the data
# from the survey to ensure consistency with our projection grid

### cancer crab larvae
dplyr::select(jsoes_bongo_cancer_crab_larvae, station, km_from_shore, X, Y) -> jsoes_bongo_cancer_crab_larvae_spatial
st_as_sf(jsoes_bongo_cancer_crab_larvae_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_cancer_crab_larvae_sf
st_distance(jsoes_bongo_cancer_crab_larvae_sf, US_west_coast_proj_km) -> jsoes_bongo_cancer_crab_larvae_dist_shore
jsoes_bongo_cancer_crab_larvae$sf_dist_shore <- as.numeric(jsoes_bongo_cancer_crab_larvae_dist_shore)

# rescale data
jsoes_bongo_cancer_crab_larvae %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> jsoes_bongo_cancer_crab_larvae

### hyperiid amphipods
dplyr::select(jsoes_bongo_hyperiid_amphipods, station, km_from_shore, X, Y) -> jsoes_bongo_hyperiid_amphipods_spatial
st_as_sf(jsoes_bongo_hyperiid_amphipods_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_hyperiid_amphipods_sf
st_distance(jsoes_bongo_hyperiid_amphipods_sf, US_west_coast_proj_km) -> jsoes_bongo_hyperiid_amphipods_dist_shore
jsoes_bongo_hyperiid_amphipods$sf_dist_shore <- as.numeric(jsoes_bongo_hyperiid_amphipods_dist_shore)

# rescale data
jsoes_bongo_hyperiid_amphipods %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> jsoes_bongo_hyperiid_amphipods

### non-cancer crab larvae
dplyr::select(jsoes_bongo_non_cancer_crab_larvae, station, km_from_shore, X, Y) -> jsoes_bongo_non_cancer_crab_larvae_spatial
st_as_sf(jsoes_bongo_non_cancer_crab_larvae_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_non_cancer_crab_larvae_sf
st_distance(jsoes_bongo_non_cancer_crab_larvae_sf, US_west_coast_proj_km) -> jsoes_bongo_non_cancer_crab_larvae_dist_shore
jsoes_bongo_non_cancer_crab_larvae$sf_dist_shore <- as.numeric(jsoes_bongo_non_cancer_crab_larvae_dist_shore)

# rescale data
jsoes_bongo_non_cancer_crab_larvae %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> jsoes_bongo_non_cancer_crab_larvae

### shrimp larvae
dplyr::select(jsoes_bongo_shrimp_larvae, station, km_from_shore, X, Y) -> jsoes_bongo_shrimp_larvae_spatial
st_as_sf(jsoes_bongo_shrimp_larvae_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> jsoes_bongo_shrimp_larvae_sf
st_distance(jsoes_bongo_shrimp_larvae_sf, US_west_coast_proj_km) -> jsoes_bongo_shrimp_larvae_dist_shore
jsoes_bongo_shrimp_larvae$sf_dist_shore <- as.numeric(jsoes_bongo_shrimp_larvae_dist_shore)

# rescale data
jsoes_bongo_shrimp_larvae %>% 
  mutate(SST_scaled = as.numeric(scale(remote_SST)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> jsoes_bongo_shrimp_larvae
