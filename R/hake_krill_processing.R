# ATM Hake and Krill data processing

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
library(viridis)

#### Load and reformat hake and krill data ####

# load the hake data
# Note that our request was only for hake from 42 to 49 degrees north
hake <- read.csv(here::here("Data", "Hake_NASC_WC_1995-2023.csv"))

# load the krill data
krill <- read.csv(here::here("Data", "Krill_NASC_WC_0.5nmi_2007-2023.csv"))
# add an index column for krill
krill$ID_new <- 1:nrow(krill)

# make the density binned
# figure out what our ranges are
# range(hake$NASC)
# range(krill$NASC)

# There are some negative values here which is really odd - should reach
# out to CPS team if we decide to include this data

# Define bins and labels
breaks <- c(0, 1, 100, 500, 1000, 10000, 1000000)
labels <- c("0-1", "1-100", "100-500", "500-1000", "1000-10000", "10000+")

hake %>% 
  mutate(`NASC (m^2 nmi^-2)` = cut(NASC, 
                                   breaks = breaks,
                                   labels = labels,
                                   include.lowest = TRUE,
                                   right = FALSE)) -> hake

krill %>% 
  mutate(`NASC (m^2 nmi^-2)` = cut(NASC, 
                                   breaks = breaks,
                                   labels = labels,
                                   include.lowest = TRUE,
                                   right = FALSE)) -> krill


#### Spatial data ####
usa_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/USA_adm0.shp", quiet = TRUE)
# load BC
CAN_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/canada/lpr_000b16a_e.shp", quiet = TRUE)
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             c(xmin = -126, ymin = 42, xmax = -123, ymax = 48.5))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -126, ymin = 42, xmax = -123, ymax = 48.5))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)

#### Visualize hake and krill data ####
survey_area_basemap <- ggplot(US_west_coast) +
  geom_sf() +
  geom_sf(data = BC_coast) +
  # coord_sf(ylim = c(44,48.5),  xlim = c(-126, -123)) +
  # scale_x_continuous(breaks = c(124,125,126)) +
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

survey_area_basemap +
  geom_point(data = hake, aes(x = Lon, y = Lat, color = `NASC (m^2 nmi^-2)`), size = 0.5) +
  facet_wrap(~year, nrow = 2) +
  ggtitle("Pacific Hake") +
  scale_color_viridis_d() + 
  guides(colour = guide_legend()) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.spacing.y = unit(0.05, 'cm'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 24)) -> hake_plot

ggsave(here::here( "figures", "hake_dist_plot.png"), hake_plot,  height = 8, width = 12)

survey_area_basemap +
  geom_point(data = krill, aes(x = Lon, y = Lat, color = `NASC (m^2 nmi^-2)`), size = 0.5) +
  facet_wrap(~Year, nrow = 2) +
  ggtitle("Krill") +
  scale_color_viridis_d() + 
  guides(colour = guide_legend()) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.spacing.y = unit(0.05, 'cm'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 24)) -> krill_plot

ggsave(here::here( "figures", "krill_dist_plot.png"), krill_plot,  height = 8, width = 12)


#### Process hake data to get zeros into the data: 2007-2023 ####

# apparently, these are binned at a 0.5 nautical mile resolution

# for 2007-2023, use the krill samples

# convert each to sf and use st_nearest_feature
st_as_sf(krill, coords = c("Lon", "Lat"), crs = 4326) -> krill_sf
st_as_sf(hake, coords = c("Lon", "Lat"), crs = 4326) -> hake_sf



krill_sf_split <- split(krill_sf, krill_sf$Year) 
hake_sf_split <- split(hake_sf, hake_sf$year)

# Find the nearest krill sample for each hake sample
krill_years <- as.character(unique(krill$Year))
hake_krillsamples_list <- map(krill_years, function(year) {
  hake_one_year <- hake_sf_split[[year]]
  krill_one_year <- krill_sf_split[[year]]
  
  # Find nearest features
  krill_index_for_hake <- st_nearest_feature(hake_one_year, krill_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  hake_one_year %>%
    mutate(matched_id = krill_index_for_hake) %>%
    bind_cols(
      krill_one_year[krill_index_for_hake, c("geometry", "ID_new")] %>% 
        dplyr::rename(krill_geometry = geometry)) %>% 
    dplyr::select(-matched_id) -> output
  
  # move the original geometry back to lat/lon columns
  output %>% 
    st_set_geometry("krill_geometry") -> output
  
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) %>% 
    dplyr::rename(hake_geometry = geometry) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
hake_krillsamples <- list_rbind(hake_krillsamples_list)

# now, add zeros for all of the krill samples that aren't found in the matched list

krill %>% 
  filter(!(ID_new %in% hake_krillsamples$ID_new)) %>% 
  dplyr::select(-c(NASC, `NASC (m^2 nmi^-2)`)) %>% 
  mutate(NASC = 0) %>% 
  dplyr::rename(year = Year) -> hake_zeros

hake_zeros %>% 
  bind_rows(dplyr::rename(hake_krillsamples, Lat = Y, Lon = X)) -> hake_new

# these aren't quite the same... but they're very close. Perhaps there are a few
# cases where hake samples were matched to the same krill sample based on distance?
# > nrow(hake_zeros) + nrow(hake_krillsamples)
# [1] 34332
# > nrow(krill)
# [1] 34219

survey_area_basemap +
  geom_point(data = hake_new, aes(x = Lon, y = Lat, color = NASC), size = 0.5) +
  facet_wrap(~year, nrow = 2) +
  ggtitle("Pacific Hake") +
  scale_color_viridis_c() +
  guides(colour = guide_legend()) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.spacing.y = unit(0.05, 'cm'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 24)) -> hake_new_plot

ggsave(here::here( "figures", "hake_new_dist_plot.png"), hake_new_plot,  height = 8, width = 12)

#### Process hake data to get zeros into the data: 2001, 2003, 2005 ####

# Here, we don't have krill data to rely on. So we'll need to go back to the old cruise reports.

# load transect info from old cruise reports
hake_2003_transects <- read.csv(here::here("Data", "hake_2003_transects.csv"))
hake_2003_transects$Start_longitude <- -1 * hake_2003_transects$Start_longitude
hake_2003_transects$End_longitude <- -1 * hake_2003_transects$End_longitude
hake_2005_transects <- read.csv(here::here("Data", "hake_2005_transects.csv"))
hake_2005_transects$Start_longitude <- -1 * hake_2005_transects$Start_longitude
hake_2005_transects$End_longitude <- -1 * hake_2005_transects$End_longitude

# 2001 doesn't have transect data in the cruise report. However, the 2003 report
# states the following:

# We elected to use transects in 2003 identical to those covered in the 2001 survey, starting
# from south of Monterey Bay, California, and covering the area to the most northern extent at
# Dixon Entrance. Seafloor depth at the nearshore end of individual transects was typically 50 m.
# The offshore extent of individual transects typically ranged to depths of about 1500 m. Transects
# were extended deeper if Pacific hake aggregations were detected at or near the predetermined
# endpoints.

# use 2003 transects as 2001 transects
hake_2003_transects %>% 
  mutate(year = 2001) -> hake_2001_transects

# visual comparison to the figure in the 2001 cruise report suggests that they are very similar

survey_area_basemap +
  geom_segment(data = hake_2003_transects, aes(x = Start_longitude, xend = End_longitude,
                                               y = Start_latitude, yend = End_latitude))

survey_area_basemap +
  geom_segment(data = hake_2005_transects, aes(x = Start_longitude, xend = End_longitude,
                                               y = Start_latitude, yend = End_latitude))


# combine all of the other transects
mutate(hake_2001_transects, year = 2001) %>% 
  bind_rows(mutate(hake_2003_transects, year = 2003)) %>% 
  bind_rows(mutate(hake_2005_transects, year = 2005)) -> hake_transects

# split this into 0.5 nm bins - same spacing as other data
# Let's look at what the typical distance is when samples are binned
test_krill_transect <- subset(krill, Transect == "x39hake69")

transect_lat_diffs <- rep(NA, nrow(test_krill_transect)-1)
transect_lon_diffs <- rep(NA, nrow(test_krill_transect)-1)
for (i in 1:(nrow(test_krill_transect)-1)){
  transect_lat_diffs[i] <- test_krill_transect$Lat[i+1] - test_krill_transect$Lat[i]
  transect_lon_diffs[i] <- test_krill_transect$Lon[i+1] - test_krill_transect$Lon[i]
}

hist(transect_lon_diffs, breaks = 30)
median_lon_spacing <- median(transect_lon_diffs)

hake_transects_binned <- data.frame()

for (i in 1:nrow(hake_transects)){
  start_lat <- hake_transects[i,"Start_latitude"]
  end_lat <- hake_transects[i,"End_latitude"]
  start_lon <- hake_transects[i,"Start_longitude"]
  end_lon <- hake_transects[i,"End_longitude"]
  
  # create sequence of longitudes
  seq_direction <- ifelse(start_lon > end_lon, -1, 1)
  lon_seq <- seq(start_lon, end_lon, by = seq_direction * median_lon_spacing)
  lat_seq <- seq(start_lat, end_lat, length = length(lon_seq))
  
  # create a dataframe that stores this transect expanded as binned points
  transect_binned <- data.frame(ID = paste0("new", i),
                                Transect = paste0(hake_transects$year[i],"_", hake_transects$Transect[i]),
                                year = hake_transects$year[i],
                                Survey = as.character(hake_transects$year[i]),
                                Date = NA,
                                Lat = lat_seq,
                                Lon = lon_seq)
  
  hake_transects_binned %>% 
    bind_rows(transect_binned) -> hake_transects_binned
 
}

# add an transects binned
hake_transects_binned$ID_new <- 1:nrow(hake_transects_binned)

# add positive data  

# convert each to sf and use st_nearest_feature
st_as_sf(hake_transects_binned, coords = c("Lon", "Lat"), crs = 4326) -> hake_transects_binned_sf
st_as_sf(hake, coords = c("Lon", "Lat"), crs = 4326) -> hake_sf



hake_transects_binned_sf_split <- split(hake_transects_binned_sf, hake_transects_binned_sf$year) 
hake_sf_split <- split(hake_sf, hake_sf$year)

# Find the nearest hake_transects_binned sample for each hake sample
hake_transects_binned_years <- as.character(unique(hake_transects_binned$year))
hake_transect_bins_matched_list <- map(hake_transects_binned_years, function(year) {
  hake_one_year <- hake_sf_split[[year]]
  hake_transects_binned_one_year <- hake_transects_binned_sf_split[[year]]
  
  # Find nearest features
  hake_transects_binned_index_for_hake <- st_nearest_feature(hake_one_year, hake_transects_binned_one_year)
  
  # Combine with matched points (optional: include matched geometry or attributes)
  hake_one_year %>%
    mutate(matched_id = hake_transects_binned_index_for_hake) %>%
    bind_cols(
      hake_transects_binned_one_year[hake_transects_binned_index_for_hake, c("geometry", "ID_new")] %>% 
        dplyr::rename(hake_transects_binned_geometry = geometry)) %>% 
    dplyr::select(-matched_id) -> output
  
  # move the original geometry back to lat/lon columns
  output %>% 
    st_set_geometry("hake_transects_binned_geometry") -> output
  
  
  # transform this to a df
  st_drop_geometry(output) %>% 
    bind_cols(st_coordinates(output)) %>% 
    dplyr::rename(hake_geometry = geometry) -> output
  
  return(output)
})

# take the list and turn it back into an sf object
hake_transect_bins_matched <- list_rbind(hake_transect_bins_matched_list)

# now, add zeros for all of the transect binned samples that aren't found in the matched list

hake_transects_binned %>%
  filter(!(ID_new %in% hake_transect_bins_matched$ID_new)) %>%
  mutate(NASC = 0) -> hake_2001_2003_2005_zeros

# get rid of any zeros below 42 degrees - that's the lower limit of our positive catches
# get rid of any zeros above 49 degrees - that's the upper limit of our positive catchces
hake_2001_2003_2005_zeros %>% 
  filter(Lat >= 42 & Lat <= 49) -> hake_2001_2003_2005_zeros


hake_2001_2003_2005_zeros %>%
  bind_rows(dplyr::rename(hake_transect_bins_matched, Lat = Y, Lon = X)) -> hake_2001_2003_2005

# visualize this

survey_area_basemap +
  geom_point(data = hake_2001_2003_2005, aes(x = Lon, y = Lat, color = NASC), size = 0.5) +
  facet_wrap(~year, nrow = 1) +
  ggtitle("Pacific Hake") +
  scale_color_viridis_c() +
  guides(colour = guide_legend()) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = "right",
        legend.spacing.y = unit(0.05, 'cm'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        plot.title = element_text(size = 24)) -> hake_2001_2003_2005_plot

ggsave(here::here( "figures", "hake_2001_2003_2005_dist_plot.png"), hake_2001_2003_2005_plot,  height = 8, width = 12)

#### combine all processed hake data and export ####
hake_2001_2003_2005 %>% 
  bind_rows(hake_new) %>% 
  dplyr::select(Transect, year, Date, Lat, Lon, NASC) -> hake_for_export

write.csv(hake_for_export, here::here("model_inputs", "hake_processed.csv"), row.names = FALSE)








