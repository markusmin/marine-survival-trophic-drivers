### 09_Fig1_map

# This script generates Figure 1, which shows a map of the study region.

## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)

## source the make mesh scripts to load the spatial data and survey data

# source the prep scripts for each of the surveys
# source(here::here("R", "CCES_make_mesh.R"))
source(here::here("R", "PRS_PWCC_make_mesh.R"))
source(here::here("R", "hake_survey_make_mesh.R"))
source(here::here("R", "JSOES_seabirds_make_mesh.R"))
source(here::here("R", "JSOES_make_mesh.R"))

#### Load Columbia River Basin data and map freshwater map ####

# Record two CRSs: WGS 84 (for visualization) and UTM zone 10 (this is the CRS that we use for modeling)
WGS_CRS <- 4326
UTM_zone_10_crs <- 32610

# record coordinates of Bonneville Dam
BON_coords <- c(-121.940792, 45.644344)

## load USA states and territories
us_states <- st_read(here("map_files", "us_states_territories", "s_05mr24.shp"))
us_states_proj <- sf::st_transform(us_states, crs = WGS_CRS)

## load CAN outline (for Vancouver Island)

BC_spdf <- st_read(here("map_files", "canada", "lpr_000b16a_e.shp"))
BC_proj <- sf::st_transform(BC_spdf, crs = WGS_CRS)

## combine CAN into one so that we don't have the province boundaries
one_CAN_spdf <- st_union(BC_spdf)


### load Major rivers
rivers_spdf <- st_read(here::here("map_files", "NA_Lakes_and_Rivers","hydrography_l_rivers_v2.shp"))
rivers_proj <- sf::st_transform(rivers_spdf, crs = WGS_CRS)

### load the Columbia river basin boundary
CRB_boundary_spdf <- st_read(here::here("map_files", "Columbia_Basin_Watershed_Boundary","Columbia_Basin_Watershed_Boundary.shp"))
CRB_boundary_proj <- sf::st_transform(CRB_boundary_spdf, crs = WGS_CRS)

# keep only the largest rivers for the inset map
rivers_proj
major_rivers_proj <- filter(rivers_proj, Shape_Leng >= 90092)

# keep only the columbia river basin area
CRB_rivers_proj <- st_crop(rivers_proj, xmin = -125.15, xmax = -113.23, ymin = 44.3, ymax = 48.45)

# subset only the columbia river
# The columbia River estuary and the McNary reservoir aren't considered parts of the Columbia or Snake
# in the shapefile, so we will also need to subset those individually

# let's make those manually, since they have this weird splitting thing
mcnary_reservoir_points <- data.frame(x = c(-119.324995, -119.253927, -119.166036, -119.064413,
                                            -119.014974, -118.946653, -118.937040, -118.985105,
                                            -119.031287, -119.104582, -119.214445, -119.242597),
                                      y = c(45.932081, 45.938767, 45.929216, 45.957864,
                                            45.976955, 46.029897, 46.093741, 46.137532,
                                            46.205528, 46.217883, 46.242583, 46.264898),
                                      NAMEEN = rep("McNary Reservoir", 12))

mcnary_reservoir_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> mcnary_reservoir_linestring

estuary_points <- data.frame(x = c(-122.876, -123.184092, -123.228724, -123.257906, -123.326914,
                                   -123.405878, -123.486216, -123.506472),
                             y = c(46.0707, 46.183222, 46.162063, 46.146128, 46.149458,
                                   46.187501, 46.250451, 46.258522),
                             NAMEEN = rep("Columbia River Estuary", 8))

estuary_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> estuary_linestring

# john day mouth
john_day_mouth_points <- data.frame(x = c(-120.513195, -120.530791, -120.547356, -120.553192, -120.564179,
                                          -120.573448, -120.602459, -120.606751, -120.622200,
                                          -120.643658, -120.652584),
                                    y = c(45.665883, 45.671761, 45.669362, 45.677398, 45.683395,
                                          45.698024, 45.703059, 45.714447, 45.721878,
                                          45.722357, 45.736976),
                                    NAMEEN = rep("John Day River Mouth", 11))

john_day_mouth_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> john_day_mouth_linestring

# walla walla mouth
walla_walla_mouth_points <- data.frame(x = c(-118.929555, -118.865462),
                                       y = c(46.063718, 46.058101),
                                       NAMEEN = rep("Walla Walla River Mouth", 2))

walla_walla_mouth_points %>% 
  st_as_sf(coords = c("x", "y"), crs = WGS_CRS) %>% 
  group_by(NAMEEN) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING") -> walla_walla_mouth_linestring



# create bbox objects with which to crop
estuary_bbox <- st_bbox(c(xmin=-124, xmax=-123, ymin=46, ymax=46.25), crs = "WGS84")
mcnary_reservoir_bbox <- st_bbox(c(xmin=-119.35, xmax=-118.79, ymin=45.87, ymax=46.27), crs = "WGS84")
rivers_proj %>% 
  st_crop(estuary_bbox) -> estuary_proj
rivers_proj %>% 
  st_crop(mcnary_reservoir_bbox) %>% 
  filter(!(NAMEEN == "Yakima River")) -> mcnary_reservoir_proj

st_geometry(mcnary_reservoir_proj)
rivers_proj %>% 
  filter(NAMEEN %in% c("Columbia River", "Snake River")) -> CS_rivers_proj

# get other major rivers from this DF
subset(rivers_proj, COUNTRY == "USA") -> USA_rivers_proj
unique(USA_rivers_proj$NAMEEN)

rivers_proj %>% 
  filter(NAMEEN %in% c("Fifteenmile Creek", "Deschutes River", 
                       "John Day River", "Umatilla River", "Yakima River", "Wenatchee River", "Entiat River",
                       "Okanogan River", "Methow River",
                       "Walla Walla River", "Tucannon River", "Clearwater River", "Grande Ronde River",
                       "Imnaha River", "Snake River", "Salmon River", "Asotin Creek")) -> major_tribs_proj

rivers_proj %>% 
  filter(NAMEEN %in% c("Clearwater River", "Salmon River",
                       "Deschutes River")) -> major_tribs_proj

unique(major_tribs_proj$NAMEEN)


### load Columbia River Basin streams
CRB_streams_spdf <- st_read(here::here("map_files", "crb_streams_over8m_100k","crb_streams_over8m_100k.shp"))
CRB_streams_proj <- sf::st_transform(CRB_streams_spdf, crs = WGS_CRS)

# # ok, so there are just missing stream segments here (at least segments that don't connect to mainstem Columbia or Snake)
# # we'll have to manually code those

# trim the streams spdf
CRB_stream_names <-  data.frame(unique(CRB_streams_proj$GNIS_NAME))
colnames(CRB_stream_names) <- c("GNIS_NAME")
# Streams we need: Columbia River, Fifteenmile Creek, Deschutes River, 
# John Day River, Umatilla River, Yakima River, Wenatchee River, Entiat River,
# Okanogan River, Methow River,
# Walla Walla River, Tucannon River, Clearwater River, Grande Ronde River,
# Imnaha River, Snake River, Salmon River, Asotin Creek

subset_streams_GNIS_name <- c(CRB_stream_names$GNIS_NAME[grep("Columbia", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Fifteenmile", CRB_stream_names$GNIS_NAME)],
                              # CRB_stream_names$GNIS_NAME[grep("Deschutes", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("John Day", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Umatilla", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Yakima", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Wenatchee", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Entiat", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Walla Walla", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Tucannon", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Clearwater", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Grande Ronde", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Imnaha", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Snake", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Salmon", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Bridge", CRB_stream_names$GNIS_NAME)],
                              # CRB_stream_names$GNIS_NAME[grep("White Creek", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Catherine", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Looking", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Lapwai", CRB_stream_names$GNIS_NAME)],
                              # CRB_stream_names$GNIS_NAME[grep("Klickitat", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Asotin", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Okanogan", CRB_stream_names$GNIS_NAME)],
                              CRB_stream_names$GNIS_NAME[grep("Methow", CRB_stream_names$GNIS_NAME)])
# CRB_stream_names$GNIS_NAME[grep("Trout", CRB_stream_names$GNIS_NAME)])

CRB_stream_names <-  data.frame(unique(CRB_streams_proj$HUCName))
colnames(CRB_stream_names) <- c("HUCName")
subset_streams <- c(
  # CRB_stream_names$HUCName[grep("Columbia", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Fifteenmile", CRB_stream_names$HUCName)],
  # CRB_stream_names$HUCName[grep("Deschutes", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("John Day", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Umatilla", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Yakima", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Wenatchee", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Entiat", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Walla Walla", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Tucannon", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Clearwater", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Grande Ronde", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Imnaha", CRB_stream_names$HUCName)],
  # CRB_stream_names$HUCName[grep("Snake", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Salmon River", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Bridge", CRB_stream_names$HUCName)],
  # CRB_stream_names$HUCName[grep("White Creek", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Catherine", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Looking", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Lapwai", CRB_stream_names$HUCName)],
  # CRB_stream_names$HUCName[grep("Klickitat", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Asotin", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Okanogan", CRB_stream_names$HUCName)],
  CRB_stream_names$HUCName[grep("Methow", CRB_stream_names$HUCName)])

# Remove the incorrectly selected ones
incorrect_streams <- c("Hoodoo Creek", "Little Wenatchee River", 
                       "Clearwater Creek", "Little Clearwater River", "Snake Creek", 
                       "No Snake Creek", "White Salmon River", "Salmon Falls Creek",
                       "Salmon Creek", "Little White Salmon River", "Little Salmon River",
                       "South Fork Salmon Falls Creek", "Little Salmon Creek",
                       "Big Salmon Creek", "Salmon la Sac Creek", "Trout Lake Creek")

subset_streams <- subset_streams[!(subset_streams %in% incorrect_streams)]

# Drop the ones that were incorrectly selected based on HUCName
# Also drop the north fork clearwater, since it's blocked by dworshak
HUCName_drop_streams <- c(CRB_stream_names$HUCName[grep("Frontal Columbia River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Willamette River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Lake Umatilla", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Little Wenatchee River", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Lake Entiat", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("North Fork Clearwater", CRB_stream_names$HUCName)],
                          CRB_stream_names$HUCName[grep("Hoodo", CRB_stream_names$HUCName)])

subset_streams <- subset_streams[!(subset_streams %in% HUCName_drop_streams)]


# CRB_streams_proj_subset <- subset(CRB_streams_proj, GNIS_NAME %in% subset_streams)
CRB_streams_proj_subset <- subset(CRB_streams_proj, HUCName %in% subset_streams)

# drop the Salmon River in Oregon
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(GNIS_NAME == "Salmon River" & ecoregion %in% c("Blue Mountains", "Cascades")))

# drop the Clearwater River in the Canadian Rockies
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Clearwater River" & ecoregion == "Canadian Rockies"))

# drop the weird Clearwater river section that's way off in Montana (?)
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Clearwater River" & UID > 2116900))

# drop weird fragments
# IDs keep changing so drop based on geography
CRB_streams_proj_subset <- subset(CRB_streams_proj_subset, !(HUCName == "Kachess River-Yakima River" & y_coord > 270000))

# create the missing river segments
# john day, connect to Columbia
# clearwater, connect to snake

# # Walla Walla missing segment
# wawa_stbbox <- st_bbox(c(xmin=-118.463298, xmax=-118.316699, ymin=45.937293, ymax=46.063931), crs = "WGS84")
# 
# CRB_streams_proj %>% 
#   st_crop(wawa_stbbox) -> wawa_missing

#### Make an inset map

inset_map <- ggplot(usa_spdf) +
  geom_sf(fill = "gray96") +
  geom_sf(data = one_CAN_spdf, fill = "gray96", color = "gray70") +
  # add the CRB layer
  geom_sf(data = CRB_boundary_proj, color = "gray85", linewidth = 0.2, fill = "gray85") +
  # add the state outlines
  geom_sf(data = us_states_proj, fill = "transparent", color = "gray70") + 
  # add the major rivers layer
  # geom_sf(data = major_rivers_proj, color = "gray70", linewidth = 0.2, fill = "gray70") +
  coord_sf(ylim = c(36,52.5),  xlim = c(-132, -95), expand = FALSE) +
  annotate(geom = "text", x = -105, y = 39, label = "United States\nof America",
           fontface = "italic", hjust = 0.5) + 
  annotate(geom = "text", x = -110, y = 51, label = "Canada",
           fontface = "italic") + 
  # remove axis text and labels
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # box for where our larger map is
  # annotate(geom = "rect", ymin = 44, ymax = 48.45,  
  #          xmin = -126, xmax = -113.23, fill = NA, color = "black")
  annotate(geom = "rect", ymin = 44.3, ymax = 48.55,  
           xmin = -126.6, xmax = -121, fill = NA, color = "black")





#### Load and reformat bathymetry data ####

bathy <- st_read(here("map_files", "bathymetry_shapefile", "NA_Bathymetry", "data", "bathy_line", "bathymetry_l_v2.shp"))
bathy_200 <- subset(bathy, DEPTH == 200)
bathy_200_proj <- sf::st_transform(bathy_200, crs = WGS_CRS)

##### MAKE BASE MAP #####
CRB_map <- ggplot(usa_spdf) +
  geom_sf(fill = "gray96") +
  geom_sf(data = one_CAN_spdf, fill = "gray96") +
  # add the CRB layer
  geom_sf(data = CRB_boundary_proj, color = "gray90", linewidth = 0.2, fill = "gray90") +
  # add the state outlines
  geom_sf(data = us_states_proj, fill = "transparent") + 
  # add the bathy 200 line
  geom_sf(data = bathy_200_proj, fill = "transparent", color = "gray60", linewidth = 0.2) + 
  ylab("Latitude")+
  xlab("Longitude")+
  # Lines for CRB streams
  # geom_sf(data = CRB_streams_proj_subset,  color = "gray30", linewidth = 0.5, fill = "gray30") +
  # Lines for North American Rivers (includes Columbia and Snake, and major tribs)
  # geom_sf(data = CS_rivers_proj, color = "gray30", linewidth = 2, fill = "gray30") +
  # geom_sf(data = estuary_linestring, color = "gray30", linewidth = 2, fill = "gray30") +
  # geom_sf(data = mcnary_reservoir_linestring, color = "gray30", linewidth = 2, fill = "gray30") +
  # geom_sf(data = major_tribs_proj, color = "gray30", linewidth = 0.5, fill = "gray30") +
  # geom_sf(data = john_day_mouth_linestring, color = "gray30", linewidth = 0.5, fill = "gray30") +
  # geom_sf(data = walla_walla_mouth_linestring, color = "gray30", linewidth = 0.5, fill = "gray30") +
  # coord_sf(ylim = c(44.3,48.45),  xlim = c(-125.15,-113.23)) +
  # coord_sf(ylim = c(44,48.45),  xlim = c(-126,-113.23)) +
  coord_sf(ylim = c(44.3,48.55),  xlim = c(-126.6,-121), expand = FALSE) +
  # ADD LABELS FOR RIVERS
  # annotate("segment", x = -122.84, y = 47.68, xend = -123.18, yend = 47.68, size = 0.5, lty = 1) + # Columbia River
  # annotate("text", x = -122.83, y = 46.15, label = "Columbia\nRiver", size = 4, fontface = 'italic', hjust = 0) + # Columbia River
  # annotate("text", x = -121.2, y = 45.44, label = "Fifteenmile Cr.", size = 3, fontface = 'italic', hjust = 1, angle = 45) + # Fifteenmile Cr.
  # annotate("text", x = -121.3, y = 44.8, label = "Deschutes R.", size = 3, fontface = 'italic', hjust = 1) + # Deschutes R.
  # annotate("text", x = -119.63, y = 44.92, label = "John Day R.", size = 3, fontface = 'italic', hjust = 1) + # John Day R.
  # annotate("text", x = -118.98, y = 45.60, label = "Umatilla R.", size = 3, fontface = 'italic', hjust = 1) + # Umatilla R.
  # annotate("text", x = -117.95, y = 46.15, label = "Walla Walla R.", size = 3, fontface = 'italic', hjust = 1) + # Walla Walla R.
  # annotate("text", x = -117.8, y = 46.366, label = "Tucannon R.", size = 3, fontface = 'italic', hjust = 1) + # Tucannon R.
  # annotate("text", x = -118, y = 45.50, label = "Grande Ronde R.", size = 3, fontface = 'italic', hjust = 1) + # Grande Ronde R.
  # annotate("text", x = -115.34, y = 46.475, label = "Clearwater R.", size = 3, fontface = 'italic', hjust = 1) + # Clearwater R.
  # annotate("text", x = -116.92, y = 45.52, label = "Imnaha R.", size = 3, fontface = 'italic', hjust = 1) + # Imnaha R.
  # annotate("text", x = -117.35, y = 44.6, label = "Snake\nRiver", size = 4, fontface = 'italic', hjust = 1) + # Snake R.
  # annotate("text", x = -115.6, y = 45.62, label = "Salmon R.", size = 3, fontface = 'italic', hjust = 1) + # Salmon R.
  # annotate("text", x = -120.44, y = 46.40, label = "Yakima R.", size = 3, fontface = 'italic', hjust = 1) + # Yakima R.
  # annotate("text", x = -120.5, y = 47.42, label = "Wenatchee R.", size = 3, fontface = 'italic', hjust = 1) + # Wenatchee R.
  # annotate("text", x = -120.6, y = 47.96, label = "Entiat R.", size = 3, fontface = 'italic', hjust = 1) + # Entiat R.
  # annotate("text", x = -120.263715, y = 48.341111, label = "Methow R.", size = 3, fontface = 'italic', hjust = 1) + # Methow R.
  # annotate("text", x = -119.434335, y = 48.468745, label = "Okanogan R.", size = 3, fontface = 'italic', hjust = 0) + # Okanogan R.
  # annotate("text", x = -117.038870, y = 46.159808, label = "Asotin Cr.", size = 3, fontface = 'italic', hjust = 1) + # Asotin Creek
  
  annotate("text", x = BON_coords[1], y = BON_coords[2], label = "|", size = 12) + # Bonneville Dam
  annotate("text", x = BON_coords[1], y = BON_coords[2]+0.35, label = "Bonneville\nDam", size = 6, hjust = 0.5) + # Bonneville Dam
  annotate("text", x = -123.33, y = 45.60, label = "Columbia River", size = 6, fontface = 'italic', hjust = 0, vjust = 0, angle = -18) + # Columbia River
  
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12))+
  guides(fill = guide_legend(title = "Legend")) +
  # add a north arrow
  # annotate("text", x = -124.9, y = 45.14, label = "N", fontface = "bold", size = 8) +
  # annotate(geom = "polygon", x = c(-125.1, -124.9, -124.9), y = c(44.6, 45, 44.7), 
  #          color = "black", fill = "gray50") +
  # annotate(geom = "polygon", x = c(-124.9, -124.9, -124.7), y = c(44.7, 45, 44.6), 
  #          color = "black", fill = "white")
  annotate("text", x = -124.9+3.3, y = 45.14+3.06,, label = "N", fontface = "bold", size = 8) +
  annotate(geom = "polygon", x = c(-125.1+3.3, -124.9+3.3, -124.9+3.3), y = c(44.6+3.06, 45+3.06, 44.7+3.06), 
           color = "black", fill = "gray50") +
  annotate(geom = "polygon", x = c(-124.9+3.3, -124.9+3.3, -124.7+3.3), y = c(44.7+3.06, 45+3.06, 44.6+3.06), 
           color = "black", fill = "white")

# save this as the base map 
# ggsave(here::here("figures", "CRB_base_map.pdf"), plot = CRB_map, height = 7.5, width  = 12.5)

#### Create the base map, add the inset

CRB_map +
  annotation_custom(grob = ggplotGrob(inset_map),
                    xmin = -123.2, xmax = -121,
                    ymin = 44.35, ymax = 45.2) -> CRB_map_plus_inset

  

# ggsave(here::here("figures", "CRB_base_map.pdf"), plot = CRB_map_plus_inset, height = 7.5, width  = 12.5)

#### process survey data ####

# extract Yearling Interior Chinook
csyif <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_fa")

csyif %>% 
  # filter(!(duplicated(station)))
  filter(year == 2013) %>% 
  mutate(survey = "JSOES") -> jsoes_stations

birds_long %>% 
  filter(species == "common_murre") %>% 
  filter(year == 2017) %>% 
  mutate(survey = "JSOES Seabirds")  -> bird_stations

rf %>% 
  filter(year == 2013) %>% 
  filter(Y >= 4905)%>% 
  mutate(survey = "PRS") -> prs_stations

hake %>% 
  filter(year == 2013) %>% 
  # filter(Lat >= 44) %>% 
  filter(Y > 4905) %>% 
  mutate(survey = "Hake") -> hake_stations

# group the hake data into transects
max(hake_stations$Y) -> hake_max_Y
min(hake_stations$Y) -> hake_min_Y

hake_stations %>% 
  mutate(transect = cut(Y, breaks = seq(hake_min_Y-1, hake_max_Y+1, length.out = 50), labels = 1:49)) -> hake_stations_new

# manually remove some stragglers
hake_stations_new %>% 
  filter(!(transect %in% c(4,6,10,42))) %>% 
  filter(!(transect == 7 & Y == max(subset(hake_stations_new, transect == 7)$Y))) %>% 
  filter(!(transect == 15 & Y == max(subset(hake_stations_new, transect == 15)$Y))) %>% 
  filter(!(transect == 17 & Y == min(subset(hake_stations_new, transect == 17)$Y))) %>% 
  filter(!(transect == 41 & Y == max(subset(hake_stations_new, transect == 41)$Y))) %>% 
  filter(!(transect == 49 & Y == min(subset(hake_stations_new, transect == 49)$Y))) -> hake_stations_new


dplyr::select(jsoes_stations, survey, X, Y) %>% 
  # bind_rows(dplyr::select(bird_stations, survey, X, Y)) %>% 
  bind_rows(dplyr::select(prs_stations, survey, X, Y)) %>% 
  bind_rows(dplyr::select(hake_stations_new, survey, X, Y, transect)) -> survey_data

survey_data$survey <- factor(survey_data$survey, levels = c("Hake", "PRS", "JSOES"))

survey_data <- arrange(survey_data, survey)

survey_data$survey <- factor(survey_data$survey, levels = c("JSOES", "PRS", "Hake"))

# convert this to WGS 84
survey_data %>% 
  mutate(X = X * 1000, Y = Y * 1000) %>% 
  st_as_sf(., coords = c("X", "Y"), crs = UTM_zone_10_crs) -> survey_data_sf
survey_data_sf_proj <- sf::st_transform(survey_data_sf, crs = WGS_CRS)

# turn this into a data frame
as.data.frame(st_coordinates(survey_data_sf_proj)) -> survey_data_coords

survey_data_for_map <- data.frame(Survey = survey_data_sf_proj$survey,
                                  transect = survey_data_sf_proj$transect,
                                  lon = survey_data_coords[,1],
                                  lat = survey_data_coords[,2])

survey_data_for_map <- arrange(survey_data_for_map, Survey)


#### Generate figure ####
survey_shapes <- c("JSOES" = 18,
                   "PRS" = 15,
                   "Hake" = 16)

survey_sizes <- c("JSOES" = 3.1,
                  "PRS" = 2.3,
                  "Hake" = 0.01)

# survey_colors <- c("JSOES" = "#440154",
#                    "PRS" = "#21918c",
#                    "Hake" = "#5ec962")

survey_colors <- c("JSOES" = "#a6cee3",
                   "PRS" = "#1f78b4",
                   "Hake" = "#b2df8a")


fig1 <- CRB_map_plus_inset +
  geom_line(data = subset(survey_data_for_map, Survey == "Hake"), aes(x = lon, y = lat, group = transect), color = "#b2df8a", linewidth = 0.8, show.legend = FALSE) +
  # geom_point(data = subset(survey_data_for_map, Survey != "Hake"), aes(x = lon, y = lat, color = Survey, shape = Survey, size = Survey)) +
  geom_point(data = survey_data_for_map, aes(x = lon, y = lat, color = Survey, shape = Survey, size = Survey)) +
  # scale_color_manual(values = survey_colors, name = "Survey", labels = c("JSOES", "PRS", "Hake")) +
  # scale_shape_manual(values = survey_shapes, name = "Survey", labels = c("JSOES", "PRS", "Hake")) +
  scale_color_manual(values = survey_colors) +
  scale_shape_manual(values = survey_shapes) +
  scale_size_manual(values = survey_sizes) +
  theme(legend.position.inside = c(0.1, 0.1)) +
  guides(size = guide_legend(override.aes = list(size = 3)))

ggsave(here::here("figures", "paper_figures", "fig1.pdf"), plot = fig1, height = 7.5, width  = 7.5)
