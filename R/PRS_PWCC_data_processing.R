# PRS + PWCC data processing

## load libraries
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
library(RODBC)

#### load and reformat the PWCC data ####
## Load "Haul" sheet
PWCC_MWT_haul <- clean_names(read_xlsx(here::here("Data", "PWCC.dataforMarcus.xlsx"), sheet = "PWCC.haul"))

# we're missing one observation for longitude
# PWCC_MWT_haul %>% 
#   filter(is.na(net_fishing_long) | is.na(net_fishing_lat))

# however, we do have a net back longitude for that tow; let's use that as the fishing observation,
# since that'll be close enough
PWCC_MWT_haul %>% 
  mutate(net_fishing_long = ifelse(is.na(net_fishing_long), net_back_long, net_fishing_long)) -> PWCC_MWT_haul


# lat and lon are in  ddmm.mm where d=degrees and m=minutes (e.g., 3642.5 is 36 degrees 42.5 minutes)

PWCC_MWT_haul %>% 
  mutate(net_fishing_lat_deg = substr(net_fishing_lat, 1, 2),
         net_fishing_lat_min = substr(net_fishing_lat, 3, 7)) %>% 
  mutate(net_fishing_long_deg = substr(net_fishing_long, 1, 3),
         net_fishing_long_min = substr(net_fishing_long, 4, 8)) %>% 
  mutate(net_fishing_lat_dd = as.numeric(as.numeric(net_fishing_lat_deg) + as.numeric(net_fishing_lat_min)/60)) %>% 
  mutate(net_fishing_long_dd = -1 * as.numeric(as.numeric(net_fishing_long_deg) + as.numeric(net_fishing_long_min)/60)) -> PWCC_MWT_haul

## Load "Catch" sheet
PWCC_MWT_catch <- clean_names(read_xlsx(here::here("Data", "PWCC.dataforMarcus.xlsx"), sheet = "PWCC.catch"))

## Load "species.codes" sheet
PWCC_MWT_species_codes <- clean_names(read_xlsx(here::here("Data", "PWCC.dataforMarcus.xlsx"), sheet = "species.codes"))

## expand the catch data so that there are zeros for every taxon
# we want different life stages of the same species to be classified as different species
# the final output should have this many rows:
# PWCC re-used haul numbers by cruise
# length(unique(paste0(PWCC_MWT_catch$species, " - ", PWCC_MWT_catch$maturity))) * length(unique(paste0(PWCC_MWT_catch$cruise, "-", PWCC_MWT_catch$haul_no)))
# table(PWCC_MWT_catch$cruise, PWCC_MWT_catch$haul_no)

PWCC_MWT_catch %>% 
  complete(nesting(cruise, haul_no), nesting(species, maturity),
           fill = list(total_no = 0)) -> PWCC


# join catch and species codes
PWCC %>% 
  left_join(PWCC_MWT_species_codes, by = "species") -> PWCC

# They're off by 35
# some hauls just have more data
# sort(-table(paste0(PWCC$cruise, PWCC$haul_no)))
# 
# subset(PWCC, cruise == "PWCC02" & haul_no == "5")$common_name == subset(PWCC, cruise == "PWCC01" & haul_no == "35")$common_name
# subset(PWCC, cruise == "PWCC02" & haul_no == "5")
# subset(PWCC, cruise == "PWCC01" & haul_no == "35")
# # in this case, it looks like maybe there are duplicate observations in cruise PWCC02 haul_no 5?
# # see observations for barracudina and dover sole
# 
# # let's look across all tows for duplicate observations
# PWCC$possible_duplicate <- duplicated(PWCC)
# subset(PWCC, possible_duplicate == TRUE)
# 
# # check one
# subset(PWCC, cruise == "PWCC01" & haul_no == "76")
# # yep, there are some duplicates. But this doesn't explain all of the extra observations.
# 
# # Let's remove the duplicates and look again at which are duplicated
# PWCC_nodup <- subset(PWCC, possible_duplicate == FALSE)
# sort(-table(paste0(PWCC_nodup$cruise, "-", PWCC_nodup$haul_no)))
# subset(PWCC, cruise == "PWCC01" & haul_no == "77")
# subset(PWCC, cruise == "PWCC01" & haul_no == "77")$common_name == subset(PWCC, cruise == "PWCC01" & haul_no == "35")$common_name
# subset(PWCC, cruise == "PWCC01" & haul_no == "77")
# subset(PWCC, cruise == "PWCC01" & haul_no == "35")
# # ok here it looks like we have the same taxon entered twice, but different numbers.
# # I think that maybe then we don't actually have duplicates, but just double entries for
# # some taxa. Perhaps someone missed one and then entered it at the end.

# Let's collapse by species + maturity within a haul
PWCC %>% 
  group_by(cruise, haul_no, species, maturity, common_name, sci_name, maturity_codes, species_group, notes) %>% 
  summarise_if(is.numeric, sum, ) %>% 
  ungroup() -> PWCC

# join catch and haul info
PWCC %>% 
  left_join(PWCC_MWT_haul, by = c("cruise", "haul_no")) -> PWCC_full

##### load and reformat the PRS data #####
## Load "Haul" sheet
PRS_MWT_haul <- clean_names(read_xls(here::here("Data", "PreRecruit MWT Master.xls"), sheet = "Haul"))

## Load "Catch" sheet
# for the catch data, we need to manually set one column to text
PRS_MWT_catch_column_types <- c("text", "numeric", "numeric", "numeric", "numeric", "numeric",
                                "text", "text", "numeric", "numeric", "numeric", "numeric",
                                "numeric", "text",  "text", "numeric", "text")
PRS_MWT_catch <- clean_names(read_xls(here::here("Data", "PreRecruit MWT Master.xls"), 
                                      sheet = "Catch", col_types = PRS_MWT_catch_column_types))

## Load "Krill Catch" sheet
PRS_MWT_krill_catch <- clean_names(read_xls(here::here("Data", "PreRecruit MWT Master.xls"), sheet = "Krill Catch"))

## Load "Krill Number"
PRS_MWT_krill_number <- clean_names(read_xls(here::here("Data", "PreRecruit MWT Master.xls"), sheet = "Krill Number"))

## Transform data
# For both haul and catch, a couple of transect_new IDs got converted to character weird, so let's manually fix them
PRS_MWT_haul %>% 
  mutate(transect_new = ifelse(transect_new == "46.200000000000003", "46.2", transect_new)) -> PRS_MWT_haul

PRS_MWT_catch %>% 
  mutate(transect_new = ifelse(transect_new == "46.200000000000003", "46.2", transect_new)) -> PRS_MWT_catch


### create a unique identifier for each haul in both the haul and the catch data
PRS_MWT_haul %>% 
  mutate(haul_id = paste0(year, "_", month, "_", day, "-", transect_new, "_", station_new)) %>% 
  relocate(haul_id) -> PRS_MWT_haul

# add a unique identifier for those deep/middle/shallow tows in 2016
PRS_MWT_haul %>% 
  mutate(haul_id = ifelse(grepl("Deep|Shallow|Middle", original_station),
                          paste0(year, "_", month, "_", day, "-", transect_new, "_", station_new, "_", original_station), haul_id)) -> PRS_MWT_haul

PRS_MWT_catch %>% 
  mutate(haul_id = paste0(year, "_", month, "_", day, "-", transect_new, "_", station_new)) %>% 
  relocate(haul_id) -> PRS_MWT_catch

# add a unique identifier for those deep/middle/shallow tows in 2016
PRS_MWT_catch %>% 
  mutate(haul_id = ifelse(grepl("Deep|Shallow|Middle", original_station),
                          paste0(year, "_", month, "_", day, "-", transect_new, "_", station_new, "_", original_station), haul_id)) -> PRS_MWT_catch

## expand the catch data so that there are zeros for every taxon
# we want different life stages of the same species to be classified as different species
# the final output should have this many rows:
length(unique(paste0(PRS_MWT_catch$taxon, " - ", PRS_MWT_catch$maturity))) * length(unique(PRS_MWT_catch$haul_id))
# Our output actually has one extra row; this is because in one tow, there are two entries for Hydrozoa
# with different comments for different species


PRS_MWT_catch %>% 
  complete(nesting(haul_id, cruise, year, month, day, original_haul_number, 
                   new_haul_number, original_station, transect_new, station_new, 
                   swfcs_station, swfcs_transect, distance_from_shore_km,
                   start_depth_m), nesting(taxon, maturity),
           fill = list(number = 0)) -> PRS_MWT_catch

# join the haul information with the catch data
PRS_MWT_catch %>% 
  left_join(dplyr::select(PRS_MWT_haul, c(setdiff(colnames(PRS_MWT_haul),
                                                  colnames(PRS_MWT_catch)), haul_id)), by = "haul_id") -> PRS_full

## deal with comments
# there are numerous comments throughout the data sheets.
# for our purposes, we don't need to act on most of these.




#### Trim both PRS and PWCC to only include north of 44 N (JSOES overlap) ####

# PRS_full %>% 
#   filter(start_latitude >= 44) -> PRS
# 
# PWCC_full %>% 
#   filter(net_fishing_lat >= 44) -> PWCC

# let's not do this - we can use the full extent to fit an SDM, and then project it north to include the JSOES domain.
# we can also always trim it if we decide to

PRS_full -> PRS
PWCC_full -> PWCC


#### Convert PRS to sf object ####
# turn this into an sf object
st_as_sf(PRS, coords = c("start_longitude", "start_latitude"), crs = 4326) -> PRS_sf

# change CRS to UTM zone 10 (to work in meters)
UTM_zone_10_crs <- 32610
PRS_sf_proj <- sf::st_transform(PRS_sf, crs = UTM_zone_10_crs)

# make this projection into kilometers to help with interpretability
PRS_sf_proj_km <- st_as_sf(PRS_sf_proj$geometry/1000, crs = UTM_zone_10_crs)

# extract geometry
as.data.frame(st_coordinates(PRS_sf_proj_km)) -> PRS_km
# add this back to jsoes_long (X and Y now represent eastings and northings in km)
bind_cols(PRS, PRS_km) -> PRS

#### Convert PWCC to sf object ####

# turn this into an sf object
st_as_sf(PWCC, coords = c("net_fishing_long_dd", "net_fishing_lat_dd"), crs = 4326) -> PWCC_sf

# change CRS to UTM zone 10 (to work in meters)
UTM_zone_10_crs <- 32610
PWCC_sf_proj <- sf::st_transform(PWCC_sf, crs = UTM_zone_10_crs)

# make this projection into kilometers to help with interpretability
PWCC_sf_proj_km <- st_as_sf(PWCC_sf_proj$geometry/1000, crs = UTM_zone_10_crs)

# extract geometry
as.data.frame(st_coordinates(PWCC_sf_proj_km)) -> PWCC_km
# add this back to jsoes_long (X and Y now represent eastings and northings in km)
bind_cols(PWCC, PWCC_km) -> PWCC

#### Explore common PRS taxa ####
## total density
PRS %>% 
  group_by(taxon) %>% 
  summarise(mean_number = mean(number),
            sd_number = sd(number)) %>% 
  arrange(desc(mean_number)) -> PRS_density_sums

## frequency of occurrence

# get count of each genus_species
PRS %>% 
  filter(number > 0) %>% 
  count(taxon) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_samples = n/length(unique(PRS$haul_id))) -> PRS_freq_occur


#### PWCC ####
## total density
PWCC %>% 
  group_by(common_name) %>% 
  summarise(mean_number = mean(total_no),
            sd_number = sd(total_no)) %>% 
  arrange(desc(mean_number)) -> PWCC_density_sums


## frequency of occurrence

# get count of each genus_species
PWCC %>% 
  filter(total_no > 0) %>% 
  count(common_name) %>% 
  arrange(desc(n)) %>% 
  mutate(prop_samples = n/length(unique(paste0(PWCC$cruise, PWCC$haul_no)))) -> PWCC_freq_occur

#### Spatial data ####
usa_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/USA_adm0.shp", quiet = TRUE)
# load BC
CAN_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/canada/lpr_000b16a_e.shp", quiet = TRUE)
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             c(xmin = -126, ymin = 40, xmax = -123, ymax = 48.5))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -126, ymin = 40, xmax = -123, ymax = 48.5))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)

# create base map
survey_area_basemap_km <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = BC_coast_proj_km) +
  # coord_sf(ylim = c(4922.052, 5342.052), xlim = c(334.8638, 404.8638)) +
  # coord_sf(ylim = c(44,48.5),  xlim = c(-126, -123)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.text = element_blank())
# temporary fix for lat/long on axes to just get rid of them - start here to actually fix: https://forum.posit.co/t/converting-axes-to-lat-lon/27181

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

# compare PRS and PWCC stations
ggplot(US_west_coast) +
  geom_sf() +
  geom_sf(data = BC_coast) +
  # coord_sf(ylim = c(33,48.5),  xlim = c(-126, -118)) +
  # # scale_x_continuous(breaks = c(124,125,126)) +
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
        legend.text = element_text(size = 12)) +
  # add PRS stations
  geom_point(data = PRS, aes(x = start_longitude, y = start_latitude), color = "blue") + 
  # add PWCC stations
  geom_point(data = PWCC, aes(x = net_fishing_long_dd, y = net_fishing_lat_dd), color = "red")

# ok yep, this makes sense. PWCC is closer to shore though? Some samples look basically onshore?
# If you look at the map generated in Field et al. 2021, this matches - PWCC is generally closer to shore
# than PRS, though perhaps not by this much.







## Plotting distributions of some common taxa

plot_distribution_PRS <- function(data, taxon_name){
  # keep only this taxon
  taxon_data <- subset(data, taxon == taxon_name)
  
  # create facet_wrap plot for distribution across all years
  taxon_data %>% 
    mutate(encounter = ifelse(number == 0, "zero", "non-zero")) -> taxon_data
  
  survey_area_basemap +
    geom_point(data = taxon_data, aes(x = start_longitude, y = start_latitude, size = number, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    facet_wrap(~year, nrow = 2) +
    ggtitle(taxon_name) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}

sanddab_dist_plot <- plot_distribution_PRS(data = PRS, taxon_name = "Citharichthys sordidus")

euphausiidae_dist_plot <- plot_distribution_PRS(data = PRS, taxon_name = "Euphausiidae")

## In situ environmental data

# Except for 2013, the following environmental data was collected in situ by the PRS:
#   
#   - Surface Temp (C)
# - Surface temp was also collected in 2013
# - Temperature (30m C)	
# - Fluorescence (30m V)	
# - Salinity (30m PSU)	
# - Density (30m sigma-theta, kg/m3)	
# - Dissolved Oxygen (30m ml/l)

# examine SST vs temperature at 30 m

ggplot(PRS_MWT_haul, aes(x = surface_temp_o_c, y = temperature_30m_o_c)) +
  geom_point() +
  ggtitle("Temperature at 30 m vs. SST")


## Krill investigation

# Krill get special treatment in the PRS data, with them being recorded in three separate sheets. These three sheets all sum to the same total number.

ggplot(subset(PRS, taxon == "Euphausiidae"), aes(x = year, y = number)) +
  geom_bar(stat = "identity") +
  ggtitle("Krill from catch sheet")

ggplot(PRS_MWT_krill_catch, aes(x = year, y = number)) +
  geom_bar(stat = "identity") +
  ggtitle("Krill from krill catch sheet")

PRS_MWT_krill_number$total_krill <- rowSums(PRS_MWT_krill_number[,14:20])

ggplot(PRS_MWT_krill_number, aes(x = year, y = total_krill)) +
  geom_bar(stat = "identity") +
  ggtitle("Krill from krill number sheet")

krill_sums <- data.frame(year = unique(PRS$year),
                         catch_sum = rep(NA, length(unique(PRS$year))),
                         krill_catch_sum = rep(NA, length(unique(PRS$year))),
                         krill_number_sum = rep(NA, length(unique(PRS$year))))

for (i in 1:length(krill_sums$year)){
  krill_sums$catch_sum[i] <- sum(subset(PRS, taxon == "Euphausiidae" & year == krill_sums$year[i])$number)
  
  krill_sums$krill_catch_sum[i] <- sum(subset(PRS_MWT_krill_catch, year == krill_sums$year[i])$number)
  
  krill_sums$krill_number_sum[i] <- sum(subset(PRS_MWT_krill_number, year == krill_sums$year[i])$total_krill)
}

#### Process species into taxa groups ####

## Selecting taxa to model


# Processing taxa - there are many notes from the spreadsheet that we will use to help us interpret taxon IDs.
# "*Did not speciate Gonatus after 2019."

# - There are four members of family Gonatidae in the dataset: "Gonatidae", "Gonatopsis borealis", "Gonatus onyx", and "Gonatus spp."
# - There is a clear break in species IDs before and after 2020. What is odd though is that only one member of the genus Gonatus was ever identified before 2019, so it's confusing that they would note that they stopped identifying to species level?

subset(PRS, taxon %in% c("Gonatidae", "Gonatopsis borealis", "Gonatus onyx", "Gonatus spp.")) -> PRS_gonatidae

ggplot(PRS_gonatidae, aes(x = year, y = number, fill = taxon)) +
  geom_bar(stat = "identity")



# Our criteria for generating SDMs are 1) data availability and 2) theoretical interactions with juvenile salmon.

# Friedman et al. (2018) identified 31 salmon forage species. These taxa were Armhook squid, Arrowtooth flounder, Barracudina, Blacktip squid, Combfish, Dover sole, Lingcod, Market squid, Myctophids, Northern anchovy, Octopus, Pacific hake, Pacific sand lance, Pacific sardine, Pacific tomcod, Painted greenling, Pandalus shrimp, Poacher, Rex sole, Rockfish, Ronquil / prickleback, Sand sole, Sanddab, Sculpin, Sergestid, Shrimp, Slender sole, Smelt, Snailfish, Turbot, and Krill.

# Daly et al. (2017) identified common diets in yearling Chinook salmon caught in the JSOES survey. The main taxa were Northern Anchovy, Clupeids, Cottids, Flatfish, Osmerids, Ronquils, Pacific s and lance, rockfishes, *Cancer* larvae, Euphausiids, Hyperiid amphipods, Non-*Cancer* crab larvae, and shrimp larvae.



friedman_2018_taxa <- c(
  "Gonatus spp.",
  "Atheresthes stomias",
  "Paralepididae",
  "Abraliopsis felis",
  "Zaniolepididae",
  "Microstomus pacificus",
  "Ophiodon elongatus",
  "Doryteuthis opalescens",
  "Myctophidae",
  "Engraulis mordax",
  "Octopoda",
  "Merluccius productus",
  "Ammodytes hexapterus",
  "Sardinops sagax",
  "Microgadus proximus",
  "Oxylebius pictus",
  "Pandalus jordani",
  "Agonidae",
  "Glyptocephalus zachirus",
  "Sebastes spp.",
  "Ronquilus / Stichaeidae",
  "Psettichthys melanostictus",
  "Citharichthys spp.",
  "Cottidae",
  "Sergestidae",
  "Natantia",
  "Lyopsetta exilis",
  "Osmeridae",
  "Liparidae",
  "Pleuronichthys spp.",
  "Euphausiidae")

# which taxa don't have exact matches?
  setdiff(friedman_2018_taxa, PRS_freq_occur$taxon)
sort(PRS_freq_occur$taxon)

subset(PRS_freq_occur, taxon %in% friedman_2018_taxa)

PRS_friedman_2018_taxa <- c(friedman_2018_taxa, "Lestidiops ringens", "Oxylebius pictus",
                            "Ammodytes spp.", "Ronquilus jordani", "Plectobranchus evides",
                            "Citharichthys sordidus", "Citharichthys stigmaeus",
                            "Liparis fucensis", "Liparis spp.",
                            # all the rockfishes
                            "Sebastes (Black/Yellowtail group)", "Sebastes (Copper group)", "Sebastes (Rosy group)", 
                            "Sebastes alutus", "Sebastes auriculatus", "Sebastes babcocki", "Sebastes crameri", 
                            "Sebastes diploproa", "Sebastes entomelas", "Sebastes flavidus", "Sebastes goodei", 
                            "Sebastes jordani", "Sebastes melanops", "Sebastes mystinus", "Sebastes paucispinis", 
                            "Sebastes pinniger", "Sebastes reedi", "Sebastes ruberrimus", "Sebastes rufus", 
                            "Sebastes saxicola", "Sebastes spp.", "Sebastes wilsoni", "Sebastes zacentrus", 
                            "Sebastolobus spp.")

# add on the ones that are members of higher order groups in the Friedman taxa
PRS_friedman_2018_taxa <- c(PRS_friedman_2018_taxa, "Tarletonbeania crenularis", "Stenobrachius leucopsarus", "Diaphus theta", "Nannobrachium regale", # myctophids!
                            "Lipolagus ochotensis", "Bathylagidae", "Leuroglossus stilbius", "Leuroglossus schmidti", # Bathylagidae - not technically part of the myctophid family, but functionally similar
                            "Gonatus onyx", "Gonatopsis borealis", "Gonatidae", # squid
                            "Hemilepidotus spp.", "Agonopsis vulsa", "Scorpaenichthys marmoratus", "Bathyagonus nigripinnis", "Odontopyxis trispinosa", "Xeneretmus latifrons", "Hemilepidotus hemilepidotus", "Hemilepidotus spinosus", # Agonidae/cottids
                            "Crangon spp.", "Streetsia challengeri", "Mysidacea", # shrimps
                            "Cancer spp." # crabs
)

# add on the ones that aren't included in Friedman but are listed in Daly
PRS_Daly_2017_taxa <- c("Phronima spp.", "Clupea pallasii", "Isopsetta isolepis", "Cancer magister", "Amphipoda",
                        "Allosmerus elongatus", "Thaleichthys pacificus", "Parophrys vetulus", "Hippoglossoides elassodon",
                        "Pleuronichthys decurrens", "Pleuronectidae", "Lepidopsetta bilineata", # lots of flatfishes
                        "Alosa sapidissima", # Clupeids
                        "Embassichthys bathybius", # more flatfishes
                        "Argentinoidei", "Hypomesus pretiosus") # Osmerids

# all important taxa

salmon_forage <- c(PRS_friedman_2018_taxa, PRS_Daly_2017_taxa)

# note what groups each species belongs to
salmon_forage_groups <- data.frame(taxon = salmon_forage,
                                   group = c("armhook squid", "flatfishes", "barracudina",
                                             "blacktip squid", "combfish", "flatfishes",
                                             "lingcod", "market squid", "myctophids",
                                             "northern anchovy", "octpus", "hake",
                                             "Pacific sand lance", "Pacific sardine", "Pacific tomcod",
                                             "combfish", "Pandalus shrimp", "poachers",
                                             "flatfishes", "rockfishes", "ronquil/prickleback",
                                             "flatfishes", "flatfishes", "sculpins",
                                             "sergestid shrimp", "shrimp", "flatfishes",
                                             "osmerids", "snailfishes", "flatfishes",
                                             "krill", "barracudina", "combfish",
                                             "Pacific sand lance", "ronquil/prickleback", "ronquil/prickleback",
                                             "flatfishes", "flatfishes", "snailfishes",
                                             "snailfishes", rep("rockfishes", 24), "myctophids", "myctophids", 
                                             "myctophids", "myctophids", "bathylagids", 
                                             "bathylagids", "bathylagids", "bathylagids",
                                             "armhook squid", "armhook squid", "armhook squid",
                                             "poachers", "poachers", "sculpins",
                                             "poachers", "poachers", "poachers",
                                             "poachers", "poachers", "crangon shrimp",
                                             "amphipods", "mysid shrimp", "Cancer crabs",
                                             "amphipods", "clupeids", "flatfishes",
                                             "Cancer crabs", "amphipods", "osmerids",
                                             "osmerids", "flatfishes", "flatfishes",
                                             "flatfishes", "flatfishes", "flatfishes",
                                             "clupeids", "flatfishes", "osmerids",
                                             "osmerids")
)



PRS_freq_occur %>% 
  left_join(salmon_forage_groups, by = "taxon") -> PRS_freq_occur

#### Common taxa PWCC ####
PWCC_freq_occur$common_name
PRS_freq_occur$taxon

# find matches for taxa between PWCC and PRS
PWCC_PRS_key <- data.frame(
  common_name = c(
    "NORTH PACIFIC HAKE", "SQUID", "LANTERNFISH", "NORTHERN ANCHOVY", "PACIFIC SARDINE",
    "REX SOLE", "BUTTER SOLE", "ARROWTOOTH FLOUNDER", "CANARY ROCKFISH", "SLENDER SOLE",
    "BARRACUDINA", "SMELT", "ROCKFISH", "WIDOW ROCKFISH", "DOVER SOLE", "CABEZON",
    "DARKBLOTCHED ROCKFISH", "YELLOWTAIL ROCKFISH", "GREENLING", "KING-OF-THE-SALMON",
    "SCULPIN", "NIGHT SMELT", "PACIFIC HERRING", "BLUE ROCKFISH", "PACIFIC SANDLANCE",
    "POACHER", "BLACK ROCKFISH", "MARKET SQUID", "PACIFIC SPINY DOGFISH", "PACIFIC SANDDAB",
    "PACIFIC TOMCOD", "SNAILFISH", "WHITEBAIT SMELT", "SHORTBELLY ROCKFISH", "COD",
    "HUMBOLDT SQUID", "SALMON", "SANDDAB", "BOCACCIO", "LINGCOD", "MEDUSAFISH",
    "LAMPREY", "TURBOT", "OCTOPUS", "YELLOWMOUTH ROCKFISH", "CHILIPEPPER", "CHUM SALMON",
    "EELPOUT", "PACIFIC COD", "PACIFIC SAND SOLE", "SABLEFISH", "STRIPETAIL ROCKFISH",
    "AMERICAN SHAD", "CURLFIN SOLE", "DEEP-SEA SMELT", "FISH", "KING SALMON",
    "PACIFIC MACKEREL", "PACIFIC SAURY", "PAINTED GREENLING", "PETRALE SOLE",
    "PLAINFIN MIDSHIPMAN", "RONQUIL", "SHRIMP", "SKATE"
  ),
  taxon = c(
    "Merluccius productus", "Cephalopoda", "Myctophidae", "Engraulis mordax", "Sardinops sagax",
    "Glyptocephalus zachirus", "Isopsetta isolepis", "Atheresthes stomias", "Sebastes pinniger", "Lyopsetta exilis",
    "Lestidiops ringens", "Osmeridae", "Sebastes spp.", "Sebastes entomelas", "Microstomus pacificus", "Scorpaenichthys marmoratus",
    "Sebastes crameri", "Sebastes flavidus", "Hexagrammidae", "Trachipterus altivelis",
    "Cottidae", "Spirinchus starksi", "Clupea pallasii", "Sebastes mystinus", "Ammodytes spp.",
    "Agonidae", "Sebastes melanops", "Doryteuthis opalescens", "Squalus suckleyi", "Citharichthys sordidus",
    "Microgadus proximus", "Liparis spp.", "Allosmerus elongatus", "Sebastes jordani", "Gadidae",
    "Dosidicus gigas", "Oncorhynchus spp.", "Citharichthys sordidus", "Sebastes paucispinis", "Ophiodon elongatus", "Icichthys lockingtoni",
    "Petromyzontidae", "Pleuronichthys decurrens", "Octopoda", "Sebastes reedi", "Sebastes goodei", "Oncorhynchus keta",
    "Liparidae", "Gadus macrocephalus", "Psettichthys melanostictus", "Anoplopoma fimbria", "Sebastes zacentrus",
    "Alosa sapidissima", "Pleuronichthys decurrens", "Bathylagidae", "Actinopterygii", "Oncorhynchus tshawytscha",
    "Scomber japonicus", "Cololabis saira", "Oxylebius pictus", "Eopsetta jordani",
    "Porichthys notatus", "Ronquilus jordani", "Natantia", "Rajidae"
  ),
  group = c("hake", "squid", "myctophids", "northern anchovy", "Pacific sardine",
            "flatfishes", "flatfishes", "flatfishes", "rockfishes", "flatfishes",
            "barracudina", "osmerids", "rockfishes", "rockfishes", "flatfishes", "sculpins",
            "rockfishes", "rockfishes", NA, NA,
            "sculpins", "osmerids", "clupeids", "rockfishes", "Pacific sandlance",
            "poachers", "rockfishes", "market squid", NA, "flatfishes",
            "Pacific tomcod", "snailfishes", "osmerids", "rockfishes", "cod",
            NA, NA, "flatfishes", "rockfishes", "lingcod", NA, 
            NA, "flatfishes", "octopus", "rockfishes", "rockfishes", NA,
            NA, "cod", "flatfishes", NA, "rockfishes",
            "osmerids", "flatfishes", "bathylagids", NA, NA,
            NA, NA, NA, "flatfishes",
            NA, "ronquil/prickleback", "shrimp", NA)
)

PWCC_freq_occur %>% 
  dplyr::rename(PWCC_prop_samples = prop_samples) %>% 
  left_join(PWCC_PRS_key, by = "common_name") %>% 
  left_join(dplyr::rename(PRS_freq_occur, PRS_prop_samples = prop_samples), by = "taxon") -> PWCC_PRS_freq_occur

salmon_forage_groups

# Observations of PWCC species identifications:
#   
#   - There aren't any gelatinous zooplankton IDs. There are also aren't any krill (which we know they must have caught). Probably were focused
# on identifying squids and fishes
# - John Field noted that they were hake fishermen, so the IDs for many taxa are spotty. But if what we care about is characterizing
# the overall forage base, and if "flatfishes" is good enough for our purposes, then these IDs should be fine.
# - There are some IDs that are very coarse in the PWCC data, for example "squid". They also ID "market squid" separately. So perhaps
# it would be safe to assume "squid" are all squid that aren't market squid. But maybe there are temporal patterns in species IDs?
#   - In the code below, there's nothing that stands out drastically


# plot usage of species IDs throughout time for PWCC data

plot_PWCC_time_series <- function(data, common_name_select){
  plot_data <- subset(data, common_name == common_name_select)
  plot_data %>% 
    mutate(year = year(net_fishing_time)) -> plot_data
  
  plot <- ggplot(plot_data, aes(x = year, y = total_no)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("\"", common_name_select, "\"")) +
    scale_x_continuous(breaks = 2001:2009)
  
  return(plot)
}


PWCC_taxa <- unique(PWCC$common_name)

for (i in 1:length(PWCC_taxa)){
  
  out_plot <- plot_PWCC_time_series(data = PWCC, common_name_select = PWCC_taxa[i])
  print(out_plot)
}



## Exploration of broad groups

# Here, we will look at patterns in functional groups of salmon forage taxa (e.g., rockfishes, flatfishes, etc.)

PRS %>% 
  left_join(salmon_forage_groups, by = "taxon") -> PRS

PWCC %>% 
  left_join(dplyr::select(PWCC_PRS_key, common_name, group), by = "common_name") -> PWCC

# make a function to plot by group
plot_PRS_group_distribution <- function(data, group_name){
  # keep only this taxon
  group_data <- subset(data, group == group_name)
  
  # collapse by this group - also collapse all maturity stages.
  # we may not want to do this in the future but for now it's okay
  group_data %>%
    group_by(across(c(-taxon, -maturity, -number))) %>%
    summarise(total = sum(number)) -> group_data
  
  # create facet_wrap plot for distribution across all years
  group_data %>% 
    mutate(encounter = ifelse(total == 0, "zero", "non-zero")) -> group_data
  
  survey_area_basemap +
    geom_point(data = group_data, aes(x = start_longitude, y = start_latitude, size = total, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    facet_wrap(~year, nrow = 2) +
    ggtitle(group_name) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}

plot_PWCC_group_distribution <- function(data, group_name){
  # keep only this taxon
  group_data <- subset(data, group == group_name)
  
  # collapse by this group - also collapse all maturity stages.
  # we may not want to do this in the future but for now it's okay
  group_data %>%
    group_by(across(c(-common_name, -sci_name, -species, -maturity, -total_no))) %>%
    summarise(total = sum(total_no)) -> group_data
  
  # create facet_wrap plot for distribution across all years
  group_data %>% 
    mutate(encounter = ifelse(total == 0, "zero", "non-zero")) %>% 
    mutate(year = year(net_fishing_time)) -> group_data
  
  survey_area_basemap +
    geom_point(data = group_data, aes(x = net_fishing_long_dd, y = net_fishing_lat_dd, size = total, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    facet_wrap(~year, nrow = 2) +
    ggtitle(group_name) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}

# plot the distributions of all of our core forage taxa
sort(table(salmon_forage_groups$group))
# let's plot rockfishes, flatfishes, poachers, osmerids, myctophids, bathylagids, squids, sand lance, clupeids, and cancer crabs

PRS_groups <- c("rockfishes", "flatfishes", "poachers", "osmerids", "myctophids", "bathylagids", "armhook squid", "market squid", "Pacific sand lance", "clupeids", "Cancer crabs")

for (i in 1:length(PRS_groups)){
  print(plot_PRS_group_distribution(data = PRS, group_name = PRS_groups[i]))
}

sort(table(PWCC_PRS_key$group))

PWCC_groups <- c("rockfishes", "flatfishes", "poachers", "osmerids", "myctophids", "bathylagids", "squid", "market squid", "Pacific sandlance", "clupeids")

for (i in 1:length(PWCC_groups)){
  print(plot_PWCC_group_distribution(data = PWCC, group_name = PWCC_groups[i]))
}

# Export PRS group data
export_PRS_group_data <- function(data, group_name, export_path){
  # keep only this taxon
  group_data <- subset(data, group == group_name)
  
  # collapse by this group - also collapse all maturity stages.
  # we may not want to do this in the future but for now it's okay
  group_data %>%
    group_by(across(c(-taxon, -maturity, -number))) %>%
    summarise(total = sum(number)) -> group_data
  
  write.csv(group_data, export_path, row.names = FALSE)
}

export_PWCC_group_data <- function(data, group_name, export_path){
  # keep only this taxon
  group_data <- subset(data, group == group_name)
  
  # collapse by this group - also collapse all maturity stages.
  # we may not want to do this in the future but for now it's okay
  group_data %>%
    group_by(across(c(-common_name, -sci_name, -species, -maturity, -total_no))) %>%
    summarise(total = sum(total_no)) -> group_data
  
  write.csv(group_data, export_path, row.names = FALSE)
  
}

for (i in 1:length(PRS_groups)){
  export_PRS_group_data(data = PRS, group_name = PRS_groups[i], export_path = here::here("model_inputs", paste0("PRS_", PRS_groups[i], ".csv")))
}

for (i in 1:length(PWCC_groups)){
  export_PWCC_group_data(data = PWCC, group_name = PWCC_groups[i], export_path = here::here("model_inputs", paste0("PWCC_", PWCC_groups[i], ".csv")))
}