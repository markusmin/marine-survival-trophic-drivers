## 11_Fig3_covariate_time_series

# Description: This script generates Figure 3, which illustrates the 
# time series of the different SDM-derived marine survival covariates.

# This version shows the bongo biomass index and the YOY rockfish index separately.
# It also does not show the index of abundance for salmon.


## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)

# source the prep scripts for each of the surveys
# source(here::here("R", "CCES_make_mesh.R"))
source(here::here("R", "PRS_PWCC_make_mesh.R"))
source(here::here("R", "hake_survey_make_mesh.R"))
source(here::here("R", "JSOES_seabirds_make_mesh.R"))
source(here::here("R", "JSOES_make_mesh.R"))


# load the SDM outputs + estimated uncertainty from the stage 1 models
# 05.1.1_jsoes_bongo_biomass_SDM
load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "jsoes_bongo_biomass_SDM_output.rda"))
jsoes_bongo_biomass_SDM_Obj <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_Obj
jsoes_bongo_biomass_SDM_Opt <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_Opt
jsoes_bongo_biomass_SDM_report <- jsoes_bongo_biomass_SDM_output$jsoes_bongo_biomass_SDM_report

load(here::here("R", "05_stage1_SDM", "05.1.1_jsoes_bongo_biomass_SDM", "estimated_SE_prey_model.rda"))
SE_pianka_o_csyif_prey_field_t_prey_model <- SE_prey_model$SE_pianka_o_csyif_prey_field_t_prey_model
SE_pianka_o_cssif_prey_field_t_prey_model <- SE_prey_model$SE_pianka_o_cssif_prey_field_t_prey_model

# 05.1.2_rf_SDM
load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "rf_SDM_output.rda"))
rf_SDM_Obj <- rf_SDM_output$rf_SDM_Obj
rf_SDM_Opt <- rf_SDM_output$rf_SDM_Opt
rf_SDM_report <- rf_SDM_output$rf_SDM_report

load(here::here("R", "05_stage1_SDM", "05.1.2_rf_SDM", "estimated_SE_rf_model.rda"))
SE_pianka_o_csyif_rf_t_rf_model <- SE_rf_model$SE_pianka_o_csyif_rf_t_rf_model
SE_pianka_o_cssif_rf_t_rf_model <- SE_rf_model$SE_pianka_o_cssif_rf_t_rf_model
SE_pianka_o_csyif_rf_t_rf_model

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

## load data on some taxa
jsoes_bongo_biomass_cancer_crab_larvae <- read.csv(here::here("model_inputs", "jsoes_bongo_biomass_cancer_crab_larvae.csv"))
jsoes_bongo_biomass_cancer_crab_larvae <- subset(jsoes_bongo_biomass_cancer_crab_larvae, !(year %in% c(1998, 2022:2025)))
sosh <- subset(birds_long, species == "sooty_shearwater")
comu <- subset(birds_long, species == "common_murre")
# extract Yearling Interior Chinook
csyif <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_fa")
# extract subyearling Interior Chinook
cssif <- subset(jsoes_long, species == "chinook_salmon_subyearling_interior_fa")

#### identify where in our covariance matrix our ADREPORTed variables are ####

### Prey field model

common_years_prey <- intersect(
  unique(jsoes_bongo_biomass_cancer_crab_larvae$year), # jsoes bongo
  unique(rf$year) # PRS/PWCC
)

# first extract the summary of the SD
prey_field_SDM_SD_summary <- summary(jsoes_bongo_biomass_SDM_Opt$SD)
rf_SDM_SD_summary <- summary(rf_SDM_Opt$SD)

# use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
# csyif index of abundance
csyif_index_of_abundance_indices <- which(grepl("csyif_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
csyif_index_of_abundance_cov_indices <- which(sqrt(diag(jsoes_bongo_biomass_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[csyif_index_of_abundance_indices, "Std. Error"])

# cssif index of abundance
cssif_index_of_abundance_indices <- which(grepl("cssif_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
cssif_index_of_abundance_cov_indices <- which(sqrt(diag(jsoes_bongo_biomass_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[cssif_index_of_abundance_indices, "Std. Error"])

# prey field
prey_field_index_of_abundance_indices <- which(grepl("prey_field_index_of_abundance", rownames(prey_field_SDM_SD_summary)))
prey_field_index_of_abundance_cov_indices <- which(sqrt(diag(jsoes_bongo_biomass_SDM_Opt$SD$cov)) %in% prey_field_SDM_SD_summary[prey_field_index_of_abundance_indices, "Std. Error"])

# rf
rf_index_of_abundance_indices <- which(grepl("rf_index_of_abundance", rownames(rf_SDM_SD_summary)))
rf_index_of_abundance_cov_indices <- which(sqrt(diag(rf_SDM_Opt$SD$cov)) %in% rf_SDM_SD_summary[rf_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for each using the indices
# csyif
csyif_index_of_abundance_cov_matrix <- as.matrix(jsoes_bongo_biomass_SDM_Opt$SD$cov[csyif_index_of_abundance_cov_indices, csyif_index_of_abundance_cov_indices])
colnames(csyif_index_of_abundance_cov_matrix) <- sort(unique(csyif$year))
rownames(csyif_index_of_abundance_cov_matrix) <- sort(unique(csyif$year))

# cssif
cssif_index_of_abundance_cov_matrix <- as.matrix(jsoes_bongo_biomass_SDM_Opt$SD$cov[cssif_index_of_abundance_cov_indices, cssif_index_of_abundance_cov_indices])
colnames(cssif_index_of_abundance_cov_matrix) <- sort(unique(cssif$year))
rownames(cssif_index_of_abundance_cov_matrix) <- sort(unique(cssif$year))

# prey field
prey_field_index_of_abundance_cov_matrix <- as.matrix(jsoes_bongo_biomass_SDM_Opt$SD$cov[prey_field_index_of_abundance_cov_indices, prey_field_index_of_abundance_cov_indices])
colnames(prey_field_index_of_abundance_cov_matrix) <- unique(jsoes_bongo_biomass_cancer_crab_larvae$year)
rownames(prey_field_index_of_abundance_cov_matrix) <- unique(jsoes_bongo_biomass_cancer_crab_larvae$year)

# rf
rf_index_of_abundance_cov_matrix <- as.matrix(rf_SDM_Opt$SD$cov[rf_index_of_abundance_cov_indices, rf_index_of_abundance_cov_indices])
colnames(rf_index_of_abundance_cov_matrix) <- unique(rf$year)
rownames(rf_index_of_abundance_cov_matrix) <- unique(rf$year)

### Seabird model

# first extract the summary of the SD
seabird_SDM_SD_summary <- summary(seabird_SDM_Opt$SD)

# use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
# comu index of abundance
comu_index_of_abundance_indices <- which(grepl("comu_index_of_abundance", rownames(seabird_SDM_SD_summary)))
comu_index_of_abundance_cov_indices <- which(sqrt(diag(seabird_SDM_Opt$SD$cov)) %in% seabird_SDM_SD_summary[comu_index_of_abundance_indices, "Std. Error"])

# sosh index of abundance
sosh_index_of_abundance_indices <- which(grepl("sosh_index_of_abundance", rownames(seabird_SDM_SD_summary)))
sosh_index_of_abundance_cov_indices <- which(sqrt(diag(seabird_SDM_Opt$SD$cov)) %in% seabird_SDM_SD_summary[sosh_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for each using the indices
# comu
comu_index_of_abundance_cov_matrix <- as.matrix(seabird_SDM_Opt$SD$cov[comu_index_of_abundance_cov_indices, comu_index_of_abundance_cov_indices])
colnames(comu_index_of_abundance_cov_matrix) <- sort(unique(comu$year))
rownames(comu_index_of_abundance_cov_matrix) <- sort(unique(comu$year))

# sosh
sosh_index_of_abundance_cov_matrix <- as.matrix(seabird_SDM_Opt$SD$cov[sosh_index_of_abundance_cov_indices, sosh_index_of_abundance_cov_indices])
colnames(sosh_index_of_abundance_cov_matrix) <- min(sosh$year):max(sosh$year)
rownames(sosh_index_of_abundance_cov_matrix) <- min(sosh$year):max(sosh$year)


### Hake model

# # first extract the summary of the SD
hake_SDM_SD_summary <- summary(hake_SDM_Opt$SD)

# # use this object to figure out what the standard errors are from our variables of interest, and then match them to the cov object
hake_index_of_abundance_indices <- which(grepl("hake_index_of_abundance", rownames(hake_SDM_SD_summary)))
hake_index_of_abundance_cov_indices <- which(sqrt(diag(hake_SDM_Opt$SD$cov)) %in% hake_SDM_SD_summary[hake_index_of_abundance_indices, "Std. Error"])

# extract the covariance matrices for hake
hake_index_of_abundance_cov_matrix <- as.matrix(hake_SDM_Opt$SD$cov[hake_index_of_abundance_cov_indices, hake_index_of_abundance_cov_indices])
colnames(hake_index_of_abundance_cov_matrix) <- min(hake$year):max(hake$year)
rownames(hake_index_of_abundance_cov_matrix) <- min(hake$year):max(hake$year)

#### Collate marine survival covariates and their uncertainties ####

### Indices of abundance
prey_field_index_df <- data.frame(year = unique(jsoes_bongo_biomass_cancer_crab_larvae$year),
                            prey_field_index = jsoes_bongo_biomass_SDM_report$prey_field_index_of_abundance/max(jsoes_bongo_biomass_SDM_report$prey_field_index_of_abundance),
                            prey_field_index_SE = sqrt(diag(prey_field_index_of_abundance_cov_matrix))/(max(jsoes_bongo_biomass_SDM_report$prey_field_index_of_abundance)))

rf_index_df <- data.frame(year = unique(rf$year),
                                  rf_index = rf_SDM_report$rf_index_of_abundance/max(rf_SDM_report$rf_index_of_abundance),
                                  rf_index_SE = sqrt(diag(rf_index_of_abundance_cov_matrix))/(max(rf_SDM_report$rf_index_of_abundance)))

csyif_index_df <- data.frame(year = sort(unique(csyif$year)),
                             csyif_index = jsoes_bongo_biomass_SDM_report$csyif_index_of_abundance/max(jsoes_bongo_biomass_SDM_report$csyif_index_of_abundance),
                             csyif_index_SE = sqrt(diag(csyif_index_of_abundance_cov_matrix))/max(jsoes_bongo_biomass_SDM_report$csyif_index_of_abundance))

cssif_index_df <- data.frame(year = sort(unique(cssif$year)),
                             cssif_index = jsoes_bongo_biomass_SDM_report$cssif_index_of_abundance/max(jsoes_bongo_biomass_SDM_report$cssif_index_of_abundance),
                             cssif_index_SE = sqrt(diag(cssif_index_of_abundance_cov_matrix))/max(jsoes_bongo_biomass_SDM_report$cssif_index_of_abundance))

comu_index_df <- data.frame(year = sort(unique(comu$year)),
                            comu_index = seabird_SDM_report$comu_index_of_abundance/max(seabird_SDM_report$comu_index_of_abundance),
                            comu_index_SE = sqrt(diag(comu_index_of_abundance_cov_matrix))/max(seabird_SDM_report$comu_index_of_abundance))

sosh_index_df <- data.frame(year = min(sosh$year):max(sosh$year),
                            sosh_index = seabird_SDM_report$sosh_index_of_abundance/max(seabird_SDM_report$sosh_index_of_abundance),
                            sosh_index_SE = sqrt(diag(sosh_index_of_abundance_cov_matrix))/max(seabird_SDM_report$sosh_index_of_abundance))
sosh_index_df <- subset(sosh_index_df, year %in% unique(sosh$year))


hake_index_df <- data.frame(year = min(hake$year):max(hake$year),
                            hake_index = hake_SDM_report$hake_index_of_abundance/max(hake_SDM_report$hake_index_of_abundance),
                            hake_index_SE = sqrt(diag(hake_index_of_abundance_cov_matrix))/max(hake_SDM_report$hake_index_of_abundance))
hake_index_df <- subset(hake_index_df, year %in% unique(hake$year))


# combine these together

csyif_index_df %>% 
  left_join(cssif_index_df, by = "year") %>% 
  left_join(prey_field_index_df, by = "year") %>% 
  left_join(rf_index_df, by = "year") %>% 
  left_join(comu_index_df, by = "year") %>% 
  left_join(sosh_index_df, by = "year") %>% 
  left_join(hake_index_df, by = "year") -> indices_df

indices_df %>% 
  pivot_longer(., cols = -c("year")) %>% 
  mutate(variable = ifelse(grepl("SE", name), "SE", "estimate")) %>% 
  mutate(index = gsub("_SE", "", name)) -> indices_long

indices_long %>% 
  dplyr::select(-name) %>% 
  pivot_wider(names_from = variable, values_from = value) -> indices_long


index_plot_titles <- data.frame(index = c("cssif_index", "csyif_index",
                                          "prey_field_index", "rf_index", "comu_index",
                                          "sosh_index", "hake_index"),
                                name = c("Subyearlings", "Yearlings", "Prey Field", "YOY Rockfishes",
                                         "Common Murres", "Sooty Shearwaters", "Hake"))
indices_long %>% 
  left_join(index_plot_titles, by = "index") -> indices_long

indices_long$name <- factor(indices_long$name, levels = c("Subyearlings", "Yearlings", "Prey Field", "YOY Rockfishes",
                                                          "Common Murres", "Sooty Shearwaters", "Hake"))


### Overlap metrics

prey_field_csyif_overlap_df <- data.frame(year = unique(jsoes_bongo_biomass_cancer_crab_larvae$year),
                                          pianka_o_csyif_prey_field = jsoes_bongo_biomass_SDM_report$pianka_o_csyif_prey_field_t,
                                          pianka_o_csyif_prey_field_SE = SE_pianka_o_csyif_prey_field_t_prey_model)

rf_csyif_overlap_df <- data.frame(year = unique(rf$year),
                                          pianka_o_csyif_rf = rf_SDM_report$pianka_o_csyif_rf_t,
                                          pianka_o_csyif_rf_SE = SE_pianka_o_csyif_rf_t_rf_model)

comu_csyif_overlap_df <- data.frame(year = sort(unique(comu$year)),
                                          pianka_o_csyif_comu = seabird_SDM_report$pianka_o_csyif_comu_t,
                                          pianka_o_csyif_comu_SE = SE_pianka_o_csyif_comu_t_seabird_model)

sosh_csyif_overlap_df <- data.frame(year = sort(unique(sosh$year)),
                                    pianka_o_csyif_sosh = seabird_SDM_report$pianka_o_csyif_sosh_t,
                                    pianka_o_csyif_sosh_SE = SE_pianka_o_csyif_sosh_t_seabird_model)

prey_field_cssif_overlap_df <- data.frame(year = unique(jsoes_bongo_biomass_cancer_crab_larvae$year),
                                          pianka_o_cssif_prey_field = jsoes_bongo_biomass_SDM_report$pianka_o_cssif_prey_field_t,
                                          pianka_o_cssif_prey_field_SE = SE_pianka_o_cssif_prey_field_t_prey_model)

rf_cssif_overlap_df <- data.frame(year = unique(rf$year),
                                          pianka_o_cssif_rf = rf_SDM_report$pianka_o_cssif_rf_t,
                                          pianka_o_cssif_rf_SE = SE_pianka_o_cssif_rf_t_rf_model)

comu_cssif_overlap_df <- data.frame(year = sort(unique(comu$year)),
                                    pianka_o_cssif_comu = seabird_SDM_report$pianka_o_cssif_comu_t,
                                    pianka_o_cssif_comu_SE = SE_pianka_o_cssif_comu_t_seabird_model)

sosh_cssif_overlap_df <- data.frame(year = sort(unique(sosh$year)),
                                    pianka_o_cssif_sosh = seabird_SDM_report$pianka_o_cssif_sosh_t,
                                    pianka_o_cssif_sosh_SE = SE_pianka_o_cssif_sosh_t_seabird_model)

hake_csyif_overlap_df <- data.frame(year = sort(unique(hake$year)),
                                    pianka_o_csyif_hake = hake_SDM_report$pianka_o_csyif_hake_t,
                                    pianka_o_csyif_hake_SE = SE_pianka_o_csyif_hake_t_hake_model)

hake_cssif_overlap_df <- data.frame(year = sort(unique(hake$year)),
                                    pianka_o_cssif_hake = hake_SDM_report$pianka_o_cssif_hake_t,
                                    pianka_o_cssif_hake_SE = SE_pianka_o_cssif_hake_t_hake_model)

# combine these together

prey_field_csyif_overlap_df %>% 
  left_join(rf_csyif_overlap_df, by = "year") %>% 
  left_join(comu_csyif_overlap_df, by = "year") %>% 
  left_join(sosh_csyif_overlap_df, by = "year") %>% 
  left_join(hake_csyif_overlap_df, by = "year") -> csyif_overlap_indices

prey_field_cssif_overlap_df %>% 
  left_join(rf_cssif_overlap_df, by = "year") %>% 
  left_join(comu_cssif_overlap_df, by = "year") %>% 
  left_join(sosh_cssif_overlap_df, by = "year") %>% 
  left_join(hake_cssif_overlap_df, by = "year")  -> cssif_overlap_indices

csyif_overlap_indices %>% 
  left_join(cssif_overlap_indices, by = "year") -> overlap_df

overlap_df %>% 
  pivot_longer(., cols = -c("year")) %>% 
  mutate(variable = ifelse(grepl("SE", name), "SE", "estimate")) %>% 
  mutate(index = gsub("_SE", "", name)) -> overlap_long

overlap_long %>% 
  dplyr::select(-name) %>% 
  pivot_wider(names_from = variable, values_from = value) -> overlap_long


overlap_plot_titles <- data.frame(index = c("pianka_o_csyif_prey_field",
                                            "pianka_o_csyif_rf",
                                            "pianka_o_csyif_comu",
                                            "pianka_o_csyif_sosh",
                                            "pianka_o_csyif_hake",
                                            "pianka_o_cssif_prey_field",
                                            "pianka_o_cssif_rf",
                                            "pianka_o_cssif_comu",
                                            "pianka_o_cssif_sosh",
                                            "pianka_o_cssif_hake"),
                                name = c("Yearlings x Prey Field",
                                         "Yearlings x YOY Rockfishes",
                                         "Yearlings x Common Murres",
                                         "Yearlings x Sooty Shearwaters",
                                         "Yearlings x Hake",
                                         "Subyearlings x Prey Field",
                                         "Subyearlings x YOY Rockfishes",
                                         "Subyearlings x Common Murres",
                                         "Subyearlings x Sooty Shearwaters",
                                         "Subyearlings x Hake"))
overlap_long %>% 
  left_join(overlap_plot_titles, by = "index") -> overlap_long

overlap_long$name <- factor(overlap_long$name, levels = c("Subyearlings x Prey Field",
                                                          "Yearlings x Prey Field",
                                                          "Subyearlings x YOY Rockfishes",
                                                          "Yearlings x YOY Rockfishes",
                                                          "Subyearlings x Common Murres",
                                                          "Yearlings x Common Murres",
                                                          "Subyearlings x Sooty Shearwaters",
                                                          "Yearlings x Sooty Shearwaters",
                                                          "Subyearlings x Hake",
                                                          "Yearlings x Hake"))

# complete this df so that the plot appropriately connects points
overlap_long %>% 
  complete(year = 2001:2021, name) -> overlap_long


#### Generate figure for presentation ####

hake_index_minimal_plot <- ggplot(subset(indices_long, index == "hake_index" & !(is.na(estimate))), aes(x = year, y = estimate)) +
  geom_line(linewidth = 0.8) +
  # geom_point() +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(here::here("figures", "presentation_figures", "indices", "hake_index_minimal_plot.png"), hake_index_minimal_plot,
       height = 2, width = 5)

sosh_index_minimal_plot <- ggplot(subset(indices_long, index == "sosh_index" & !(is.na(estimate))), aes(x = year, y = estimate)) +
  geom_line(linewidth = 0.8) +
  # geom_point() +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(here::here("figures", "presentation_figures", "indices", "sosh_index_minimal_plot.png"), sosh_index_minimal_plot,
       height = 2, width = 5)

# sosh index plot but for SDM visualization

sosh_index_plot <- ggplot(subset(indices_long, index == "sosh_index" & !(is.na(estimate)) & year %in% 2005:2011), aes(x = year, y = estimate)) +
  # geom_line(linewidth = 0.8) +
  geom_point() +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_continuous(breaks = 2005:2011)

ggsave(here::here("figures", "presentation_figures", "indices", "sosh_index_plot.png"), sosh_index_plot,
       height = 2, width = 5)

# sosh x yearlings overlap plot for presentation

sosh_yearling_overlap_plot <- ggplot(subset(overlap_long, index == "pianka_o_csyif_sosh" & !(is.na(estimate)) & year %in% 2005:2011), aes(x = year, y = estimate)) +
  # geom_line(linewidth = 0.8) +
  geom_point() +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_continuous(breaks = 2005:2011)

ggsave(here::here("figures", "presentation_figures", "indices", "sosh_yearling_overlap_plot.png"), sosh_yearling_overlap_plot,
       height = 2, width = 5)

comu_index_minimal_plot <- ggplot(subset(indices_long, index == "comu_index" & !(is.na(estimate))), aes(x = year, y = estimate)) +
  geom_line(linewidth = 0.8) +
  # geom_point() +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

ggsave(here::here("figures", "presentation_figures", "indices", "comu_index_minimal_plot.png"), comu_index_minimal_plot,
       height = 2, width = 5)




#### Generate figures ####

# SANITY CHECK
hake %>% group_by(year) %>% summarise(total = sum(NASC))
sosh %>% group_by(year) %>% summarise(total = sum(n_per_km2))
comu %>% group_by(year) %>% summarise(total = sum(n_per_km2))
csyif %>% group_by(year) %>% summarise(total = sum(n_per_km))
cssif %>% group_by(year) %>% summarise(total = sum(n_per_km))


abundance_ts_plot <- ggplot(indices_long, aes(x = year, y = estimate, 
                         ymin = estimate - 1.96*SE,
                         ymax = estimate + 1.96*SE)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  geom_line() +
  facet_wrap(~name, ncol = 1) +
  ylab("Relative Abundance Index") +
  xlab("Year") +
  theme(plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

ggsave(here::here("figures", "paper_figures", "biomass_models", "fig3a_abundance_ts_plot.png"), abundance_ts_plot,  
       height = 8, width = 4)

overlap_ts_plot <- ggplot(overlap_long, aes(x = year, y = estimate, 
                                              ymin = estimate - 1.96*SE,
                                              ymax = estimate + 1.96*SE)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  geom_line() +
  facet_wrap(~name, ncol = 1) +
  ylab("Local Index of Collocation") +
  xlab("Year") +
  theme(plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        # panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

ggsave(here::here("figures", "paper_figures", "biomass_models", "fig3b_overlap_ts_plot.png"), overlap_ts_plot,  
       height = 8, width = 4)


fig3_covariates_ts <- ggarrange(abundance_ts_plot,
                               overlap_ts_plot, 
                               ncol = 2, nrow = 1)

ggsave(here::here("figures", "paper_figures", "biomass_models", "fig3_covariates_ts.png"), fig3_covariates_ts,  
       height = 8, width = 8)

#### Version 2 - new layout ####

# subyearlings and yearlings together on one plot
chinook_abundance_ts_plot <- ggplot(subset(indices_long, name %in% c("Yearlings", "Subyearlings")), aes(x = year, y = estimate, 
                                              ymin = estimate - 1.96*SE,
                                              ymax = estimate + 1.96*SE)) +
  geom_errorbar(width = 0.2, color = "gray50") +
  geom_point() +
  geom_line() +
  facet_wrap(~name, nrow = 1) +
  scale_y_continuous(labels = c("0.00", "0.50", "1.00", "1.50", "2.00")) +
  ylab("Relative Abundance Index") +
  xlab("Year") +
  theme(plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # these plot margins are to leave space for the panel label
        plot.margin = unit(c(1.0, 0.2, 0.2, 0.2),"cm"))

# prey field, common murres, sooty shearwaters, and hake on one plot

non_chinook_abundance_ts_plot <- ggplot(subset(indices_long, !(name %in% c("Yearlings", "Subyearlings"))), aes(x = year, y = estimate, 
                                                                                                        ymin = estimate - 1.96*SE,
                                                                                                        ymax = estimate + 1.96*SE)) +
  geom_errorbar(width = 0.2, color = "gray50") +
  geom_point() +
  geom_line() +
  facet_wrap(~name, ncol = 1) +
  ylab("Relative Abundance Index") +
  scale_y_continuous(labels = c("0.00", "0.50", "1.00", "1.50", "2.00")) +
  xlab("Year") +
  theme(plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # these plot margins are to leave space for the panel label
        plot.margin = unit(c(1.0, 0.2, 0.2, 0.2),"cm"))


# overlap on one plot

interior_overlap_ts_plot <- ggplot(overlap_long, aes(x = year, y = estimate, 
                                            ymin = estimate - 1.96*SE,
                                            ymax = estimate + 1.96*SE)) +
  geom_errorbar(width = 0.2, color = "gray50") +
  geom_point() +
  geom_line() +
  facet_wrap(~name, nrow = 5) +
  ylab("Local Index of Collocation") +
  scale_x_continuous(lim = c(min(indices_long$year), max(indices_long$year))) +
  xlab("Year") +
  theme(plot.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        # these plot margins are to leave space for the panel label
        plot.margin = unit(c(1.0, 0.2, 0.2, 0.2),"cm"))

fig3_v2_covariates_ts <- ggarrange(NULL,
                                   chinook_abundance_ts_plot,
                                   non_chinook_abundance_ts_plot,
                                   interior_overlap_ts_plot, 
                                   widths = c(3.5, 6.5, 3.5, 6.5),
                                   heights = c(2.3, 8),
                                   labels = c("", "(A)", "(B)", "(C)"),
                                   font.label = list(size = 20, face = "plain"),
                                   label.x = 0.025, label.y = 0.985,
                                ncol = 2, nrow = 2)

ggsave(here::here("figures", "paper_figures", "biomass_models", "fig3_v2_covariates_ts.png"), fig3_v2_covariates_ts,  
       height = 10, width = 12)


fig3_no_salmon_index_covariates_ts <- ggarrange(non_chinook_abundance_ts_plot,
                                   interior_overlap_ts_plot, 
                                   widths = c(3.5, 6.5),
                                   # heights = c(2.3, 8),
                                   labels = c("(A)", "(B)                                             (C)"),
                                   font.label = list(size = 20, face = "plain"),
                                   label.x = 0.025, label.y = 0.985,
                                   hjust = 0,
                                   ncol = 2, nrow = 1)

# Create a version without the salmon index of abundance

ggsave(here::here("figures", "paper_figures", "no_salmon_index_models", "fig3_covariates_ts.png"), fig3_no_salmon_index_covariates_ts,  
       height = 8, width = 12)


#### get summary values for the paper ####

overlap_long %>% 
  group_by(name) %>% 
  summarise(mean_overlap = mean(estimate, na.rm = T)) %>% 
  arrange(desc(mean_overlap))

max(subset(overlap_long, name == "Yearlings x Sooty Shearwaters")$estimate, na.rm = T)/min(subset(overlap_long, name == "Yearlings x Sooty Shearwaters")$estimate, na.rm = T)
max(subset(overlap_long, name == "Subyearlings x Sooty Shearwaters")$estimate, na.rm = T)/min(subset(overlap_long, name == "Subyearlings x Sooty Shearwaters")$estimate, na.rm = T)
max(subset(overlap_long, name == "Yearlings x Common Murres")$estimate, na.rm = T)/min(subset(overlap_long, name == "Yearlings x Common Murres")$estimate, na.rm = T)
max(subset(overlap_long, name == "Subyearlings x Common Murres")$estimate, na.rm = T)/min(subset(overlap_long, name == "Subyearlings x Common Murres")$estimate, na.rm = T)
max(subset(overlap_long, name == "Yearlings x Hake")$estimate, na.rm = T)/min(subset(overlap_long, name == "Yearlings x Hake")$estimate, na.rm = T)
max(subset(overlap_long, name == "Subyearlings x Hake")$estimate, na.rm = T)/min(subset(overlap_long, name == "Subyearlings x Hake")$estimate, na.rm = T)
max(subset(overlap_long, name == "Yearlings x Prey Field")$estimate, na.rm = T)/min(subset(overlap_long, name == "Yearlings x Prey Field")$estimate, na.rm = T)
max(subset(overlap_long, name == "Subyearlings x Prey Field")$estimate, na.rm = T)/min(subset(overlap_long, name == "Subyearlings x Prey Field")$estimate, na.rm = T)
max(subset(overlap_long, name == "Yearlings x YOY Rockfishes")$estimate, na.rm = T)/min(subset(overlap_long, name == "Yearlings x YOY Rockfishes")$estimate, na.rm = T)
max(subset(overlap_long, name == "Subyearlings x YOY Rockfishes")$estimate, na.rm = T)/min(subset(overlap_long, name == "Subyearlings x YOY Rockfishes")$estimate, na.rm = T)









