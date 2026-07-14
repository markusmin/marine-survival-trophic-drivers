## 13_EiV_fit_fig

# Description: This script generates figures showing the “fit” of the observed value 
# for covariates based on a latent covariate being estimated

common_years_seabirds_prey <- c(2003, 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2015, 2016, 2017, 2018, 2019)
common_years_hake_prey <- seq(2001, 2021, 2)
jsoes_years <- c(1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 
                 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 
                 2020, 2021)

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
library(ggrepel)

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
                                  prey_field_index = jsoes_bongo_biomass_SDM_report$prey_field_index_of_abundance,
                                  prey_field_index_SE = sqrt(diag(prey_field_index_of_abundance_cov_matrix)))

rf_index_df <- data.frame(year = unique(rf$year),
                          rf_index = rf_SDM_report$rf_index_of_abundance,
                          rf_index_SE = sqrt(diag(rf_index_of_abundance_cov_matrix)))

csyif_index_df <- data.frame(year = sort(unique(csyif$year)),
                             csyif_index = jsoes_bongo_biomass_SDM_report$csyif_index_of_abundance,
                             csyif_index_SE = sqrt(diag(csyif_index_of_abundance_cov_matrix)))

cssif_index_df <- data.frame(year = sort(unique(cssif$year)),
                             cssif_index = jsoes_bongo_biomass_SDM_report$cssif_index_of_abundance,
                             cssif_index_SE = sqrt(diag(cssif_index_of_abundance_cov_matrix)))

comu_index_df <- data.frame(year = sort(unique(comu$year)),
                            comu_index = seabird_SDM_report$comu_index_of_abundance,
                            comu_index_SE = sqrt(diag(comu_index_of_abundance_cov_matrix)))

sosh_index_df <- data.frame(year = min(sosh$year):max(sosh$year),
                            sosh_index = seabird_SDM_report$sosh_index_of_abundance,
                            sosh_index_SE = sqrt(diag(sosh_index_of_abundance_cov_matrix)))
sosh_index_df <- subset(sosh_index_df, year %in% unique(sosh$year))


hake_index_df <- data.frame(year = min(hake$year):max(hake$year),
                            hake_index = hake_SDM_report$hake_index_of_abundance,
                            hake_index_SE = sqrt(diag(hake_index_of_abundance_cov_matrix)))
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

# get all the SDM covariate estimates together
overlap_long %>% 
  bind_rows(indices_long) %>% 
  mutate(stage = "SDM") -> SDM_estimates_long


# get SDM estimates relevant for all four SAR models separately
# z-score each to make sure they match the SAR estimates

SDM_estimates_long %>% 
  filter(!(grepl("sosh|comu|csyif_index|cssif_index", index))) %>% 
  filter(!(grepl("cssif", index))) %>% 
  filter(year %in% common_years_hake_prey) %>% 
  group_by(index) %>% 
  mutate(SE = sqrt(SE^2*(1/(sd(estimate))^2))) %>% 
  # mutate(SE = SE * (1/sd(estimate))^2) %>% 
  mutate(estimate = (estimate - mean(estimate)) / sd(estimate)) %>%
  ungroup() -> SDM_estimates_hake_csyif_long

SDM_estimates_long %>% 
  filter(!(grepl("sosh|comu|csyif_index|cssif_index", index))) %>% 
  filter(!(grepl("csyif", index))) %>% 
  filter(year %in% common_years_hake_prey) %>% 
  group_by(index) %>% 
  mutate(SE = sqrt(SE^2*(1/(sd(estimate))^2))) %>% 
  # mutate(SE = SE * (1/sd(estimate))^2) %>% 
  mutate(estimate = (estimate - mean(estimate)) / sd(estimate)) %>%
  ungroup() -> SDM_estimates_hake_cssif_long

SDM_estimates_long %>% 
  filter(!(grepl("hake|csyif_index|cssif_index", index))) %>% 
  filter(!(grepl("cssif", index))) %>% 
  filter(year %in% common_years_seabirds_prey) %>% 
  group_by(index) %>% 
  mutate(SE = sqrt(SE^2*(1/(sd(estimate))^2))) %>% 
  # mutate(SE = SE * (1/sd(estimate))^2) %>% 
  mutate(estimate = (estimate - mean(estimate)) / sd(estimate)) %>%
  ungroup() -> SDM_estimates_seabirds_csyif_long

SDM_estimates_long %>% 
  filter(!(grepl("hake|csyif_index|cssif_index", index))) %>% 
  filter(!(grepl("csyif", index))) %>% 
  filter(year %in% common_years_seabirds_prey) %>% 
  group_by(index) %>% 
  mutate(SE = sqrt(SE^2*(1/(sd(estimate))^2))) %>% 
  # mutate(SE = SE * (1/sd(estimate))^2) %>% 
  mutate(estimate = (estimate - mean(estimate)) / sd(estimate)) %>%
  ungroup() -> SDM_estimates_seabirds_cssif_long




#### Load the models for SAR for combined interior ####

### seabird model

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.2.6_seabird_two_prey_no_salmon_index_one_model_SAR", "interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Obj <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Obj
interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt
interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report

as.data.frame(summary(interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt$SD)) %>% 
  rownames_to_column("parameter") %>% 
  janitor::clean_names() -> seabird_csyif_SAR_SD_summary

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2.6_seabird_two_prey_no_salmon_index_one_model_SAR", "interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Obj <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Obj
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report

as.data.frame(summary(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt$SD)) %>% 
  rownames_to_column("parameter") %>% 
  janitor::clean_names() -> seabird_cssif_SAR_SD_summary

### hake model

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3.4_hake_two_prey_no_salmon_index_one_model_SAR", "interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Obj <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Obj
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report

as.data.frame(summary(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt$SD)) %>% 
  rownames_to_column("parameter") %>% 
  janitor::clean_names() -> hake_csyif_SAR_SD_summary

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.3.4_hake_two_prey_no_salmon_index_one_model_SAR", "interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Obj <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Obj
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report

as.data.frame(summary(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt$SD)) %>% 
  rownames_to_column("parameter") %>% 
  janitor::clean_names() -> hake_cssif_SAR_SD_summary




#### Reformat the SAR outputs ####

### Seabirds x subyearlings ###

# prey field index
seabird_cssif_prey_field_index_df <- data.frame(year = common_years_seabirds_prey,
                                  prey_field_index = subset(seabird_cssif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$estimate,
                                  prey_field_index_SE = subset(seabird_cssif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$std_error)

# prey field overlap
seabird_cssif_prey_field_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                                  pianka_o_cssif_prey_field = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_prey_field_t", parameter))$estimate,
                                                  pianka_o_cssif_prey_field_SE = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_prey_field_t", parameter))$std_error)

# rf index
seabird_cssif_rf_index_df <- data.frame(year = common_years_seabirds_prey,
                                          rf_index = subset(seabird_cssif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$estimate,
                                          rf_index_SE = subset(seabird_cssif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$std_error)

# rf overlap
seabird_cssif_rf_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                          pianka_o_cssif_rf = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_rf_t", parameter))$estimate,
                                          pianka_o_cssif_rf_SE = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_rf_t", parameter))$std_error)

# sosh index
seabird_cssif_sosh_index_df <- data.frame(year = common_years_seabirds_prey,
                                        sosh_index = subset(seabird_cssif_SAR_SD_summary, grepl("sosh_index_t_latent", parameter))$estimate,
                                        sosh_index_SE = subset(seabird_cssif_SAR_SD_summary, grepl("sosh_index_t_latent", parameter))$std_error)

# sosh overlap
seabird_cssif_sosh_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                            pianka_o_cssif_sosh = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_sosh_t", parameter))$estimate,
                                            pianka_o_cssif_sosh_SE = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_sosh_t", parameter))$std_error)

# comu index
seabird_cssif_comu_index_df <- data.frame(year = common_years_seabirds_prey,
                                        comu_index = subset(seabird_cssif_SAR_SD_summary, grepl("comu_index_t_latent", parameter))$estimate,
                                        comu_index_SE = subset(seabird_cssif_SAR_SD_summary, grepl("comu_index_t_latent", parameter))$std_error)

# comu overlap
seabird_cssif_comu_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                            pianka_o_cssif_comu = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_comu_t", parameter))$estimate,
                                            pianka_o_cssif_comu_SE = subset(seabird_cssif_SAR_SD_summary, grepl("pianka_o_cssif_comu_t", parameter))$std_error)

## join all together

seabird_cssif_prey_field_index_df  %>% 
  left_join(seabird_cssif_prey_field_overlap_df, by = "year") %>% 
  left_join(seabird_cssif_rf_index_df, by = "year") %>% 
  left_join(seabird_cssif_rf_overlap_df, by = "year") %>% 
  left_join(seabird_cssif_sosh_index_df, by = "year") %>% 
  left_join(seabird_cssif_sosh_overlap_df, by = "year") %>% 
  left_join(seabird_cssif_comu_index_df, by = "year") %>% 
  left_join(seabird_cssif_comu_overlap_df, by = "year") %>% 
  mutate(model = "Seabirds x Subyearlings") -> seabirds_cssif_SAR_estimates


### Seabirds x yearlings ###

# prey field index
seabird_csyif_prey_field_index_df <- data.frame(year = common_years_seabirds_prey,
                                                prey_field_index = subset(seabird_csyif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$estimate,
                                                prey_field_index_SE = subset(seabird_csyif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$std_error)

# prey field overlap
seabird_csyif_prey_field_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                                  pianka_o_csyif_prey_field = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_prey_field_t", parameter))$estimate,
                                                  pianka_o_csyif_prey_field_SE = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_prey_field_t", parameter))$std_error)

# rf index
seabird_csyif_rf_index_df <- data.frame(year = common_years_seabirds_prey,
                                        rf_index = subset(seabird_csyif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$estimate,
                                        rf_index_SE = subset(seabird_csyif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$std_error)

# rf overlap
seabird_csyif_rf_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                          pianka_o_csyif_rf = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_rf_t", parameter))$estimate,
                                          pianka_o_csyif_rf_SE = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_rf_t", parameter))$std_error)

# sosh index
seabird_csyif_sosh_index_df <- data.frame(year = common_years_seabirds_prey,
                                          sosh_index = subset(seabird_csyif_SAR_SD_summary, grepl("sosh_index_t_latent", parameter))$estimate,
                                          sosh_index_SE = subset(seabird_csyif_SAR_SD_summary, grepl("sosh_index_t_latent", parameter))$std_error)

# sosh overlap
seabird_csyif_sosh_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                            pianka_o_csyif_sosh = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_sosh_t", parameter))$estimate,
                                            pianka_o_csyif_sosh_SE = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_sosh_t", parameter))$std_error)

# comu index
seabird_csyif_comu_index_df <- data.frame(year = common_years_seabirds_prey,
                                          comu_index = subset(seabird_csyif_SAR_SD_summary, grepl("comu_index_t_latent", parameter))$estimate,
                                          comu_index_SE = subset(seabird_csyif_SAR_SD_summary, grepl("comu_index_t_latent", parameter))$std_error)

# comu overlap
seabird_csyif_comu_overlap_df <- data.frame(year = common_years_seabirds_prey,
                                            pianka_o_csyif_comu = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_comu_t", parameter))$estimate,
                                            pianka_o_csyif_comu_SE = subset(seabird_csyif_SAR_SD_summary, grepl("pianka_o_csyif_comu_t", parameter))$std_error)

## join all together

seabird_csyif_prey_field_index_df  %>% 
  left_join(seabird_csyif_prey_field_overlap_df, by = "year") %>% 
  left_join(seabird_csyif_rf_index_df, by = "year") %>% 
  left_join(seabird_csyif_rf_overlap_df, by = "year") %>% 
  left_join(seabird_csyif_sosh_index_df, by = "year") %>% 
  left_join(seabird_csyif_sosh_overlap_df, by = "year") %>% 
  left_join(seabird_csyif_comu_index_df, by = "year") %>% 
  left_join(seabird_csyif_comu_overlap_df, by = "year") %>% 
  mutate(model = "Seabirds x Yearlings") -> seabirds_csyif_SAR_estimates

### Hake x subyearlings ###

# prey field index
hake_cssif_prey_field_index_df <- data.frame(year = common_years_hake_prey,
                                                prey_field_index = subset(hake_cssif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$estimate,
                                                prey_field_index_SE = subset(hake_cssif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$std_error)

# prey field overlap
hake_cssif_prey_field_overlap_df <- data.frame(year = common_years_hake_prey,
                                               pianka_o_cssif_prey_field = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_prey_field_t", parameter))$estimate,
                                               pianka_o_cssif_prey_field_SE = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_prey_field_t", parameter))$std_error)

# rf index
hake_cssif_rf_index_df <- data.frame(year = common_years_hake_prey,
                                        rf_index = subset(hake_cssif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$estimate,
                                        rf_index_SE = subset(hake_cssif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$std_error)

# rf overlap
hake_cssif_rf_overlap_df <- data.frame(year = common_years_hake_prey,
                                       pianka_o_cssif_rf = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_rf_t", parameter))$estimate,
                                       pianka_o_cssif_rf_SE = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_rf_t", parameter))$std_error)

# hake index
hake_cssif_hake_index_df <- data.frame(year = common_years_hake_prey,
                                          hake_index = subset(hake_cssif_SAR_SD_summary, grepl("hake_index_t_latent", parameter))$estimate,
                                          hake_index_SE = subset(hake_cssif_SAR_SD_summary, grepl("hake_index_t_latent", parameter))$std_error)

# hake overlap
hake_cssif_hake_overlap_df <- data.frame(year = common_years_hake_prey,
                                         pianka_o_cssif_hake = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_hake_t", parameter))$estimate,
                                         pianka_o_cssif_hake_SE = subset(hake_cssif_SAR_SD_summary, grepl("pianka_o_cssif_hake_t", parameter))$std_error)

## join all together

hake_cssif_prey_field_index_df  %>% 
  left_join(hake_cssif_prey_field_overlap_df, by = "year") %>% 
  left_join(hake_cssif_rf_index_df, by = "year") %>% 
  left_join(hake_cssif_rf_overlap_df, by = "year") %>% 
  left_join(hake_cssif_hake_index_df, by = "year") %>% 
  left_join(hake_cssif_hake_overlap_df, by = "year") %>% 
  mutate(model = "Hake x Subyearlings") -> hake_cssif_SAR_estimates

### hake x yearlings ###

# prey field index
hake_csyif_prey_field_index_df <- data.frame(year = common_years_hake_prey,
                                                prey_field_index = subset(hake_csyif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$estimate,
                                                prey_field_index_SE = subset(hake_csyif_SAR_SD_summary, grepl("prey_field_index_t_latent", parameter))$std_error)

# prey field overlap
hake_csyif_prey_field_overlap_df <- data.frame(year = common_years_hake_prey,
                                               pianka_o_csyif_prey_field = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_prey_field_t", parameter))$estimate,
                                               pianka_o_csyif_prey_field_SE = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_prey_field_t", parameter))$std_error)

# rf index
hake_csyif_rf_index_df <- data.frame(year = common_years_hake_prey,
                                        rf_index = subset(hake_csyif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$estimate,
                                        rf_index_SE = subset(hake_csyif_SAR_SD_summary, grepl("rf_index_t_latent", parameter))$std_error)

# rf overlap
hake_csyif_rf_overlap_df <- data.frame(year = common_years_hake_prey,
                                       pianka_o_csyif_rf = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_rf_t", parameter))$estimate,
                                       pianka_o_csyif_rf_SE = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_rf_t", parameter))$std_error)

# hake index
hake_csyif_hake_index_df <- data.frame(year = common_years_hake_prey,
                                          hake_index = subset(hake_csyif_SAR_SD_summary, grepl("hake_index_t_latent", parameter))$estimate,
                                          hake_index_SE = subset(hake_csyif_SAR_SD_summary, grepl("hake_index_t_latent", parameter))$std_error)

# hake overlap
hake_csyif_hake_overlap_df <- data.frame(year = common_years_hake_prey,
                                         pianka_o_csyif_hake = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_hake_t", parameter))$estimate,
                                         pianka_o_csyif_hake_SE = subset(hake_csyif_SAR_SD_summary, grepl("pianka_o_csyif_hake_t", parameter))$std_error)

## join all together

hake_csyif_prey_field_index_df  %>% 
  left_join(hake_csyif_prey_field_overlap_df, by = "year") %>% 
  left_join(hake_csyif_rf_index_df, by = "year") %>% 
  left_join(hake_csyif_rf_overlap_df, by = "year") %>% 
  left_join(hake_csyif_hake_index_df, by = "year") %>% 
  left_join(hake_csyif_hake_overlap_df, by = "year") %>% 
  mutate(model = "Hake x Yearlings") -> hake_csyif_SAR_estimates

#### reformat and join together SDM + SAR estimates ####

## function to join + reformat ##

join_SDM_SAR_estimates <- function(SAR_estimates, SDM_estimates){
  SAR_estimates %>% 
    pivot_longer(cols = -c("year", "model"),
                 names_to = "name",
                 values_to = "value") %>% 
    mutate(variable = ifelse(grepl("SE", name), "SE", "estimate")) %>% 
    mutate(index = gsub("_SE", "", name)) %>% 
    dplyr::select(-name) %>% 
    pivot_wider(names_from = variable, values_from = value) %>% 
    mutate(stage = "SAR") -> SAR_estimates_long
  
  SAR_estimates_long %>% 
    bind_rows(SDM_estimates) %>% 
    dplyr::select(-c(model, name)) %>% 
    pivot_wider(names_from = "stage", 
                values_from = c("estimate", "SE")) -> SAR_SDM_comp
  
  return(SAR_SDM_comp)
}

## function to plot ##
plot_SAR_SDM_EiV_fit <- function(SAR_SDM_comp, index_name,
                                 plot_title){
  plot <- ggplot(subset(SAR_SDM_comp, index == index_name), 
         aes(x = estimate_SDM, y = estimate_SAR, label = year)) +
    geom_errorbar(aes(ymin = estimate_SAR - 1.96*SE_SAR, ymax = estimate_SAR + 1.96*SE_SAR), width = 0) +
    geom_errorbar(aes(xmin = estimate_SDM - 1.96*SE_SDM, xmax = estimate_SDM + 1.96*SE_SDM), width = 0) +
    geom_point() +
    # geom_text_repel(size = 2) +
    ylab("SAR estimate (stage 2)") +
    xlab("SDM estimate (stage 1)") +
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    theme(plot.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = "white"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          panel.background = element_rect(fill="white", color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          # these plot margins are to leave space for the panel label
          plot.margin = unit(c(1.0, 0.2, 0.2, 0.2),"cm")) +
    ggtitle(plot_title)
  
  return(plot)
}

#### Hake x Yearlings ####

# create comp df
hake_csyif_SAR_SDM_comp <- join_SDM_SAR_estimates(SAR_estimates = hake_csyif_SAR_estimates, 
                                                  SDM_estimates = SDM_estimates_hake_csyif_long)

# generate plots
## Prey Field Index (Hake x Yearlings)
hake_csyif_prey_field_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                     index_name = "prey_field_index",
                     plot_title = "Prey Field Index (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_prey_field_index_eiv_plot.png"), 
       hake_csyif_prey_field_index_eiv_plot, height = 6, width = 6)

## Prey Field Overlap (Hake x Yearlings)
hake_csyif_prey_field_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                                                             index_name = "pianka_o_csyif_prey_field",
                                                             plot_title = "Prey Field Overlap (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_prey_field_overlap_eiv_plot.png"), 
       hake_csyif_prey_field_overlap_eiv_plot, height = 6, width = 6)

## YOY Rockfish Index (Hake x Yearlings)
hake_csyif_rf_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                                                             index_name = "rf_index",
                                                             plot_title = "YOY Rockfish Index (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_rf_index_eiv_plot.png"), 
       hake_csyif_rf_index_eiv_plot, height = 6, width = 6)

## YOY Rockfish Overlap (Hake x Yearlings)
hake_csyif_rf_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                                                               index_name = "pianka_o_csyif_rf",
                                                               plot_title = "YOY Rockfish Overlap (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_rf_overlap_eiv_plot.png"), 
       hake_csyif_rf_overlap_eiv_plot, height = 6, width = 6)

## Hake Index (Hake x Yearlings)
hake_csyif_hake_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                                                             index_name = "hake_index",
                                                             plot_title = "Hake Index (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_hake_index_eiv_plot.png"), 
       hake_csyif_hake_index_eiv_plot, height = 6, width = 6)

## Hake Overlap (Hake x Yearlings)
hake_csyif_hake_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_csyif_SAR_SDM_comp, 
                                                               index_name = "pianka_o_csyif_hake",
                                                               plot_title = "Hake Overlap (Hake x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_yearlings", "hake_csyif_hake_overlap_eiv_plot.png"), 
       hake_csyif_hake_overlap_eiv_plot, height = 6, width = 6)

#### Hake x Subyearlings ####

# create comp df
hake_cssif_SAR_SDM_comp <- join_SDM_SAR_estimates(SAR_estimates = hake_cssif_SAR_estimates, 
                                                  SDM_estimates = SDM_estimates_hake_cssif_long)

# generate plots
## Prey Field Index (Hake x Subyearlings)
hake_cssif_prey_field_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                             index_name = "prey_field_index",
                                                             plot_title = "Prey Field Index (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_prey_field_index_eiv_plot.png"), 
       hake_cssif_prey_field_index_eiv_plot, height = 6, width = 6)

## Prey Field Overlap (Hake x Subyearlings)
hake_cssif_prey_field_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                               index_name = "pianka_o_cssif_prey_field",
                                                               plot_title = "Prey Field Overlap (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_prey_field_overlap_eiv_plot.png"), 
       hake_cssif_prey_field_overlap_eiv_plot, height = 6, width = 6)

## YOY Rockfish Index (Hake x Subyearlings)
hake_cssif_rf_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                     index_name = "rf_index",
                                                     plot_title = "YOY Rockfish Index (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_rf_index_eiv_plot.png"), 
       hake_cssif_rf_index_eiv_plot, height = 6, width = 6)

## YOY Rockfish Overlap (Hake x Subyearlings)
hake_cssif_rf_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                       index_name = "pianka_o_cssif_rf",
                                                       plot_title = "YOY Rockfish Overlap (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_rf_overlap_eiv_plot.png"), 
       hake_cssif_rf_overlap_eiv_plot, height = 6, width = 6)

## Hake Index (Hake x Subyearlings)
hake_cssif_hake_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                       index_name = "hake_index",
                                                       plot_title = "Hake Index (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_hake_index_eiv_plot.png"), 
       hake_cssif_hake_index_eiv_plot, height = 6, width = 6)

## Hake Overlap (Hake x Subyearlings)
hake_cssif_hake_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = hake_cssif_SAR_SDM_comp, 
                                                         index_name = "pianka_o_cssif_hake",
                                                         plot_title = "Hake Overlap (Hake x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "hake_subyearlings", "hake_cssif_hake_overlap_eiv_plot.png"), 
       hake_cssif_hake_overlap_eiv_plot, height = 6, width = 6)


#### Seabirds x Yearlings ####

# create comp df
seabirds_csyif_SAR_SDM_comp <- join_SDM_SAR_estimates(SAR_estimates = seabirds_csyif_SAR_estimates, 
                                                      SDM_estimates = SDM_estimates_seabirds_csyif_long)

# generate plots
## Prey Field Index (Seabirds x Yearlings)
seabirds_csyif_prey_field_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                                 index_name = "prey_field_index",
                                                                 plot_title = "Prey Field Index (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_prey_field_index_eiv_plot.png"), 
       seabirds_csyif_prey_field_index_eiv_plot, height = 6, width = 6)

## Prey Field Overlap (Seabirds x Yearlings)
seabirds_csyif_prey_field_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                                   index_name = "pianka_o_csyif_prey_field",
                                                                   plot_title = "Prey Field Overlap (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_prey_field_overlap_eiv_plot.png"), 
       seabirds_csyif_prey_field_overlap_eiv_plot, height = 6, width = 6)

## YOY Rockfish Index (Seabirds x Yearlings)
seabirds_csyif_rf_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                         index_name = "rf_index",
                                                         plot_title = "YOY Rockfish Index (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_rf_index_eiv_plot.png"), 
       seabirds_csyif_rf_index_eiv_plot, height = 6, width = 6)

## YOY Rockfish Overlap (Seabirds x Yearlings)
seabirds_csyif_rf_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                           index_name = "pianka_o_csyif_rf",
                                                           plot_title = "YOY Rockfish Overlap (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_rf_overlap_eiv_plot.png"), 
       seabirds_csyif_rf_overlap_eiv_plot, height = 6, width = 6)

## COMU Index (Seabirds x Yearlings)
seabirds_csyif_comu_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                           index_name = "comu_index",
                                                           plot_title = "COMU Index (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_comu_index_eiv_plot.png"), 
       seabirds_csyif_comu_index_eiv_plot, height = 6, width = 6)

## COMU Overlap (Seabirds x Yearlings)
seabirds_csyif_comu_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                             index_name = "pianka_o_csyif_comu",
                                                             plot_title = "COMU Overlap (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_comu_overlap_eiv_plot.png"), 
       seabirds_csyif_comu_overlap_eiv_plot, height = 6, width = 6)

## SOSH Index (Seabirds x Yearlings)
seabirds_csyif_sosh_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                           index_name = "sosh_index",
                                                           plot_title = "SOSH Index (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_sosh_index_eiv_plot.png"), 
       seabirds_csyif_sosh_index_eiv_plot, height = 6, width = 6)

## SOSH Overlap (Seabirds x Yearlings)
seabirds_csyif_sosh_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_csyif_SAR_SDM_comp, 
                                                             index_name = "pianka_o_csyif_sosh",
                                                             plot_title = "SOSH Overlap (Seabirds x Yearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_yearlings", "seabirds_csyif_sosh_overlap_eiv_plot.png"), 
       seabirds_csyif_sosh_overlap_eiv_plot, height = 6, width = 6)

#### Seabirds x Subyearlings ####

# create comp df
seabirds_cssif_SAR_SDM_comp <- join_SDM_SAR_estimates(SAR_estimates = seabirds_cssif_SAR_estimates, 
                                                      SDM_estimates = SDM_estimates_seabirds_cssif_long)

# generate plots
## Prey Field Index (Seabirds x Subyearlings)
seabirds_cssif_prey_field_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                                 index_name = "prey_field_index",
                                                                 plot_title = "Prey Field Index (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_prey_field_index_eiv_plot.png"), 
       seabirds_cssif_prey_field_index_eiv_plot, height = 6, width = 6)

## Prey Field Overlap (Seabirds x Subyearlings)
seabirds_cssif_prey_field_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                                   index_name = "pianka_o_cssif_prey_field",
                                                                   plot_title = "Prey Field Overlap (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_prey_field_overlap_eiv_plot.png"), 
       seabirds_cssif_prey_field_overlap_eiv_plot, height = 6, width = 6)

## YOY Rockfish Index (Seabirds x Subyearlings)
seabirds_cssif_rf_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                         index_name = "rf_index",
                                                         plot_title = "YOY Rockfish Index (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_rf_index_eiv_plot.png"), 
       seabirds_cssif_rf_index_eiv_plot, height = 6, width = 6)

## YOY Rockfish Overlap (Seabirds x Subyearlings)
seabirds_cssif_rf_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                           index_name = "pianka_o_cssif_rf",
                                                           plot_title = "YOY Rockfish Overlap (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_rf_overlap_eiv_plot.png"), 
       seabirds_cssif_rf_overlap_eiv_plot, height = 6, width = 6)

## COMU Index (Seabirds x Subyearlings)
seabirds_cssif_comu_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                           index_name = "comu_index",
                                                           plot_title = "COMU Index (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_comu_index_eiv_plot.png"), 
       seabirds_cssif_comu_index_eiv_plot, height = 6, width = 6)

## COMU Overlap (Seabirds x Subyearlings)
seabirds_cssif_comu_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                             index_name = "pianka_o_cssif_comu",
                                                             plot_title = "COMU Overlap (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_comu_overlap_eiv_plot.png"), 
       seabirds_cssif_comu_overlap_eiv_plot, height = 6, width = 6)

## SOSH Index (Seabirds x Subyearlings)
seabirds_cssif_sosh_index_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                           index_name = "sosh_index",
                                                           plot_title = "SOSH Index (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_sosh_index_eiv_plot.png"), 
       seabirds_cssif_sosh_index_eiv_plot, height = 6, width = 6)

## SOSH Overlap (Seabirds x Subyearlings)
seabirds_cssif_sosh_overlap_eiv_plot <- plot_SAR_SDM_EiV_fit(SAR_SDM_comp = seabirds_cssif_SAR_SDM_comp, 
                                                             index_name = "pianka_o_cssif_sosh",
                                                             plot_title = "SOSH Overlap (Seabirds x Subyearlings)")

ggsave(here::here("figures", "paper_figures", "eiv_fit_plots", "seabirds_subyearlings", "seabirds_cssif_sosh_overlap_eiv_plot.png"), 
       seabirds_cssif_sosh_overlap_eiv_plot, height = 6, width = 6)




