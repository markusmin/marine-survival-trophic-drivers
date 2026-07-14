## 14 - re-generate plot of parameter estimates but assuming covariates are known
# with certainty (no EiV)

# Description: This script generates Figure 4, which shows the parameter estimates
# for the various SAR models
# The model inputs are the new two prey (JSOES bongo biomass) models with the salmon index of abundance removed.

# need to run 06.2.6.1 and 06.3.4.1 scripts before running this!

## Load libraries
library(tidyverse)
library(readxl)
library(here)
library(viridis)
library(broom)
library(ggpubr)
library(sf)
library(lubridate)
library(grid)

# source the function scripts
source("R/functions/cAIC.R")
source("R/functions/rmvnorm_prec.R")
source("R/functions/sample_var.R")
source("R/functions/plot_anisotropy_MM.R")
source("R/functions/make_anisotropy_spde.R")
source("R/functions/H_matrix_prior_predictive_check_for_ggplot.R")
source("R/functions/plot_distribution.R")
source("R/functions/plot_distribution_PRS_PWCC.R")
source("R/functions/plot_distribution_jsoes_bongo.R")
# source("R/functions/plot_distribution_hake_survey.R")
source("R/functions/generate_prediction_maps.R")
source("R/functions/generate_prediction_maps_PRS_PWCC.R")
source("R/functions/generate_prediction_maps_jsoes_bongo.R")
source("R/functions/generate_prediction_maps_hake_survey.R")
source("R/functions/generate_prediction_maps_cces.R")
source("R/functions/stage2_helper_functions.R")
# source the prep scripts for each of the surveys
# source(here::here("R", "PRS_PWCC_make_mesh.R"))
# source(here::here("R", "JSOES_seabirds_make_mesh.R"))
# source(here::here("R", "JSOES_make_mesh.R"))

# # # load taxa data
# # # extract Yearling Interior Chinook
# # csyif <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_fa")
# # # extract subyearling Interior Chinook
# # cssif <- subset(jsoes_long, species == "chinook_salmon_subyearling_interior_fa")
# # # extract seabirds
# # sosh <- subset(birds_long, species == "sooty_shearwater")
# # comu <- subset(birds_long, species == "common_murre")
# # # extract bongo data
# # jsoes_bongo_biomass_cancer_crab_larvae <- read.csv(here::here("model_inputs", "jsoes_bongo_biomass_cancer_crab_larvae.csv"))
# # 
# # common_years_seabird_prey <- intersect(intersect(unique(csyif$year), # jsoes trawl
# #                                                  unique(jsoes_bongo_biomass_cancer_crab_larvae$year)), # jsoes bongo
# #                                        intersect(unique(rf$year), # PRS/PWCC
# #                                                  unique(sosh$year)))
# # 
# # common_years_prey <- intersect(
# #   unique(jsoes_bongo_biomass_cancer_crab_larvae$year), # jsoes bongo
# #   unique(rf$year) # PRS/PWCC
# # )
# 
# load the SAR data
SAR_data <- read.csv(here::here("model_inputs", "chinook_det_hist.csv"))

# inspect distribution of outmigration timing

SAR_data %>%
  mutate(BON_juv_det_date = substr(BON_juv_det_time, 1, 10)) %>%
  mutate(BON_juv_det_date = ymd(BON_juv_det_date)) %>%
  mutate(outmigration_date = yday(BON_juv_det_date)) -> SAR_data

### 06.2 seabird model

## Interior Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.2.6.1_seabird_two_prey_no_salmon_index_one_model_SAR_no_EiV", "interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_output_zscored.rda"))

interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_Obj <- interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_Obj
interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_Opt = interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_Opt
interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_report = interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2.6.1_seabird_two_prey_no_salmon_index_one_model_SAR_no_EiV", "interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_output_zscored.rda"))

interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_Obj <- interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_Obj
interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_Opt <- interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_Opt
interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_report <- interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_output$interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_report


### 06.3 hake model

## Interior Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3.4.1_hake_two_prey_no_salmon_index_one_model_SAR_no_EiV", "interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output_zscored.rda"))

interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Obj <- interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Obj
interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt <- interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt
interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_report <- interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.3.4.1_hake_two_prey_no_salmon_index_one_model_SAR_no_EiV", "interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output_zscored.rda"))

interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Obj <- interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Obj
interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt = interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt
interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_report = interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_output$interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_report


#### Compare fixed effects across models ####

# add full names and abbreviations
effect_name_df <- data.frame(effect = c("hake_index", "ov_hake", "ov_prey_field", "prey_field_index",
                                        "ov_rf", "rf_index",
                                        "sosh_index", "ov_sosh", "comu_index", "ov_comu"),
                             name = c("   Hake Abundance", "Hake Overlap", "Prey Overlap", "Prey Abundance",
                                      "YOY RF Overlap", "YOY RF Abundance", "SOSH Abundance", "SOSH Overlap",
                                      "COMU Abundance", "COMU Overlap"),
                             abbrev = c("A[h]", "O[h]", "O[p]", "A[p]",
                                        "O[r]", "A[r]", "A[s]", "O[s]",
                                        "A[c]", "O[c]"))

effect_name_df$name <- fct_rev(factor(effect_name_df$name, 
                              levels = c("Prey Abundance", "Prey Overlap", 
                                         "YOY RF Abundance", "YOY RF Overlap", 
                                         "   Hake Abundance", "Hake Overlap", 
                                         "SOSH Abundance", "SOSH Overlap",
                                         "COMU Abundance", "COMU Overlap")))

effect_name_df$abbrev <- fct_rev(factor(effect_name_df$abbrev, 
                              levels = c("A[p]", "O[p]",
                                         "A[r]", "O[r]", 
                                         "A[h]", "O[h]", 
                                         "A[s]", "O[s]",
                                         "A[c]", "O[c]")))


# function to extract fixed effects
extract_fixed_effects_uncertainty <- function(model_opt_sd){
  as.data.frame(summary(model_opt_sd)) %>% 
    rownames_to_column("parameter") %>% 
    janitor::clean_names() -> SD_summary
  
  # extract the fixed effects
  fixed_effects_SD_summary <- subset(SD_summary, grepl("beta", parameter))
  
  fixed_effects_SD_summary %>% 
    mutate(upper = estimate + 1.96 * std_error,
           lower = estimate - 1.96 * std_error) %>% 
    mutate(significance = ifelse(lower > 0 | upper < 0, 
                                 "significant", "not significant")) %>% 
    mutate(effect = gsub("beta_", "", parameter)) %>% 
    left_join(effect_name_df, by = "effect") -> fixed_effects_SD_summary
  
  return(fixed_effects_SD_summary)
}

# function to visualize fixed effects estimates

plot_FE_estimates <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
                             drop_outmigration = TRUE){
  if(drop_intercept == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")))
  }
  if(drop_outmigration == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                    "beta_outmigration_UCSF", "beta_outmigration_UCSF2")))
  }
  
  significance_colors = c("significant" = "black",
                          "not significant" = "gray80")
  
  plot <- ggplot(fixed_effects_SD_summary, aes(x = name, y = estimate,
                                               ymax = upper, ymin = lower)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 5) +
    geom_errorbar(width = 0.2) +
    # scale_color_manual(values = significance_colors) +
    guides(color = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.1, 0.1)) +
    xlab("Marine Survival Covariate") +
    ylab("Parameter Estimate") +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 


### Extract the estimates from each of the model objects

## SRF models

interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_Opt$SD)
interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt$SD)
interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_Opt$SD)
interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_Opt$SD)

### Visualize the parameter estimates

interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE_plot <- plot_FE_estimates(interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE)
interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE_plot <- plot_FE_estimates(interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE)
interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE_plot <- plot_FE_estimates(interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE)
interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE_plot <- plot_FE_estimates(interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE)

seabird_FE_comparison_plot <- ggarrange(interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE_plot,
                                        interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE_plot, 
                                        labels = c("CSYIF", "CSSIF"))

hake_FE_comparison_plot <- ggarrange(interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE_plot, 
                                     interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE_plot,
                                        labels = c("CSYIF", "CSSIF"))


#### Version for paper figure

plot_FE_estimates_forfig <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
                                     drop_outmigration = TRUE){
  if(drop_intercept == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")))
  }
  if(drop_outmigration == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                    "beta_outmigration_UCSF", "beta_outmigration_UCSF2")))
  }
  
  significance_colors = c("significant" = "black",
                          "not significant" = "gray80")
  
  plot <- ggplot(fixed_effects_SD_summary, aes(x = abbrev, y = estimate, color = significance,
                                               ymax = upper, ymin = lower)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 5, show.legend = FALSE) +
    geom_errorbar(width = 0.2, show.legend = FALSE) +
    scale_color_manual(values = significance_colors, guide = "none") +
    guides(color = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.1, 0.1)) +
    xlab("Marine Survival Covariate") +
    ylab("Parameter Estimate") +
    facet_wrap(~model, ncol = 1) +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 18),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 

plot_FE_estimates_forfig_combined <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
                                              drop_outmigration = TRUE){
  if(drop_intercept == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")))
  }
  if(drop_outmigration == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                    "beta_outmigration_UCSF", "beta_outmigration_UCSF2")))
  }
  
  significance_colors = c("significant" = "black",
                          "not significant" = "gray80")
  
  plot <- ggplot(fixed_effects_SD_summary, aes(x = abbrev, y = estimate, color = significance,
                                               ymax = upper, ymin = lower)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 5, show.legend = FALSE) +
    geom_errorbar(width = 0.2, show.legend = FALSE) +
    scale_color_manual(values = significance_colors, guide = "none") +
    guides(color = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.1, 0.1)) +
    xlab("Marine Survival Covariate") +
    ylab("Parameter Estimate") +
    facet_wrap(~model, ncol = 2, scales = "free_x") +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    # facet_wrap(~model, ncol = 2) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 


plot_FE_estimates_forfig_v2 <- function(fixed_effects_SD_summary, drop_intercept = TRUE,
                                     drop_outmigration = TRUE){
  if(drop_intercept == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")))
  }
  if(drop_outmigration == TRUE){
    fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                    "beta_outmigration_UCSF", "beta_outmigration_UCSF2")))
  }
  
  significance_colors = c("significant" = "black",
                          "not significant" = "gray65")
  
  plot <- ggplot(fixed_effects_SD_summary, aes(x = name, y = estimate, color = significance,
                                               ymax = upper, ymin = lower)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(size = 5, show.legend = FALSE) +
    geom_errorbar(width = 0.2, show.legend = FALSE) +
    scale_color_manual(values = significance_colors, guide = "none") +
    guides(color = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.1, 0.1)) +
    xlab("Marine Survival Covariate") +
    ylab("Parameter Estimate") +
    facet_wrap(~model, ncol = 2) +
    # scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_y_continuous(limits = c(-0.85, 0.85), breaks = seq(-0.75, 0.75, 0.25)) +
    coord_flip() +
    theme(panel.grid.major = element_line(color = "gray95"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          legend.key.height = unit(1.25, "cm"),
          legend.key.width = unit(1.25, "cm"),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12),
          # axis.title.x = element_text(size = 20, margin = margin(t = 10)),
          axis.title.x = element_blank(),
          # axis.title.y = element_text(size = 20, margin = margin(r = 10)),
          axis.title.y = element_blank(),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
} 


interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE$model  <- "Yearlings (1)"
interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE$model <- "Yearlings (1)"
interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE$model <- "Subyearlings (0)"
interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE$model  <- "Subyearlings (0)"


interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE %>% 
  bind_rows(interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE) %>% 
  bind_rows(interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE) %>% 
  bind_rows(interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE) -> SAR_models_FE

interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE %>% 
  bind_rows(interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE) -> SAR_hake_models_FE

interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE %>% 
  bind_rows(interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE) -> SAR_seabird_models_FE

SAR_seabird_models_FE_plot <- plot_FE_estimates_forfig_v2(SAR_seabird_models_FE)
SAR_hake_models_FE_plot <- plot_FE_estimates_forfig_v2(SAR_hake_models_FE)

fig4_SAR_parameter_estimates <- annotate_figure(ggarrange(SAR_seabird_models_FE_plot,
                                          SAR_hake_models_FE_plot,
                                          labels = c("                           (A)                                                                           (B)", 
                                                     "                           (C)                                                                           (D)"),
                                          font.label = list(size = 16, face = "plain"),
                                          label.x = 0.025, label.y = 0.985, hjust = 0,
                                   heights = c(1.2, 1),
                                   ncol = 1),
                                   left = textGrob("Marine Survival Covariate", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                   bottom = textGrob("Parameter Estimate", gp = gpar(cex = 1.3), hjust = 0.15))

ggsave(here::here("figures", "paper_figures", "no_EiV_sensitivity","fig4_SAR_parameter_estimates_v2.png"), fig4_SAR_parameter_estimates,  
       height = 6, width = 12)


#### Plot the other parameters: outmigration timing, transport ####

# extract the effect of the abundance indices from each of these
#### Plot effect of group/transport/intercept ####
subset(interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Seabirds x Yearlings") -> interior_seabirds_csyif_transport
subset(interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Hake x Yearlings") -> interior_hake_csyif_transport
subset(interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Seabirds x Subyearlings") -> interior_seabirds_cssif_transport
subset(interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Hake x Subyearlings") -> interior_hake_cssif_transport


interior_seabirds_csyif_transport %>% 
  bind_rows(interior_hake_csyif_transport) %>% 
  bind_rows(interior_seabirds_cssif_transport) %>% 
  bind_rows(interior_hake_cssif_transport) -> transport_param_estimates

significance_colors = c("significant" = "black",
                        "not significant" = "gray80")

# add model names for plotting
effect_names_for_plot <- data.frame(effect = c("UCSF", "SRF_transported", "SRF_non_transported"),
                                    effect_name = c("Upper Columbia Summer/Fall",
                                                    "Snake River Fall (transported)",
                                                    "Snake River Fall (in-river)"))

transport_param_estimates <- left_join(transport_param_estimates, effect_names_for_plot, by = "effect")

beta_group_param_estimates_plot <- ggplot(transport_param_estimates, aes(x = effect_name, y = estimate, color = significance,
                                                                           ymax = upper, ymin = lower)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size = 5, show.legend = FALSE) +
  geom_errorbar(width = 0.2, show.legend = FALSE) +
  scale_color_manual(values = significance_colors) +
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.1, 0.1)) +
  xlab("Effect of Stock Group/Transport Status") +
  ylab("Parameter Estimate") +
  facet_wrap(~model, ncol = 1, scales = "free_x") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm")) +
  coord_flip()

ggsave(here::here("figures", "paper_figures", "no_EiV_sensitivity","beta_group_param_estimates_plot.png"), beta_group_param_estimates_plot,  
       height = 8, width = 8)

#### Plot effect of outmigration date ####

subset(interior_06_2_6_1_seabird_prey_csyif_only_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                        "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>%
  mutate(model = "Seabirds x Yearlings") -> SRF_seabirds_csyif_outmigration

subset(interior_06_3_4_1_hake_prey_csyif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                                     "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>% 
  mutate(model = "Hake x Yearlings") -> SRF_hake_csyif_outmigration

subset(interior_06_2_6_1_seabird_prey_cssif_only_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                        "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>% 
  mutate(model = "Seabirds x Subyearlings") -> SRF_seabirds_cssif_outmigration

subset(interior_06_3_4_1_hake_prey_cssif_only_no_salmon_index_one_model_SAR_no_EiV_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                                     "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>% 
  mutate(model = "Hake x Subyearlings") -> SRF_hake_cssif_outmigration

outmigration_param_names_for_fig <- data.frame(parameter = c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                 "beta_outmigration_UCSF", "beta_outmigration_UCSF2"),
                                               param_name = c("Outmigration Date (SRF)", "Outmigration Date^2 (SRF)",
                                                        "Outmigration Date (UCSF)", "Outmigration Date^2 (UCSF)"))


SRF_seabirds_csyif_outmigration %>% 
  bind_rows(SRF_hake_csyif_outmigration) %>% 
  bind_rows(SRF_seabirds_cssif_outmigration) %>% 
  bind_rows(SRF_hake_cssif_outmigration) %>% 
  left_join(outmigration_param_names_for_fig, by = "parameter") -> outmigration_param_estimates


beta_outmigration_param_estimates_plot <- ggplot(outmigration_param_estimates, aes(x = model, y = estimate, color = significance,
                                                                             ymax = upper, ymin = lower)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size = 5, show.legend = FALSE) +
  geom_errorbar(width = 0.2, show.legend = FALSE) +
  scale_color_manual(values = significance_colors) +
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.1, 0.1)) +
  xlab("Model") +
  ylab("Parameter Estimate") +
  facet_wrap(~param_name, ncol = 2, scales = "free_x") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm")) +
  coord_flip()

ggsave(here::here("figures", "paper_figures", "no_EiV_sensitivity","beta_outmigration_param_estimates_plot.png"), beta_outmigration_param_estimates_plot,  
       height = 6, width = 12)

