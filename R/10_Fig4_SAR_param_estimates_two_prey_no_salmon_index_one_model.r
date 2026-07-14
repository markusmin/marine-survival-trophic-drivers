## 12_Fig4_SAR_param_estimates

# Description: This script generates Figure 4, which shows the parameter estimates
# for the various SAR models
# The model inputs are the new two prey (JSOES bongo biomass) models with the salmon index of abundance removed.

# need to run 06.2.6 and 06.3.4 scripts before running this!

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
load(here::here("R", "06_stage2_SAR", "06.2.6_seabird_two_prey_no_salmon_index_one_model_SAR", "interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Obj <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Obj
interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt
interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.2.6_seabird_two_prey_no_salmon_index_one_model_SAR", "interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Obj <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Obj
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_output$interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report


### 06.3 hake model

## Interior Models

# CSYIF output
load(here::here("R", "06_stage2_SAR", "06.3.4_hake_two_prey_no_salmon_index_one_model_SAR", "interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Obj <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Obj
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report

# CSSIF output
load(here::here("R", "06_stage2_SAR", "06.3.4_hake_two_prey_no_salmon_index_one_model_SAR", "interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output_zscored.rda"))

interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Obj <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Obj
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_output$interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report


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

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_Opt$SD)
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_Opt$SD)
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_Opt$SD)
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE <- extract_fixed_effects_uncertainty(model_opt_sd = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_Opt$SD)

### Visualize the parameter estimates

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE_plot <- plot_FE_estimates(interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE)
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE_plot <- plot_FE_estimates(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE)
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE_plot <- plot_FE_estimates(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE)
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE_plot <- plot_FE_estimates(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE)

seabird_FE_comparison_plot <- ggarrange(interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE_plot,
                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE_plot, 
                                        labels = c("CSYIF", "CSSIF"))

ggsave(here::here("figures", "paper_figures", "one_model","parameter_estimate_plots", "seabird_FE_comparison_plot.png"), seabird_FE_comparison_plot,  
       height = 12, width = 16)

hake_FE_comparison_plot <- ggarrange(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE_plot, 
                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE_plot,
                                        labels = c("CSYIF", "CSSIF"))

ggsave(here::here("figures", "paper_figures", "one_model","parameter_estimate_plots", "hake_FE_comparison_plot.png"), hake_FE_comparison_plot,  
       height = 12, width = 16)


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


interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$model  <- "Yearlings (1)"
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$model <- "Yearlings (1)"
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$model <- "Subyearlings (0)"
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$model  <- "Subyearlings (0)"


interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE %>% 
  bind_rows(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE) %>% 
  bind_rows(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE) %>% 
  bind_rows(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE) -> SAR_models_FE

interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE %>% 
  bind_rows(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE) -> SAR_hake_models_FE

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE %>% 
  bind_rows(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE) -> SAR_seabird_models_FE

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


ggsave(here::here("figures", "paper_figures", "one_model","fig4_SAR_parameter_estimates_v2.png"), fig4_SAR_parameter_estimates,  
       height = 6, width = 12)


#### Plot the other parameters: outmigration timing, transport ####

# extract the effect of the abundance indices from each of these
#### Plot effect of group/transport/intercept ####
subset(interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Seabirds x Yearlings") -> interior_seabirds_csyif_transport
subset(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Hake x Yearlings") -> interior_hake_csyif_transport
subset(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
  mutate(model = "Seabirds x Subyearlings") -> interior_seabirds_cssif_transport
subset(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE, parameter %in% c("beta_UCSF", "beta_SRF_transported", "beta_SRF_non_transported")) %>% 
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

ggsave(here::here("figures", "paper_figures", "one_model","beta_group_param_estimates_plot.png"), beta_group_param_estimates_plot,  
       height = 8, width = 8)

#### Plot effect of outmigration date ####

subset(interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                        "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>%
  mutate(model = "Seabirds x Yearlings") -> SRF_seabirds_csyif_outmigration

subset(interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                                     "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>% 
  mutate(model = "Hake x Yearlings") -> SRF_hake_csyif_outmigration

subset(interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
                                                                                        "beta_outmigration_UCSF", "beta_outmigration_UCSF2")) %>% 
  mutate(model = "Seabirds x Subyearlings") -> SRF_seabirds_cssif_outmigration

subset(interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE, parameter %in% c("beta_outmigration_SRF", "beta_outmigration_SRF2",
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

ggsave(here::here("figures", "paper_figures", "one_model","beta_outmigration_param_estimates_plot.png"), beta_outmigration_param_estimates_plot,  
       height = 6, width = 12)

#### Plot model fit to data ####

# function to plot fits to data
plot_model_predictions_comp <- function(model_predictions_comp_data){
  model_predictions_comp_plot <- ggplot(model_predictions_comp_data, aes(x = run_year, y = SAR)) +
    # geom_errorbar(data = model_predictions_comp_data, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "#2166ac") +
    geom_errorbar(data = model_predictions_comp_data, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "black") +
    geom_point(aes(x = run_year, shape = "Empirical", color = "Empirical"), alpha = 0.5, size = 5) +
    geom_point(data = model_predictions_comp_data, aes(x = run_year, y = predicted_prob, shape = "Predicted", color = "Predicted"), size = 5, stroke = 0.7) +
    scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 13)) +
    # scale_color_manual(name = NULL, values = c("Empirical" = "#2166ac", "Predicted" = "#b2182b")) +
    scale_color_manual(name = NULL, values = c("Empirical" = "black", "Predicted" = "#b2182b")) +
    coord_cartesian(ylim = c(0, 0.1), clip = "off") +
    xlab("Run Year") +
    ylab("Marine Survival") +
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
  
  return(model_predictions_comp_plot)
}



# prep SAR data

#### Snake River Fall Chinook

# extract Snake River Fall Chinook from SAR data
SRF <- subset(SAR_data, run_name == "Fall" & group == "Snake River")

#### subset by outmigration date
hist(SRF$outmigration_date)
summary(SRF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# SRF_subset <- filter(SRF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
SRF_subset <- filter(SRF, outmigration_date <= 167)
nrow(SRF_subset)/nrow(SRF)
# this leaves us with 54% of our initial dataset


#### Join with the transport data

transport_data_export <- read.csv(here::here("model_inputs", "chinook_transport.csv"))

# join the juvenile transport data with the SAR data
SRF_subset %>% 
  left_join(transport_data_export, join_by(tag_code == tag_id)) -> SRF_subset

SRF_subset$transport[is.na(SRF_subset$transport)] <- 0

# inspect number of transported fish
table(SRF_subset$transport)


#### Drop the few natural origin fish

# fish with rear type U or H are hatchery fish
SRF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> SRF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
SRF_subset <- subset(SRF_subset, rear_type_code != "W")


#### Upper Columbia Summer/Fall Chinook

# extract Snake River Fall Chinook from SAR data
UCSF <- subset(SAR_data, run_name %in% c("Summer", "Fall") & group == "Upper Columbia")

#### subset by outmigration date
hist(UCSF$outmigration_date)
summary(UCSF$outmigration_date)

# restrict the analysis to fish that could have theoretically been caught in the June JSOES survey.
# we will call this April 15 - June 15 (for now!)
# day 105 is April 15 in a non-leap year and day 167 is June 15 in a leap year
# UCSF_subset <- filter(UCSF, outmigration_date >= 105 & outmigration_date <= 167)

# let's just have a final cutoff - based on the fact that fish can hang out
# off the coast for a while, it's less justified to have a front-end cutoff date
UCSF_subset <- filter(UCSF, outmigration_date <= 167)
nrow(UCSF_subset)/nrow(UCSF)
# this leaves us with 72% of our initial dataset

#### Drop the few natural origin fish (28 of ~97,500)

table(UCSF_subset$rear_type_code)

# fish with rear type U or H are hatchery fish
UCSF_subset %>% 
  mutate(rear_numeric = ifelse(rear_type_code %in% c("H", "U"), 1, 0)) -> UCSF_subset

# ok, so they're all hatchery. That's perhaps not surprising, given what we see in the JSOES survey is basically all hatchery
# Let's only keep hatchery fish
UCSF_subset <- subset(UCSF_subset, rear_type_code != "W")

# JOIN SRF AND UCSF TOGETHER
bind_rows(SRF_subset, UCSF_subset) -> interior_subset


#### function for analytical variance calculation
var_2_rv <- function(est1, var1, est2, var2){
  total_var <- var1*var2 + var1*est2^2 + var2*est1^2
  return(total_var)
}


# look at our latent estimates vs. the input data
# ok good, it's recovering the input well
# plot(x = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$csyif_index_t_latent, y = (SRF_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(SRF_v1_seabird_prey_predictors$csyif_index_of_abundance))


#### interior - Yearlings x Seabirds ####
interior_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> interior_subset_SAR

UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_non_transported_subset_SAR

SRF_subset %>% 
  filter(transport == 1) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_transported_subset_SAR

interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$variance <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$std_error^2

## variance of beta_UCSF parameter
var_beta_UCSF <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_UCSF", "variance"]

## variance of beta_SRF_transported parameter
var_beta_SRF_transported <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported", "variance"]

## variance of beta_SRF_non_transported parameter
var_beta_SRF_non_transported <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported", "variance"]

# ## csyif index of abundance
# # expected value of beta_csyif index parameter
# est_beta_csyif_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index", "estimate"]
# 
# # variance of beta_csyif index parameter
# var_beta_csyif_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index", "variance"]
# 
# # expected value of csyif_index of abundance predictor
# est_csyif_index_of_abundance <- (interior_v1_seabird_prey_predictors$csyif_index_of_abundance - mean(interior_v1_seabird_prey_predictors$csyif_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$csyif_index_of_abundance)
# 
# # variance of csyif index of abundance predictor
# var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (interior_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t - mean(interior_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_prey_field_t))^2))

## rf index
# expected value of beta_rf index parameter
est_beta_rf_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "estimate"]

# variance of beta_rf index parameter
var_beta_rf_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "variance"]

# expected value of rf_index of abundance predictor
est_rf_index_of_abundance <- (interior_v1_seabird_prey_predictors$rf_index_of_abundance - mean(interior_v1_seabird_prey_predictors$rf_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$rf_index_of_abundance)

# variance of rf index of abundance predictor
var_rf_index_of_abundance <- sqrt(diag(rf_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$rf_index_of_abundance))^2))

## rf overlap
# expected value of beta_rf parameter
est_beta_ov_rf <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "estimate"]

# variance of beta_rf parameter
var_beta_ov_rf <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "variance"]

# expected value of rf overlap predictor
est_rf_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_csyif_rf_t - mean(interior_v1_seabird_prey_predictors$pianka_o_csyif_rf_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_rf_t)

# variance of rf index of abundance predictor
var_rf_overlap <- sqrt(diag(pianka_o_csyif_rf_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_rf_t))^2))

## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (interior_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(interior_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t - mean(interior_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_csyif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (interior_v1_seabird_prey_predictors$comu_index_of_abundance - mean(interior_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_csyif_comu_t - mean(interior_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_csyif_comu_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_csyif_comu_t))^2))


#### Calculate and plot total variance by year
# get predictors in a df
interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                                        # csyif_index_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$csyif_index_t_latent,
                                                                                        pianka_o_csyif_sosh_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$pianka_o_csyif_sosh_t_latent,
                                                                                        sosh_index_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$sosh_index_t_latent,
                                                                                        pianka_o_csyif_comu_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$pianka_o_csyif_comu_t_latent,
                                                                                        comu_index_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$comu_index_t_latent,
                                                                                        pianka_o_csyif_prey_field_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                                                                                        prey_field_index_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$prey_field_index_t_latent,
                                                                                        pianka_o_csyif_rf_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$pianka_o_csyif_rf_t_latent,
                                                                                        rf_index_t_latent = interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_report$rf_index_t_latent)

# Get the total variance by group/transport (just change the beta term)

### UCSF ###
UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var <- var_beta_UCSF + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                      linear_predictor = 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_UCSF","estimate"] +
                                                                                        # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                        # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent)

UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$var <- UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var

UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions

UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp


UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot.png"), UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_non_transported ###
SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var <- var_beta_SRF_non_transported + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                                     linear_predictor = 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported","estimate"] +
                                                                                                       # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                                       # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent)

SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$var <- SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var

SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions

SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions %>% 
  left_join(SRF_non_transported_subset_SAR, by = "run_year") -> SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp


SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot.png"), SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_transported ###
SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var <- var_beta_SRF_transported + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                                 linear_predictor = 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported","estimate"] +
                                                                                                   # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                                   # interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_sosh_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_comu_t_latent)

SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions$var <- SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_SEcov_pred_var

SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions

SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions %>% 
  left_join(SRF_transported_subset_SAR, by = "run_year") -> SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp


SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot.png"), SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)








#### interior - Subyearlings x Seabirds ####
interior_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> interior_subset_SAR

UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_non_transported_subset_SAR

SRF_subset %>% 
  filter(transport == 1) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_transported_subset_SAR

interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$variance <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$std_error^2

## variance of beta_UCSF parameter
var_beta_UCSF <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_UCSF", "variance"]

## variance of beta_SRF_transported parameter
var_beta_SRF_transported <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported", "variance"]

## variance of beta_SRF_non_transported parameter
var_beta_SRF_non_transported <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported", "variance"]

# ## cssif index of abundance
# # expected value of beta_cssif index parameter
# est_beta_cssif_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index", "estimate"]
# 
# # variance of beta_cssif index parameter
# var_beta_cssif_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index", "variance"]
# 
# # expected value of cssif_index of abundance predictor
# est_cssif_index_of_abundance <- (interior_v1_seabird_prey_predictors$cssif_index_of_abundance - mean(interior_v1_seabird_prey_predictors$cssif_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$cssif_index_of_abundance)
# 
# # variance of cssif index of abundance predictor
# var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (interior_v1_seabird_prey_predictors$prey_field_index_of_abundance - mean(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t - mean(interior_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_prey_field_t))^2))

## rf index
# expected value of beta_rf index parameter
est_beta_rf_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "estimate"]

# variance of beta_rf index parameter
var_beta_rf_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "variance"]

# expected value of rf_index of abundance predictor
est_rf_index_of_abundance <- (interior_v1_seabird_prey_predictors$rf_index_of_abundance - mean(interior_v1_seabird_prey_predictors$rf_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$rf_index_of_abundance)

# variance of rf index of abundance predictor
var_rf_index_of_abundance <- sqrt(diag(rf_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$rf_index_of_abundance))^2))

## rf overlap
# expected value of beta_rf parameter
est_beta_ov_rf <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "estimate"]

# variance of beta_rf parameter
var_beta_ov_rf <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "variance"]

# expected value of rf overlap predictor
est_rf_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_cssif_rf_t - mean(interior_v1_seabird_prey_predictors$pianka_o_cssif_rf_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_rf_t)

# variance of rf index of abundance predictor
var_rf_overlap <- sqrt(diag(pianka_o_cssif_rf_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_rf_t))^2))

## sosh index
# expected value of beta_sosh index parameter
est_beta_sosh_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index", "estimate"]

# variance of beta_sosh index parameter
var_beta_sosh_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index", "variance"]

# expected value of sosh_index of abundance predictor
est_sosh_index_of_abundance <- (interior_v1_seabird_prey_predictors$sosh_index_of_abundance - mean(interior_v1_seabird_prey_predictors$sosh_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$sosh_index_of_abundance)

# variance of sosh index of abundance predictor
var_sosh_index_of_abundance <- sqrt(diag(sosh_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$sosh_index_of_abundance))^2))

## sosh overlap
# expected value of beta_sosh parameter
est_beta_ov_sosh <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh", "estimate"]

# variance of beta_sosh parameter
var_beta_ov_sosh <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh", "variance"]

# expected value of sosh overlap predictor
est_sosh_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t - mean(interior_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t)

# variance of sosh index of abundance predictor
var_sosh_overlap <- sqrt(diag(pianka_o_cssif_sosh_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_sosh_t))^2))

## comu index
# expected value of beta_comu index parameter
est_beta_comu_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index", "estimate"]

# variance of beta_comu index parameter
var_beta_comu_index <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index", "variance"]

# expected value of comu_index of abundance predictor
est_comu_index_of_abundance <- (interior_v1_seabird_prey_predictors$comu_index_of_abundance - mean(interior_v1_seabird_prey_predictors$comu_index_of_abundance))/sd(interior_v1_seabird_prey_predictors$comu_index_of_abundance)

# variance of comu index of abundance predictor
var_comu_index_of_abundance <- sqrt(diag(comu_index_of_abundance_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$comu_index_of_abundance))^2))

## comu overlap
# expected value of beta_comu parameter
est_beta_ov_comu <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu", "estimate"]

# variance of beta_comu parameter
var_beta_ov_comu <- interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu", "variance"]

# expected value of comu overlap predictor
est_comu_overlap <- (interior_v1_seabird_prey_predictors$pianka_o_cssif_comu_t - mean(interior_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_comu_t)

# variance of comu index of abundance predictor
var_comu_overlap <- sqrt(diag(pianka_o_cssif_comu_cov_matrix_common_years_seabird_prey * (1/sd(interior_v1_seabird_prey_predictors$pianka_o_cssif_comu_t))^2))


#### Calculate and plot total variance by year
# get predictors in a df
interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df <- data.frame(year = common_years_seabird_prey,
                                                                                        # cssif_index_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$cssif_index_t_latent,
                                                                                        pianka_o_cssif_sosh_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$pianka_o_cssif_sosh_t_latent,
                                                                                        sosh_index_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$sosh_index_t_latent,
                                                                                        pianka_o_cssif_comu_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$pianka_o_cssif_comu_t_latent,
                                                                                        comu_index_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$comu_index_t_latent,
                                                                                        pianka_o_cssif_prey_field_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$pianka_o_cssif_prey_field_t_latent,
                                                                                        prey_field_index_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$prey_field_index_t_latent,
                                                                                        pianka_o_cssif_rf_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$pianka_o_cssif_rf_t_latent,
                                                                                        rf_index_t_latent = interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_report$rf_index_t_latent)

# Get the total variance by group/transport (just change the beta term)

### UCSF ###
UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var <- var_beta_UCSF + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                      linear_predictor = 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_UCSF","estimate"] +
                                                                                        # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                        # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_sosh_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                        interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_comu_t_latent)

UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$var <- UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var

UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions

UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp


UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot.png"), UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_non_transported ###
SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var <- var_beta_SRF_non_transported + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                                     linear_predictor = 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported","estimate"] +
                                                                                                       # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                                       # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_sosh_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                                       interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_comu_t_latent)

SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$var <- SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var

SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions

SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions %>% 
  left_join(SRF_non_transported_subset_SAR, by = "run_year") -> SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp


SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot.png"), SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_transported ###
SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var <- var_beta_SRF_transported + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_sosh_index, var1 = var_beta_sosh_index,
           est2 = est_sosh_index_of_abundance, var2 = var_sosh_index_of_abundance) +
  var_2_rv(est1 = est_beta_comu_index, var1 = var_beta_comu_index,
           est2 = est_comu_index_of_abundance, var2 = var_comu_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_sosh, var1 = var_beta_ov_sosh,
           est2 = est_sosh_overlap, var2 = var_sosh_overlap) +
  var_2_rv(est1 = est_beta_ov_comu, var1 = var_beta_ov_comu,
           est2 = est_comu_overlap, var2 = var_comu_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_seabird_prey,
                                                                                                 linear_predictor = 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported","estimate"] +
                                                                                                   # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                                   # interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_sosh_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$sosh_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_sosh","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_sosh_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_comu_index","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$comu_index_t_latent +
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE[interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_FE$parameter == "beta_ov_comu","estimate"] * 
                                                                                                   interior_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_comu_t_latent)

SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions$var <- SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_SEcov_pred_var

SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions

SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions %>% 
  left_join(SRF_transported_subset_SAR, by = "run_year") -> SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp


SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot.png"), SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### interior - Yearlings x hake ####
interior_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> interior_subset_SAR

UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_non_transported_subset_SAR

SRF_subset %>% 
  filter(transport == 1) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_transported_subset_SAR

interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$variance <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$std_error^2

## variance of beta_UCSF parameter
var_beta_UCSF <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_UCSF", "variance"]

## variance of beta_SRF_transported parameter
var_beta_SRF_transported <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported", "variance"]

## variance of beta_SRF_non_transported parameter
var_beta_SRF_non_transported <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported", "variance"]

# ## csyif index of abundance
# # expected value of beta_csyif index parameter
# est_beta_csyif_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index", "estimate"]
# 
# # variance of beta_csyif index parameter
# var_beta_csyif_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index", "variance"]
# 
# # expected value of csyif_index of abundance predictor
# est_csyif_index_of_abundance <- (interior_v1_hake_prey_predictors$csyif_index_of_abundance - mean(interior_v1_hake_prey_predictors$csyif_index_of_abundance))/sd(interior_v1_hake_prey_predictors$csyif_index_of_abundance)
# 
# # variance of csyif index of abundance predictor
# var_csyif_index_of_abundance <- sqrt(diag(csyif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$csyif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (interior_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(interior_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(interior_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (interior_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t - mean(interior_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_csyif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_prey_field_t))^2))

## rf index
# expected value of beta_rf index parameter
est_beta_rf_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "estimate"]

# variance of beta_rf index parameter
var_beta_rf_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "variance"]

# expected value of rf_index of abundance predictor
est_rf_index_of_abundance <- (interior_v1_hake_prey_predictors$rf_index_of_abundance - mean(interior_v1_hake_prey_predictors$rf_index_of_abundance))/sd(interior_v1_hake_prey_predictors$rf_index_of_abundance)

# variance of rf index of abundance predictor
var_rf_index_of_abundance <- sqrt(diag(rf_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$rf_index_of_abundance))^2))

## rf overlap
# expected value of beta_rf parameter
est_beta_ov_rf <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "estimate"]

# variance of beta_rf parameter
var_beta_ov_rf <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "variance"]

# expected value of rf overlap predictor
est_rf_overlap <- (interior_v1_hake_prey_predictors$pianka_o_csyif_rf_t - mean(interior_v1_hake_prey_predictors$pianka_o_csyif_rf_t))/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_rf_t)

# variance of rf index of abundance predictor
var_rf_overlap <- sqrt(diag(pianka_o_csyif_rf_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_rf_t))^2))

## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (interior_v1_hake_prey_predictors$hake_index_of_abundance - mean(interior_v1_hake_prey_predictors$hake_index_of_abundance))/sd(interior_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (interior_v1_hake_prey_predictors$pianka_o_csyif_hake_t - mean(interior_v1_hake_prey_predictors$pianka_o_csyif_hake_t))/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_csyif_hake_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_csyif_hake_t))^2))


#### Calculate and plot total variance by year
# get predictors in a df
interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                                     # csyif_index_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$csyif_index_t_latent,
                                                                                     pianka_o_csyif_hake_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_csyif_hake_t_latent,
                                                                                     hake_index_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$hake_index_t_latent,
                                                                                     pianka_o_csyif_prey_field_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_csyif_prey_field_t_latent,
                                                                                     prey_field_index_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$prey_field_index_t_latent,
                                                                                     pianka_o_csyif_rf_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_csyif_rf_t_latent,
                                                                                     rf_index_t_latent = interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_report$rf_index_t_latent)

# Get the total variance by group/transport (just change the beta term)

### UCSF ###
UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_UCSF + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                   linear_predictor = 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_UCSF","estimate"] +
                                                                                     # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                     # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_hake_t_latent)

UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions

UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp


UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_non_transported ###
SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_SRF_non_transported + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                                  linear_predictor = 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported","estimate"] +
                                                                                                    # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                                    # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_hake_t_latent)

SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions

SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(SRF_non_transported_subset_SAR, by = "run_year") -> SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp


SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_transported ###
SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_SRF_transported + 
  # var_2_rv(est1 = est_beta_csyif_index, var1 = var_beta_csyif_index,
  #          est2 = est_csyif_index_of_abundance, var2 = var_csyif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                              linear_predictor = 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported","estimate"] +
                                                                                                # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_csyif_index","estimate"] * 
                                                                                                # interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$csyif_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_prey_field_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_rf_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_csyif_hake_t_latent)

SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions

SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(SRF_transported_subset_SAR, by = "run_year") -> SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp


SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)








#### interior - Subyearlings x hake ####
interior_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> interior_subset_SAR

UCSF_subset %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> UCSF_subset_SAR

SRF_subset %>% 
  filter(transport == 0) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_non_transported_subset_SAR

SRF_subset %>% 
  filter(transport == 1) %>% 
  group_by(run_year) %>% 
  summarise(SAR = mean(adult_det),
            N = n()) -> SRF_transported_subset_SAR

interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$variance <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$std_error^2

## variance of beta_UCSF parameter
var_beta_UCSF <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_UCSF", "variance"]

## variance of beta_SRF_transported parameter
var_beta_SRF_transported <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported", "variance"]

## variance of beta_SRF_non_transported parameter
var_beta_SRF_non_transported <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported", "variance"]

# ## cssif index of abundance
# # expected value of beta_cssif index parameter
# est_beta_cssif_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index", "estimate"]
# 
# # variance of beta_cssif index parameter
# var_beta_cssif_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index", "variance"]
# 
# # expected value of cssif_index of abundance predictor
# est_cssif_index_of_abundance <- (interior_v1_hake_prey_predictors$cssif_index_of_abundance - mean(interior_v1_hake_prey_predictors$cssif_index_of_abundance))/sd(interior_v1_hake_prey_predictors$cssif_index_of_abundance)
# 
# # variance of cssif index of abundance predictor
# var_cssif_index_of_abundance <- sqrt(diag(cssif_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$cssif_index_of_abundance))^2))

## prey_field index
# expected value of beta_prey_field index parameter
est_beta_prey_field_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "estimate"]

# variance of beta_prey_field index parameter
var_beta_prey_field_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index", "variance"]

# expected value of prey_field_index of abundance predictor
est_prey_field_index_of_abundance <- (interior_v1_hake_prey_predictors$prey_field_index_of_abundance - mean(interior_v1_hake_prey_predictors$prey_field_index_of_abundance))/sd(interior_v1_hake_prey_predictors$prey_field_index_of_abundance)

# variance of prey_field index of abundance predictor
var_prey_field_index_of_abundance <- sqrt(diag(prey_field_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$prey_field_index_of_abundance))^2))

## prey_field overlap
# expected value of beta_prey_field parameter
est_beta_ov_prey_field <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "estimate"]

# variance of beta_prey_field parameter
var_beta_ov_prey_field <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field", "variance"]

# expected value of prey_field overlap predictor
est_prey_field_overlap <- (interior_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t - mean(interior_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t)

# variance of prey_field index of abundance predictor
var_prey_field_overlap <- sqrt(diag(pianka_o_cssif_prey_field_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_prey_field_t))^2))

## rf index
# expected value of beta_rf index parameter
est_beta_rf_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "estimate"]

# variance of beta_rf index parameter
var_beta_rf_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index", "variance"]

# expected value of rf_index of abundance predictor
est_rf_index_of_abundance <- (interior_v1_hake_prey_predictors$rf_index_of_abundance - mean(interior_v1_hake_prey_predictors$rf_index_of_abundance))/sd(interior_v1_hake_prey_predictors$rf_index_of_abundance)

# variance of rf index of abundance predictor
var_rf_index_of_abundance <- sqrt(diag(rf_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$rf_index_of_abundance))^2))

## rf overlap
# expected value of beta_rf parameter
est_beta_ov_rf <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "estimate"]

# variance of beta_rf parameter
var_beta_ov_rf <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf", "variance"]

# expected value of rf overlap predictor
est_rf_overlap <- (interior_v1_hake_prey_predictors$pianka_o_cssif_rf_t - mean(interior_v1_hake_prey_predictors$pianka_o_cssif_rf_t))/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_rf_t)

# variance of rf index of abundance predictor
var_rf_overlap <- sqrt(diag(pianka_o_cssif_rf_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_rf_t))^2))

## hake index
# expected value of beta_hake index parameter
est_beta_hake_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index", "estimate"]

# variance of beta_hake index parameter
var_beta_hake_index <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index", "variance"]

# expected value of hake_index of abundance predictor
est_hake_index_of_abundance <- (interior_v1_hake_prey_predictors$hake_index_of_abundance - mean(interior_v1_hake_prey_predictors$hake_index_of_abundance))/sd(interior_v1_hake_prey_predictors$hake_index_of_abundance)

# variance of hake index of abundance predictor
var_hake_index_of_abundance <- sqrt(diag(hake_index_of_abundance_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$hake_index_of_abundance))^2))

## hake overlap
# expected value of beta_hake parameter
est_beta_ov_hake <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake", "estimate"]

# variance of beta_hake parameter
var_beta_ov_hake <- interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake", "variance"]

# expected value of hake overlap predictor
est_hake_overlap <- (interior_v1_hake_prey_predictors$pianka_o_cssif_hake_t - mean(interior_v1_hake_prey_predictors$pianka_o_cssif_hake_t))/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_hake_t)

# variance of hake index of abundance predictor
var_hake_overlap <- sqrt(diag(pianka_o_cssif_hake_cov_matrix_common_years_hake_prey * (1/sd(interior_v1_hake_prey_predictors$pianka_o_cssif_hake_t))^2))


#### Calculate and plot total variance by year
# get predictors in a df
interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df <- data.frame(year = common_years_hake_prey,
                                                                                     # cssif_index_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$cssif_index_t_latent,
                                                                                     pianka_o_cssif_hake_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_cssif_hake_t_latent,
                                                                                     hake_index_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$hake_index_t_latent,
                                                                                     pianka_o_cssif_prey_field_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_cssif_prey_field_t_latent,
                                                                                     prey_field_index_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$prey_field_index_t_latent,
                                                                                     pianka_o_cssif_rf_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$pianka_o_cssif_rf_t_latent,
                                                                                     rf_index_t_latent = interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_report$rf_index_t_latent)

# Get the total variance by group/transport (just change the beta term)

### UCSF ###
UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_UCSF + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)

# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                   linear_predictor = 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_UCSF","estimate"] +
                                                                                     # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                     # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                     interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_hake_t_latent)

UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions

UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(UCSF_subset_SAR, by = "run_year") -> UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp


UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_non_transported ###
SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_SRF_non_transported + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                                  linear_predictor = 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_non_transported","estimate"] +
                                                                                                    # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                                    # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                                    interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_hake_t_latent)

SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions

SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(SRF_non_transported_subset_SAR, by = "run_year") -> SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp


SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)

### SRF_transported ###
SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var <- var_beta_SRF_transported + 
  # var_2_rv(est1 = est_beta_cssif_index, var1 = var_beta_cssif_index,
  #          est2 = est_cssif_index_of_abundance, var2 = var_cssif_index_of_abundance) +
  var_2_rv(est1 = est_beta_prey_field_index, var1 = var_beta_prey_field_index,
           est2 = est_prey_field_index_of_abundance, var2 = var_prey_field_index_of_abundance) +
  var_2_rv(est1 = est_beta_rf_index, var1 = var_beta_rf_index,
           est2 = est_rf_index_of_abundance, var2 = var_rf_index_of_abundance) +
  var_2_rv(est1 = est_beta_hake_index, var1 = var_beta_hake_index,
           est2 = est_hake_index_of_abundance, var2 = var_hake_index_of_abundance) +
  var_2_rv(est1 = est_beta_ov_prey_field, var1 = var_beta_ov_prey_field,
           est2 = est_prey_field_overlap, var2 = var_prey_field_overlap) +
  var_2_rv(est1 = est_beta_ov_rf, var1 = var_beta_ov_rf,
           est2 = est_rf_overlap, var2 = var_rf_overlap) +
  var_2_rv(est1 = est_beta_ov_hake, var1 = var_beta_ov_hake,
           est2 = est_hake_overlap, var2 = var_hake_overlap)


# Calculate MLE predictions, then calculate upper and lower bounds based on the analytical variance estimate
SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions <- data.frame(run_year = common_years_hake_prey,
                                                                                              linear_predictor = 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_SRF_transported","estimate"] +
                                                                                                # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_cssif_index","estimate"] * 
                                                                                                # interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$cssif_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_prey_field_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$prey_field_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_prey_field","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_prey_field_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_rf_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$rf_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_rf","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_rf_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_hake_index","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$hake_index_t_latent +
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE[interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_FE$parameter == "beta_ov_hake","estimate"] * 
                                                                                                interior_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_predictors_df$pianka_o_cssif_hake_t_latent)

SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$predicted_prob <- inv.logit(SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$linear_predictor)

SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions$var <- SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_SEcov_pred_var

SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_MLE_predictions %>% 
  mutate(linear_predictor_upper = linear_predictor + sqrt(var)*1.96,
         linear_predictor_lower = linear_predictor - sqrt(var)*1.96) %>% 
  mutate(predicted_prob_upper = inv.logit(linear_predictor_upper),
         predicted_prob_lower = inv.logit(linear_predictor_lower)) -> SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions

SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions %>% 
  left_join(SRF_transported_subset_SAR, by = "run_year") -> SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp


SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot <- plot_model_predictions_comp(SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp)

ggsave(here::here("figures", "paper_figures", "one_model","SAR_prediction_plots", "SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot.png"), SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp_plot,  
       height = 6, width = 8)


#### Combine figure ####

### UCSF ###
UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp$model <- "Yearlings x Seabirds"
UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model <- "Yearlings x Hake"
UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Seabirds"
UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Hake"


UCSF_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp %>% 
  bind_rows(UCSF_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp) -> UCSF_SAR_seabird_model_predictions_comp

UCSF_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp %>%
  bind_rows(UCSF_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp) -> UCSF_SAR_hake_model_predictions_comp

UCSF_SAR_seabird_model_predictions_comp %>% 
  bind_rows(UCSF_SAR_hake_model_predictions_comp) -> UCSF_SAR_model_predictions_comp

UCSF_SAR_model_predictions_comp$model <- factor(UCSF_SAR_model_predictions_comp$model, 
                                                levels = c("Subyearlings x Seabirds",
                                                           "Yearlings x Seabirds",
                                                           "Subyearlings x Hake",
                                                           "Yearlings x Hake"))
labels_df <- data.frame(model = c("Subyearlings x Seabirds",
                                  "Yearlings x Seabirds",
                                  "Subyearlings x Hake",
                                  "Yearlings x Hake"),
                        label = c("(A) Subyearlings x Seabirds",
                                  "(B) Yearlings x Seabirds",
                                  "(C) Subyearlings x Hake",
                                  "(D) Yearlings x Hake"))

UCSF_SAR_model_predictions_comp %>% 
  left_join(labels_df, by = "model") -> UCSF_SAR_model_predictions_comp


UCSF_SAR_model_predictions_comp_plot <- ggplot(UCSF_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  # geom_errorbar(data = UCSF_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "#2166ac") +
  geom_errorbar(data = UCSF_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "blue") +
  geom_point(aes(x = run_year, shape = "Empirical", color = "Empirical", size = N), alpha = 0.5) +
  geom_point(data = UCSF_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted", color = "Predicted"), size = 5, stroke = 0.7) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 13)) +
  # scale_color_manual(name = NULL, values = c("Empirical" = "#2166ac", "Predicted" = "#b2182b")) +
  scale_color_manual(name = NULL, values = c("Empirical" = "black", "Predicted" = "blue")) +
  scale_size_continuous(breaks = c(10, 1000, 10000)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  guides(size = guide_legend(order = 2)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.position = "bottom",
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, hjust = 0),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))+
  facet_wrap(~label, ncol = 2)

ggsave(here::here("figures", "paper_figures", "one_model","UCSF_SAR_model_predictions_comp_plot.png"), 
       UCSF_SAR_model_predictions_comp_plot,  
       height = 6, width = 12)

### SRF transported ###
SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp$model <- "Yearlings x Seabirds"
SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model <- "Yearlings x Hake"
SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Seabirds"
SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Hake"


SRF_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp %>% 
  bind_rows(SRF_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp) -> SRF_transported_SAR_seabird_model_predictions_comp

SRF_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp %>%
  bind_rows(SRF_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp) -> SRF_transported_SAR_hake_model_predictions_comp

SRF_transported_SAR_seabird_model_predictions_comp %>% 
  bind_rows(SRF_transported_SAR_hake_model_predictions_comp) -> SRF_transported_SAR_model_predictions_comp

SRF_transported_SAR_model_predictions_comp$model <- factor(SRF_transported_SAR_model_predictions_comp$model, 
                                                           levels = c("Subyearlings x Seabirds",
                                                                      "Yearlings x Seabirds",
                                                                      "Subyearlings x Hake",
                                                                      "Yearlings x Hake"))
labels_df <- data.frame(model = c("Subyearlings x Seabirds",
                                  "Yearlings x Seabirds",
                                  "Subyearlings x Hake",
                                  "Yearlings x Hake"),
                        label = c("(A) Subyearlings x Seabirds",
                                  "(B) Yearlings x Seabirds",
                                  "(C) Subyearlings x Hake",
                                  "(D) Yearlings x Hake"))

SRF_transported_SAR_model_predictions_comp %>% 
  left_join(labels_df, by = "model") -> SRF_transported_SAR_model_predictions_comp


SRF_transported_SAR_model_predictions_comp_plot <- ggplot(SRF_transported_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  # geom_errorbar(data = SRF_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "#2166ac") +
  geom_errorbar(data = SRF_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "blue") +
  geom_point(aes(x = run_year, shape = "Empirical", color = "Empirical", size = N), alpha = 0.5) +
  geom_point(data = SRF_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted", color = "Predicted"), size = 5, stroke = 0.7) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 13)) +
  # scale_color_manual(name = NULL, values = c("Empirical" = "#2166ac", "Predicted" = "#b2182b")) +
  scale_color_manual(name = NULL, values = c("Empirical" = "black", "Predicted" = "blue")) +
  scale_size_continuous(breaks = c(10, 1000, 10000)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  guides(size = guide_legend(order = 2)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.position = "bottom",
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, hjust = 0),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))+
  facet_wrap(~label, ncol = 2)

ggsave(here::here("figures", "paper_figures", "one_model","SRF_transported_SAR_model_predictions_comp_plot.png"), 
       SRF_transported_SAR_model_predictions_comp_plot,  
       height = 6, width = 12)

### SRF non transported ###
SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp$model <- "Yearlings x Seabirds"
SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model <- "Yearlings x Hake"
SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Seabirds"
SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp$model  <- "Subyearlings x Hake"


SRF_non_transported_06_2_6_seabird_prey_csyif_only_one_model_SAR_model_predictions_comp %>% 
  bind_rows(SRF_non_transported_06_2_6_seabird_prey_cssif_only_one_model_SAR_model_predictions_comp) -> SRF_non_transported_SAR_seabird_model_predictions_comp

SRF_non_transported_06_3_4_hake_prey_csyif_only_no_salmon_index_one_model_SAR_model_predictions_comp %>%
  bind_rows(SRF_non_transported_06_3_4_hake_prey_cssif_only_no_salmon_index_one_model_SAR_model_predictions_comp) -> SRF_non_transported_SAR_hake_model_predictions_comp

SRF_non_transported_SAR_seabird_model_predictions_comp %>% 
  bind_rows(SRF_non_transported_SAR_hake_model_predictions_comp) -> SRF_non_transported_SAR_model_predictions_comp

SRF_non_transported_SAR_model_predictions_comp$model <- factor(SRF_non_transported_SAR_model_predictions_comp$model, 
                                                               levels = c("Subyearlings x Seabirds",
                                                                          "Yearlings x Seabirds",
                                                                          "Subyearlings x Hake",
                                                                          "Yearlings x Hake"))
labels_df <- data.frame(model = c("Subyearlings x Seabirds",
                                  "Yearlings x Seabirds",
                                  "Subyearlings x Hake",
                                  "Yearlings x Hake"),
                        label = c("(A) Subyearlings x Seabirds",
                                  "(B) Yearlings x Seabirds",
                                  "(C) Subyearlings x Hake",
                                  "(D) Yearlings x Hake"))

SRF_non_transported_SAR_model_predictions_comp %>% 
  left_join(labels_df, by = "model") -> SRF_non_transported_SAR_model_predictions_comp


SRF_non_transported_SAR_model_predictions_comp_plot <- ggplot(SRF_non_transported_SAR_model_predictions_comp, aes(x = run_year, y = SAR)) +
  # geom_errorbar(data = SRF_non_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "#2166ac") +
  geom_errorbar(data = SRF_non_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, ymax = predicted_prob_upper, ymin = predicted_prob_lower), width = 0.1, color = "blue") +
  geom_point(aes(x = run_year, shape = "Empirical", color = "Empirical", size = N), alpha = 0.5) +
  geom_point(data = SRF_non_transported_SAR_model_predictions_comp, aes(x = run_year, y = predicted_prob, shape = "Predicted", color = "Predicted"), size = 5, stroke = 0.7) +
  scale_shape_manual(name = NULL, values = c("Empirical" = 19, "Predicted" = 13)) +
  # scale_color_manual(name = NULL, values = c("Empirical" = "#2166ac", "Predicted" = "#b2182b")) +
  scale_color_manual(name = NULL, values = c("Empirical" = "black", "Predicted" = "blue")) +
  scale_size_continuous(breaks = c(10, 1000, 10000)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  guides(size = guide_legend(order = 2)) +
  xlab("Run Year") +
  ylab("Marine Survival") +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.position = "bottom",
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, hjust = 0),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))+
  facet_wrap(~label, ncol = 2)

ggsave(here::here("figures", "paper_figures", "one_model","SRF_non_transported_SAR_model_predictions_comp_plot.png"), 
       SRF_non_transported_SAR_model_predictions_comp_plot,  
       height = 6, width = 12)
