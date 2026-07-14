## stage2_helper_functions

# These scripts help visualize stage 2 model outputs

# function to get the linear predictor into an interpretable space
inv.logit <- function(linear_predictor){
  response <- 1/(1 + exp(-linear_predictor))
  return(response)
}

generate_posterior_predictive <- function(SD_summary,
                                          fixed_effects_SD_summary,
                                          parameter_name, 
                                          covariate_values){
  param_estimate <- subset(SD_summary, parameter == parameter_name)$estimate
  param_SE <- subset(SD_summary, parameter == parameter_name)$std_error
  
  # estimate posteriors using the param estimate, and the 95% CI
  param_estimate_low <- param_estimate - 1.96 * param_SE
  param_estimate_high <- param_estimate + 1.96 * param_SE
  
  # Marginalize over the parameter of interest - set every other parameter to the mean value
  other_fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, parameter != parameter_name)
  
  mean_predictor <- param_estimate * covariate_values +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  low_predictor <- param_estimate_low * covariate_values +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  high_predictor <- param_estimate_high * covariate_values +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  data.frame(covariate_value = covariate_values,
             lower = low_predictor,
             mean = mean_predictor,
             upper = high_predictor) %>% 
    pivot_longer(cols = c("lower", "mean", "upper"), names_to = "estimate", values_to = "effect") %>% 
    mutate(prob_survival = inv.logit(effect)) %>% 
    dplyr::select(-effect) %>% 
    pivot_wider(names_from = "estimate",
                values_from = "prob_survival") -> posterior_predictive
  
  
  
  return(posterior_predictive)
}

plot_posterior_predictive <- function(posterior_predictive, covariate_name){
  plot <- ggplot(posterior_predictive, aes(x = covariate_value, y = mean, ymin = lower, ymax = upper)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    xlab(covariate_name) +
    ylab("Probability of marine survival") +
    scale_y_continuous(expand = c(0,0), breaks = seq(0, 0.05, 0.01),
                       labels = c("0", "0.01", "0.02", "0.03", "0.04", "0.05")) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
    coord_cartesian(ylim = c(0, 0.055), clip="off") +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 15),
          panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
          legend.position = "none",
          plot.margin = unit(c(0, 0.5, 0.2, 0.2),"cm"))
  
  return(plot)
}


# Outmigration needs separate code, because it is composed of two parameters
generate_posterior_predictive_outmigration <- function(SD_summary,
                                                       fixed_effects_SD_summary,
                                                       covariate_values){
  outmigration_param_estimate <- subset(SD_summary, parameter == "beta_outmigration")$estimate
  outmigration_param_SE <- subset(SD_summary, parameter == "beta_outmigration")$std_error
  
  # estimate posteriors using the param estimate, and the 95% CI
  outmigration_param_estimate_low <- outmigration_param_estimate - 1.96 * outmigration_param_SE
  outmigration_param_estimate_high <- outmigration_param_estimate + 1.96 * outmigration_param_SE
  
  outmigration2_param_estimate <- subset(SD_summary, parameter == "beta_outmigration2")$estimate
  outmigration2_param_SE <- subset(SD_summary, parameter == "beta_outmigration2")$std_error
  
  # estimate posteriors using the param estimate, and the 95% CI
  outmigration2_param_estimate_low <- outmigration2_param_estimate - 1.96 * outmigration2_param_SE
  outmigration2_param_estimate_high <- outmigration2_param_estimate + 1.96 * outmigration2_param_SE
  
  # Marginalize over the parameter of interest - set every other parameter to the mean value
  other_fixed_effects_SD_summary <- subset(fixed_effects_SD_summary, !(parameter %in% c("beta_outmigration", "beta_outmigration2")))
  
  mean_predictor <- outmigration_param_estimate * covariate_values + outmigration2_param_estimate * covariate_values^2 +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  low_predictor <- outmigration_param_estimate_low * covariate_values + outmigration2_param_estimate_low * covariate_values^2 +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  high_predictor <- outmigration_param_estimate_high * covariate_values + outmigration2_param_estimate_high * covariate_values^2 +
    sum(other_fixed_effects_SD_summary$mean_effect)
  
  data.frame(covariate_value = covariate_values + mean(SRF_subset$outmigration_date),
             lower = low_predictor,
             mean = mean_predictor,
             upper = high_predictor) %>% 
    pivot_longer(cols = c("lower", "mean", "upper"), names_to = "estimate", values_to = "effect") %>% 
    mutate(prob_survival = inv.logit(effect)) %>% 
    dplyr::select(-effect) %>% 
    pivot_wider(names_from = "estimate",
                values_from = "prob_survival") -> posterior_predictive
  
  return(posterior_predictive)
}

# EiV visualization
generate_eiv_comp_data <- function(SDM_estimate, SDM_SE, SAR_estimate, SAR_SE, common_years){
  
  data.frame(year = common_years,
             SDM_estimate = SDM_estimate,
             SDM_SE = SDM_SE) %>% 
    mutate(SDM_lower = SDM_estimate - 1.96*SDM_SE,
           SDM_upper = SDM_estimate + 1.96*SDM_SE) %>% 
    dplyr::select(-SDM_SE) %>% 
    dplyr::rename(lower = SDM_lower, upper = SDM_upper, estimate = SDM_estimate) %>% 
    mutate(model = "SDM") -> SDM_estimate_df
  
  data.frame(year = common_years,
             SAR_estimate = SAR_estimate,
             SAR_SE = SAR_SE) %>% 
    mutate(SAR_lower = SAR_estimate - 1.96*SAR_SE,
           SAR_upper = SAR_estimate + 1.96*SAR_SE) %>% 
    dplyr::select(-SAR_SE) %>% 
    dplyr::rename(lower = SAR_lower, upper = SAR_upper, estimate = SAR_estimate) %>% 
    mutate(model = "SAR") -> SAR_estimate_df
  
  SDM_estimate_df %>% 
    bind_rows(SAR_estimate_df) %>% 
    mutate(model = factor(model, levels = c("SDM", "SAR"))) -> SDM_SAR_comp
  
  return(SDM_SAR_comp)
  
}

compare_EiV_plot <- function(eiv_comp_data, covariate_name){
  plot <- ggplot(eiv_comp_data, aes(x = year, y = estimate, ymin = lower, ymax = upper, color = model)) +
    geom_line() +
    geom_point() + 
    geom_ribbon(alpha = 0.2, lty = 2, fill = NA) +
    xlab("Year") +
    ylab(covariate_name)
  
  return(plot)
}