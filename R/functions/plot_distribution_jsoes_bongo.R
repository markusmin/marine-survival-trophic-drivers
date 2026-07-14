# plot_distribution_jsoes_bongo - function to plot the density of a taxon from the JSOES Bongo data

plot_distribution_jsoes_bongo <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(total_sum_of_density_number_m3 == 0, "zero", "non-zero")) -> data
  
  # Drop 1998 (no water jellies recorded in that year)
  data %>% 
    filter(year != 1998) -> data
  
  survey_area_basemap +
    geom_point(data = data, aes(x = dec_long, y = dec_lat, size = total_sum_of_density_number_m3, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    guides(size = guide_legend(title = "Density (N/m^3)"),
           color = guide_legend(title = "Encounter")) +
    facet_wrap(~year, nrow = 3) +
    ggtitle(taxon_name) +
    theme(legend.position = c(0.925, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}

# plot_distribution_jsoes_bongo - function to plot the density of a taxon from the JSOES Bongo data

plot_distribution_jsoes_bongo_biomass <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(total_sum_of_final_carbon_mg_m3 == 0, "zero", "non-zero")) -> data
  
  # Drop 1998 (no water jellies recorded in that year)
  data %>% 
    filter(year != 1998) -> data
  
  survey_area_basemap +
    geom_point(data = data, aes(x = dec_long, y = dec_lat, size = total_sum_of_final_carbon_mg_m3, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    guides(size = guide_legend(title = "Biomass (mg/m^3)"),
           color = guide_legend(title = "Encounter")) +
    facet_wrap(~year, nrow = 3) +
    ggtitle(taxon_name) +
    theme(legend.position = c(0.925, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}
