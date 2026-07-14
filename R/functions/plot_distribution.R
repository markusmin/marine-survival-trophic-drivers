# plot_distribution - function to plot the density of a taxon

plot_distribution <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(n_per_km == 0, "zero", "non-zero")) -> data
  
  # Drop 1998 (no water jellies recorded in that year)
  data %>% 
    filter(year != 1998) -> data
  
  survey_area_basemap +
    geom_point(data = data, aes(x = mid_long, y = mid_lat, size = n_per_km, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    guides(size = guide_legend(title = "Density\n(N/km trawled)"),
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
