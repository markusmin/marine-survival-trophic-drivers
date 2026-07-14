# plot_distribution_jsoes_seabirds - function to plot the density of a taxon from the JSOES seabirds data

plot_distribution_jsoes_seabirds <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(n_per_km2 == 0, "zero", "non-zero")) -> data
  
  survey_area_basemap +
    geom_point(data = data, aes(x = dec_long, y = dec_lat, size = n_per_km2, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    guides(size = guide_legend(title = "Density (N/km^2)"),
           color = guide_legend(title = "Encounter")) +
    facet_wrap(~year, nrow = 3) +
    ggtitle(taxon_name) +
    theme(legend.position = "right",
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}
