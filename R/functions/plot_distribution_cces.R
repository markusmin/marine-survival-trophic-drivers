# plot_distribution_cces - function to plot the density of a taxon from the JSOES seabirds data

plot_distribution_cces <- function(data, taxon_name){
  
  survey_area_basemap +
    geom_point(data = data, aes(x = long, y = lat, size = density_bins, color = density_bins),
               alpha = 0.5) +
    scale_color_manual(values = c("gray", viridis(5))) +
    facet_wrap(~year, nrow = 2) +
    ggtitle(taxon_name) +
    scale_size_manual(values = c(0.4, 2, 4, 6, 8, 10)) +
    guides(colour = guide_legend(), size = guide_legend()) +
    theme(legend.key.height = unit(1, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.position = "right",
          legend.spacing.y = unit(0.05, 'cm'),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 24)) -> species_distribution_plot
  
  species_distribution_plot
  
  return(species_distribution_plot)
}
