# plot_distribution_hake - function to plot the density of hake in the hake survey data

plot_distribution_hake <- function(data){
  survey_area_basemap_km_PRS_PWCC +
    geom_point(data = data, aes(x = X, y = Y, color = NASC, size = NASC)) +
    facet_wrap(~year, nrow = 2) +
    ggtitle("Pacific Hake") +
    scale_color_viridis_c(limits = c(0,100000)) +
    scale_size_area(limits = c(0,100000)) +
    guides(colour = guide_legend(title = "NASC (m^2 nmi^-2)"),
           size = guide_legend(title = "NASC (m^2 nmi^-2)")) +
    theme(legend.key.height = unit(1, "cm"),
          legend.key.width = unit(1, "cm"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.position = "right",
          legend.spacing.y = unit(0.05, 'cm'),
          # axis.text = element_text(size = 12),
          axis.title = element_blank(),
          axis.text = element_blank(),
          # axis.title = element_text(size = 16),
          strip.text = element_text(size = 20),
          plot.title = element_text(size = 24)) -> hake_distribution_plot
  
  return(hake_distribution_plot)
}
