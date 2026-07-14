# plot_distribution - function to plot the density of a taxon

plot_distribution_PRS_PWCC <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(total == 0, "zero", "non-zero")) -> data
  
  # Drop 1998 (no water jellies recorded in that year)
  data %>% 
    filter(year != 1998) -> data
  # 
  # # crop them to our desired area
  # US_west_coast <- sf::st_crop(usa_spdf,
  #                              c(xmin = -126, ymin = 40.42, xmax = -123, ymax = 48.5))
  # 
  # BC_coast <- sf::st_crop(BC_proj,
  #                         c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))
  # 
  # 
  # 
  # # convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
  # UTM_zone_10_crs <- 32610
  # 
  # US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
  # BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)
  # 
  # # make this projection into kilometers
  # US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
  # BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
  # 
  # 
  # #### create base map for visualizing data
  # 
  # survey_area_basemap_km <- ggplot(US_west_coast_proj_km) +
  #   geom_sf() +
  #   geom_sf(data = BC_coast_proj_km) + 
  #   ylab("Latitude")+
  #   xlab("Longitude")+
  #   theme(plot.background = element_rect(fill = "white"),
  #         panel.background = element_rect(fill="white", color = "black"),
  #         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         legend.position = c(0.14, 0.2),
  #         legend.title = element_text(size = 14),
  #         legend.text = element_text(size = 12),
  #         axis.ticks = element_blank(),
  #         axis.text = element_blank(),
  #         axis.title = element_blank())
  

  
  survey_area_basemap_km_PRS_PWCC +
    geom_point(data = data, aes(x = X, y = Y, size = total, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    guides(size = guide_legend(title = "Density (N/km trawled)"),
           color = guide_legend(title = "Encounter")) +
    facet_wrap(~year, nrow = 2) +
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
