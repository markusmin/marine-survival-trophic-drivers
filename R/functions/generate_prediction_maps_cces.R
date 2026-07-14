# generate figures for density and SE density for one SDM

generate_prediction_maps_cces <- function(report, obj, opt, model_name){

  if (!(exists("cces_inla_mesh_cutoff15_jsoes_domain"))){
    source(here::here("R", "CCES_make_mesh.R"))
  }
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid_jsoes, year %in% 2001:2021)
  SDM_predict_grid$ln_d_gt <- as.vector(report$ln_d_gt)
  
  # get the standard errors
  SE_ln_d_gt = sample_var( obj=obj, var_name="ln_d_gt", mu=obj$env$last.par.best, prec=opt$SD$jointPrecision )
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  ### create a basemap that's tailored to cces
  usa_spdf <- st_read(here::here("Data", "map_files", "USA_adm0.shp"))
  # load BC
  CAN_spdf <- st_read(here::here("Data", "map_files", "canada", "lpr_000b16a_e.shp"))
  BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
  BC_proj <- st_transform(BC_spdf, crs = 4326)


  # crop them to our desired area
  US_west_coast <- sf::st_crop(usa_spdf,
                               c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))

  BC_coast <- sf::st_crop(BC_proj,
                          c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))



  # convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
  UTM_zone_10_crs <- 32610

  US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
  BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

  # make this projection into kilometers
  US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
  BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)


  #### create base map for visualizing data
  cces_survey_area_basemap_km <- ggplot(US_west_coast_proj_km) +
    geom_sf() +
    geom_sf(data = BC_coast_proj_km) +
    ylab("Latitude")+
    xlab("Longitude")+
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill="white", color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.14, 0.2),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  
  
  ### Visualize predicted density
  
  density_predict_map <- cces_survey_area_basemap_km +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = ln_d_gt),
              width = 7, height = 7) +
    # scale_fill_viridis_c(
    #                       # trim extreme high values to make spatial variation more visible
    #                       na.value = "yellow", limits = c(0, quantile(SDM_predict_grid$ln_d_gt, 0.999)),
    #                       name = "Predicted\nlog-density\n(total)") +
    scale_fill_viridis_c( name = "Predicted\nlog-density\n(total)") +
    facet_wrap(~year, nrow = 2) +
    scale_y_continuous(labels = c(expression(~"44°N"), 
                                  expression(~"45°N"),
                                  expression(~"46°N"),
                                  expression(~"47°N"),
                                  expression(~"48°N"))) +
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
          plot.title = element_text(size = 24))
  
  ggsave(here::here("two_stage_models", "SDM_stage1", "cces", "figures", paste0(model_name, "_density_predict_map.png")), density_predict_map,  
         height = 10, width = 10)
  
  ### Visualize SE in predicted density
  
  SE_density_predict_map <- cces_survey_area_basemap_km +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = "SE [log-\ndensity]\n(total)") +
    facet_wrap(~year, nrow = 2) +
    scale_y_continuous(labels = c(expression(~"44°N"), 
                                  expression(~"45°N"),
                                  expression(~"46°N"),
                                  expression(~"47°N"),
                                  expression(~"48°N"))) +
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
          plot.title = element_text(size = 24))
  
  ggsave(here::here("two_stage_models", "SDM_stage1", "cces", "figures",  paste0(model_name, "_density_SE_predict_map.png")), SE_density_predict_map,  
         height = 10, width = 10)
  
  return(list(density_predict_map, SE_density_predict_map))
}
