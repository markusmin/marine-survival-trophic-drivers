# generate figures for density and SE density for one SDM

generate_prediction_maps <- function(report, obj, opt, model_name){
  
  source(here::here("R", "JSOES_make_mesh.R"))
  
  # survey_predict_grid is a grid for making predictions that is created in the sourced scripts
  # trim survey_predict_grid to only the years of interest
  SDM_predict_grid <- subset(survey_predict_grid, year %in% 1999:2021)
  SDM_predict_grid$ln_d_gt <- as.vector(report$ln_d_gt)
  
  # get the standard errors
  SE_ln_d_gt = sample_var( obj=obj, var_name="ln_d_gt", mu=obj$env$last.par.best, prec=opt$SD$jointPrecision )
  SDM_predict_grid$SE_ln_d_gt <- as.vector(SE_ln_d_gt)
  
  ### Visualize predicted density
  
  
  density_predict_map <- survey_area_basemap_km +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = exp(ln_d_gt)),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow", limits = c(0, quantile(exp(SDM_predict_grid$ln_d_gt), 0.999)),
                          name = "Predicted\ndensity\n(N per km)") +
    # scale_fill_viridis_c(
    #   # trim extreme low values to make spatial variation more visible
    #   na.value = "gray", limits = c(0.05, quantile(exp(SDM_predict_grid$ln_d_gt), 0.999)),
    #   name = "Predicted\ndensity\n(N per km)") +
    facet_wrap(~year, nrow = 3) +
    theme(axis.text = element_blank(),
          legend.position = c(0.93, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'))
  
  ggsave(here::here("two_stage_models", "SDM_stage1", "jsoes_trawl", "csyif", "figures", paste0(model_name, "_density_predict_map.png")), density_predict_map,  
         height = 8, width = 10)
  
  ### Visualize SE in predicted density
  
  SE_density_predict_map <- survey_area_basemap_km +
    geom_tile(data = SDM_predict_grid, aes(x = X, y = Y, fill = SE_ln_d_gt),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow",  limits = c(0, quantile(SDM_predict_grid$SE_ln_d_gt, 0.995)),
                          name = "SE [log-\ndensity]\n(N per km)") +
    facet_wrap(~year, nrow = 3) +
    theme(axis.text = element_blank(),
          legend.position = c(0.93, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm'))
  
  ggsave(here::here("two_stage_models", "SDM_stage1", "jsoes_trawl", "csyif", "figures",  paste0(model_name, "_density_SE_predict_map.png")), SE_density_predict_map,  
         height = 8, width = 10)
  
  return(list(density_predict_map, SE_density_predict_map))
}