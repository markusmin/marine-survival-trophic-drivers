# marine-survival-trophic-drivers
Code and data for the manuscript "Trophic drivers of salmon marine survival revealed by spatiotemporal data integration"

The directories in this repository are set us as follows
```
-- Data ------- all input data files for project
-- R/ ------- all scripts necessary to reproduce analysis
-- docs ------- includes supplemental materials HTML output, so that it can be easily viewed on GitHub
-- figures ------- various figures
-- model_inputs/ ------- data processing output from scripts in the R/ folder that are inputs to later scripts.
```

Please note that some files are not uploaded to GitHub because of their large size. These include:
- the shapefiles used to generate the maps
- the raw PIT tag data (this can be queried from PTAGIS by following the steps in the `01_PIT_tag_data.Rmd` script)