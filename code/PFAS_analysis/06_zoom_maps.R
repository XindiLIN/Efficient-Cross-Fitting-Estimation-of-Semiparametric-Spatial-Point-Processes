library(ggplot2)
library(sf)
library(gstat)

load("output/PFAS/prepared_data.RData")
load("output/PFAS/national_grid.RData")

##### Zoom in Alabama

us_Alabama <- usa_states[usa_states$name == "Alabama", ]
us_Alabama <- st_transform(us_Alabama, crs = target_epsg_code)

grid_Alabama      <- st_make_grid(us_Alabama, cellsize = 3000, square = TRUE)
grid_Alabama_sf   <- st_as_sf(grid_Alabama)
grid_Alabama_sf   <- grid_Alabama_sf[us_Alabama, ]
grid_Alabama_sf_centroid <- st_centroid(grid_Alabama_sf)

sf_pred_projected_Alabama      <- st_filter(sf_pred_projected,      grid_Alabama_sf)
sf_pred_para_projected_Alabama <- st_filter(sf_pred_para_projected, grid_Alabama_sf)

idw_result_Alabama      <- gstat::idw(intensity ~ 1, locations = sf_pred_projected_Alabama,      newdata = grid_Alabama_sf_centroid, idp = 2)
idw_result_para_Alabama <- gstat::idw(intensity ~ 1, locations = sf_pred_para_projected_Alabama, newdata = grid_Alabama_sf_centroid, idp = 2)

grid_Alabama_sf$pred_intensity      <- idw_result_Alabama$var1.pred
grid_Alabama_sf$pred_intensity_para <- idw_result_para_Alabama$var1.pred

grid_Alabama_sf$pred_intensity_qt      <- cut(grid_Alabama_sf$pred_intensity,      breaks = quartile_breaks, include.lowest = TRUE)
grid_Alabama_sf$pred_intensity_qt_para <- cut(grid_Alabama_sf$pred_intensity_para, breaks = quartile_breaks, include.lowest = TRUE)

grid_Alabama_sf_ll <- st_transform(grid_Alabama_sf, crs = 4326)
us_Alabama_ll      <- st_transform(us_Alabama,      crs = 4326)

ggplot() +
  geom_sf(data = grid_Alabama_sf_ll, aes(fill = pred_intensity_qt_para), color = NA, lwd = 0.1) +
  geom_sf(data = us_Alabama_ll, color = 'grey', fill = NA) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Parametric Model — Alabama") +
  theme_void() +
  theme(legend.position = "right")

ggplot() +
  geom_sf(data = grid_Alabama_sf_ll, aes(fill = pred_intensity_qt), color = NA, lwd = 0.1) +
  geom_sf(data = us_Alabama_ll, color = 'grey', fill = NA) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Semi-Parametric Model — Alabama") +
  theme_void() +
  theme(legend.position = "right")


##### Zoom in California

us_California <- usa_states[usa_states$name == "California", ]
us_California <- st_transform(us_California, crs = target_epsg_code)

grid_California      <- st_make_grid(us_California, cellsize = 3000, square = TRUE)
grid_California_sf   <- st_as_sf(grid_California)
grid_California_sf   <- grid_California_sf[us_California, ]
grid_California_sf_centroid <- st_centroid(grid_California_sf)

sf_pred_projected_California      <- st_filter(sf_pred_projected,      grid_California_sf)
sf_pred_para_projected_California <- st_filter(sf_pred_para_projected, grid_California_sf)

idw_result_California      <- gstat::idw(intensity ~ 1, locations = sf_pred_projected_California,      newdata = grid_California_sf_centroid, idp = 2)
idw_result_para_California <- gstat::idw(intensity ~ 1, locations = sf_pred_para_projected_California, newdata = grid_California_sf_centroid, idp = 2)

grid_California_sf$pred_intensity      <- idw_result_California$var1.pred
grid_California_sf$pred_intensity_para <- idw_result_para_California$var1.pred

grid_California_sf$pred_intensity_qt      <- cut(grid_California_sf$pred_intensity,      breaks = quartile_breaks, include.lowest = TRUE)
grid_California_sf$pred_intensity_qt_para <- cut(grid_California_sf$pred_intensity_para, breaks = quartile_breaks, include.lowest = TRUE)

grid_California_sf_ll <- st_transform(grid_California_sf, crs = 4326)
us_California_ll      <- st_transform(us_California,      crs = 4326)

sf_All_ab_California <- st_filter(sf_All_ab, usa_states[usa_states$name == "California", ])

ggplot() +
  geom_sf(data = grid_California_sf_ll, aes(fill = pred_intensity_qt_para), color = NA, lwd = 0.1) +
  geom_sf(data = us_California_ll, color = 'grey', fill = NA) +
  geom_sf(data = sf_All_ab_California, color = 'lightgreen', size = 0.5) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Parametric Model — California") +
  theme_void() +
  theme(legend.position = "right")

ggplot() +
  geom_sf(data = grid_California_sf_ll, aes(fill = pred_intensity_qt), color = NA, lwd = 0.1) +
  geom_sf(data = us_California_ll, color = 'grey', fill = NA) +
  geom_sf(data = sf_All_ab_California, color = 'lightgreen', size = 0.5) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Semi-Parametric Model — California") +
  theme_void() +
  theme(legend.position = "right")
