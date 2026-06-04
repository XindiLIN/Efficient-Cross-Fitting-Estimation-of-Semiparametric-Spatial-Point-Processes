library(ggplot2)
library(sf)
library(gstat)

source("code/PFAS_analysis/03_fit_ipp.R")

##### get fitted intensity
pred_intensity      = exp(predict(gamfit_semi, df_ext))
pred_intensity_para = exp(predict(gamfit_para, df_ext))

df_pred      = data.frame(Latitude = df_ext$Latitude, Longitude = df_ext$Longitude, intensity = pred_intensity)
df_pred_para = data.frame(Latitude = df_ext$Latitude, Longitude = df_ext$Longitude, intensity = pred_intensity_para)

sf_pred      = st_as_sf(df_pred,      coords = c("Longitude", "Latitude"), crs = 4326)
sf_pred_para = st_as_sf(df_pred_para, coords = c("Longitude", "Latitude"), crs = 4326)

##### build continental US grid and interpolate fitted intensity
us_mainland <- usa_states[!usa_states$name %in% c("Alaska", "Hawaii"), ]
us_mainland <- st_transform(us_mainland, crs = target_epsg_code)

sf_pred_projected      <- st_transform(sf_pred,      crs = target_epsg_code)
sf_pred_para_projected <- st_transform(sf_pred_para, crs = target_epsg_code)

grid    <- st_make_grid(us_mainland, cellsize = 5000, square = TRUE)
grid_sf <- st_as_sf(grid)
grid_sf <- grid_sf[us_mainland, ]
grid_sf_centroid <- st_centroid(grid_sf)

idw_result      <- gstat::idw(intensity ~ 1, locations = sf_pred_projected,      newdata = grid_sf_centroid, idp = 2)
idw_result_para <- gstat::idw(intensity ~ 1, locations = sf_pred_para_projected, newdata = grid_sf_centroid, idp = 2)

grid_sf$pred_intensity      <- idw_result$var1.pred
grid_sf$pred_intensity_para <- idw_result_para$var1.pred

##### quantile breaks shared across both models for a fair visual comparison
quartile_breaks <- quantile(
  c(grid_sf$pred_intensity, grid_sf$pred_intensity_para),
  probs = c(0, 0.40, 0.60, 0.70, 0.80, 0.90, 0.95, 0.97, 0.99, 1)
)

grid_sf$pred_intensity_qt      <- cut(grid_sf$pred_intensity,      breaks = quartile_breaks, include.lowest = TRUE)
grid_sf$pred_intensity_qt_para <- cut(grid_sf$pred_intensity_para, breaks = quartile_breaks, include.lowest = TRUE)

num_intervals      <- length(quartile_breaks) - 1
color_ramp         <- colorRampPalette(c("beige", "orangered", "darkred", "black"))
ramp_colors        <- color_ramp(num_intervals - 1)
my_custom_palette  <- c("white", ramp_colors)
qt_labels          <- c('<40%', '<60%', '<70%', '<80%', '<90%', '<95%', '<97%', '<99%', '<100%')

grid_sf_ll      <- st_transform(grid_sf, crs = 4326)
us_mainland_ll  <- st_transform(us_mainland, crs = 4326)

##### semiparametric model map
p_semi <- ggplot() +
  geom_sf(data = grid_sf_ll, aes(fill = pred_intensity_qt), color = NA, lwd = 0.1) +
  geom_sf(data = us_mainland_ll, color = "grey90", fill = NA) +
  scale_fill_manual(values = my_custom_palette, name = 'Quantiles', labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Map of Fitted Intensity of Semi-Parametric Model") +
  theme_void() +
  theme(legend.position = "right")

print(p_semi)
ggsave("output/PFAS/intensity_map_semiparametric.png", plot = p_semi,
       width = 10, height = 6, dpi = 300)

##### parametric model map
p_para <- ggplot() +
  geom_sf(data = grid_sf_ll, aes(fill = pred_intensity_qt_para), color = NA, lwd = 0.1) +
  geom_sf(data = us_mainland_ll, color = 'grey', fill = NA) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Map of Fitted Intensity of Parametric Model") +
  theme_void() +
  theme(legend.position = "right")

print(p_para)
ggsave("output/PFAS/intensity_map_parametric.png", plot = p_para,
       width = 10, height = 6, dpi = 300)

save(grid_sf, sf_pred_projected, sf_pred_para_projected,
     quartile_breaks, my_custom_palette, qt_labels,
     file = "output/PFAS/national_grid.RData")
