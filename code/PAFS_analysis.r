library(spatstat)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(mgcv)
library(stats)
library(tidyverse)
library(gratia)
source('code/semi_spp_functions.R')

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install the high-res data package
remotes::install_github("ropensci/rnaturalearthhires") # download package for ne_states()


#### load data
# df_IPP is the information of PFAS
df_IPP = read.csv("data/PFAS/IPP Ready Model Data.csv")
# df_ext is the information of PFAS
df_ext = read.csv("data/PFAS/IPP Ready Extrap Data.csv")

###### draw locations of PFAS


df_IPP_thresholded = df_IPP[df_IPP$Concentration>=8,]
sf_PFAS = st_as_sf(df_IPP, coords = c("Longitude","Latitude"),crs = 4326)


us_mainland <- usa_states[!usa_states$name %in% c("Alaska", "Hawaii"),]
us_mainland_ll <- st_transform(us_mainland, crs = 4326)

ggplot() +
  geom_sf(data = us_mainland_ll, fill = 'ivory',color = 'grey80') + # color=NA removes
  geom_sf(data = PA_states_geometry, fill = NA, color = "green",linewidth = 0.6)+
  geom_sf(data = sf_PFAS[sf_PFAS$Concentration<8,], color = 'grey',size = 0.05)+
  geom_sf(data = sf_PFAS[sf_PFAS$Concentration>=8,], color = 'darkred',size = 0.05)+
  theme_minimal()+
  theme_void() 


##### get data in regions with enough observations

# PA means the states have enough PFAS measurement

PA_states <- c("California", "Colorado", "Maryland", "Massachusetts", 
               "Michigan", "New Hampshire", "New Jersey", "North Carolina", 
               "Ohio", "Rhode Island", "South Carolina", "West Virginia", "Wisconsin")
usa_states <- ne_states(country = "united states of america", returnclass = "sf")
PA_states_geometry <- usa_states[usa_states$name %in% PA_states, ]
df_PA <- df_IPP %>%
  filter(STUSPS %in% state.abb[match(PA_states, state.name)])


# Find locations where PFAS concentration >= 8 ppt
PFAS_threshold = 8
df_PA_ab = df_PA[df_PA$Concentration >= 8,]
df_All_ab = df_IPP[df_IPP$Concentration >= 8,]

sf_PA_ab = st_as_sf(df_PA_ab,coords = c("Longitude","Latitude"),crs = 4326)
sf_All_ab = st_as_sf(df_All_ab,coords = c("Longitude","Latitude"),crs = 4326)

# transfer to ppp class in spatstat

# find the observational window
target_epsg_code = 26910 # This is NAD83 for California
PA_states_geometry_proj = st_transform(PA_states_geometry, crs = target_epsg_code)
window_PA = as.owin(PA_states_geometry_proj$geometry)
pp_PFAS_PA =as.ppp(st_transform(sf_PA_ab, crs = target_epsg_code)$geometry, W = window_PA)





###### get covariates on observation

load("data/PFAS/cov_extrapolation_im.RData")



# build covariates
## combine covariates in the ext_CA and IPP_CA
## create a pixel image in spatstat
selected_covariates = colnames(df_ext)
cov_df = rbind(df_ext[,selected_covariates],df_IPP[,selected_covariates])
cov_df = cov_df[cov_df$STUSPS %in% state.abb[match(PA_states, state.name)],]
cov_df$STUSPS = NULL



cov_sf = st_as_sf(cov_df, coords = c("Longitude","Latitude"), crs = 4326)
cov_sf = st_transform(cov_sf, crs = target_epsg_code)
# further filter the points outside the boundary observational window
cov_sf <- st_filter(cov_sf, PA_states_geometry_proj, .predicate = st_intersects)

covariate_names <- c("DistanceToAirport", 
                     "DistanceToMilitaryBase", 
                     "DistanceToLandfill",
                     "DistanceToApparel",
                     "DistanceToFurniture",
                     "PopDen",
                     "Median.Earnings.Inflation.Adj.2020",
                     "Bachelor.Degree.Percent.2020")

# create covariate image to fit IPP using package


covariate_images <- list()
for (cov_name in covariate_names){
  cat("Processing covariate:", cov_name, "\n")
  ppp.base = as.ppp(cov_sf$geometry, W = window_PA)
  marks(ppp.base) = cov_sf[[cov_name]]
  # eps is width and height of pixels, unit is meter
  im_cov = idw(ppp.base, at="pixels", eps = 1000) # this is very slow if eps is 500 and 1000, okay for 2000
  # im_cov = nnmark(ppp.base, at="pixels", eps = 1000)
  covariate_images[[cov_name]] <- im_cov
}

save(covariate_images, file = "data/PFAS/cov_extrapolation_im.RData")

###### Fit IPP

covariate_names <- c("DistanceToAirport", 
                     "DistanceToMilitaryBase", 
                     "DistanceToLandfill",
                     "DistanceToApparel",
                     "DistanceToFurniture",
                     "PopDen",
                     "Median.Earnings.Inflation.Adj.2020",
                     "Bachelor.Degree.Percent.2020")

#### Fit IPP

parametric_formula = as.formula(paste("pp_PFAS_PA~", paste(covariate_names, collapse = "+")))

semiparametric_formula = as.formula(paste0("pp_PFAS_PA~",paste(covariate_names[1:5], collapse = "+"),'+', paste('s(',covariate_names[-(1:5)],')', collapse = "+")))


ppmfit_para = ppm(parametric_formula,data = covariate_images, use.gam = TRUE, method="mpl")

ppmfit_semi = ppm(semiparametric_formula,data = covariate_images, use.gam = TRUE, method="mpl")



##### covariance estimation for semiparametric model assuming Poisson point process

target_covariate_names = covariate_names[1:5]
nuisance_covariate_names = covariate_names[6:8]

ppmfit =  ppmfit_semi
gamfit = ppmfit$internal$glmfit

cov.df = gamfit$model
lambda = GAMpredict(gamfit = gamfit,ppmfit = ppmfit)
wt = ppmfit$Q$w

lfd_projection = list()
for (target_cov in target_covariate_names){
  formula <- as.formula(paste0(target_cov,'~',paste(nuisance_covariate_names,collapse = '+')))
  weighted_loess <-loess(formula, weights = lambda*wt, data = cov.df)
  lfd = predict(weighted_loess,newdata = cov.df)
  lfd_projection[[target_cov]] = lfd
}

lfd_projection = as.data.frame(lfd_projection)

lfd_projection = cov.df[,target_covariate_names] - lfd_projection 

lfd_projection = as.matrix(lfd_projection)

sensitivity = t(lfd_projection) %*% (lfd_projection * lambda * wt)

# calculate S.E. of the target estimation 
sqrt(diag(solve(sensitivity)))

###### plot nuisance effect nonlinear term

draw(ppmfit_semi$internal$glmfit, select = "s(Median.Earnings.Inflation.Adj.2020)", rug = FALSE, ci_alpha = 0) + 
  geom_line(color = "#0479A8", linewidth = 1.0) +
  coord_cartesian(xlim = c(15000, 66000)) +  # Set the visible range for the x-axis
  labs(
    title = "",
    x = "Median Earnings Per Year ($)", # Set custom x-axis label
    y = "Estimated Nuisance Function"         # Set custom y-axis label
  )+
  # theme_minimal() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank()
  )



#### national intensity


##### get fitted intensity

pred_intensity = exp(predict(ppmfit_semi$internal$glmfit, df_ext))
df_pred = data.frame(Latitude = df_ext$Latitude, Longitude = df_ext$Longitude, intensity = pred_intensity)

pred_intensity_parametric = exp(predict(ppmfit_para$internal$glmfit, df_ext))
df_pred_para = data.frame(Latitude = df_ext$Latitude, Longitude = df_ext$Longitude, intensity = pred_intensity_parametric)



write.csv(df_pred,'pred_intensity.csv',row.names = FALSE)

df_pred <- read.csv('output/pred_intensity.csv')

sf_pred = st_as_sf(df_pred, coords = c("Longitude","Latitude"), crs = 4326)
sf_pred_para = st_as_sf(df_pred_para, coords = c("Longitude","Latitude"), crs = 4326)


##### interpolate fitted intensity

us_mainland <- usa_states[!usa_states$name %in% c("Alaska", "Hawaii"),]
us_mainland <- st_transform(us_mainland, crs = target_epsg_code)

sf_pred_projected <- st_transform(sf_pred, crs = target_epsg_code)
sf_pred_para_projected <- st_transform(sf_pred_para, crs = target_epsg_code)

grid <- st_make_grid(us_mainland, cellsize = 5000, square = TRUE) # 2 km resolution
grid_sf <- st_as_sf(grid)
grid_sf <- grid_sf[us_mainland,]

grid_sf_centroid <- st_centroid(grid_sf) # use center points

idw_result <- gstat::idw(intensity ~ 1, locations = sf_pred_projected, newdata = grid_sf_centroid, idp = 2)
idw_result_para <- gstat::idw(intensity ~ 1, locations = sf_pred_para_projected, newdata = grid_sf_centroid, idp = 2)

grid_sf$pred_intensity <- idw_result$var1.pred
grid_sf$pred_intensity_para <- idw_result_para$var1.pred

##### draw quantile map for fitted/interpolated intensity

quartile_breaks <- quantile(c(grid_sf$pred_intensity,grid_sf$pred_intensity) , probs = c(0, 0.40,0.60, 0.70,0.80,0.90,0.95,0.97,0.99, 1))

quantile_intervals <- cut(grid_sf$pred_intensity,
                          breaks = quartile_breaks,
                          include.lowest = TRUE) 
grid_sf$pred_intensity_qt <- quantile_intervals

quantile_intervals_para <- cut(grid_sf$pred_intensity_para,
                               breaks = quartile_breaks,
                               include.lowest = TRUE) 

grid_sf$pred_intensity_qt_para <- quantile_intervals_para


num_intervals <- length(quartile_breaks) - 1
color_ramp <- colorRampPalette(c("beige","orangered", "darkred","black"))
ramp_colors <- color_ramp(num_intervals - 1)
my_custom_palette <- c("white", ramp_colors)

grid_sf_ll <- st_transform(grid_sf, crs = 4326)
us_mainland_ll <- st_transform(us_mainland, crs = 4326)


ggplot() +
  geom_sf(data = grid_sf_ll,
          aes(fill = pred_intensity_qt), # Map the fill color to your new quantile factor
          color = NA,              # Set a static border color for all polygons
          lwd = 0.1                      # Line width for borders
  ) +
  geom_sf(data = us_mainland_ll, color = "grey90", fill = NA )+
  scale_fill_manual(
    values = my_custom_palette,         
    name = 'Quantiles',
    labels = c('<40%','<60%', '<70%','<80%','<90%',"<95%","<97%",'<99%', '<100%'), 
    na.value = "transparent"
  ) +
  theme(
    legend.position = "right"
  ) +
  labs(title = "Quantile Map of Fitted Intensity of Semi-Parametric Model") +
  theme_void() 

ggplot() +
  geom_sf(data = grid_sf_ll,
          aes(fill = pred_intensity_qt_para), # Map the fill color to your new quantile factor
          color = NA,              # Set a static border color for all polygons
          lwd = 0.1                      # Line width for borders
  ) +
  geom_sf(data = us_mainland_ll, color = 'grey', fill = NA) +
  # geom_sf(data = PA_states_geometry, fill = NA, color = "green",linewidth = 0.6)+
  scale_fill_manual(
    values = my_custom_palette,             # Apply your custom color vector
    # name = "Predicted Intensity\n(Quantiles)" # Set the legend title
    name = NULL, 
    labels = c('<40%','<60%', '<70%','<80%','<90%',"<95%","<97%",'<99%', '<100%'),
    na.value = "transparent"
  ) +
  # theme_minimal() + # Use a clean theme
  theme(
    # plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Quantile Map of Fitted Intensity of Parametric Model") +
  theme_void()


##### Zoom in Alabama

qt_labels <- c('<40%','<60%', '<70%','<80%','<90%',"<95%","<97%",'<99%', '<100%')

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
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Parametric Model — California") +
  theme_void() +
  theme(legend.position = "right")

ggplot() +
  geom_sf(data = grid_California_sf_ll, aes(fill = pred_intensity_qt), color = NA, lwd = 0.1) +
  geom_sf(data = us_California_ll, color = 'grey', fill = NA) +
  scale_fill_manual(values = my_custom_palette, name = NULL, labels = qt_labels, na.value = "transparent") +
  labs(title = "Quantile Intensity Map using Semi-Parametric Model — California") +
  theme_void() +
  theme(legend.position = "right")
