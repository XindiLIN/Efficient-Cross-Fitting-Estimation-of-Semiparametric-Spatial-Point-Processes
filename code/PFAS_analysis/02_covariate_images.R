library(spatstat)
library(sf)
library(tidyverse)

load("output/PFAS/prepared_data.RData")

covariate_names <- c("DistanceToAirport",
                     "DistanceToMilitaryBase",
                     "DistanceToLandfill",
                     "DistanceToApparel",
                     "DistanceToFurniture",
                     "PopDen",
                     "Median.Earnings.Inflation.Adj.2020",
                     "Bachelor.Degree.Percent.2020")

selected_covariates = colnames(df_ext)
cov_df = rbind(df_ext[, selected_covariates], df_IPP[, selected_covariates])
cov_df = cov_df[cov_df$STUSPS %in% state.abb[match(PA_states, state.name)], ]
cov_df$STUSPS = NULL

cov_sf = st_as_sf(cov_df, coords = c("Longitude", "Latitude"), crs = 4326)
cov_sf = st_transform(cov_sf, crs = target_epsg_code)
cov_sf <- st_filter(cov_sf, PA_states_geometry_proj, .predicate = st_intersects)

# Build pixel images for each covariate (slow; run once and reload via RData)
covariate_images <- list()
for (cov_name in covariate_names) {
  cat("Processing covariate:", cov_name, "\n")
  ppp.base = as.ppp(cov_sf$geometry, W = window_PA)
  marks(ppp.base) = cov_sf[[cov_name]]
  # eps controls pixel resolution in meters; smaller = slower
  im_cov = idw(ppp.base, at = "pixels", eps = 1000)
  covariate_images[[cov_name]] <- im_cov
}

save(covariate_images, covariate_names, file = "output/PFAS/cov_extrapolation_im.RData")
