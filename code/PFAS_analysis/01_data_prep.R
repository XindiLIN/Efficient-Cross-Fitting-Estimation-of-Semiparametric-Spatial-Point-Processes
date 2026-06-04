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
source('code/simulation/semi_spp_functions.R')

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  remotes::install_github("ropensci/rnaturalearthhires")
}

#### load data
df_IPP = read.csv("data/PFAS/IPP Ready Model Data.csv")
df_ext = read.csv("data/PFAS/IPP Ready Extrap Data.csv")

##### get data in regions with enough observations
PA_states <- c("California", "Colorado", "Maryland", "Massachusetts",
               "Michigan", "New Hampshire", "New Jersey", "North Carolina",
               "Ohio", "Rhode Island", "South Carolina", "West Virginia", "Wisconsin")
usa_states <- ne_states(country = "united states of america", returnclass = "sf")
PA_states_geometry <- usa_states[usa_states$name %in% PA_states, ]
df_PA <- df_IPP %>%
  filter(STUSPS %in% state.abb[match(PA_states, state.name)])

###### draw locations of PFAS
sf_PFAS = st_as_sf(df_IPP, coords = c("Longitude", "Latitude"), crs = 4326)

us_mainland <- usa_states[!usa_states$name %in% c("Alaska", "Hawaii"), ]
us_mainland_ll <- st_transform(us_mainland, crs = 4326)
ggplot() +
  geom_sf(data = us_mainland_ll, fill = 'ivory', color = 'grey80') +
  geom_sf(data = PA_states_geometry, fill = NA, color = "green", linewidth = 0.6) +
  geom_sf(data = sf_PFAS[sf_PFAS$Concentration < 8, ], color = 'grey', size = 0.05) +
  geom_sf(data = sf_PFAS[sf_PFAS$Concentration >= 8, ], color = 'darkred', size = 0.05) +
  theme_void()

# Find locations where PFAS concentration >= 8 ppt
PFAS_threshold = 8
df_PA_ab  = df_PA[df_PA$Concentration >= PFAS_threshold, ]
df_All_ab = df_IPP[df_IPP$Concentration >= PFAS_threshold, ]

sf_PA_ab  = st_as_sf(df_PA_ab,  coords = c("Longitude", "Latitude"), crs = 4326)
sf_All_ab = st_as_sf(df_All_ab, coords = c("Longitude", "Latitude"), crs = 4326)

# transfer to ppp class in spatstat
target_epsg_code = 26910  # NAD83 / UTM Zone 10N (California-centric; consider EPSG 5070 for multi-state)
PA_states_geometry_proj = st_transform(PA_states_geometry, crs = target_epsg_code)
window_PA = as.owin(PA_states_geometry_proj$geometry)
pp_PFAS_PA = as.ppp(st_transform(sf_PA_ab, crs = target_epsg_code)$geometry, W = window_PA)

save(df_IPP, df_ext, df_PA,
     PA_states, PA_states_geometry, PA_states_geometry_proj,
     usa_states,
     sf_PFAS, sf_PA_ab, sf_All_ab,
     pp_PFAS_PA, window_PA, target_epsg_code,
     PFAS_threshold,
     file = "output/PFAS/prepared_data.RData")
