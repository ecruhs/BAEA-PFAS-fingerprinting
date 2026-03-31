# ---- Load libraries ----
library(sf)
library(dplyr)
library(units)
library(tigris)

options(tigris_use_cache = TRUE)

# ---- Step 1: Set working directory and load nest data ----
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# CSV with Nest_no, long, lat
nests <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert to sf points
nests_sf <- st_as_sf(nests, coords = c("long", "lat"), crs = 4326)  # WGS84

# ---- Step 2: Download all Wisconsin water polygons ----
wi_counties <- tigris::counties(state = "55")
county_fips <- wi_counties$COUNTYFP

# Lakes
water_lakes_list <- lapply(county_fips, function(fip){
  area_water(state = "55", county = fip)
})

# Rivers
water_rivers_list <- lapply(county_fips, function(fip){
  linear_water(state = "55", county = fip)
})

# Combine
wi_lakes <- do.call(rbind, water_lakes_list)
wi_rivers <- do.call(rbind, water_rivers_list)
waterbodies <- bind_rows(wi_lakes, wi_rivers)

# ---- Step 3: Use waterbody centroids to avoid zero distances ----
waterbodies_centroid <- st_centroid(waterbodies)

# ---- Step 4: Project both layers to planar CRS (meters) ----
nests_sf_proj <- st_transform(nests_sf, 3071)
waterbodies_centroid_proj <- st_transform(waterbodies_centroid, 3071)

# ---- Step 5: Calculate nearest distance and store nearest waterbody ID ----
nests_sf_proj <- nests_sf_proj %>%
  rowwise() %>%
  mutate(
    # Compute distances to all waterbody centroids
    distances = list(st_distance(geometry, waterbodies_centroid_proj)),
    nearest_index = which.min(distances[[1]]),
    nearest_water_dist_m = distances[[1]][nearest_index],
    nearest_water_ID = waterbodies_centroid_proj$HYDROID[nearest_index]
  ) %>%
  ungroup() %>%
  mutate(nearest_water_dist_km = set_units(nearest_water_dist_m, "km"))

# ---- Step 6: Save results to CSV ----
results_df <- nests_sf_proj %>%
  select(nest_no, nearest_water_ID, nearest_water_dist_m, nearest_water_dist_km) %>%
  as.data.frame()

write.csv(results_df, "output-data/nest_distances_to_water.csv", row.names = FALSE)
