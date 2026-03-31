# Load libraries
library(sf)
library(dplyr)
library(units)
library(tigris)

options(tigris_use_cache = TRUE)

# ---- Step 1: Load your nest coordinates ----
# Example: assuming you have a CSV with columns Nest_no, lon, lat
# ---- set working directory
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load data ----
nests <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

nests_sf <- st_as_sf(nests, coords = c("long", "lat"), crs = 4326)  # WGS84


# ---- 2. Get all Wisconsin counties ----
wi_counties <- tigris::counties(state = "55")
county_fips <- wi_counties$COUNTYFP

# ---- 3. Download water for all counties ----
water_lakes_list <- lapply(county_fips, function(fip){
  area_water(state = "55", county = fip)
})
water_rivers_list <- lapply(county_fips, function(fip){
  linear_water(state = "55", county = fip)
})

# Combine into single layers
wi_lakes <- do.call(rbind, water_lakes_list)
wi_rivers <- do.call(rbind, water_rivers_list)


# ---- 4. Filter to your waterbodies of interest ----
# names(wi_lakes)
# 
# water_lakes <- wi_lakes %>%
#   filter(NAME %in% c("Lake Superior", "Lake Michigan", "Lake Winnebago"))
# 
# water_rivers <- wi_rivers %>%
#   filter(NAME %in% c("Wisconsin River", "Fox River"))

# waterbodies <- bind_rows(water_lakes, water_rivers)
waterbodies <- bind_rows(wi_lakes, wi_rivers)

# ---- 5. Project to planar CRS ----
nests_sf_proj <- st_transform(nests_sf, 3071)
waterbodies_proj <- st_transform(waterbodies, 3071)

# ---- 6. Calculate distance to nearest waterbody ----
nests_sf_proj <- nests_sf_proj %>%
  rowwise() %>%
  mutate(nearest_water_dist_m = min(st_distance(geometry, waterbodies_proj))) %>%
  ungroup() %>%
  mutate(nearest_water_dist_km = set_units(nearest_water_dist_m, "km"))

# ---- 7. View results ----
nests_sf_proj %>%
  select(nest_no, nearest_water_dist_m, nearest_water_dist_km) %>%
  as.data.frame()

# ---- 8. Save results -----
# Convert to a data frame
results_df <- nests_sf_proj %>%
  select(nest_no, nearest_water_dist_m, nearest_water_dist_km) %>%
  as.data.frame()

# Save to CSV
write.csv(results_df, "output-data/nest_distances_to_water.csv", row.names = FALSE)

