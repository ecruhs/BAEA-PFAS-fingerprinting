library(sf)
library(dplyr)
library(units)

# 1. Read your point data
my_points <- read.csv("your_points.csv", stringsAsFactors = FALSE)
# Must have columns like lat, lon
points_sf <- st_as_sf(my_points, coords = c("lon","lat"), crs = 4326)

# 2. Get coastline or water polygons/lines
# For example, if you have a shapefile of Wisconsin lakes + Great Lakes coastline, lake polygons etc:
water <- st_read("path/to/wisconsin_waterbodies.shp")   # polygons of lakes, rivers, etc.
coastline <- st_read("path/to/wisconsin_coastline.shp") # lines

# 3. (Optional) Project both to a projected CRS (units in meters), this makes distance calculations accurate
points_proj <- st_transform(points_sf, crs = 3857)    # e.g. Web Mercator, or better local UTM / equal area
water_proj <- st_transform(water, crs = st_crs(points_proj))
coast_proj <- st_transform(coastline, crs = st_crs(points_proj))

# 4. Compute distance to nearest water polygon (if you want “any water”)
dist_to_water <- st_distance(points_proj, water_proj)  # this gives a matrix: each point × each water feature
# pick the minimum for each point
min_dist_water <- apply(dist_to_water, 1, min)

# 5. Also compute distance to coastline if that’s a separate layer
dist_to_coast <- st_distance(points_proj, coast_proj)
min_dist_coast <- apply(dist_to_coast, 1, min)

# 6. Add distances back to your data
results <- my_points %>%
  mutate(
    dist_to_water_m = set_units(min_dist_water, "m"),
    dist_to_coast_m = set_units(min_dist_coast, "m")
  )

# 7. Inspect or export results
head(results)
write.csv(results, "points_with_distances.csv", row.names = FALSE)
