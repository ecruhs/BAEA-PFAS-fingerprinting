library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(viridis)
library(colorspace)

homewd <- "/Users/emilyruhs/Desktop/"
setwd(homewd)

options(tigris_use_cache = TRUE)

# Get counties in Wisconsin
wi_counties <- tigris::counties(state = "55", cb = TRUE, class = "sf")

# Read county PFAS data
county_data <- read.csv(file = paste0(homewd, "averages_20232024_PFASmap.csv"), 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

# Join county geometries with PFAS data
wi_data <- wi_counties %>%
  left_join(county_data, by = "NAME")

# Compute centroids for pie chart placement
wi_data <- wi_data %>%
  mutate(
    lon = st_coordinates(st_centroid(geometry))[,1],
    lat = st_coordinates(st_centroid(geometry))[,2]
  )

# PFAS compound columns
compound_cols <- c("TOTAL_PFAS", "PFUNDA", "PFTRDA", "PFTEDA", "PFOS", "PFOA", "PFNS",
                   "PFNA", "PFHXS", "PFHPS", "PFHPA", "PFDS", "PFDOA", "PFDA",
                   "PFBA", "N_MEFOSAA", "N_ETFOSAA", "FOSA", "FTSA10_2", "FTSA8_2")

# Ensure PFAS columns are numeric
wi_data <- wi_data %>%
  mutate(across(all_of(compound_cols), ~ as.numeric(as.character(.))))

# Keep only counties with at least one non-NA, non-zero PFAS value
wi_sampled <- wi_data %>%
  filter(if_any(all_of(compound_cols), ~ !is.na(.) & . > 0))

# Pie columns = all compounds except TOTAL_PFAS
pie_cols <- compound_cols[compound_cols != "TOTAL_PFAS"]

# Create a separate data frame for scatterpie (only numeric + lon/lat)
scatter_data <- wi_sampled %>%
  st_drop_geometry() %>%    
  select(lon, lat, all_of(compound_cols))

# Add jittered coordinates for pies
set.seed(123)  # reproducibility
scatter_data <- scatter_data %>%
  mutate(
    lon_jitter = jitter(lon, amount = 0.15),
    lat_jitter = jitter(lat, amount = 0.15)
  )

# Colors
bright_colors <- viridisLite::turbo(length(pie_cols))

# Dummy data for size legend
scatter_data_size <- scatter_data %>%
  mutate(size_scale = TOTAL_PFAS)

# ---- Final Plot ----
ggplot() +
  # Map
  geom_sf(data = wi_data, fill = "gray95", color = "white") +
  # Pies for counties
  geom_scatterpie(
    aes(x = lon_jitter, y = lat_jitter, r = sqrt(TOTAL_PFAS)/80),
    data = scatter_data,
    cols = pie_cols,
    color = NA
  ) +
  # Dummy points for size legend only
  geom_point(
    data = scatter_data_size,
    aes(x = lon_jitter, y = lat_jitter, size = size_scale),
    shape = 21, fill = "white", color = "black",
    alpha = 0  # invisible on map
  ) +
  # Map coordinates
  coord_sf(crs = st_crs(4269)) +
  # Color legend
  scale_fill_manual(
    values = bright_colors,
    name = "PFAS Compounds"
  ) +
  # Size legend
  scale_size_area(
    name = "Total PFAS (ng/g)",
    max_size = 10
  ) +
  # Layout & theme
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(5, 20, 50, 5), # extra space for bottom legends
    axis.text = element_text(),
    axis.title = element_text()
  ) +
  guides(
    fill = guide_legend(order = 1),
    size = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, fill = "white", color = "black")
    )
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "\nLongitude",
    y = "Latitude\n"
  )

