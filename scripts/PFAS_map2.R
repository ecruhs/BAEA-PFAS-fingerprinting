# ---- Libraries ----
library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(scatterpie)
library(viridisLite)

homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)
options(tigris_use_cache = TRUE)

# ---- Load WI counties ----
wi_counties <- tigris::counties(state = "55", cb = TRUE, class = "sf")

# ---- Load PFAS data ----
county_data <- read.csv(
  file = paste0(homewd, "data/averages_202320242025_PFASmap.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# ---- Join geometry with data ----
wi_data <- wi_counties %>%
  left_join(county_data, by = "NAME")

# ---- Compound columns ----
compound_cols <- c(
  "TOTAL_PFAS", "PFUNDA", "PFTRDA", "PFTEDA", "PFOS", "PFOA", "PFNS",
  "PFNA", "PFHXS", "PFHPS", "PFHPA", "PFDS", "PFDOA", "PFDA",
  "PFBA", "N_MEFOSAA", "N_ETFOSAA", "FOSA", "FTSA10_2", "FTSA8_2",
  "X5:3FTCA", "X7:3FTCA", "X8:2FTS", "PFHxA"
)
compound_cols <- intersect(compound_cols, names(wi_data))


# ---- Force numeric ----
wi_data <- wi_data %>%
  mutate(across(all_of(compound_cols), ~ suppressWarnings(as.numeric(.))))

# ---- Keep counties with PFAS detected ----
wi_sampled <- wi_data %>%
  filter(if_any(all_of(compound_cols), ~ !is.na(.) & . > 0))

# ---- Extract centroids ----
wi_sampled <- wi_sampled %>%
  mutate(centroid = st_centroid(geometry)) %>%
  mutate(
    lon = st_coordinates(centroid)[,1],
    lat = st_coordinates(centroid)[,2]
  )

# ---- Reproject to Wisconsin Albers for circular pies ----
scatter_data_sf <- st_as_sf(
  wi_sampled %>%
    st_drop_geometry() %>%
    select(lon, lat, all_of(compound_cols)),
  coords = c("lon", "lat"),
  crs = 4269
) %>%
  st_transform(3071) %>%
  mutate(
    lon_proj = st_coordinates(.)[,1],
    lat_proj = st_coordinates(.)[,2]
  )

# ---- Remove unwanted columns ----
cols_to_remove <- c("FTSA10_2", "FTSA8_2", "PFHxA")
scatter_data_sf <- scatter_data_sf %>%
  select(-all_of(cols_to_remove))

# ---- Define pie columns ----
pie_cols <- compound_cols[compound_cols != "TOTAL_PFAS"]
pie_cols_valid <- pie_cols[pie_cols %in% names(scatter_data_sf)]
pie_cols_valid <- pie_cols_valid[
  colSums(!is.na(st_drop_geometry(scatter_data_sf)[pie_cols_valid])) > 0
]

# ---- Colors ----
bright_colors <- turbo(length(pie_cols_valid))



# ---- Compute pie radii ----
scatter_data_sf <- scatter_data_sf %>%
  mutate(pie_r = sqrt(TOTAL_PFAS) * 1200)  # radius in meters

# ---- Repel pies to avoid overlaps ----
repel_positions <- function(lon, lat, r) {
  n <- length(lon)
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      dist <- sqrt((lon[i]-lon[j])^2 + (lat[i]-lat[j])^2)
      min_dist <- r[i] + r[j]
      if(dist < min_dist) {
        dx <- lon[j]-lon[i]
        dy <- lat[j]-lat[i]
        angle <- atan2(dy, dx)
        move <- (min_dist - dist)/2
        lon[i] <- lon[i] - cos(angle)*move
        lat[i] <- lat[i] - sin(angle)*move
        lon[j] <- lon[j] + cos(angle)*move
        lat[j] <- lat[j] + sin(angle)*move
      }
    }
  }
  return(data.frame(lon = lon, lat = lat))
}

repelled_coords <- repel_positions(scatter_data_sf$lon_proj, scatter_data_sf$lat_proj, scatter_data_sf$pie_r)
scatter_data_sf$lon_rep <- repelled_coords$lon
scatter_data_sf$lat_rep <- repelled_coords$lat

# ---- Dummy for size legend ----
scatter_data_size <- scatter_data_sf %>%
  mutate(size_scale = TOTAL_PFAS)


# ---- rename pfas
pfas_labels <- c(
  "TOTAL_PFAS" = "Total PFAS",
  "PFUNDA"     = "PFUnA",
  "PFTRDA"     = "PFTriA",
  "PFTEDA"     = "PFTeA",
  "PFOS"       = "PFOS",
  "PFOA"       = "PFOA",
  "PFNS"       = "PFNS",
  "PFNA"       = "PFNA",
  "PFHXS"      = "PFHxS",
  "PFHPS"      = "PFHPS",
  "PFHPA"      = "PFHPA",
  "PFDS"       = "PFDS",
  "PFDOA"      = "PFDoA",
  "PFDA"       = "PFDA",
  "PFBA"       = "PFBA",
  "N_MEFOSAA"  = "NMeFOSAA",
  "N_ETFOSAA"  = "NEtFOSAA",
  "FOSA"       = "FOSA",
  "FTSA10_2"   = "10:2 FTS",
  "FTSA8_2"    = "8:2 FTS",
  "X5:3FTCA"   = "5:3 FTCA",
  "X7:3FTCA"   = "7:3 FTCA",
  "X8:2FTS"    = "8:2 FTS",
  "PFHxA"      = "PFHxA"
)


# Example: assign FOSA a standout color
names(bright_colors) <- pie_cols_valid

# Override FOSA color
bright_colors["FOSA"] <- "violet" 


# ---- Final Plot ----
ggplot() +
  geom_sf(data = st_transform(wi_data, 3071), fill = "gray80", color = "white") +
  
  geom_scatterpie(
    aes(x = lon_rep, y = lat_rep, r = pie_r),
    data = st_drop_geometry(scatter_data_sf),
    cols = pie_cols_valid,
    color = NA
  ) +
  
  geom_point(
    data = st_drop_geometry(scatter_data_size),
    aes(x = lon_rep, y = lat_rep, size = size_scale),
    shape = 21, fill = "white", color = "black", alpha = 0
  ) +
  
  coord_sf(crs = st_crs(3071), datum = st_crs(4269), expand = FALSE) +
  scale_fill_manual(values = bright_colors, name = "PFAS Compounds", 
                   labels = pfas_labels[pie_cols_valid]) +
  scale_size_area(name = "Total PFAS (ng/g)", max_size = 10) +
  
  theme_minimal() +
  theme(
    legend.text = element_text(size=12),
    legend.position = c(1.3, 0.5),
    legend.box = "vertical",
    plot.margin = margin(5, 20, 50, 5),
    axis.text = element_text(size=14),
    axis.title = element_text(size=14)
  ) +
  guides(
    fill = guide_legend(order = 1),
    size = guide_legend(order = 2, override.aes = list(alpha = 1, fill = "white", color = "black"))
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "\nLongitude",
    y = "Latitude\n"
  )

# ---- Save map ----
ggsave(
  filename = paste0(homewd, "PFAS_Wisconsin_map.png"),
  width = 18,
  height = 12,
  dpi = 1000,
  units = "in",
  bg = "white",
  limitsize = FALSE
)
