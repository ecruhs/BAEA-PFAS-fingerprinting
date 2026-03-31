library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

homewd <- "/Users/gavindehnert/PFAS_2025/"
setwd(homewd)
options(tigris_use_cache = TRUE)

# ---- Load WI counties ----
wi_counties <- tigris::counties(state = "55", cb = TRUE, class = "sf") %>%
  st_transform(3071)  # Wisconsin Albers projection

# ---- Load PFAS sampling data ----
sample_data <- read.csv(
  file = paste0(homewd, "data/Allyears_PFAS.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

sample_data <- sample_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))


# ---- Filter out points outside Wisconsin bounds ----
sample_data_clean <- sample_data %>%
  filter(long >= -93.5, long <= -86.5, lat >= 42.5, lat <= 47)

# ---- Keep only rows with PFAS detected ----
compound_cols <- c(
  "TOTAL_PFAS", "PFUNDA", "PFTRDA", "PFTEDA", "PFOS", "PFOA", "PFNS",
  "PFNA", "PFHXS", "PFHPS", "PFHPA", "PFDS", "PFDOA", "PFDA",
  "PFBA", "N_MEFOSAA", "N_ETFOSAA", "FOSA", "FTSA10_2", "FTSA8_2",
  "X5:3FTCA", "X7:3FTCA", "X8:2FTS", "PFHxA"
)
compound_cols <- intersect(compound_cols, names(sample_data_clean))

sample_data_clean <- sample_data_clean %>%
  mutate(across(all_of(compound_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(if_any(all_of(compound_cols), ~ !is.na(.) & . > 0))

# ---- Convert to sf points ----
sample_points <- st_as_sf(
  sample_data_clean,
  coords = c("long", "lat"),
  crs = 4269
) %>%
  st_transform(3071)

### set color scale
cb_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#D55E00",  # vermillion/red-orange
  "#0072B2",  # blue
  "#CC79A7"   # reddish purple
)


# ---- Plot map ----
sampling_locations <- 
  ggplot() +
  geom_sf(data = wi_counties, fill = "white", color = "gray9") +
  geom_sf(data = sample_points, aes(color = region), size = 2, alpha = 0.8) +
  scale_color_manual(values = cb_colors) +  # use colorblind-friendly palette
  coord_sf(datum = NA) +  # fixes aspect ratio
  theme_minimal() +
  theme(
    legend.text = element_text(size=13),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(color = "Region")



##### combine the map and the boxplot
library(cowplot)
library(ggplot2)

# Example: assuming these are your plot objects
# p1 <- ggplot(...) + geom_sf(...) + ...   # Sampling locations map
# p2 <- ggplot(...) + geom_boxplot(...) + ... # PFAS boxplot

# Combine plots
# combined_plot <- plot_grid(
#   sampling_locations, boxplot,
#   labels = c("A", "B"),        # add labels
#   label_size = 14,             # size of labels
#   ncol = 2,                    # stack vertically
#   rel_widths = c(1, 1.5)        # adjust relative heights: map smaller, boxplots larger
# )

# Save or display
# combined_plot
# To save:
# ggsave("combined_plot.png", combined_plot, width = 17, height = 9, bg="white",
#        dpi=600)

ggsave("sampling_plot_simple.png", sampling_locations, width = 10, height = 9, bg="white",
       dpi=600)
