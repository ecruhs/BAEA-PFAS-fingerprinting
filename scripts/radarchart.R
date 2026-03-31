# ---- Load libraries ----
library(dplyr)
library(fmsb)
library(scales)  # for alpha()

# ---- Set working directory ----
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load data ----
wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

# ---- Factorize region to preserve order ----
wi_data <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# ---- Define PFAS variables ----
vars <- c("FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",
          "PFBA","PFDA","PFDOA","PFDS","PFHPA","PFHPS",
          "PFHXA","PFHXS","PFNA","PFNS","PFOA","PFOS",
          "PFTEDA","PFTRDA","PFUNDA", "TOTAL_PFAS")

# ---- Color-blind-friendly palette ----
cb_colors <- c(
  "#FF5733",  # AIS - bright orange/red
  "#CCCCCC",  # SSS - light gray
  "#DDDDDD",  # FRGB - light gray
  "#BBBBBB",  # GBLM - muted gray
  "#AAAAAA",  # LWR - muted gray
  "#1F77B4",  # MWR - bright blue
  "#999999"   # UPWR - light gray
)

# ---- Compute mean PFAS per region in correct order ----
pfas_means <- wi_data %>%
  group_by(region) %>%
  summarize(across(all_of(vars), mean, na.rm = TRUE)) %>%
  arrange(region)  # ensures rows follow factor levels

# ---- Prepare matrix for radar chart ----
radar_data <- pfas_means[,-1]  # remove region column
radar_data <- rbind(
  apply(radar_data, 2, max),  # max values (first row)
  apply(radar_data, 2, min),  # min values (second row)
  radar_data                   # actual data
)
rownames(radar_data) <- c("max", "min", pfas_means$region)

# ---- Adjust margins for legend ----
op <- par(mar = c(1, 1, 1, 5))  # extra right margin

# ---- Plot radar chart ----
radar <- radarchart(
  radar_data,
  axistype = 1,
  pcol = cb_colors[1:nrow(pfas_means)],
  pfcol = alpha(cb_colors[1:nrow(pfas_means)], 0.3),
  plwd = 2,
  plty = 1,
  cglcol = "gray80",
  cglty = 1,
  cglwd = 0.8,
  vlcex = 0.8
)

# ---- Add legend outside chart ----
# legend(
#   title = "Region",
#   x = 1.35,        # move far right
#   y = 1,
#   legend = pfas_means$region,
#   bty = "n",
#   pch = 20,
#   col = cb_colors[1:nrow(pfas_means)],
#   pt.cex = 2,
#   cex = 0.8
# )

# ---- Reset margins ----
par(op)
