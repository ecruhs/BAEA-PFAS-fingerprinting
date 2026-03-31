library(dplyr)
library(ggplot2)

# ---- set working directory
homewd <- "/Users/gavindehnert/PFAS_2025/"
setwd(homewd)

# ---- Load data ----
wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

wi_data <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# --- Select just the PCA variables (Z through BT) ---
# Adjust column numbers if needed (Z = 26th col, BT = 72nd col)
vars <- (c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA", "N_MEFOSAA",
           "PFBA", "PFBS", "PFDA",
           "PFDOA","PFDS","PFHPA","PFHPS",
           "PFHXA","PFHXS","PFNA","PFNS","PFOA",
           "PFOS","PFPES", "PFTEDA","PFTRDA","PFUNDA"))


pca_vars <- wi_data %>%
  select(all_of(vars))



# --- Remove rows with NAs in those columns ---
pca_vars <- na.omit(pca_vars)

# --- Run PCA (scale = TRUE standardizes the variables) ---
pca_res <- prcomp(pca_vars, scale. = TRUE)

# --- Create a dataframe for plotting ---
# ensure we only use rows with complete data for the PCA variables
complete_idx <- complete.cases(wi_data[, vars])

# run PCA
pca_res <- prcomp(wi_data[complete_idx, vars], scale. = TRUE)

# build PCA dataframe
pca_df <- data.frame(
  pca_res$x,
  region = wi_data$region[complete_idx]
)

cb_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#D55E00",  # vermillion/red-orange
  "#0072B2",  # blue
  "#CC79A7"   # reddish purple
)


# --- Plot PCA, colored by county ---
ggplot(pca_df, aes(x = PC1, y = PC2, color = region)) +  # replace 'region' with your categorical variable
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cb_colors) +  # use your colorblind-friendly palette
  labs(
    color = "Region",
    x = paste0("PC1 (", round(100*summary(pca_res)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca_res)$importance[2,2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = region, fill = region)) +
  # Add ellipses first so points are on top
  stat_ellipse(aes(fill = region), geom = "polygon", alpha = 0.1, color = NA) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = cb_colors) +
  scale_fill_manual(values = cb_colors) +
  labs(
    color = "Region",
    fill = "Region",
    x = paste0("PC1 (", round(100*summary(pca_res)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100*summary(pca_res)$importance[2,2], 1), "%)")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

pca_plot









# Loadings for the first few PCs
loadings <- pca_res$rotation   # each column = PC, each row = original PFAS

# Check top contributors for PC1 and PC2
top_PC1 <- sort(abs(loadings[, "PC1"]), decreasing = TRUE)
top_PC2 <- sort(abs(loadings[, "PC2"]), decreasing = TRUE)

top_PC1
top_PC2



#### loading plot
library(ggplot2)
library(tidyr)

# ---- Convert loadings to long format ----
loadings_df <- as.data.frame(loadings) %>%
  mutate(PFAS = rownames(.)) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "Loading") %>%
  filter(PC %in% c("PC1", "PC2", "PC3", "PC4"))   # limit to PC1 and PC2

# ---- Plot ----
S1 <- ggplot(loadings_df, aes(x = PFAS, y = Loading, fill = PC)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("black", "gray30", "gray50", "gray75")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),
        axis.text.y = element_text(size=12)) +
  labs(
    x = "",
    y = "PCA Loading\n",
    fill = "Principal Component"
  ) +
  theme(legend.position = "right")

S1


ggsave("plotS1.png", S1, width = 12, height = 10, bg="white",
       dpi=600)





# ---- Scale the loadings for plotting arrows ----
# Multiply by a constant to fit arrows nicely on the PCA plot

# ---- Scale the loadings for plotting arrows ----
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- PCA variables ----
vars <- c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",
          "PFBA","PFBS","PFDA","PFDOA","PFDS","PFHPA","PFHPS",
          "PFHXA","PFHXS","PFNA","PFNS","PFOA","PFOS","PFPES",
          "PFTEDA","PFTRDA","PFUNDA")

# ---- PFAS groups ----
pfas_groups <- data.frame(
  PFAS = vars,
  Group = c(
    rep("Precursors",5),
    rep("Short-chain",3),
    "Legacy Acids", "Legacy Sulfonates", "Legacy Acids",
    rep("Legacy Acids",5),
    rep("Legacy Acids",4),
    rep("Legacy Sulfonates",2)
  ),
  stringsAsFactors = FALSE
)

# ---- Load data ----
wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

# ---- Prepare PCA ----
pca_vars <- wi_data %>%
  select(all_of(vars)) %>%
  na.omit()  # remove rows with NAs

pca_res <- prcomp(pca_vars, scale. = TRUE)

# ---- Extract loadings for PC1 and PC2 ----
loadings <- pca_res$rotation[,1:2]  # PC1 and PC2
loadings_df <- as.data.frame(loadings) %>%
  mutate(PFAS = rownames(.)) %>%
  left_join(pfas_groups, by = "PFAS")

# ---- Scale arrows for visibility ----
arrow_scale <- 1.5
loadings_df <- loadings_df %>%
  mutate(PC1 = PC1 * arrow_scale, PC2 = PC2 * arrow_scale)



# ---- Manual label replacements ----
label_map <- c(
  FTSA6_2 = "6:2 FTSA",
  FTSA8_2 = "8:2 FTSA",
  FOSA = "FOSA",
  N_ETFOSAA = "NEtFOSAA",
  N_MEFOSAA = "NMeFOSAA",
  PFBA = "PFBA",
  PFBS = "PFBS",
  PFDA = "PFDA",
  PFDOA = "PFDoA",
  PFDS = "PFDS",
  PFHPA = "PFHpA",
  PFHPS = "PFHpS",
  PFHXA = "PFHxA",
  PFHXS = "PFHxS",
  PFNA = "PFNA",
  PFNS = "PFNS",
  PFOA = "PFOA",
  PFOS = "PFOS",
  PFPES = "PFPeS",
  PFTEDA = "PFTeA",
  PFTRDA = "PFTriA",
  PFUNDA = "PFUnA"
)

loadings_df <- loadings_df %>%
  mutate(PlotLabel = recode(PFAS, !!!label_map))

# ---- Plot ----
S2 <- ggplot(loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group)) +
  geom_segment(arrow = arrow(length = unit(0.25, "cm")), size = 1) +
  geom_text_repel(
    aes(x = PC1, y = PC2, label = PlotLabel),
    size = 4,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Precursors" = "#E69F00",
    "Short-chain" = "#56B4E9",
    "Legacy Acids" = "#0072B2",
    "Legacy Sulfonates" = "#D55E00"
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(
    x = "\nPC1",
    y = "PC2\n",
    color = "",
    title = ""
  ) +
  coord_fixed()

print(S2)


ggsave("plotS2_renamed.png", S2, width = 10, height = 10, bg="white",
       dpi=600)


### ---- put them together
library(cowplot)
combined_plot <- plot_grid(
  pca_plot, S2,
  labels = c("A", "B"),        # add labels
  label_size = 14,             # size of labels
  ncol = 2,                    # stack vertically
  rel_widths = c(1, 1) ,
  rel_heights = c(0.9,1.7) # adjust relative heights: map smaller, boxplots larger
)

# Save or display
combined_plot
# To save:
ggsave("combined_pca_plot.png", combined_plot, width = 17, height = 9, bg="white",
       dpi=600)




#### ---- EXTRACT Q score per bird
# PCA Residual Distance (Q-residuals)
# If you want a measure of how unusual a bird is given the PCA structure, you can use the reconstruction error (residuals I showed earlier). This is sometimes called the Q statistic: 👉 Q_residual = how much PFAS variation in that bird is not captured by the first k PCs.
#👉 Large values = more "outlier-like" birds.
k <- 2  # number of PCs you want to retain
scaled_data <- scale(wi_data[complete_idx, vars], 
                     center = pca_res$center, 
                     scale = pca_res$scale)

reconstructed <- pca_res$x[, 1:k] %*% t(pca_res$rotation[, 1:k])
residuals_k <- scaled_data - reconstructed

# one value per bird = squared residual distance
Q_resid <- apply(residuals_k^2, 1, sum)

bird_df <- data.frame(
  BirdID = wi_data$bird_id[complete_idx],
  region = wi_data$region[complete_idx],
  Q_residual = Q_resid
)



##### ----
# PCA Scores (Projection into PC space)
# Each bird already has PCA scores (pca_res$x), one per PC.
# If you want one value per bird, you can take the distance in PC space, e.g.:
# keep the first k PCs (explaining most variance). PCA_score is a scaled summary of where each bird lies in the reduced PC space (higher values = farther from the average PFAS profile).
k <- 2
scores <- pca_res$x[, 1:k]

# compute Euclidean distance from the origin (the "center bird")
bird_score <- apply(scores, 1, function(row) sqrt(sum(row^2)))

# attach metadata
bird_df <- data.frame(
  BirdID = wi_data$bird_id[complete_idx],
  region = wi_data$region[complete_idx],
  PCA_score = bird_score,
  Q_residual = Q_resid
)

write.csv(bird_df, 
          "output-data/PFAS_PCA_values.csv", 
          row.names = FALSE)
