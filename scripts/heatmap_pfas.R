library(pheatmap)
library(dplyr)

homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load PFAS data ----
wi_data <- read.csv(
  file = paste0(homewd, "data/Allyears_PFAS.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# --- Define PFAS variables only (no region here)
pfas_vars <- c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",
               "PFBA","PFBS","PFDA","PFDOA","PFDS","PFHPA","PFHPS",
               "PFHXA","PFHXS","PFNA","PFNS","PFOA","PFOS","PFPES",
               "PFTEDA","PFTRDA","PFUNDA")

# --- Summarize PFAS by region
pfas_matrix <- wi_data %>%
  group_by(region) %>%
  summarize(across(all_of(pfas_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

# --- Set rownames to region
rownames(pfas_matrix) <- pfas_matrix$region
pfas_matrix <- pfas_matrix[, -1]

# --- Heatmap (scaled by column)
pheatmap(scale(pfas_matrix),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2")





################################
#### making it more complex
library(pheatmap)
library(dplyr)

# --- Define PFAS variables ---
pfas_vars <- c(
  "FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",   # Precursors
  "PFBA","PFBS","PFHPA",                                 # Short-chain
  "PFHPS","PFDS","PFHXS","PFNS","PFHXA","PFDA","PFDOA",
  "PFNA","PFOA","PFOS","PFPES","PFTEDA","PFTRDA","PFUNDA" # Legacy
)

# --- Summarize PFAS by region ---
pfas_matrix <- wi_data %>%
  group_by(region) %>%
  summarize(across(all_of(pfas_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()

# Set rownames to region
rownames(pfas_matrix) <- pfas_matrix$region
pfas_matrix <- pfas_matrix[, -1]

# --- Assign PFAS groups ---
pfas_groups <- data.frame(
  Group = c(
    rep("Precursors", 5),
    rep("Short-chain", 3),
    rep("Legacy", 14)
  )
)
rownames(pfas_groups) <- colnames(pfas_matrix)

# --- Reorder columns by PFAS group so colors are contiguous ---
pfas_matrix <- pfas_matrix[, rownames(pfas_groups)]

# --- Color scheme for groups ---
ann_colors <- list(
  Group = c(
    "Legacy" = "#222222",       # dark gray
    "Short-chain" = "#888888",  # medium gray
    "Precursors" = "#CCCCCC"    # light gray
  )
)

# --- Heatmap with column annotation ---
pheatmap(
  scale(pfas_matrix),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_col = pfas_groups,
  annotation_colors = ann_colors
)


