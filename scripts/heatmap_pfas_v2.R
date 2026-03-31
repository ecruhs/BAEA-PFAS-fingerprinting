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




# ---- Define PFAS variables ----
pfas_vars <- c(
  "FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",   # Precursors
  "PFBA","PFBS","PFHPA",                                 # Short-chain
  "PFHPS","PFDS","PFHXS","PFNS","PFHXA","PFDA","PFDOA", # Legacy acids
  "PFNA","PFOA","PFOS","PFPES","PFTEDA","PFTRDA","PFUNDA" # Legacy continued
)

# ---- Summarize PFAS by region ----
pfas_matrix <- wi_data %>%
  dplyr::group_by(region) %>%
  dplyr::summarise_at(vars(pfas_vars), ~ mean(.x, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  as.data.frame()




# 1. PFAS matrix with regions as rows
rownames(pfas_matrix) <- pfas_matrix$region
pfas_matrix <- pfas_matrix[, pfas_vars]

# 2. Column labels for display
colnames(pfas_matrix) <- pfas_labels[colnames(pfas_matrix)]

# 3. PFAS groups
pfas_groups <- data.frame(
  Group = c(rep("Precursors",5), rep("Short-chain",3), rep("Legacy",14)),
  row.names = pfas_vars,
  stringsAsFactors = FALSE
)

# 4. Rename annotation row names to match column labels
rownames(pfas_groups) <- pfas_labels[rownames(pfas_groups)]

# 5. Annotation colors
ann_colors <- list(
  Group = c(
    "Precursors" = "#CCCCCC",
    "Short-chain" = "#888888",
    "Legacy" = "#222222"
  )
)

# 6. Heatmap
pheatmap(
  scale(pfas_matrix),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_col = pfas_groups,
  annotation_colors = ann_colors,
  fontsize_col = 10,
  fontsize_row = 10
)
