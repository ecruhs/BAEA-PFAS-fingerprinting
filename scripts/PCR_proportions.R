library(dplyr)
library(ggplot2)
library(compositions)

# ---- set working directory
homewd <- "/Users/gavindehnert/PFAS_2025/"
setwd(homewd)

# ---- Load data ----
wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)
wi_data <- subset(wi_data, bird_id!="LI30-001")

wi_data <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# --- Select just the PCA variables (Z through BT) ---
# Adjust column numbers if needed (Z = 26th col, BT = 72nd col)
vars <- (c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA", "N_MEFOSAA",
           "PFBA", "PFBS", "PFDA",
           "PFDOA","PFDS","PFHPA","PFHPS",
           "PFHXA","PFHXS","PFNA","PFNS","PFOA",
           "PFOS","PFPES", "PFTEDA","PFTRDA","PFUNDA", "TOTAL_PFAS"))

pca_vars <- wi_data %>%
  select(all_of(vars)) 



# --- Create proportions out of TOTAL_PFAS ---
prop_data <- pca_vars %>%
  mutate(across(-TOTAL_PFAS, ~ .x / TOTAL_PFAS))

# --- Remove TOTAL_PFAS column ---
prop_data <- prop_data %>%
  select(-TOTAL_PFAS)

# --- Replace zeros (required for CLR) ---
prop_data[prop_data == 0] <- 1e-6

# --- Apply centered log-ratio transformation ---
clr_data <- clr(as.matrix(prop_data))

# Optional: convert back to dataframe
clr_data <- as.data.frame(clr_data)

pca_res <- prcomp(clr_data, scale. = TRUE)
summary(pca_res)

pca_scores <- as.data.frame(pca_res$x)
# limit to first 4
pca_use <- pca_scores[, 1:4] 

library(mclust)

# start with pca scores
# assume you've already run
# pca_res <- prcomp(my_data, scale. = TRUE)


pca_scores <- as.data.frame(pca_res$x)

# use the top PCs to reduce noise (typical choice: 3–5 PCs)
pca_use <- pca_scores[, 1:4]

# set row names
# pca_use <- pca_use %>% 
#   mutate(row_id = row_number())
# 
# rownames(pca_use) <- rownames(wi_data[complete.cases(wi_data[,vars_used_in_PCA]), ])

# run model based clustering
set.seed(123)  # for reproducibility

mclust_res <- Mclust(pca_use)  # automatic model and cluster number selection
summary(mclust_res)
# Mclust VEV (ellipsoidal, equal shape) model with 6 components: 
#   log-likelihood   n df       BIC       ICL
# -550.4622 114 74 -1451.403 -1454.728

mclust_res$modelName   # "VEV"
mclust_res$G #6

uncertainty <- mean(1 - apply(mclust_res$z, 1, max))
uncertainty #0.0131
mean_max_prob <- mean(apply(mclust_res$z, 1, max))
mean_max_prob #0.99
# 
# Clustering table:
#   1  2  3  4  5  6 
# 40 15 25 18  8  8 

# inspect results
# Number of clusters selected
mclust_res$G # this selects for 6 clusters

# Cluster assignment for each observation
head(mclust_res$classification) 

# BIC values for different models
plot(mclust_res, what = "BIC")

plot(mclust_res, what = "classification")  # nice auto plot



ggplot(data.frame(pca_use, cluster = factor(mclust_res$classification)),
       aes(x = PC1, y = PC2, color = cluster)) +
  #stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.1, color = NA) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(color = "Cluster", title = "Model-based clustering (Mclust)")

# add cluster assignment back to original dataset
wi_data$cluster <- factor(mclust_res$classification)

vars <- (c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA", "N_MEFOSAA",
           "PFBA", "PFBS", "PFDA",
           "PFDOA","PFDS","PFHPA","PFHPS",
           "PFHXA","PFHXS","PFNA","PFNS","PFOA",
           "PFOS","PFPES", "PFTEDA","PFTRDA","PFUNDA"))

pfas_summary <- wi_data %>%
  group_by(cluster) %>%
  summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)))

pfas_summary

# heat map of contribution by cluster
library(pheatmap)

pfas_mat <- as.matrix(pfas_summary[,-1])        # drop cluster column
rownames(pfas_mat) <- paste0("Cluster_", pfas_summary$cluster)

pheatmap(
  pfas_mat,
  scale = "row",                                # standardize within PFAS
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "PFAS profiles by Mclust cluster"
)

# barplot of contribution by cluster
library(tidyr)
library(ggplot2)

pfas_long <- pfas_summary %>%
  pivot_longer(-cluster, names_to = "PFAS", values_to = "mean_conc")

pfas_long2 <- subset(pfas_long, PFAS!="PFOS")

library(ggbreak)

pfas_long2 <- subset(pfas_long, PFAS!="PFOS")
ggplot(pfas_long2, aes(x = PFAS, y = mean_conc, fill = cluster)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    #title = "PFAS concentration by cluster",
    y = "Mean PFAS concentration",
    x = "PFAS compound"
  )


ggplot(pfas_long, aes(x = PFAS, y = mean_conc, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_y_break(c(60, 200)) +   # adjust upper value as needed
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    y = "Mean PFAS concentration",
    x = "PFAS compound"
  )



pfas_long %>%
  filter(PFAS == "PFOS") %>%
  group_by(cluster) %>%
  summarise(mean_PFOS = round(mean(mean_conc, na.rm = TRUE), 4))

# # A tibble: 6 × 2
# cluster mean_PFOS
# <fct>       <dbl>
# 1 1            140.
# 2 2            151.
# 3 3            210.
# 4 4            305.
# 5 5            176.
# 6 6            155.

# rank pfas compounds by cluster separation (Kruskal-Wallis)
library(purrr)

pfas_tests <- map_dfr(
  vars,
  function(var) {
    test <- kruskal.test(as.formula(paste(var, "~ cluster")), data = wi_data)
    data.frame(PFAS = var, p_value = test$p.value)
  }
) %>%
  arrange(p_value) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

# The PFAS compounds with the smallest p-values are the strongest drivers of cluster structure.
pfas_tests
# ✅ Pro tip: combine the statistical ranking (pfas_tests) with the heatmap to interpret which PFAS define each cluster (e.g., “Cluster 1 is high in PFOS/PFOA, Cluster 2 dominated by short-chain PFAS”).





#### ------- PCA with cluster as points and color as regions ####
cb_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#D55E00",  # vermillion/red-orange
  "#0072B2",  # blue
  "#CC79A7"   # reddish purple
)

plot <- ggplot(
  data.frame(
    pca_use,
    cluster = factor(mclust_res$classification),
    region = wi_data$region
  ),
  aes(x = PC1, y = PC2)
) +
  # Ellipses by region
  stat_ellipse(
    aes(fill = region),
    geom = "polygon",
    alpha = 0.15,
    color = NA
  ) +
  # Points: color by region, shape by cluster
  geom_point(aes(color = region, shape = cluster), size = 2.5) +
  scale_color_manual(values = cb_colors) +
  scale_fill_manual(values = cb_colors) +
  theme_minimal() +
  labs(
    color = "Region",
    shape = "Cluster",
    fill = "Region"
    #title = "Mclust clusters (shape) and region grouping (color/ellipse)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size=14)
  )

print(plot)

ggsave("pca_plot_cluster.png", plot, width = 8, height = 8, bg="white",
       dpi=600)
