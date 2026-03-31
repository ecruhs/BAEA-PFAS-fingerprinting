library(dplyr)
library(ggplot2)
library(reshape2)

# ---- set working directory
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load data
wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

wi_data <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# ---- PFAS variables (numeric only)
pfas_vars <- c(
  "FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",   # Precursors
  "PFBA","PFBS","PFHPA",                                 # Short-chain
  "PFHXA","PFOA","PFNA","PFDA","PFDOA","PFTRDA","PFTEDA","PFUNDA",
  "PFHXS","PFOS","PFPES","PFHPS","PFDS"         # Legacy
)

# ---- Pretty names for PFAS
vars_renamed <- c(
  "FOSA"       = "FOSA",
  "FTSA6_2"    = "6:2 FTS",
  "FTSA8_2"    = "8:2 FTS",
  "PFBS"       = "PFBS",
  "N_ETFOSAA"  = "NEtFOSAA",
  "N_MEFOSAA"  = "NMeFOSAA",
  "PFBA"       = "PFBA",
  "PFDA"       = "PFDA",
  "PFDOA"      = "PFDoA",
  "PFDS"       = "PFDS",
  "PFHPA"      = "PFHPA",
  "PFHPS"      = "PFHPS",
  "PFHXA"      = "PFHxA",
  "PFHXS"      = "PFHxS",
  "PFNA"       = "PFNA",
  "PFNS"       = "PFNS",
  "PFOA"       = "PFOA",
  "PFOS"       = "PFOS",
  "PFTEDA"     = "PFTeA",
  "PFTRDA"     = "PFTriA",
  "PFUNDA"     = "PFUnA",
  "PFPES"      = "PFPeS"
)

# ---- Select PFAS columns
df <- wi_data %>%
  select(region, all_of(pfas_vars))

# ---- Compute correlations by region and melt
cor_long <- df %>%
  group_by(region) %>%
  group_modify(~ {
    mat <- cor(.x[, pfas_vars], use = "pairwise.complete.obs")
    
    # Keep only lower triangle
    mat[upper.tri(mat)] <- NA
    
    melt(mat) %>%
      as.data.frame() %>%
      mutate(region = unique(.x$region))
  }) %>%
  ungroup()

# ---- Filter out NAs and unused PFAS
cor_long <- cor_long %>%
  filter(!is.na(value)) %>%
  filter(Var1 %in% names(vars_renamed) & Var2 %in% names(vars_renamed))

# ---- Factor PFAS in desired order for axes
cor_long$Var1 <- factor(cor_long$Var1, levels = pfas_vars, labels = vars_renamed[pfas_vars])
cor_long$Var2 <- factor(cor_long$Var2, levels = pfas_vars, labels = vars_renamed[pfas_vars])

# ---- Plot
ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "black") +  # outline each heatmap cell
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#A44400",
    midpoint = 0, limit = c(-1,1)
  ) +
  facet_wrap(~ region, scales = "free", ncol = 3, strip.position = "top") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold", size = 10)  # bold region titles
  ) +
  labs(fill = "Correlation", x = "", y = "")







# Get top correlations across all regions
library(dplyr)
library(knitr)


# Your PFAS variables
pfas_vars <- c("FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",
               "PFBA","PFBS","PFHPA",
               "PFHXA","PFOA","PFNA","PFDA","PFDOA","PFTRDA","PFUNDA",
               "PFHXS","PFOS","PFPES","PFHPS","PFDS","PFTEDA")

# PFAS groups
pfas_groups <- data.frame(
  PFAS = pfas_vars,
  Group = c(rep("Precursor",5),
            rep("Short-chain",3),
            rep("Legacy Acids",7),
            rep("Legacy Sulfonates",3),
            rep("Legacy Sulfonamides",3)),
  stringsAsFactors = FALSE
)

cor_all <- data.frame()

# Loop over regions
for(reg in unique(wi_data$region)){
  df_region <- wi_data %>% filter(region == reg) %>% select(all_of(pfas_vars))
  mat <- cor(df_region, use = "pairwise.complete.obs")
  
  # Convert to long format
  cor_long <- as.data.frame(as.table(mat))
  colnames(cor_long) <- c("PFAS_1","PFAS_2","Correlation")
  cor_long$Region <- reg
  
  # Remove self-correlations
  cor_long <- cor_long %>% filter(PFAS_1 != PFAS_2)
  
  # Bind to main table
  cor_all <- rbind(cor_all, cor_long)
}

# add PFAS group labels
top_cor_table <- cor_all %>%
  left_join(pfas_groups, by = c("PFAS_1" = "PFAS")) %>%
  setNames(c("PFAS_1","PFAS_2","Correlation","Region","Group_1")) %>%
  left_join(pfas_groups, by = c("PFAS_2" = "PFAS")) %>%
  setNames(c("PFAS_1","PFAS_2","Correlation","Region","Group_1","Group_2")) %>%
  arrange(Region, desc(Correlation))

# keep top 5 correlations per region
top_cor_table <- top_cor_table %>%
  group_by(Region) %>%
  slice_max(order_by = Correlation, n = 5) %>%
  ungroup()

# Print table
kable(top_cor_table, digits = 2, caption = "Top 5 correlations between different PFAS per region with group labels")

# Export to CSV
write.csv(top_cor_table, "Top_PFAS_Correlations_by_Region_no_self.csv", row.names = FALSE)

hist(data$PFNA)
