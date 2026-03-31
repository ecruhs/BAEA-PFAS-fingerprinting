# ==============================
# Libraries
# ==============================
library(dplyr)
library(ggplot2)
library(reshape2)
library(knitr)
library(purrr)

# ==============================
# Set working directory & load data
# ==============================
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

wi_data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

wi_data <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# ==============================
# PFAS variables
# ==============================
pfas_vars <- c(
  "FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",   # Precursors
  "PFBA","PFBS","PFHPA",                                # Short-chain
  "PFHXA","PFOA","PFNA","PFDA","PFDOA","PFTRDA","PFTEDA","PFUNDA",
  "PFHXS","PFOS","PFPES","PFHPS","PFDS"                 # Legacy
)

# Group labels
pfas_groups <- data.frame(
  PFAS = pfas_vars,
  Group = c(rep("Precursor",5),
            rep("Short-chain",3),
            rep("Legacy Acids",7),
            rep("Legacy Sulfonates",3),
            rep("Legacy Sulfonamides",3)),
  stringsAsFactors = FALSE
)

# Pretty names for plotting
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

# ==============================
# Heatmap correlations (Spearman)
# ==============================
df <- wi_data %>%
  select(region, all_of(pfas_vars))

cor_long <- df %>%
  group_by(region) %>%
  group_modify(~ {
    mat <- cor(.x[, pfas_vars], use = "pairwise.complete.obs", method = "spearman")
    mat[upper.tri(mat)] <- NA
    melt(mat) %>%
      as.data.frame() %>%
      mutate(region = unique(.x$region))
  }) %>%
  ungroup()

# Filter out NAs and unused PFAS
cor_long <- cor_long %>%
  filter(!is.na(value)) %>%
  filter(Var1 %in% names(vars_renamed) & Var2 %in% names(vars_renamed))

# Reorder axes
cor_long$Var1 <- factor(cor_long$Var1, levels = pfas_vars, labels = vars_renamed[pfas_vars])
cor_long$Var2 <- factor(cor_long$Var2, levels = pfas_vars, labels = vars_renamed[pfas_vars])

# Plot heatmap
ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "black") +
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
    strip.text = element_text(face = "bold", size = 10)
  ) +
  labs(fill = "Spearman ρ", x = "", y = "")

# ==============================
# Spearman correlation tests with p-values
# ==============================
library(dplyr)
library(purrr)

spearman_corr_test <- function(data, vars, reg) {
  # Generate all unique PFAS pairs
  pairs <- combn(vars, 2, simplify = FALSE)
  
  # Map over pairs safely
  results <- map_dfr(pairs, function(pair) {
    x <- data[[pair[1]]]
    y <- data[[pair[2]]]
    
    # Only compute if at least 3 non-missing pairs
    if(sum(complete.cases(x, y)) >= 3){
      test <- cor.test(x, y, method = "spearman", use = "pairwise.complete.obs")
      tibble(
        PFAS_1 = pair[1],
        PFAS_2 = pair[2],
        Correlation = unname(test$estimate),
        p_value = test$p.value
      )
    } else {
      tibble(
        PFAS_1 = pair[1],
        PFAS_2 = pair[2],
        Correlation = NA_real_,
        p_value = NA_real_
      )
    }
  })
  
  results <- results %>%
    mutate(Region = reg)
  
  return(results)
}

# Apply to all regions
cor_all_tests <- map_dfr(unique(wi_data$region), function(reg) {
  df_region <- wi_data %>%
    filter(region == reg) %>%
    select(all_of(pfas_vars))
  spearman_corr_test(df_region, pfas_vars, reg)
})

# Add PFAS group labels + significance stars
cor_all_tests <- cor_all_tests %>%
  # Join for PFAS_1 group
  left_join(setNames(pfas_groups$Group, pfas_groups$PFAS) %>% 
              tibble::enframe(name = "PFAS_1", value = "Group_1"),
            by = "PFAS_1") %>%
  # Join for PFAS_2 group
  left_join(setNames(pfas_groups$Group, pfas_groups$PFAS) %>% 
              tibble::enframe(name = "PFAS_2", value = "Group_2"),
            by = "PFAS_2") %>%
  # Add significance stars
  mutate(
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.1   ~ ".",
      TRUE            ~ ""
    )
  )

write.csv(cor_all_tests, 
          "output-data/PFAS_SpearmanCorrelations_withGroups_Significance.csv", 
          row.names = FALSE)


