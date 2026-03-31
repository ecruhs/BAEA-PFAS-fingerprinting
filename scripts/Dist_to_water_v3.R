# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(tidytext)  # for reorder_within

### --- Load data
homewd <- "/Users/gavindehnert/PFAS_2025/"
setwd(homewd)

df <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

# --- Function to calculate correlations ---
calc_corr <- function(data, vars, dist_var = "distance_water_m") {
  results <- lapply(vars, function(v) {
    tmp <- data[, c(v, dist_var)]
    tmp <- tmp[complete.cases(tmp), ]
    
    if (nrow(tmp) > 2) {
      test <- cor.test(tmp[[v]], tmp[[dist_var]], method = "spearman")
      data.frame(
        PFAS = v,
        rho = as.numeric(test$estimate),
        p_value = test$p.value,
        n = nrow(tmp)
      )
    } else {
      data.frame(PFAS = v, rho = NA, p_value = NA, n = nrow(tmp))
    }
  })
  do.call(rbind, results)
}

# --- PFAS variables of interest ---
pfas_vars <- c(
  "FTSA6_2","FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",   # Precursors
  "PFBA","PFBS","PFHPA",                                # Short-chain
  "PFHXA","PFOA","PFNA","PFDA","PFDOA","PFTRDA","PFTEDA","PFUNDA",
  "PFHXS","PFOS","PFPES","PFHPS","PFDS"                 # Legacy
)

# --- Renaming dictionary ---
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

# --- Correlations by region ---
region_corr <- df %>%
  group_by(region) %>%
  group_modify(~ calc_corr(.x, pfas_vars)) %>%
  ungroup()

region_corr <- region_corr %>%
  group_by(region) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))


# --- Apply renaming ---
region_corr <- region_corr %>%
  mutate(PFAS_clean = recode(PFAS, !!!vars_renamed),
         signif = ifelse(p_adj < 0.05, "p < 0.05", "ns"),
         PFAS_clean = reorder_within(PFAS_clean, abs(rho), region))

# --- Plot ---
p <- ggplot(region_corr, aes(x = PFAS_clean, y = rho, fill = signif)) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~ region, scales = "free_y") +
  scale_x_reordered() +
  scale_fill_manual(values = c("p < 0.05" = "red", "ns" = "grey80")) +
  labs(
    #title = "PFAS correlations with distance to nearest water, by region",
    x = "PFAS compound",
    y = "Spearman correlation (rho)",
    fill = "Significance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(p)


#### ---- to add in correlations across all of Wisconsin
# all_corr <- calc_corr(df, pfas_vars)
# 
# # --- Apply renaming and significance labeling ---
# all_corr <- all_corr %>%
#   mutate(
#     PFAS_clean = recode(PFAS, !!!vars_renamed),
#     signif = ifelse(p_value < 0.05, "p < 0.05", "ns"),
#     PFAS_clean = reorder(PFAS_clean, abs(rho))
#   )
# 
# p <- ggplot(all_corr, aes(x = PFAS_clean, y = rho, fill = signif)) +
#   geom_col(color = "black") +
#   coord_flip() +
#   scale_fill_manual(values = c("p < 0.05" = "red", "ns" = "grey80")) +
#   labs(
#     title = "",
#     x = "",
#     y = "",
#     fill = "Significance"
#   ) +
#   theme_minimal(base_size = 11) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     axis.title.y = element_text(face = "bold"),
#     axis.title.x = element_text(face = "bold"),
#     legend.position = "none"
#   )
# 
# print(p)


# --- Save outputs with clean PFAS names ---
region_corr_out <- region_corr %>%
  mutate(PFAS = recode(PFAS, !!!vars_renamed)) %>%
  select(region, PFAS, rho, p_value, p_adj, n, signif)

write.csv(region_corr_out, "output-data/PFAS_distance_correlations_byregion.csv", row.names = FALSE)



# Spearman’s rho (ρ) tells you both the strength and the direction of a monotonic relationship:
#   Sign (positive or negative):
#   Positive ρ → as distance to water increases, PFAS levels also tend to increase.
# Negative ρ → as distance to water increases, PFAS levels tend to decrease.
# Magnitude (absolute value):
#   ρ close to 1 or –1 → strong monotonic relationship.
# ρ close to 0 → little or no monotonic relationship.

### test this with MWR
data_MWR <- subset(df, region == "MWR")
ggplot(data = data_MWR, aes(x = distance_water_m, y = FOSA)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")

#### test this with SSS
data_SSS <- subset(df, region == "SSS")
ggplot(data = data_SSS, aes(x = distance_water_m, y = TOTAL_PFAS)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")

ggplot(data = data_SSS, aes(x = distance_water_m, y = PFUNDA)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")
