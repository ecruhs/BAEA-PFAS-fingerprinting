library(ggplot2)
library(dplyr)
library(tidyr)  # for pivot_longer

homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load PFAS data ----
wi_data <- read.csv(
  file = paste0(homewd, "data/Allyears_PFAS.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# --- Define PFAS variables ---
vars <- c("FTSA8_2","FOSA","N_ETFOSAA","N_MEFOSAA",
          "PFBA","PFDA","PFDOA","PFDS","PFHPA","PFHPS",
          "PFHXA","PFHXS","PFNA","PFNS","PFOA","PFOS",
          "PFTEDA","PFTRDA","PFUNDA", "TOTAL_PFAS")

# Manual renaming for facet labels
vars_renamed <- c(
  "FOSA"       = "FOSA",
  "FTSA8_2"    = "8:2 FTS",
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
  "TOTAL_PFAS" = "Total PFAS"
)

# --- Colorblind-friendly palette ---
cb_colors <- c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#D55E00",  # vermillion/red-orange
  "#0072B2",  # blue
  "#CC79A7"   # reddish purple
)

# --- Ensure region is a factor with desired order ---
wi_data2 <- wi_data %>%
  mutate(region = factor(region, levels = c("AIS", "SSS", "FRGB", "GBLM", 
                                            "LWR", "MWR", "UPWR")))

# --- Pivot longer and plot ---
boxplot <- 
wi_data2 %>%
  pivot_longer(
    all_of(vars),        # use original column names here
    names_to = "PFAS",
    values_to = "conc"
  ) %>%
  ggplot(aes(x = region, y = conc, fill = region)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~PFAS, scales = "free_y", labeller = labeller(PFAS = vars_renamed)) +
  scale_fill_manual(values = cb_colors) +  # fill, not color
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) + ylab("Concentration\n") + xlab("")

boxplot
