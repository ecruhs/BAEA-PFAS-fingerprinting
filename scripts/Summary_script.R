# ---- Libraries ----
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# set working directory
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/"
setwd(homewd)

# ---- Load data ----
data <- read.csv("data/Allyears_PFAS.csv", header = TRUE, stringsAsFactors = FALSE)

# ---- Summarize to just pfas data and counties
names(data)
pfas_data <- select(data, c(1, 7, 26:73,74)) #seven is year
names(pfas_data)

# --- sample sizes ----
sample.sizes <- ddply(pfas_data, .(region, year), summarise, 
                       nTot = length(unique(bird_id)))

sample.sizes <- ddply(data, .(DNA_sex), summarise, 
                      nTot = length(unique(bird_id)))

length(unique(data$nest_no)) #65


# ---- make data long format ----
str(pfas_data)
pfas_data$year <- as.factor(pfas_data$year)
pfas_data <- pfas_data %>%
  mutate(across(where(is.integer), as.numeric))

long.data <- melt(pfas_data, by="bird_id")
names(long.data)[names(long.data)=="variable"] <- "toxicant"
names(long.data)[names(long.data)=="value"] <- "value"


region.sum <- ddply(
  long.data,
  .(region, toxicant),
  summarise,
  summary = paste0(
    round(mean(value, na.rm = TRUE), 3), " (",
    round(min(value, na.rm = TRUE), 3), "–",
    round(max(value, na.rm = TRUE), 3), ")"
  )
)

region_wide <- region.sum %>%
  pivot_wider(
    names_from = toxicant,
    values_from = summary
  )

write.csv(region_wide,
          "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/output-data/PFAS_summary_by_region_simple.csv",
          row.names = FALSE)




region.sum <- ddply(
  long.data,
  .(region, toxicant),
  summarise,
  summary = paste0(
    round(mean(value, na.rm = TRUE), 3), " ± ",
    round(sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))), 3)
  )
)

region_wide <- region.sum %>%
  pivot_wider(
    names_from = toxicant,
    values_from = summary
  )

write.csv(
  region_wide,
  "/Users/emilyruhs/Desktop/GitHub_repos/PFAS_2025/output-data/PFAS_mean_error_by_region.csv",
  row.names = FALSE
)
