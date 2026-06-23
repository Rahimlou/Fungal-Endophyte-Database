###############################################################
# Global Distribution Map of Genus Records
# Robinson Projection
# Vector-output (PDF and SVG)
###############################################################

# Install packages if needed
# install.packages(c(
#   "tidyverse",
#   "sf",
#   "rnaturalearth",
#   "rnaturalearthdata",
#   "viridis",
#   "ggplot2"
# ))

setwd("/Users/srahimlou/OneDrive_The_Pennsylvania_State_University/Fungal_Diversity_Final/Final_Revision")

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggplot2)
library(svglite)
###############################################################
# 1. Read dataset
###############################################################

# Expected columns:
# Genus, Lat, Lon, Country, Red_Type

df <- read.csv(
  "FunEndo_Map_Data.csv",
  stringsAsFactors = FALSE
)

#names(df)[names(df) == "X...Genus"] <- "Genus"
###############################################################
# 2. Count records per country
###############################################################

country_counts <- df %>%
  group_by(Country) %>%
  summarise(
    Genus_Records = n(),
    Unique_Genera = n_distinct(Genus),
    .groups = "drop"
  )

# Removing the records with unknown countries

country_counts <- country_counts %>%
  filter(!is.na(Country))
###############################################################
# 3. Download world map
###############################################################

world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

###############################################################
# 4. Match country names
###############################################################
# Depending on your dataset, some country names may need
# adjustment before joining.
#
# Examples:
# "USA" -> "United States of America"
# "UK"  -> "United Kingdom"
#
# Uncomment and edit as needed:
#
# country_counts$Country <- recode(
#    country_counts$Country,
#    "USA" = "United States of America",
#    "UK"  = "United Kingdom"
# )

###############################################################
# 5. Join country counts to world map
###############################################################

world_map <- world %>%
  left_join(
    country_counts,
    by = c("name" = "Country")
  )

# Countries without records become zero
world_map$Genus_Records[
  is.na(world_map$Genus_Records)
] <- 0

world_map$Unique_Genera[
  is.na(world_map$Unique_Genera)
] <- 0

# Copy China's values to Taiwan

china_records <- world_map$Genus_Records[
  world_map$name == "China"
]

china_genera <- world_map$Unique_Genera[
  world_map$name == "China"
]

world_map$Genus_Records[
  world_map$name == "Taiwan"
] <- china_records

world_map$Unique_Genera[
  world_map$name == "Taiwan"
] <- china_genera
###############################################################
# 6. Robinson projection
###############################################################

robin_crs <- "+proj=robin +datum=WGS84"

###############################################################
# 7. Create map
###############################################################

p <- ggplot(world_map) +
  
  geom_sf(
    aes(fill = Genus_Records),
    color = "grey50",
    linewidth = 0.15
  ) +
  
  scale_fill_viridis(
    option = "plasma",
    trans = "sqrt",
    name = "Number of\nRecords",
    na.value = "grey95"
  ) +
  
  coord_sf(crs = robin_crs) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    )
  )

###############################################################
# 8. Display map
###############################################################

print(p)

###############################################################
# 9. Save publication-quality vector outputs
###############################################################

ggsave(
  filename = "Global_Genus_Map_Robinson.pdf",
  plot = p,
  width = 14,
  height = 8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "Global_Genus_Map_Robinson.svg",
  plot = p,
  width = 14,
  height = 8,
  units = "in"
)

###############################################################
# Use Unique_Genera instead of total records
###############################################################

p_unique <- ggplot(world_map) +
  
  geom_sf(
    aes(fill = Unique_Genera),
    color = "grey50",
    linewidth = 0.15
  ) +
  
  scale_fill_viridis(
    option = "magma",
    trans = "sqrt",
    name = "Unique\nGenera",
    na.value = "grey95"
  ) +
  
  coord_sf(crs = robin_crs) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_unique)

# Save second map
ggsave(
  "Global_Unique_Genera_Map_Robinson.pdf",
  p_unique,
  width = 14,
  height = 8,
  device = cairo_pdf
)

ggsave(
  filename = "Global_Unique_Genera_Map_Robinson.svg",
  plot = p_unique,
  width = 14,
  height = 8,
  units = "in"
)
