# Doing db-RDA-based variance partitioning

# Loading required libraries

library(tidyr)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(phyloseq)
library(adespatial)
library(geodata)
library(elevatr)
library(sf)
install.packages("raster")
library(raster)
library(ggplot2)
library(ggrepel)
library(viridis)

# Setting up the working directory
setwd("/Users/srahimlou/OneDrive_The_Pennsylvania_State_University/FunEndo_R_Codes/db_RDA")

df <- read.csv("FunEndo_DB.csv", header = TRUE, fileEncoding = "latin1")

# Edit Hosts in df
df <- df %>%
  mutate(Host = recode(
    Host,
    "Agrostis hiemalis"            = "Agrostis hyemalis",
    "Sequoia sempervire"           = "Sequoia sempervirens",
    "Serjania laruotteana"         = "Serjania lauretana",
    "Quercus phillyraeoides"       = "Quercus phillyreoides",
    "Quercus semicarpifolia"       = "Quercus semecarpifolia",
    "Dendrobium offcinale"         = "Dendrobium officinale",
    "Dioscorea opposite"           = "Dioscorea opposita",
    "Ixora¬†chinensis"             = "Ixora chinensis",
    "Punctelia borrerii"           = "Punctelia borreri",
    "Tylophora ovate"              = "Tylophora ovata",
    "Cryptolepis buchanani"        = "Cryptolepis buchananii",
    "Paris polyphyllaSmith"        = "Paris polyphylla",
    "Coriospermum declinatum"      = "Cardiospermum declinatum",
    "Coriospermum tibeticum"       = "Cardiospermum tibeticum",
    "Combretum latifolium¬†"       = "Combretum latifolium"
  ))

# Remove rows where hosts are wrong
hosts_to_remove <- c(
  "Betulaceous",
  "Euphorbia",
  "Pyrenacantha",
  "Umbelliferae",
  "Viridiplantae",
  "Cladonia coniocraea",
  "lichens",
  "Melanelia sorediata",
  "Parmelia",
  "Punctelia borreri",
  "Xanthoria mandschurica"
)

df <- df %>%
  filter(!Host %in% hosts_to_remove)

write.csv(df, "FunEndo_DB_Host_Edited.csv")

# counting the number of filled cells in Lat and Lon columns
sum(!is.na(df$Lat))
# filled = 24938

sum(!is.na(df$Lon))
# filled = 24982

df$Lat <- as.numeric(gsub("[^0-9.-]", "", df$Lat))
df$Lon <- as.numeric(gsub("[^0-9.-]", "", df$Lon))

# Replace "NA" with actual NA
df[df == "NA"] <- ""

# Keep kingdoms that are Fungi
df <- subset(df, Kingdom == "Fungi")

# Remove rows where Lat are missing
df <- df %>%
  filter(!(is.na(Lat)))

# Remove rows where Lon are missing
df <- df %>%
  filter(!(is.na(Lon)))

######################### Adding climate zone variable #########################
classify_climate <- function(lat) {
  if (lat >= -23.5 && lat <= 23.5) {
    return("Tropical")
  } else if ((lat > 23.5 && lat <= 66.5) || (lat < -23.5 && lat >= -66.5)) {
    return("Temperate")
  } else {
    return("Boreal")
  }
}

df <- df %>% 
  mutate(ClimateZone = sapply(Lat, classify_climate))
############################ GET MAT and MAP ###################################

# Extract MAT and MAP
df_coords <- data.frame(
  Lon = df$Lon,
  Lat = df$Lat
)

# Specify the directory containing the extracted BIOCLIM files
bioclim_dir <- "/Users/srahimlou/OneDrive_The_Pennsylvania_State_University/FunEndo_R_Codes/wc2.1_2.5m_bio"

# List all BIOCLIM raster files in the directory
bioclim_files <- list.files(bioclim_dir, pattern = "bio.*tif", full.names = TRUE)

# Read BIOCLIM raster files into a list of raster layers
bioclim_rasters <- lapply(bioclim_files, raster)

# Define the Coordinate Reference System (CRS) for longitude/latitude (WGS84)
crs <- CRS("+proj=longlat +datum=WGS84")

# Create a SpatialPoints object from longitude and latitude
points <- SpatialPoints(coords = df_coords[, c("Lon", "Lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract climate values (MAT and MAP) at the coordinates
climate_values <- extract(stack(bioclim_rasters), points)
climate_values <- as.data.frame(climate_values)

# Extract MAT (BIO1) and MAP (BIO12) values
MAT <- climate_values$wc2.1_2.5m_bio_1   # Mean Annual Temperature (BIO1)
MAP <- climate_values$wc2.1_2.5m_bio_12  # Mean Annual Precipitation (BIO12)

# Add MAT and MAP values to the original dataframe
df$MAT <- MAT
df$MAP <- MAP

# Remove rows where MAT are missing
df <- df %>%
  filter(!(is.na(MAT)))

# Remove rows where MAP are missing
df <- df %>%
  filter(!(is.na(MAP)))
################################ Fixing Organ Data #############################
# Remove the rows where Organ is empty
df <- df[!is.na(df$Organ),]

# Correct the organ typos
df <- subset(df, Organ != "")
df$Organ <- gsub("Leaf", "leaf", df$Organ)
df$Organ <- gsub("Flower", "flower", df$Organ)
df$Organ <- gsub("Others", "other", df$Organ)
df$Organ <- gsub("Stem", "stem", df$Organ)
df$Organ <- gsub("Root", "root", df$Organ)
df$Organ <- gsub("Rhizome", "rhizome", df$Organ)
df$Organ <- gsub("Wood", "wood", df$Organ)
df$Organ <- gsub(" ", "", df$Organ)
df$Organ <- gsub("Stolon", "stolon", df$Organ)
df$Organ <- gsub("Needle", "needle", df$Organ)
df$Organ <- gsub("Bark", "bark", df$Organ)
df$Organ <- gsub("Trunk", "trunk", df$Organ)
df$Organ <- gsub("Gall", "gall", df$Organ)
df$Organ <- gsub("Seed", "seed", df$Organ)
df$Organ <- gsub("Fruit", "fruit", df$Organ)
df$Organ <- gsub("Shoot", "shoot", df$Organ)
df$Organ <- gsub("Tuber", "tuber", df$Organ)
df$Organ <- gsub("Branch\\(Twig\\)", "twig", df$Organ)
df <- subset(df, Organ != "Photosynthetictissue")

########################################Getting Elevation in Batch ##########################################
# Install and load required packages
if (!requireNamespace("geodata", quietly = TRUE)) {
  install.packages("geodata")
}
if (!requireNamespace("elevatr", quietly = TRUE)) {
  install.packages("elevatr")
}
if (!requireNamespace("sf", quietly = TRUE)) {
  install.packages("sf")
}

# ---- Batch elevation extraction ----

# Create sf object (WGS84)
locations <- st_as_sf(
  df,
  coords = c("Lon", "Lat"),
  crs = 4326,
  remove = FALSE
)

# Get elevations for ALL points at once
elevations <- get_elev_point(
  locations = locations,
  src = "aws"
)

elevations$elev_units <- NULL
elevations$geometry <- NULL
df <- elevations

# Remove rows where the elevations are missing
df <- df %>%
  filter(!(is.na(elevation)))
######################################## Getting Ecoregions #################################################
spr_df <- read.csv("Joined_Layers_Super_Ecoregions.csv", header = T)

# checking how many unique rows are there in spr_df if we IGNORE the columns spr_rgn and area?
spr_df %>%
  distinct(across(-c(spr_rgn, area))) %>%
  nrow()
# 16197 unique rows identified

# Checking whether df has any duplicate rows across all columns
nrow(distinct(df))
# 16084 records identified

# keep the unique rows in spr_df ignoring the spr_rgn and area columns. Prefer rows where spr_rgn is NOT NA
spr_df_unique <- spr_df %>%
  arrange(is.na(spr_rgn)) %>%   # non-NA first
  distinct(across(-c(spr_rgn, area)), .keep_all = TRUE)

# Remove rows where the spr_rgn are missing
spr_df_unique <- spr_df_unique %>%
  filter(!(is.na(spr_rgn)))

df_final <- spr_df_unique

# 15375 records remained
############################### Adding Growth Forms ############################
growth_forms <- read.csv("host_growth_form.csv", header = T)

df_final <- df_final %>%
  left_join(growth_forms, by = c("Host" = "Host"))

# Remove rows where the Genus are missing
df_final <- df_final %>%
  filter(!(is.na(Genus)))
# 9054 records remained

# Remove rows where the Host are missing
#df_final <- df_final %>%
#  filter(!(is.na(Host)))
# 8775 records remained
################################################################################
# Subset data frame
subset_df <- df_final[, c("Reference", "Genus", "Species", "Host", "Study_Type", "Lat", 
                           "Lon", "MAT", "MAP", "Organ", "ClimateZone", "spr_rgn", "elevation", 
                          "Growth_Form")]

sum(subset_df$Organ == "leaf")
# leaf count = 3654
sum(subset_df$Organ == "root")
# root count = 3498
sum(subset_df$Organ == "fruit")
# fruit count = 487
sum(subset_df$Organ == "flower")
# flower count = 13
sum(subset_df$Organ == "stem")
# stem count = 0
sum(subset_df$Organ == "rhizome")
# rhizome count = 53
sum(subset_df$Organ == "wood")
# wood count = 0
sum(subset_df$Organ == "stolon")
# stolon count = 1
sum(subset_df$Organ == "needle")
# needle count = 0
sum(subset_df$Organ == "bark")
# bark count = 3
sum(subset_df$Organ == "trunk")
# trunk count = 39
sum(subset_df$Organ == "gall")
# gall count = 0
sum(subset_df$Organ == "seed")
# seed count = 79
sum(subset_df$Organ == "shoot")
# shoot count = 0
sum(subset_df$Organ == "tuber")
# tuber count = 0
sum(subset_df$Organ == "twig")
# twig count = 600
sum(subset_df$Organ == "other")
# other count = 348

# keep the rows where the Organ is only leaf or root
subset_df <- subset(subset_df, Organ %in% c("leaf", "root"))

# 7152 records remained

# Creating a new column Sample that assigns a unique ID like sample1, sample2, … 
# to each unique combination of Reference, Lat, Lon, and Organ.
subset_df <- subset_df %>%
  mutate(Sample = paste0("sample", dense_rank(paste(Reference, Lat, Lon, Organ)))) 

# Count the number of unique samples created
length(unique(subset_df$Sample))

# 627 unique samples were created. 

# Create an OTU table
# Count occurrences of each Genus in each Sample
df_counts <- subset_df %>%
  group_by(Sample, Genus) %>%
  summarise(Count = n(), .groups = "drop")

# Pivot the table to make Genus as row names and Sample as columns
otu_table <- df_counts %>%
  pivot_wider(names_from = Genus, values_from = Count, values_fill = list(Count = 0))

# Convert Sample to row names
otu_table <- otu_table %>%
  column_to_rownames("Sample")

# Count the number of unique Genera
length(unique(df_counts$Genus))

# Genus = 534
# Samples = 627

# Write the OTU table
write.csv(otu_table, "OTU_table.csv")

# NOTE: Unload the raster library before running it.

# Remove columns (Genera) where the sum is less than or equal to 3
otu_table <- otu_table %>%
  # Calculate the sum for each column and filter columns with sum > 3
  select(where(~sum(.) > 3))

# The number of fungal genera remained in OTU table
length(unique(colnames(otu_table)))
# 190 genera remained

# Remove rows where the sum is 0
otu_table <- otu_table %>%
  # Filter out rows where the sum of values in the row is 0
  filter(rowSums(.) > 0)

# Samples = 594
# Genera = 190

write.csv(otu_table, "final_otu_table.csv")

# OTU table ("final_otu_table.csv") has 594 samples and 190 genera.

# Set the taxa names (column names) in the OTU table
taxa_names(otu_table) <- colnames(otu_table)

# Create a metadata frame for your climate zones
metadata <- data.frame(
  Sample = subset_df$Sample,
  Lat = subset_df$Lat,
  Lon = subset_df$Lon,
  Elevation = subset_df$elevation,
  MAT = subset_df$MAT,
  MAP = subset_df$MAP,
  spr_rgn = subset_df$spr_rgn,
  ClimateZone = subset_df$ClimateZone,
  Organ = subset_df$Organ)

# keeping the unique rows which are the unique samples.
metadata <- metadata %>%
  distinct()

# Filter metadata so it keeps only the samples that exist in otu_table.
metadata <- metadata %>%
  filter(Sample %in% rownames(otu_table))

# Calculate how many unique samples exist.
length(unique(metadata$Sample))
# 594 unique samples exist which is matching the number of samples.


######################## Creating a phylogeny tree ##########################
fungal_genera <- unique(colnames(otu_table))

# Subset the data frame to include only rows where Genus is in fungal_genera
filtered_df <- df[df$Genus %in% fungal_genera, 
                  c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus")]

filtered_df$Genus <- trimws(filtered_df$Genus)

# keeping the unique rows
filtered_df <- filtered_df %>%
  distinct()

length(unique(filtered_df$Genus))

# Matched to the number of genera in OTU table which is 190.

write.csv(filtered_df, "Genera_for_phylogeny_construction.csv")

write.table(
  filtered_df,               # your data frame
  file = "Genera_for_phylogeny_construction.txt",  # output file name
  sep = "\t",           # tab delimiter
  row.names = FALSE,    # do not write row names
  col.names = FALSE,
  quote = FALSE         # do not quote strings
)

phylo_tree <- read.tree("endophytes_tree.nwk")

# Replace underscores with spaces in the tip labels (not used for genus-level analysis)
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)

###################### Calculating UniFrac Distance ############################
# Keep only the taxa that are in OTU table

# Identify the taxa in the tree
taxa_in_tree <- phylo_tree$tip.label

# Keep only the taxa that are in your list of fungal genera
taxa_to_keep <- taxa_in_tree[taxa_in_tree %in% fungal_genera]

# Subset the tree to keep only the selected taxa
phylo_tree <- drop.tip(phylo_tree, setdiff(taxa_in_tree, taxa_to_keep))
length(phylo_tree$tip.label)

setdiff(metadata$Sample, rownames(otu_table))  # in metadata but not in OTU
setdiff(rownames(otu_table), metadata$Sample)  # in OTU but not in metadata

otu_table_sample_names <- rownames(otu_table)
write.csv(otu_table_sample_names, "otu_table_sample_names.csv")

write.csv(metadata, "metadata.csv")

tip_lables <- phylo_tree$tip.label
write.csv(tip_lables, "tip_lables.csv")

otu_table <- t(otu_table)
otu_table <- as.data.frame(otu_table)

# Convert the first column named "Sample" to row names
rownames(metadata) <- metadata$Sample
metadata$Sample <- NULL

OTU_taxa <- rownames(otu_table)
write.csv(OTU_taxa, "OTU_taxa.csv")

# Create a phyloseq object
ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), sample_data(metadata), phy_tree(phylo_tree))

# Calculate unweighted UniFrac distance
unweighted_unifrac <- phyloseq::distance(ps, method = "unifrac", weighted = FALSE)


########################## Doing db-RDA analysis ###############################
# Environmental matrix
env <- data.frame(
  MAT = sample_data(ps)$MAT,
  MAP = sample_data(ps)$MAP,
  Elevation = sample_data(ps)$Elevation,
  ClimateZone = sample_data(ps)$ClimateZone
)

env[, c("MAT", "MAP", "Elevation")] <- scale(env[, c("MAT", "MAP", "Elevation")])

# Host matrix
host <- data.frame(
  Organ = sample_data(ps)$Organ
)
host$Organ <- factor(host$Organ)

# Spatial matrix (important!)
# Simple lat/long is not enough — you need spatial eigenvectors:

## ================================
## Spatial db-RDA with UniFrac
## ================================

## 1. Extract coordinates
coords <- data.frame(
  lon = sample_data(ps)$Lon,
  lat = sample_data(ps)$Lat
)

## 2. Build spatial eigenvectors (MEM / PCNM)
pcnm <- dbmem(coords)

## 3. Convert MEMs to data frame and name them
pcnm_df <- as.data.frame(pcnm)
colnames(pcnm_df) <- paste0("MEM", seq_len(ncol(pcnm_df)))

## 4. Prepare metadata
meta <- as.data.frame(sample_data(ps))

## Combine metadata + spatial predictors
meta_pcnm <- cbind(meta, pcnm_df)


## 6. Forward selection of MEMs using UniFrac-based db-RDA
space_null <- capscale(
  unweighted_unifrac ~ 1,
  data = meta_pcnm
)

space_global <- capscale(
  unweighted_unifrac ~ .,
  data = pcnm_df
)

space_sel <- ordiR2step(
  space_null,
  scope = formula(space_global),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)

## View selected spatial eigenvectors
space_sel$anova

## 7. FINAL db-RDA model
## (Replace MEMs below with those selected above)
full_model <- capscale(
  unweighted_unifrac ~ MAT + MAP + Elevation + ClimateZone +
    Organ +
    MEM1 + MEM2 + MEM3 + MEM8 + MEM5 +
    MEM15 + MEM12 + MEM6 + MEM9 + MEM51,
  data = meta_pcnm
)

vif.cca(full_model) # check collinearity: If any MEM has VIF > 10, drop it.

## 8. Model significance
anova(full_model, permutations = 999)

## Optional: test individual terms
anova(full_model, by = "terms", permutations = 999)

################################### Plotting Data ##############################
# Create a table of growth forms per sample
gf_summary <- subset_df %>%
  group_by(Sample) %>%
  summarise(
    Growth_Form = {
      gf <- unique(na.omit(Growth_Form))
      if (length(gf) == 0) {
        NA_character_
      } else if (length(gf) == 1) {
        gf
      } else {
        "mixed"
      }
    },
    .groups = "drop"
  )


#### 1. extract site scores (constrained axes) ----
# Use vegan::scores (works for capscale objects)
site_scores <- as.data.frame(scores(full_model, display = "sites", choices = c(1,2)))
colnames(site_scores) <- c("CAP1", "CAP2")

# add metadata for plotting
site_scores$Sample <- rownames(site_scores)
site_scores <- cbind(site_scores, meta_pcnm[rownames(site_scores), , drop = FALSE])

# Add growth form
site_scores <- site_scores %>%
  left_join(gf_summary, by = "Sample")

# keep row names exactly as they are
rownames(site_scores) <- site_scores$Sample

#### 2. get percent variance explained for axis labels ----
# use constrained canonical axes importance (constrained components)
# summary(full_model)$concont may be NULL for some capscale objects; fall back to eigenvalues if necessary
cont_summary <- tryCatch(summary(full_model)$concont$importance, error = function(e) NULL)

if (!is.null(cont_summary)) {
  var1 <- cont_summary[2,1] * 100
  var2 <- cont_summary[2,2] * 100
} else {
  # fallback: use eigenvalues from capscale object
  ev <- eigenvals(full_model, model = "constrained")
  var1 <- ev[1] / sum(ev) * 100
  var2 <- ev[2] / sum(ev) * 100
}
xlab <- paste0("db-RDA1 (", round(var1, 1), "%)")
ylab <- paste0("db-RDA2 (", round(var2, 1), "%)")

#### 3. compute centroids for categorical variable(s) e.g. Organ ----
centroids <- site_scores %>%
  group_by(Organ) %>%
  summarise(CAP1 = mean(CAP1, na.rm = TRUE),
            CAP2 = mean(CAP2, na.rm = TRUE))

#### 4. envfit: fit environmental variables (numeric) to ordination ----
# choose numeric environmental variables to fit (change as needed)
env_vars <- meta_pcnm %>% dplyr::select(MAT, MAP, Elevation)

envfit_res <- envfit(full_model, env_vars, permutations = 999)

# extract arrow coordinates and p-values
env_arrows <- as.data.frame(scores(envfit_res, display = "vectors"))  # columns CAP1, CAP2
env_arrows$var <- rownames(env_arrows)
env_arrows$pval <- envfit_res$vectors$pvals[rownames(env_arrows)]

# keep only significant vectors (p < 0.05), optional — change threshold if needed
sig_arrows <- env_arrows %>% filter(pval <= 0.05)
if (nrow(sig_arrows) == 0) {
  sig_arrows <- env_arrows  # fallback: if none significant, include all
}

#### 5. scale arrows so they fit nicely on the plot ----
# compute scaling factor: make longest arrow ~ 50% of range of sites
site_range_x <- max(site_scores$CAP1, na.rm = TRUE) - min(site_scores$CAP1, na.rm = TRUE)
site_range_y <- max(site_scores$CAP2, na.rm = TRUE) - min(site_scores$CAP2, na.rm = TRUE)
site_range_max <- max(site_range_x, site_range_y)

arrow_max <- max(abs(sig_arrows$CAP1), abs(sig_arrows$CAP2))
# safety check
if(arrow_max == 0) arrow_max <- 1

# multiplier: adjust 0.6 to make arrows longer/shorter
mul <- (site_range_max * 0.6) / arrow_max

sig_arrows <- sig_arrows %>%
  mutate(arrow_x = CAP1 * mul,
         arrow_y = CAP2 * mul)

#### 6. build ggplot ----
p <- ggplot(site_scores, aes(x = CAP1, y = CAP2)) +
  geom_point(aes(color = Organ), size = 2.2, alpha = 0.8) +
  geom_point(data = centroids, aes(x = CAP1, y = CAP2),
             color = "black", shape = 4, size = 4, stroke = 1.2) +
  geom_text_repel(data = centroids,
                  aes(x = CAP1, y = CAP2, label = Organ),
                  color = "black", fontface = "bold", size = 3.5) +
  geom_segment(data = sig_arrows,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "black", size = 0.8) +
  geom_text_repel(data = sig_arrows,
                  aes(x = arrow_x, y = arrow_y, label = var),
                  color = "black", size = 3.5) +
  labs(
    x = xlab,
    y = ylab,
    color = "Host organ",   # 🔹 legend title
    title = "db-RDA (unweighted UniFrac) of fungal endophyte communities"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 5, alpha = 1)  # 🔹 bigger legend points
    )
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  coord_fixed()

print(p)
################### Alternative plot with growth form shapes ###################
site_scores <- site_scores %>%
  mutate(
    Growth_Form_plot = ifelse(
      Growth_Form %in% c(
        "Grass/Graminoid",
        "Herb/Forb",
        "Shrub",
        "Vine/Liana",
        "Tree"
      ),
      Growth_Form,
      "Other"
    )
  )

site_scores$Growth_Form_plot <- factor(site_scores$Growth_Form_plot)
p <- ggplot(site_scores, aes(x = CAP1, y = CAP2)) +
  geom_point(
    aes(color = Organ, shape = Growth_Form_plot),
    size = 2.2, alpha = 0.8
  ) +
  
  geom_point(data = centroids,
             aes(x = CAP1, y = CAP2),
             color = "black", shape = 4, size = 4, stroke = 1.2) +
  
  geom_text_repel(data = centroids,
                  aes(x = CAP1, y = CAP2, label = Organ),
                  color = "black", fontface = "bold", size = 3.5) +
  
  geom_segment(data = sig_arrows,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "black", size = 0.8) +
  
  geom_text_repel(data = sig_arrows,
                  aes(x = arrow_x, y = arrow_y, label = var),
                  color = "black", size = 3.5) +
  
  scale_shape_manual(
    values = c(
      "Grass/Graminoid" = 15,  # square
      "Herb/Forb"       = 16,  # circle
      "Shrub"           = 17,  # triangle
      "Vine/Liana"      = 18,  # diamond
      "Tree"            = 8,   # star
      "Other"           = 1    # open circle
    )
  ) +
  
  labs(
    x = xlab,
    y = ylab,
    color = "Host organ",
    shape = "Growth form",
    title = "db-RDA (unweighted UniFrac) of fungal endophyte communities"
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 4, alpha = 1))
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  coord_fixed()


print(p)
################################### Plotting Data ##############################

# Make Sample a column
coords <- meta_pcnm %>%
  tibble::rownames_to_column(var = "Sample") %>%
  select(Sample, Lon, Lat, MEM1:MEM10)  # MEMs selected

# Convert to long format
coords_long <- coords %>%
  pivot_longer(cols = starts_with("MEM"), names_to = "MEM", values_to = "Value")

# Convert to sf
coords_sf <- st_as_sf(coords_long, coords = c("Lon", "Lat"), crs = 4326)

# Plot with fixed scales
ggplot(coords_sf, aes(geometry = geometry, color = Value)) +
  geom_sf(size = 1) +  # smaller points
  facet_wrap(~MEM, ncol = 2) +  # force 2 columns
  scale_color_viridis_c(option = "C") +
  theme_minimal() +
  labs(color = "MEM value", title = "Spatial patterns of selected MEMs") +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )
############################ Test Growth Form ##################################
growth_model <- capscale(
  unweighted_unifrac ~ Growth_Form_plot +
    MAT + MAP + Elevation + ClimateZone +
    Organ +
    MEM1 + MEM2 + MEM3 + MEM8 + MEM5 +
    MEM15 + MEM12 + MEM6 + MEM9 + MEM51,
  data = site_scores
)

anova(growth_model, by = "terms", permutations = 999)









