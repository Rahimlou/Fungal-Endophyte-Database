# Performing Unifrac test
library(tidyr)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(phyloseq)

setwd("E:/Fungal_Endophyte_Database/Final_Dataset/Organ")

# Preparing the data for PERMANOVA
df <- read.csv("FunEndo_Final_Edit.csv", header = TRUE, fileEncoding = "Windows-1252")

df$Lat <- as.numeric(df$Lat)
df$Lon <- as.numeric(df$Lon)

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

################################################# Adding climate zone variable ####################################################
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

####################################################### GET MAT and MAP ###########################################################
install.packages("raster")
library(raster)

# Extract MAT and MAP
df_coords <- data.frame(
  Lon = df$Lon,
  Lat = df$Lat
)

# Specify the directory containing the extracted BIOCLIM files
bioclim_dir <- "E:/Fungal_Endophyte_Database/Final_Dataset/wc2.1_2.5m_bio/"

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

############################################################ Getting Organ Data ####################################################

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

########################################################## Getting Ecoregions and Elevation ########################################
# Subset data frame
subset_df <- df[, c("Reference","Species", "Lat", "Lon", "MAT", "MAP", "Organ", "ClimateZone")]

elev_ecoregions <- read.csv("df_ecoregions.csv", header = TRUE)

# Select specific columns
elev_ecoregions <- elev_ecoregions[, c("Lat", "Lon", "Elevation", "spr_rgn")]

elev_ecoregions_unique <- elev_ecoregions[!duplicated(elev_ecoregions[, c("Lat", "Lon")]), ]

# Perform the join to get spr_rgn based on matching Lat and Lon
df_joined <- merge(subset_df, elev_ecoregions_unique, by = c("Lat", "Lon"), all.x = TRUE)

# Remove rows with any NA values using complete.cases
df_joined <- df_joined[complete.cases(df_joined), ]

###################################################################################################################################

# Subset data frame
subset_df <- df_joined[, c("Reference", "Species", "Lat", "Lon", "MAT", "MAP", "Organ", "ClimateZone", "spr_rgn", "Elevation")]

sum(subset_df$Organ == "leaf")
# leaf count = 3802
sum(subset_df$Organ == "root")
# root count = 3544
sum(subset_df$Organ == "fruit")
# fruit count = 502
sum(subset_df$Organ == "flower")
# flower count = 13
sum(subset_df$Organ == "stem")
# stem count = 0
sum(subset_df$Organ == "rhizome")
# rhizome count = 54
sum(subset_df$Organ == "wood")
# wood count = 1
sum(subset_df$Organ == "stolon")
# stolon count = 1
sum(subset_df$Organ == "needle")
# needle count = 0
sum(subset_df$Organ == "bark")
# bark count = 3
sum(subset_df$Organ == "trunk")
# trunk count = 40
sum(subset_df$Organ == "gall")
# gall count = 0
sum(subset_df$Organ == "seed")
# seed count = 79
sum(subset_df$Organ == "shoot")
# shoot count = 0
sum(subset_df$Organ == "tuber")
# tuber count = 0
sum(subset_df$Organ == "twig")
# twig count = 645
sum(subset_df$Organ == "other")
# other count = 405

###### keep only the rows where Organ is leaf or root#######
subset_df <- subset(subset_df, Organ %in% c("leaf", "root", "fruit", "twig"))


# Create 'Sample' column using dplyr
subset_df <- subset_df %>%
  # Create a unique identifier for each combination of 'spr_rgn' and 'Organ'
  mutate(Sample = paste0("sample", dense_rank(paste(Reference, Lat, Lon, Organ)))) 

# Create an OTU table
# Count occurrences of each Genus in each Sample
df_counts <- subset_df %>%
  group_by(Sample, Species) %>%
  summarise(Count = n(), .groups = "drop")

# Pivot the table to make Genus as row names and Sample as columns
otu_table <- df_counts %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = list(Count = 0))

# Convert Sample to row names
otu_table <- otu_table %>%
  column_to_rownames("Sample")

# 836 Samples
# 611 Genus

# 592 Samples
# 972 Species


# Calculate the sum for each column (genus)
sum_row <- colSums(otu_table)

# Create a data frame for the sum row with the same column names
sum_row_df <- as.data.frame(t(sum_row))  # Transpose to make it a row

# Add a name to the row
rownames(sum_row_df) <- "sum"

# Append the sum row to the original OTU table
otu_table_with_sum <- rbind(otu_table, sum_row_df)

#write.csv(otu_table_with_sum, "otu_table_with_sum.csv")

# Remove columns where the sum is less than or equal to 3
otu_table <- otu_table %>%
  # Calculate the sum for each column and filter columns with sum > 3
  select(where(~sum(.) > 3))

# Remove rows where the sum is 0
otu_table <- otu_table %>%
  # Filter out rows where the sum of values in the row is 0
  filter(rowSums(.) > 0)

# Set the taxa names (column names) in the OTU table
taxa_names(otu_table) <- colnames(otu_table)

# Create a metadata frame for your climate zones
metadata <- data.frame(
  Sample = subset_df$Sample,
  Lat = subset_df$Lat,
  Lon = subset_df$Lon,
  Elevation = subset_df$Elevation,
  MAT = subset_df$MAT,
  MAP = subset_df$MAP,
  spr_rgn = subset_df$spr_rgn,
  ClimateZone = subset_df$ClimateZone,
  Organ = subset_df$Organ)

# Calculate the mean Lat for each sample
metadata <- metadata %>%
  group_by(Sample) %>%
  mutate(mean_lat = mean(Lat, na.rm = TRUE)) %>%
  ungroup()

# Calculate the mean Lon for each sample
metadata <- metadata %>%
  group_by(Sample) %>%
  mutate(mean_lon = mean(Lon, na.rm = TRUE)) %>%
  ungroup()

# Calculate the mean Elevation for each sample
metadata <- metadata %>%
  group_by(Sample) %>%
  mutate(mean_elev = mean(Elevation, na.rm = TRUE)) %>%
  ungroup()

# Calculate the average MAT for each sample
metadata <- metadata %>%
  group_by(Sample) %>%
  mutate(avg_MAT = mean(MAT, na.rm = TRUE)) %>%
  ungroup()

# Calculate the average MAP for each sample
metadata <- metadata %>%
  group_by(Sample) %>%
  mutate(avg_MAP = mean(MAP, na.rm = TRUE)) %>%
  ungroup()

# Remove duplicate samples
metadata <- metadata %>%
  distinct(Sample, .keep_all = TRUE)

metadata <- metadata %>%
  column_to_rownames("Sample")

# Just keep the Organ
#metadata <- subset(metadata, select = "Organ")

fungal_genera <- unique(colnames(otu_table))

# Now, let's load the phylogenetic tree
phylo_tree <- read.tree("Endophytes_Tree_ClimateZone.nwk")


######################## Creating a tree for endophytic species #############################
fungal_species <- unique(colnames(otu_table))

# Subset the data frame to include only rows where Species is in fungal_species
filtered_df <- df[df$Species %in% fungal_species, c("Kingdom", "Phylum", "Order", "Class", "Family", "Genus", "Species")]

filtered_df$Species <- trimws(filtered_df$Species)

# Remove rows with duplicate species where all other columns are the same
unique_df <- filtered_df[!duplicated(filtered_df[c("Species", "Kingdom", "Phylum", "Order", "Class", "Family", "Genus")]), ]

write.csv(unique_df, "unique_species.csv")

phylo_tree <- read.tree("Endophytes_Taxonomy_Species.nwk")

# Replace underscores with spaces in the tip labels
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)


####################################################################################################################################
# Keep only the taxa that are in OTU table

# Identify the taxa in the tree
taxa_in_tree <- phylo_tree$tip.label

# Keep only the taxa that are in your list of fungal genera
taxa_to_keep <- taxa_in_tree[taxa_in_tree %in% fungal_genera]

# Subset the tree to keep only the selected taxa
phylo_tree <- drop.tip(phylo_tree, setdiff(taxa_in_tree, taxa_to_keep))
length(phylo_tree$tip.label)

# Create a phyloseq object
ps <- phyloseq(otu_table(otu_table, taxa_are_rows = FALSE), sample_data(metadata), phy_tree(phylo_tree))

# Calculate unweighted UniFrac distance
unweighted_unifrac <- distance(ps, method = "unifrac", weighted = FALSE)

sample_data_df <- data.frame(sample_data(ps))

# PERMANOVA on unweighted UniFrac distance
adonis_unweighted <- adonis2(unweighted_unifrac ~ Organ+mean_lat+mean_lon+mean_elev+avg_MAT+avg_MAP+ClimateZone+spr_rgn, data = sample_data_df)
print(adonis_unweighted)


######################################################## Ordination plot#############################################################
# Perform PCoA ordination
ordination_unweighted <- ordinate(ps, method = "PCoA", distance = unweighted_unifrac)

# View the ordination results
ordination_unweighted

# Load ggplot2 for visualization
library(ggplot2)

# Extract the PCoA results
ordination_data <- as.data.frame(ordination_unweighted$vectors)

# Add the sample metadata (Organ column) to the ordination results
ordination_data$Organ <- sample_data_df$Organ

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = Organ)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  scale_color_manual(values = c("red", "blue", "green", "purple")) +  # Custom colors (adjust as needed)
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
#####################################################################################################################################
ordination_data$avg_MAT <- sample_data_df$avg_MAT

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = avg_MAT)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

#####################################################################################################################################
ordination_data$avg_MAP <- sample_data_df$avg_MAP

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = avg_MAP)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

#####################################################################################################################################
ordination_data$mean_elev <- sample_data_df$mean_elev

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = mean_elev)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

#####################################################################################################################################
ordination_data$mean_lat <- sample_data_df$mean_lat

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = mean_lat)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

#####################################################################################################################################
ordination_data$mean_lon <- sample_data_df$mean_lon

# Create the ordination plot
ggplot(ordination_data, aes(x = Axis.1, y = Axis.2, color = mean_lon)) +
  geom_point(size = 4) +  # Add points with a specific size
  labs(title = "PCoA of Unweighted UniFrac Distance",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +  # A minimal theme
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

