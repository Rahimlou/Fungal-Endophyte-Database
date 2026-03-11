library(tidyr)
library(tibble)
library(dplyr)
library(ape)
library(vegan)
library(phyloseq)
library(ggplot2)

setwd("/Users/srahimlou/OneDrive_The_Pennsylvania_State_University/FunEndo_R_Codes")

# Preparing the data for PERMANOVA
df <- read.csv("FunEndo_DB.csv", header = TRUE, fileEncoding = "latin1")

# Keep kingdoms that are Fungi
df <- subset(df, Kingdom == "Fungi")

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
df$Organ <- gsub("Branch\\(Twig\\)", "stem", df$Organ)
df <- subset(df, Organ != "Photosynthetictissue")

df <- subset(df, Organ != "other")
###########################################################
# Classify Organs to below-ground and above-ground
# df$Organ <- gsub("leaf", "above-ground", df$Organ)
# df$Organ <- gsub("flower", "above-ground", df$Organ)
# df$Organ <- gsub("stem", "above-ground", df$Organ)
# df$Organ <- gsub("root", "below-ground", df$Organ)
# df$Organ <- gsub("rhizome", "below-ground", df$Organ)
# df$Organ <- gsub("stolon", "above-ground", df$Organ)
# df$Organ <- gsub("needle", "above-ground", df$Organ)
# df$Organ <- gsub("bark", "above-ground", df$Organ)
# df$Organ <- gsub("trunk", "above-ground", df$Organ)
# df$Organ <- gsub("bark", "above-ground", df$Organ)
# df$Organ <- gsub("fruit", "above-ground", df$Organ)
# df$Organ <- gsub("shoot", "above-ground", df$Organ)
# df$Organ <- gsub("tuber", "below-ground", df$Organ)
# df$Organ <- gsub("twig", "above-ground", df$Organ)

# df <- subset(df, Organ != "other")
# df <- subset(df, Organ != "wood")
# df <- subset(df, Organ != "gall")
# df <- subset(df, Organ != "seed")
# 
# unique(df$Organ)


# Remove rows where Genus is missing
df <- df %>%
  filter(!(is.na(Genus)))

# Subset data frame
subset_df <- df[, c("Genus","Organ")]

# Create a contingency table
contingency_table <- table(subset_df$Genus, subset_df$Organ)
contingency_table <- as.matrix(contingency_table)
dim(contingency_table)
print(contingency_table)

# Genus = 1359
# Organ = 14

# Do not run this if you want to calculate the p-value for above and below-ground organs
# Calculate the sum of each column (organ)
# column_sums <- colSums(contingency_table)

# Identify columns where the sum is less than 10
# columns_to_remove <- column_sums < 10

# Subset the contingency table to keep only columns where the sum is >= 10
# contingency_table <- contingency_table[, !columns_to_remove]
# dim(contingency_table)

# Remove the column named "other" and "gall"
# contingency_table <- contingency_table[, !colnames(contingency_table) %in% "other"]
# contingency_table <- contingency_table[, !colnames(contingency_table) %in% "gall"]
# dim(contingency_table)

# Do not run this if you want to calculate the p-value for subset of data
# Calculate the sum of each row (Genus)
# row_sums <- rowSums(contingency_table)

# Identify rows where the sum is less than 3
# rows_to_remove <- row_sums < 3

# Subset the contingency table to keep only rows where the sum is >= 3
# contingency_table <- contingency_table[!rows_to_remove, ]
# dim(contingency_table)
# 497 genera remained

write.csv(contingency_table, "contingency_table.csv")
#########################################################

# Convert numbers to "Presence" and "Absence"
# presence_absence_table <- ifelse(contingency_table > 0, 1, 0)
# dim(presence_absence_table)
# 
# # Calculate sum for each column (organ)
# column_sums <- colSums(presence_absence_table)
# print(column_sums)
# 
# presence_absence_table <- as.data.frame(presence_absence_table)
# 
# #Identify rows where wood is present (1) and where wood is absent (0)
# present_genera <- which(contingency_table[, "wood"] > 0)
# absent_genera <- which(contingency_table[, "wood"] == 0)
# 
# # Randomly sample 3 times the number of present genera from absent genera
# num_present <- length(present_genera)
# 
# set.seed(42)  # Set seed for reproducibility
# subsampled_absent_genera <- sample(absent_genera, size = 3 * num_present, replace = FALSE)

# # Combine the present genera and the subsampled absent genera
# final_genera <- c(present_genera, subsampled_absent_genera)
# 
# # Create the final table using the selected genera
# final_table_wood <- presence_absence_table[final_genera, ]
# dim(final_table_wood)
#########################################################
# Run this if you want to calculate the p-value for a subset of data
# rownames <- rownames(final_table_wood)

# Subset the contingency table to keep only the genera in 'genera_list'
# subset_table <- contingency_table[rownames(contingency_table) %in% rownames, ]
# dim(subset_table)
# 
# # Convert to data frame
# subset_table <- as.data.frame(as.table(subset_table))
# 
# # Expand data frame to repeat combinations based on the frequency
# subset_table <- subset_table[rep(1:nrow(subset_table), subset_table$Freq), c("Var1", "Var2")]
# 
# # Rename columns for clarity
# colnames(subset_table) <- c("Genus", "Organ")
# 
# length(unique(subset_table$Genus))

#######################################################################################################
# Run this if you want to calculate p-value for all data
# Convert to data frame
contingency_df <- as.data.frame(as.table(contingency_table))

# Expand data frame to repeat combinations based on the frequency
contingency_df <- contingency_df[rep(1:nrow(contingency_df), contingency_df$Freq), c("Var1", "Var2")]

# Rename columns for clarity
colnames(contingency_df) <- c("Genus", "Organ")


#########################################################################################################
# data <- subset_table

data <- contingency_df

# data <- subset_df

# Get unique genera and organs
genera <- unique(data$Genus)
organs <- unique(data$Organ)

length(genera)
length(organs)

# Create a data frame to store the Fisher's test results (P-values)
fisher_p_values <- data.frame(
  Genus = character(),
  P_Value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# Iterate over each genus to create a contingency table and run Fisher's Exact Test
for (genus in genera) {
  
  # Create an empty 10x2 matrix for the contingency table (one row per organ)
  contingency_table <- matrix(0, nrow = length(organs), ncol = 2)
  colnames(contingency_table) <- c("Genus_Present", "Genus_Absent")
  rownames(contingency_table) <- organs
  
  # Fill the contingency table: count presence and absence of the genus for each organ
  for (i in 1:length(organs)) {
    organ <- organs[i]
    
    # Count presence of the genus in this organ
    genus_present <- sum(data$Genus == genus & data$Organ == organ)
    
    # Count absence of the genus in this organ
    genus_absent <- sum(data$Genus != genus & data$Organ == organ)
    
    # Populate the contingency table
    contingency_table[i, "Genus_Present"] <- genus_present
    contingency_table[i, "Genus_Absent"] <- genus_absent
  }
  
  # Run Fisher's Exact Test on the contingency table (for each organ) with simulation
  fisher_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 9999)
  
  # Get the p-value and assign the significance
  p_value <- fisher_result$p.value
  if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  # Store the P-value and significance for this genus
  fisher_p_values <- rbind(fisher_p_values, data.frame(
    Genus = genus,
    P_Value = p_value,
    Significance = significance
  ))
}

# View the Fisher's Exact Test results with P-values and significance annotations
write.csv(fisher_p_values, "fisher_p_values.csv")
write.csv(fisher_p_values, "fisher_p_values_subset.csv")
write.csv(fisher_p_values, "fisher_p_values_above_vs_below_ground.csv")



############################## Plotting the results ############################
# Summary Plot: How Many Genera Are Organ-Specific?

summary_stats <- data.frame(
  Significant = c("Yes","No"),
  Count = c(sum(fisher_p_values$P_Value < 0.05),
            sum(fisher_p_values$P_Value >= 0.05))
)

ggplot(summary_stats, aes(x = Significant, y = Count, fill = Significant)) +
  geom_col() +
  scale_fill_manual(values = c("#d91811", "#0d7807")) +
  labs(title = "Number of organ-specific fungal genera") +
  theme_minimal()

# Update the plot using FDR
fisher_p_values$FDR <- p.adjust(fisher_p_values$P_Value, method = "BH")
sig_genera <- fisher_p_values[fisher_p_values$FDR < 0.05, ]
sum(fisher_p_values$P_Value < 0.05)
sum(fisher_p_values$P_Value < 0.05)
summary_stats <- data.frame(
  Significant = c("Yes", "No"),
  Count = c(
    sum(fisher_p_values$FDR < 0.05),
    sum(fisher_p_values$FDR >= 0.05)
  )
)

ggplot(summary_stats, aes(x = Significant, y = Count, fill = Significant)) +
  geom_col() +
  scale_fill_manual(values = c("#d91811", "#0d7807")) +
  labs(title = "Number of organ-specific fungal genera (FDR-adjusted)") +
  theme_minimal()
################################################################################
# Heatmap of Enrichment Strength (Best Complement)

# Removed the following organs
data_filtered <- data %>%
  filter(!Organ %in% c("bark", "stolon", "tuber", "gall"))

assoc_matrix <- data_filtered %>%
  group_by(Genus, Organ) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Genus) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
# Filter to significant genera
assoc_sig <- assoc_matrix %>%
  filter(Genus %in% sig_genera$Genus)
# Plot Heatmap
ggplot(assoc_sig, aes(x = Organ, y = Genus, fill = prop)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#ffd700", "#ff4500", "#8b0000", "#2b0000"),
    name = "Normalized association",
    limits = c(0, 1)
  ) +
  labs(
    title = "Organ association strength of significant fungal genera",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

################################################################################

# Volcano-Style Significance Plot (Recommended)

fisher_p_values$logP <- -log10(fisher_p_values$P_Value)

fisher_p_values <- fisher_p_values %>%
  mutate(Significant = ifelse(logP >= -log10(0.05), "Significant", "Not significant"))

write.csv(fisher_p_values, "fisher_p_values.csv")

ggplot(fisher_p_values, aes(x = reorder(Genus, logP), y = logP, color = Significant)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  coord_flip() +
  scale_color_manual(values = c(
    "Significant" = "#1B9E77",   # green
    "Not significant" = "#D95F02" # orange-red
  )) +
  labs(
    x = "Fungal genus",
    y = expression(-log[10](P-value)),
    title = "Organ specificity of fungal genera (Fisher's Exact Test)",
    color = "FDR status"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", size = 3),  # <--- Y-axis font size reduced
    axis.text.x = element_text(color = "black", size = 10), # optional: control X-axis separately
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggplot(fisher_p_values, aes(x = reorder(Genus, logP), y = logP, color = Significant)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  coord_flip() +
  scale_color_manual(values = c(
    "Significant" = "#1B9E77",   # green
    "Not significant" = "#D95F02" # orange-red
  )) +
  labs(
    x = "Fungal Genera",  # Y-axis label after coord_flip()
    y = expression(-log[10](P-value)),
    title = "Organ specificity of fungal genera (Fisher's Exact Test)",
    color = "FDR status"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),    # remove genus names
    axis.ticks.y = element_blank(),   # remove Y-axis ticks
    axis.text.x = element_text(color = "black", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )



# -------------------------
# Step 1: Filter out unwanted organs
# -------------------------
data_filtered <- data %>%
  filter(!Organ %in% c("bark", "stolon", "tuber", "gall"))

# -------------------------
# Step 2: Create normalized association matrix
# -------------------------
assoc_matrix <- data_filtered %>%
  group_by(Genus, Organ) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Genus) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# -------------------------
# Step 3: Calculate total number of records per genus
# -------------------------
n_records_df <- assoc_matrix %>%
  group_by(Genus) %>%
  summarise(n_records = sum(n), .groups = "drop")

# -------------------------
# Step 4: Prepare Fisher's test results
# (Assuming fisher_p_values exists with Genus, P_Value, FDR)
# -------------------------

# Merge n_records into fisher_p_values
fisher_p_values <- fisher_p_values %>%
  left_join(n_records_df, by = "Genus")

# Add logP column and significance stars
fisher_p_values <- fisher_p_values %>%
  mutate(
    logP = -log10(P_Value),
    stars = cut(FDR,
                breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("***", "**", "*", ""))
  )

# -------------------------
# Step 5: Identify the most significant genera for labeling
# -------------------------

# Choose top N genera by smallest FDR (or highest logP)
top_genera <- fisher_p_values %>%
  filter(FDR < 0.05) %>%        # only significant
  arrange(FDR) %>%
  slice_head(n = 10)             # label top 5

# -------------------------
# Step 6: Plot publication-ready figure
# -------------------------

# Significance threshold line
sig_threshold <- -log10(0.05)

# Colorblind-friendly palette
color_palette <- c("Significant" = "#0072B2", "Not significant" = "#E69F00")

# Plot
ggplot(fisher_p_values, aes(x = reorder(Genus, logP), y = logP, 
                            color = Significant, size = n_records)) +
  geom_point(alpha = 0.85) +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", linewidth = 0.8, color = "black") +
  # Add significance stars
  geom_text(aes(label = stars), nudge_y = 0.15, color = "black", size = 3) +
  # Label top significant genera
  geom_text(data = top_genera,
            aes(label = Genus),
            nudge_y = 0.5,
            color = "black",
            fontface = "bold",
            size = 3) +
  coord_flip() +
  scale_color_manual(values = color_palette) +
  scale_size_continuous(range = c(2,6)) +
  labs(
    x = "Fungal Genera",
    y = expression(-log[10](P-value)),
    title = "Organ specificity of fungal genera (Fisher's Exact Test)",
    color = "Significance (FDR)",
    size = "Number of reports"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# -------------------------
# Step 7: Optional save as high-resolution PDF
# -------------------------
ggsave("Fisher_Organ_Specificity_Labeled.pdf", width = 6, height = 8, dpi = 300)
