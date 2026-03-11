# ----------------------------
# Bipartite Network Analysis (FAST VERSION)
# Using igraph Louvain Modularity
# ----------------------------

library(igraph)
library(ape)
library(picante)
library(dplyr)
library(tidyr)
library(vegan)

setwd("/Users/srahimlou/OneDrive_The_Pennsylvania_State_University/FunEndo_R_Codes/")

# ----------------------------
# Load Data
# ----------------------------
df <- read.csv("FundEndo_Final_v.02.csv", header = TRUE)
df <- df %>%
  mutate(Host_Genus = ifelse(Host_Kingdom == "" | is.na(Host_Kingdom), NA, Host_Genus))
df$Host_Kingdom <- gsub("†", "", df$Host_Kingdom)

# remove rows where the Kingdom is "Chromista" and "Stramenopiles"
df <- df[df$Kingdom != "Chromista", ]
df <- df[df$Kingdom != "Stramenopiles", ]

write.csv(df, "FundEndo_Final_v.03.csv")

# subset data to Fungi and their hosts
subset_df <- df %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, 
         Host_Kingdom, Host_Phylum, Host_Class, Host_Order, Host_Family, Host_Genus)

subset_df$Host_Kingdom <- gsub("†", "", subset_df$Host_Kingdom)

# remove rows where the Genus is empty
subset_df <- subset_df %>%
  filter(!is.na(Genus))

# remove rows where the Host_Genus is empty
subset_df <- subset_df %>%
  filter(!is.na(Host_Genus))

# Constructing presence-absence matrix
presence_absence <- subset_df %>%
  # Keep only relevant columns
  select(Genus, Host_Genus) %>%
  
  # Remove empty values (important!)
  filter(!is.na(Genus),
         !is.na(Host_Genus),
         trimws(Genus) != "",
         trimws(Host_Genus) != "") %>%
  
  # Remove duplicate fungus–host pairs
  distinct() %>%
  
  # Add presence column
  mutate(presence = 1) %>%
  
  # Convert to wide format
  pivot_wider(
    names_from = Host_Genus,
    values_from = presence,
    values_fill = 0
  )

# Make the first column as row names
pa_matrix <- presence_absence %>%
  tibble::column_to_rownames(var = colnames(presence_absence)[1])

# Extract Fungal Taxonomy Matching Matrix Rows
fungal_tax <- subset_df %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  distinct() %>%
  filter(Genus %in% rownames(pa_matrix)) %>%
  arrange(match(Genus, rownames(pa_matrix)))

# Extract Host Taxonomy Matching Matrix Columns
host_tax <- subset_df %>%
  select(Host_Kingdom, Host_Phylum, Host_Class,
         Host_Order, Host_Family, Host_Genus) %>%
  distinct() %>%
  filter(Host_Genus %in% colnames(pa_matrix)) %>%
  arrange(match(Host_Genus, colnames(pa_matrix)))

write.csv(fungal_tax, "fungal_tax.csv")
write.csv(host_tax, "host_tax.csv")
write.table(fungal_tax, "fungal_tax.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(host_tax, "host_tax.txt", sep = "\t", row.names = FALSE, quote = FALSE)

host_tree <- read.tree("host_tree.nwk")
fungal_tree <- read.tree("fungal_tree.nwk")

# ----------------------------
# Remove Rare Nodes (Speed Optimization)
# ----------------------------

host_sums <- colSums(pa_matrix)
pa_matrix_small <- pa_matrix[, host_sums >= 2]

fungi_sums <- rowSums(pa_matrix_small)
pa_matrix_small <- pa_matrix_small[fungi_sums >= 2, ]

cat("Final matrix size:",
    nrow(pa_matrix_small),
    "fungi x",
    ncol(pa_matrix_small),
    "hosts\n")

# ----------------------------
# Build Bipartite Graph
# ----------------------------

g <- graph_from_biadjacency_matrix(pa_matrix_small,
                                   weighted = TRUE,
                                   directed = FALSE)

# ----------------------------
# Compute Modularity (FAST)
# ----------------------------

set.seed(123)

louvain_mod <- cluster_louvain(g)

Q <- modularity(louvain_mod)
cat("Observed modularity (Q):", Q, "\n")

membership_vec <- membership(louvain_mod)

# Separate host memberships only
host_nodes <- colnames(pa_matrix_small)
host_modules <- membership_vec[host_nodes]

# ----------------------------
# Permutation Test for Modularity
# ----------------------------

nperm <- 999
null_Q <- numeric(nperm)

for(i in 1:nperm){
  
  shuffled_mat <- pa_matrix_small
  shuffled_mat[] <- sample(shuffled_mat)
  
  g_null <- graph_from_incidence_matrix(shuffled_mat)
  
  mod_null <- cluster_louvain(g_null)
  
  null_Q[i] <- modularity(mod_null)
}

p_mod <- mean(null_Q >= Q)

cat("Modularity permutation p-value:", p_mod, "\n")

# ----------------------------
# Host Phylogenetic Clustering Within Modules
# ----------------------------

host_tree2 <- drop.tip(host_tree,
                       setdiff(host_tree$tip.label,
                               names(host_modules)))

host_dist <- cophenetic(host_tree)

calc_module_mpds <- function(module_id){
  
  hosts <- names(host_modules[host_modules == module_id])
  
  if(length(hosts) < 3) return(NA)
  
  mean(host_dist[hosts, hosts][lower.tri(host_dist[hosts, hosts])])
}

module_ids <- unique(host_modules)

observed_mpds <- sapply(module_ids, calc_module_mpds)

# ----------------------------
# Null Model for Phylogenetic Clustering
# ----------------------------

null_mpds <- replicate(nperm, {
  
  shuffled <- sample(host_modules)
  
  sapply(module_ids, function(m){
    
    hosts <- names(shuffled[shuffled == m])
    
    if(length(hosts) < 3) return(NA)
    
    mean(host_dist[hosts, hosts][lower.tri(host_dist[hosts, hosts])])
  })
  
})

# Z-scores + p-values
z_scores <- (observed_mpds -
               rowMeans(null_mpds, na.rm = TRUE)) /
  apply(null_mpds, 1, sd, na.rm = TRUE)

p_values <- rowMeans(null_mpds <= observed_mpds, na.rm = TRUE)

phylo_module_df <- data.frame(
  Module = module_ids,
  MPD = observed_mpds,
  Z = z_scores,
  P = p_values
)

print(phylo_module_df)

# ----------------------------
# Visualization
# ----------------------------

# Color nodes by module
V(g)$color <- membership_vec

# Basic network visualization
plot(g,
     vertex.size = 3,
     vertex.label = NA,
     edge.width = 0.4,
     layout = layout_with_fr,
     main = "Host–Fungus Bipartite Network (Louvain Modules)")
