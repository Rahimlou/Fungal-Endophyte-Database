library(ape)
library(dplyr)
library(tidyr)
library(vegan)


# Permutation Test
# To test if this correlation is stronger than expected by chance, 
# we randomly shuffle the host-fungus associations and recompute the correlation multiple times.

setwd("E:/Fungal_Endophyte_Database/Final_Dataset/Organ")

# Preparing the data for Mantel test
df <- read.csv("FunEndo_Final_Edit_Host.csv", header = TRUE, fileEncoding = "Windows-1252")

subset_df <- df[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Host")]

# Remove rows where Genus are missing
subset_df <- subset_df %>%
  filter(!(is.na(Genus)))

# Remove rows where Hosts are missing
subset_df <- subset_df %>%
  filter(!(is.na(Host)))


# Split the 'Host' column
subset_df <- subset_df %>%
  separate(Host, into = c("Host_Genus", "Host_Species"), sep = " ", extra = "merge")

# Remove the Host_Species column
subset_df <- subset(subset_df, select = -Host_Species)


fungal_tree <- read.tree("endophyte_taxonomy.nwk")
host_tree <- read.tree("host_taxonomy.nwk")
fungal_host_matrix <- read.csv("presence_absence_matrix.csv", row.names = 1)
fungal_host_matrix <- as.matrix(association_matrix)
fungal_host_matrix <- t(fungal_host_matrix)

# Calculate phylogenetic distances for fungal tree
fungal_dist_matrix <- cophenetic(fungal_tree)

# Calculate phylogenetic distances for host tree
host_dist_matrix <- cophenetic(host_tree)

# Check Dimensions of Distance Matrices
dim(fungal_host_matrix)
dim(fungal_dist_matrix)
dim(host_dist_matrix)

common_hosts <- intersect(colnames(fungal_host_matrix), rownames(host_dist_matrix))
common_fungi <- intersect(rownames(fungal_host_matrix), rownames(fungal_dist_matrix))


# Compute Observed Mean Pairwise Taxonomic Distance
# We calculate the mean pairwise distance between hosts of the same fungus and compare it to random associations.

# Compute Mean Host Taxonomic Distance Per Fungus

# Initialize vector to store mean host distances per fungus
observed_host_distances <- numeric(nrow(fungal_host_matrix))

# calculating the mean pairwise taxonomic distance between hosts that are associated with a particular fungus.
for (i in 1:nrow(fungal_host_matrix)) {
  associated_hosts <- which(fungal_host_matrix[i, ] == 1)  # Get associated host indices
  if (length(associated_hosts) > 1) {
    observed_host_distances[i] <- mean(as.dist(host_dist_matrix[associated_hosts, associated_hosts]))
  } else {
    observed_host_distances[i] <- NA  # No valid comparison
  }
}

# Compute overall observed mean distance across fungi
observed_mean_distance <- mean(observed_host_distances, na.rm = TRUE)

# Permutation Test (Shuffle Fungal-Host Associations)
# Now, we randomly shuffle host-fungus associations to create a null expectation.

set.seed(123)  # Ensure reproducibility
n_perm <- 999
random_distances <- numeric(n_perm)

for (i in 1:n_perm) {
  # Shuffle host associations randomly
  shuffled_matrix <- fungal_host_matrix
  shuffled_matrix[] <- sample(as.vector(fungal_host_matrix))
  
  permuted_host_distances <- numeric(nrow(shuffled_matrix))
  
  for (j in 1:nrow(shuffled_matrix)) {
    associated_hosts <- which(shuffled_matrix[j, ] == 1)
    if (length(associated_hosts) > 1) {
      permuted_host_distances[j] <- mean(as.dist(host_dist_matrix[associated_hosts, associated_hosts]))
    } else {
      permuted_host_distances[j] <- NA
    }
  }
  
  random_distances[i] <- mean(permuted_host_distances, na.rm = TRUE)
}

# Compute p-value: Proportion of permuted distances smaller than observed
p_value <- sum(random_distances <= observed_mean_distance) / n_perm
print(paste("P-value:", p_value))

# Interpretation
# If p < 0.05, the observed mean taxonomic distance of hosts per fungus is significantly lower than expected by chance, meaning closely related fungi associate with closely related hosts.
# If p > 0.05, fungal-host associations appear random with respect to taxonomy.









