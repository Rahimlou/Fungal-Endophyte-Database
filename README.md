# FunEndo: Global Fungal Endophyte Database and Analysis Pipeline

## Overview

FunEndo is a comprehensive database and reproducible analytical framework for investigating global patterns of fungal endophyte diversity, host associations, and ecological drivers. The repository integrates curated datasets with modular R scripts for statistical, phylogenetic, and network-based analyses.

The dataset comprises ~105,000 records from >4,500 studies, covering ~4,300 fungal species across ~1,700 genera associated with ~4,100 plant species.

---

## Repository Structure


├── Bipartite_network_analysis.R  
├── Fisher's_exact_test.R  
├── Permutation_Test_Host.R  
├── UniFrac.R  
├── db_RDA_variance_partitioning.R  
└── README.md  

---

## Analysis Workflow

The analytical pipeline consists of five main components. Scripts are modular but designed to be used in the following logical order:

---

### **Step 1: Data Preparation (External to this repo)**
Before running the scripts, ensure:
- Taxonomy is standardized (e.g., GBIF backbone)
- Host plant metadata is harmonized
- Data are formatted as presence/absence or abundance matrices

---

### **Step 2: Organ Specificity Analysis**

**Script:** `Fisher's_exact_test.R`

**Purpose:**
- Identify fungal genera significantly enriched in specific plant organs

**Method:**
- Fisher’s exact test for each genus–organ combination
- Multiple testing correction (FDR)

**Output:**
- Table of enriched genera with adjusted p-values
- Identification of organ specialists vs generalists

---

### **Step 3: Host Association Structure**

**Script:** `Permutation_Test_Host.R`

**Purpose:**
- Test whether fungal–host associations differ from random expectations

**Method:**
- Permutation-based null models
- Random reshuffling of host–fungus associations
- Comparison of observed vs null distributions

**Output:**
- Z-scores (standardized effect sizes)
- P-values indicating deviation from randomness

---

### **Step 4: Network Analysis**

**Script:** `Bipartite_network_analysis.R`

**Purpose:**
- Characterize structure of host–fungus interaction networks

**Method:**
- Construction of bipartite networks
- Modularity detection
- Community clustering

**Output:**
- Network modularity metrics
- Module assignments
- Visualization-ready outputs

---

### **Step 5: Phylogenetic Community Structure**

**Script:** `UniFrac.R`

**Purpose:**
- Quantify phylogenetic similarity among fungal communities

**Method:**
- UniFrac distance (weighted/unweighted)
- Phylogenetic tree-based comparisons

**Output:**
- Distance matrices
- Community similarity metrics

---

### **Step 6: Environmental and Host Drivers**

**Script:** `db_RDA_variance_partitioning.R`

**Purpose:**
- Identify drivers of community composition

**Method:**
- Distance-based redundancy analysis (db-RDA)
- Variance partitioning across:
  - Host traits (organ, growth form)
  - Climate variables
  - Spatial structure

**Output:**
- Explained variance components
- Statistical significance of predictors
- Ordination plots

---

## Dependencies

R (≥ 4.0)

Required packages:

```r
install.packages(c(
  "tidyverse",
  "vegan",
  "ape",
  "phyloseq",
  "picante",
  "bipartite",
  "igraph"
))
