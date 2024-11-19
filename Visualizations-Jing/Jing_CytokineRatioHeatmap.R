#The following code for plotting a heatmap from already normalized data of protein expression
#Rows have Gene/Protein names, columns have sample names
library(dplyr)
library(tidyr)
library(pheatmap)
library(scales)
getwd()
# setwd("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/Visualizations/Visualizations-Jing")
data <- read.table("input/cytokine/Input data (before normalization) for cytokine heatmap_7 pts.csv", header = TRUE, sep = ",")

# Data wrangling - removed unnecessary curly braces
cytokine_cols <- setdiff(names(data), c("X", "Subject", "Group", "Timepoint"))

# Step 1: Get cytokine column names
cytokine_cols <- setdiff(names(data), c("X", "Subject", "Group", "Timepoint"))

# Step 2: Remove X column and reshape data to wide format
wide_data <- data %>%
  select(-X) %>%
  pivot_wider(
    id_cols = Subject,
    names_from = Timepoint,
    values_from = all_of(cytokine_cols),
    names_glue = "{.value}_T{Timepoint}"
  )

# Step 3: Calculate log2 ratios for each cytokine
cytokine_ratios <- data.frame(Subject = wide_data$Subject)

##Select Either Log2Fold Changes or Standardized Fold Changes across cytokine or all
# Step 4: Calculate ratios for each cytokine

# for(cytokine in cytokine_cols) {
#   t1_col <- paste0(cytokine, "_T1")
#   t2_col <- paste0(cytokine, "_T2")
#   
#   cytokine_ratios[[cytokine]] <- log2(wide_data[[t2_col]] / wide_data[[t1_col]])
# }

# Step 4: Calculate z-scores for fold changes per cytokine

# for(cytokine in cytokine_cols) {
#   t1_col <- paste0(cytokine, "_T1")
#   t2_col <- paste0(cytokine, "_T2")
#   
#   # Calculate fold change first
#   fold_changes <- wide_data[[t2_col]] / wide_data[[t1_col]]
#   
#   # Calculate z-score of fold changes
#   cytokine_ratios[[cytokine]] <- (fold_changes - mean(fold_changes, na.rm = TRUE)) / sd(fold_changes, na.rm = TRUE)
# }

# Step 4: For global z-score across all cytokines:
# First calculate all fold changes
for(cytokine in cytokine_cols) {
  t1_col <- paste0(cytokine, "_T1")
  t2_col <- paste0(cytokine, "_T2")
  
  # Calculate fold change first
  cytokine_ratios[[cytokine]] <- wide_data[[t2_col]] / wide_data[[t1_col]]
}

# Then apply z-score normalization across all values
all_values <- unlist(cytokine_ratios[,-1])  # exclude Subject column
global_mean <- mean(all_values, na.rm = TRUE)
global_sd <- sd(all_values, na.rm = TRUE)

# Apply global z-score normalization
cytokine_ratios[,cytokine_cols] <- (cytokine_ratios[,cytokine_cols] - global_mean) / global_sd







# Step 5 Convert to matrix for heatmap
ratio_matrix <- cytokine_ratios %>%
  select(-Subject) %>%
  as.matrix()

# Set row names as subject numbers
rownames(ratio_matrix) <- cytokine_ratios$Subject

# Calculate breaks for color scale
max_abs_value <- max(abs(ratio_matrix[is.finite(ratio_matrix)]), na.rm = TRUE)
breaks <- seq(-max_abs_value, max_abs_value, length.out = 100)

# Create color palette
colors <- colorRampPalette(c("blue", "white", "red"))(99)

# Plot heatmap
pheatmap(
  t(ratio_matrix),  # Transpose to match your example
  cluster_cols = TRUE,  # Cluster subjects
  cluster_rows = FALSE,  # Cluster cytokines
  clustering_method = "complete",  # Default method. #Other method: "ward.D2"
  clustering_distance_cols = "euclidean",  # Change distance method
  treeheight_col = 50,                      # Adjust dendrogram height
  na_col = "grey",     # Color for NA values
  breaks = breaks,
  color = colors,
  main = "Cytokine Expression Changes\n(Global Z-score normalization across all measurements)", #"Cytokine Expression Changes\n(Z-scores calculated per cytokine)", #"Cytokine Expression Ratio (log2(T2/T1))", #
  fontsize_row = 15,
  fontsize_col = 15,
  display_numbers = TRUE,
  angle_col = 0
)

# Different clustering methods available:
clustering_options <- list(
  "complete" = "maximum/complete linkage",
  "single" = "minimum/single linkage",
  "average" = "average linkage (UPGMA)",
  "mcquitty" = "McQuitty's method",
  "median" = "median linkage (WPGMC)",
  "centroid" = "centroid linkage (UPGMC)",
  "ward.D" = "Ward's method",
  "ward.D2" = "Ward's method version 2"
)

# Distance methods available:
distance_options <- list(
  "euclidean" = "Euclidean distance",
  "manhattan" = "Manhattan distance",
  "maximum" = "Maximum distance",
  "canberra" = "Canberra distance",
  "binary" = "Binary distance",
  "minkowski" = "Minkowski distance"
)
# "complete": Tends to produce compact clusters
# "single": Can produce long, straggly clusters
# "average": Moderate approach, often good for biological data
# "ward.D2": Tends to produce balanced, similarly-sized clusters
# "centroid": Based on cluster centers, can be unstable
# "mcquitty": Alternative that can handle unequal cluster sizes