#The following code for plotting a heatmap from already normalized data of protein expression
#Rows have Gene/Protein names, columns have sample names
library(dplyr)
library(openxlsx)
getwd()
# setwd("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/Visualizations/Visualizations-Jing")

library(dplyr)
library(tidyr)
library(pheatmap)
library(scales)

data <- read.table("input/cytokine/Input data (before normalization) for cytokine heatmap_7 pts.csv", header = TRUE, sep = ",")

# Load required libraries
library(dplyr)

# Data wrangling - removed unnecessary curly braces
cytokine_cols <- setdiff(names(data), c("X", "Subject", "Group", "Timepoint"))


# Load required libraries
library(dplyr)
library(tidyr)
library(pheatmap)
library(scales)

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

# Step 4: Calculate ratios for each cytokine
for(cytokine in cytokine_cols) {
  t1_col <- paste0(cytokine, "_T1")
  t2_col <- paste0(cytokine, "_T2")
  
  cytokine_ratios[[cytokine]] <- log2(wide_data[[t2_col]] / wide_data[[t1_col]])
}

# Convert to matrix for heatmap
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
  na_col = "grey",     # Color for NA values
  breaks = breaks,
  color = colors,
  main = "Cytokine Expression Ratio (log2(T2/T1))",
  fontsize_row = 15,
  fontsize_col = 15,
  display_numbers = TRUE,
  angle_col = 0
)
