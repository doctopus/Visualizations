#To make PCA Plots of the metabolites grouped by Clinical Characteristics, like GBM or PFS
#Input gbm_metab_info and gbm_sample_info_combined from Dahlia_Heatmap Script

# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(stats)

# Function to perform PCA and create plots
create_pca_plots <- function(metab_data, sample_info, color_var, title, custom_labels = NULL) {
  # Transpose metabolite data to have samples as rows
  metab_matrix <- t(as.matrix(metab_data))
  
  # Perform PCA ensuring rownames are preserved
  pca_result <- prcomp(metab_matrix, scale. = TRUE, center = TRUE)
  # Explicitly set the rotation rownames to be the metabolite names (For later extracting top variables)
  rownames(pca_result$rotation) <- colnames(metab_matrix)
  
  # Create a dataframe for plotting
  pca_df <- as.data.frame(pca_result$x)
  
  # Add sample information
  pca_df$SampleID <- rownames(pca_df)
  pca_df <- merge(pca_df, 
                  sample_info %>% mutate(SampleID = rownames(sample_info)), 
                  by = "SampleID")
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  # Create the plot
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(!!sym(color_var)))) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, show.legend = FALSE) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = color_var
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # Apply custom labels if provided
  if (!is.null(custom_labels)) {
    pca_plot <- pca_plot + 
      scale_color_discrete(name = names(custom_labels)[1],  # Gets the legend title
                           labels = unname(custom_labels[[1]]))  # Gets the labels
  }
  
  return(list(plot = pca_plot, pca_obj = pca_result))
}

# Create updated PCA plots
# For GBM categories with custom labels


# Create PCA plot for GBM categories
gbm_pca <- create_pca_plots(
  gbm_metab_info,
  gbm_sample_info_combined,
  "GBM",
  "Metabolite Profiles in Control vs Glioblastoma",
  custom_labels = list("Diagnosis" = c("0" = "Control", "1" = "Glioblastoma"))
)

# Create PCA plot for PFS categories
pfs_pca <- create_pca_plots(
  gbm_metab_info,
  gbm_sample_info_combined,
  "PFS",
  "Metabolite Profiles by PFS Status",
  custom_labels = list("PFS Status" = c("0" = "Low", "1" = "Medium", "2" = "High"))  # If you want to rename PFS categories
)


# Display the plots
print(gbm_pca$plot + 
        scale_color_manual(name = "Diagnosis",
                           values = c("0" = "#2E86C1", "1" = "#E74C3C"),
                           labels = c("0" = "Control", "1" = "Glioblastoma")
                           )
      )
print(pfs_pca$plot + 
        scale_color_manual(name = "PFS Status",
                           values = c("0" = "#27AE60", "1" = "#F4D03F", "2" = "#E67E22"),
                           labels = c("0" = "PFS<9Months", "1" = "PFS>9Months", "2" = "Control")
                           )
        )

# Optional: Save the plots with higher resolution
# ggsave("gbm_pca_plot.pdf", gbm_pca$plot, width = 10, height = 8, dpi = 300)
# ggsave("pfs_pca_plot.pdf", pfs_pca$plot, width = 10, height = 8, dpi = 300)

# Optional: Get loadings for the most important metabolites
get_top_loadings <- function(pca_obj, n = 10) {
  loadings <- as.data.frame(pca_obj$rotation)
  
  # Create named vectors with absolute values
  pc1_abs <- abs(loadings$PC1)
  pc2_abs <- abs(loadings$PC2)
  names(pc1_abs) <- rownames(loadings)
  names(pc2_abs) <- rownames(loadings)
  
  # Get top contributors
  pc1_contrib <- sort(pc1_abs, decreasing = TRUE)[1:n]
  pc2_contrib <- sort(pc2_abs, decreasing = TRUE)[1:n]
  
  # Create data frames
  pc1_df <- data.frame(
    metabolite = names(pc1_contrib),
    loading = pc1_contrib,
    stringsAsFactors = FALSE
  )
  
  pc2_df <- data.frame(
    metabolite = names(pc2_contrib),
    loading = pc2_contrib,
    stringsAsFactors = FALSE
  )
  
  return(list(
    PC1_top = pc1_df,
    PC2_top = pc2_df
  ))
}

# Get top contributing metabolites
gbm_top_metabolites <- get_top_loadings(gbm_pca$pca_obj)
pfs_top_metabolites <- get_top_loadings(pfs_pca$pca_obj)

# To view the results:
print("Top PC1 contributors:")
print(gbm_top_metabolites$PC1_top)
print("\nTop PC2 contributors:")
print(gbm_top_metabolites$PC2_top)


#### # Function to create loading plots ----
plot_top_loadings <- function(pca_loadings, title) {
  # Combine PC1 and PC2 data
  pc1_data <- pca_loadings$PC1_top %>%
    mutate(PC = "PC1")
  pc2_data <- pca_loadings$PC2_top %>%
    mutate(PC = "PC2")
  
  combined_data <- rbind(pc1_data, pc2_data)
  
  # Create the plot
  ggplot(combined_data, aes(x = reorder(metabolite, loading), y = abs(loading), fill = PC)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    coord_flip() +  # Flip coordinates for horizontal bars
    scale_fill_manual(values = c("PC1" = "#2E86C1", "PC2" = "#E74C3C")) +
    labs(
      title = title,
      x = "Metabolite",
      y = "Absolute Loading Value",
      fill = "Principal Component"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor.y = element_blank()
    )
}

# Create plots for both GBM and PFS analyses
gbm_loadings_plot <- plot_top_loadings(gbm_top_metabolites, 
                                       "Top Contributing Metabolites in GBM Analysis")
pfs_loadings_plot <- plot_top_loadings(pfs_top_metabolites, 
                                       "Top Contributing Metabolites in PFS Analysis")

# Display the plots
print(gbm_loadings_plot)
print(pfs_loadings_plot)

# Optionally save the plots
# ggsave("gbm_loadings_plot.pdf", gbm_loadings_plot, width = 10, height = 8)
# ggsave("pfs_loadings_plot.pdf", pfs_loadings_plot, width = 10, height = 8)


