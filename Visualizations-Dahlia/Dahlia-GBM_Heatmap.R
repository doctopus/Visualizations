library(dplyr)
library(openxlsx)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

gbm_metab <- read.csv("input/GBM + Healthy Controls for ClustVis.csv", header = TRUE, check.names = FALSE)
# Transpose the gbm_metab data
transposed_gbm_metab <- gbm_metab %>% t()
transposed_gbm_metab <- cbind(Samples = rownames(transposed_gbm_metab), transposed_gbm_metab)

# Reset the row names to NULL
# rownames(transposed_gbm_metab) <- NULL

colnames(transposed_gbm_metab) <- as.character(unlist(transposed_gbm_metab[1, ]))
transposed_gbm_metab <- transposed_gbm_metab[-1, ]

colnames(transposed_gbm_metab)[1] <- "Samples"

transposed_gbm_metab <- as.data.frame(transposed_gbm_metab) %>% 
  mutate(GBM = as.numeric(GBM))

#Keep only first two columns
gbm_sample_info <- transposed_gbm_metab %>% select(1:2) %>% mutate(GBM =as.factor(GBM))

#Make the Sample column as the rownames and remove that column
# rownames(gbm_sample_info) <- gbm_sample_info$Samples

# Remove the first column from the dataset
# gbm_sample_info <- gbm_sample_info[, -1] #Base R method does not work

# Alternative way Remove the first column from the dataset using select
gbm_sample_info <- gbm_sample_info %>% select(-Samples)

gbm_sample_info <- as.data.frame(gbm_sample_info)
metabGroups <- gbm_sample_info

#Count How many POsitive and Negative GBMs
check_counts <- metabGroups %>% 
  filter(GBM == 1) %>% 
  nrow() #100 GBM0 and 60 GBM1

gbm_metab_info <- gbm_metab[-1, ] #Remove first row GBM
rownames(gbm_metab_info) <- gbm_metab_info[, 1] #Assign Rownames
gbm_metab_info <- gbm_metab_info[, -1] # Remove the first column from the dataset

metabData <- gbm_metab_info #for GBM Data
##Now I have metabData and metabGroups

####PLOT-HEATMAP Dahlia GBM Data----
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(grid)

create_heatmap_gbm <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), , drop = FALSE]
  
  # Define the desired order for Experiment
  experiment_order <- c(1, 0) #
  # sex_order <- c("Male", "Female")
  
  # Convert Sex and Experiment to factors with specified levels
  # sample_data$Sex <- factor(sample_data$Sex, levels = sex_order)
  sample_data$GBM <- factor(sample_data$GBM, levels = experiment_order)
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#5757f9", "#606060"), experiment_order)#, "#049193", "#5B2897", pink#FF3079 brick"#EC756B" cement"#5399CA"
  # neuron_colors <- setNames(c("#007DEF", "#F08C00"), sex_order)
  
  # Create a diverging color palette for z-Score Normalized data
  # colors <- colorRampPalette(c("#377eb8", "white", "#e41a1c"))(101)
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  # max_abs <- max(abs(count_scores))
  # max_abs <- max(5, max(abs(count_scores))) #To cap the range of colors at 5. (The range is till 9, more than 5 would have only 21 extreme values)
  max_abs <- 2
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  # breaks <- seq(-3, 12, length.out = 101)
  
  ################START-For Grouping by Sex [INPUT_NEEDED]
  column_order <- order(sample_data$GBM)
  # column_order <- order(sample_data$Sex, sample_data$Experiment) #Order columns first by Sex, then by Experiment
  ################END-For Grouping by Sex  
  
  ################START-For Grouping by Experiment [INPUT_NEEDED]
  # sample_data$OrderGroup <- paste(sample_data$Experiment, sample_data$Sex, sep="_")
  # sample_data$OrderGroup <- factor(sample_data$OrderGroup,
  #                                  levels = paste(rep(experiment_order, each=2), rep(sex_order, times=4), sep="_"))
  # column_order <- order(sample_data$OrderGroup)
  ################END-For Grouping by Experiment
  
  # Ensure count_scores and sample_data are in the correct order
  count_scores <- count_scores[, column_order]
  sample_data <- sample_data[column_order, , drop = FALSE]
  
  # Create top annotation with custom colors
  ha_top <- HeatmapAnnotation(
    # Sex = sample_data$Sex,
    GBM = sample_data$GBM,
    # col = list(Sex = neuron_colors, Experiment = slide_colors),
    col = list(GBM = slide_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "right",
    gap = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 20, fontface = "bold"),
    show_legend = FALSE  # Hide default legends for annotations
  )
  
  
  
  
  
  
  # Create a two-level split for columns
  # column_split <- factor(paste(sample_data$Experiment, sample_data$Sex, sep = "_"),
  #                        levels = paste(rep(experiment_order, each = 2), rep(sex_order, times = 4), sep = "_"))
  
  # Create main heatmap
  ht <- Heatmap(count_scores,
                name = paste(pathway, "Scores"),
                col = colorRamp2(breaks, colors),
                column_title = pathway,
                column_title_gp = gpar(fontsize = 26, fontface = "bold", col = "black"),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_row_names = TRUE,
                show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                # row_title_side = c("left", "right"),
                row_names_gp = gpar(fontsize = 20, fontface = "bold"),
                top_annotation = ha_top,
                clustering_distance_rows = "spearman",
                clustering_method_rows = "average", #default is complete, ward is renamed to ward.D and there is ward.D2
                clustering_distance_columns = "spearman",
                clustering_method_columns = "average",
                # column_split = sample_data$GBM, #[INPUT_NEEDED] Change between sample_data$Sex or sample_data$Experiment
                # column_order = column_order, #Needed for matching legend order with column order
                # column_gap = unit(2, "mm"), #Option
                border = TRUE, #Option
                show_heatmap_legend = FALSE  # Hide default heatmap legend
  )
  
  # Create custom legends
  # sex_legend <- Legend(
  #   labels = sex_order,
  #   labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
  #   legend_gp = gpar(fill = neuron_colors),
  #   column_gap = unit(5, "mm"), row_gap = unit(2, "mm"),
  #   title = "Sex",
  #   title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  # )
  # 
  experiment_legend <- Legend(
    # labels = experiment_order,
    labels = rev(c("Healthy", "GBM")), #Reverse names (need to reverse legend_gp too)
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    legend_gp = gpar(fill = slide_colors),
    # legend_gp = gpar(fill = rev(slide_colors)), #Reverse associated colors too (need to have labels reversed too)
    row_gap = unit(2, "mm"),
    title = "GBM",
    title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  )
  
  expression_legend <- Legend(
    col_fun = colorRamp2(breaks, colors),
    at = c(-max_abs, 0, max_abs),
    labels = c("Low", "Medium", "High"),
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    column_gap = unit(5, "mm"), 
    row_gap = unit(5, "mm"),
    title = "Expression",
    title_gp = gpar(fontsize = 22, fontface='bold'), #Increase size of legend label
    direction = "horizontal",  # Add this line to make the legend horizontal
    legend_height = unit(2, "cm"), # Optionally adjust the legend size
    legend_width = unit(8, "cm")
  )
  
  # Combine all legends into a single column
  combined_legend <- packLegend(
    # sex_legend,
    experiment_legend,
    expression_legend,
    direction = "horizontal", #vertical or horizontal
    gap = unit(15, "mm")
  )
  
  num_rows <- nrow(count_scores)
  total_height <- unit(min(15, max(10, num_rows * 0.4)), "inch")
  
  # Draw the heatmap with only the combined legend on the left/bottom
  draw(ht, 
       annotation_legend_side = "bottom",
       annotation_legend_list = combined_legend,
       padding = unit(c(2, 20, 2, 80), "mm"),#x, left, y, right
       height = total_height)
}
plot_and_save_heatmap_gbm <- function(normalizedCountsData, sample_data_subset, pathway, output_file) {
  # Calculate height based on number of genes
  # height <- max(12, nrow(normalizedCountsData) * 0.2)  # Adjust the multiplier (0.2) as needed
  height <- 20 #19
  width <- 15 #11 or 10.5
  
  pdf(output_file, width = width, height = height)  # Increased width to accommodate legends
  create_heatmap_gbm(normalizedCountsData, pathway, sample_data_subset)
  dev.off()
  print(create_heatmap_gbm(normalizedCountsData, pathway, sample_data_subset))
  print(paste("Heatmap for", pathway, "pathway saved to", output_file))
}

#Z score Normalize the gbm_metab_info
#Per Gene Scaling of filtered Normalized count  Data
# Function to scale each row (gene) (Decide if to scale across all samples or after the filtering below)
scale_rows <- function(x) {
  (x - mean(x)) / sd(x)
}
# Apply scaling to each row
scaled_data_gbm <- t(apply(metabData, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data_gbm[is.nan(scaled_data)] <- 0

#Plot Histogram to see distribution of scaled data to decide if need to cap the color range
hist(scaled_data_gbm, breaks = 50, main="Distribution of GBM Scaled Data", xlab= "Scaled Values of Metabolites")
abline(v = c(-2, 2), col = "red")

# Plot and save the heatmap
plot_and_save_heatmap_gbm(
  scaled_data_gbm, #data_metab_final or scaled_data
  metabGroups, 
  "Metabolites in GBM vs Control", 
  "20250613_GBM Metabolites -Distance-spearman Clustering-average.pdf"
)
