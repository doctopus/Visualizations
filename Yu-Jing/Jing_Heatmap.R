#The following code for plotting a heatmap from already normalized data of protein expression
#Rows have Gene/Protein names, columns have sample names
library(dplyr)
library(openxlsx)
getwd()
setwd("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/Visualizations/Yu-Jing")

data <- read.table("input/Jenifer_Yu_Spatial_Protein_Dataset.txt", header = TRUE, sep = "\t")

#Need to remove columns that contains these terms
data1 <- data %>% select(-matches("Norm|Nom|Edge"))
# dataRemoved <- data %>% select(matches("Norm|Nom|Edge"))

# write.xlsx(data1, "Spatial_Protein_Dataset_Without_Norm_OR_Edge.xlsx", showNA = FALSE)

#All Male
# data_M <- data1 %>% select(matches("M.Tumor|Protien"))
# write.xlsx(data_M, "Male_All.xlsx", showNA = FALSE)

#All Female
# data_F <- data1 %>% select(matches("F.Tumor|Protien"))
# write.xlsx(data_F, "Female_All.xlsx", showNA = FALSE)

#Male Control
data_M_C.M <- data1 %>% select(matches("C.M.Tumor|Protien"))
# write.xlsx(data_M_C.M, "Male_Control.xlsx", showNA = FALSE)

#Male Radiation
data_M_R.M <- data1 %>% select(matches("R.M.Tumor|Protien"))
# write.xlsx(data_M_R.M, "Male_Radiation.xlsx", showNA = FALSE)

#Male Hyperthermia
data_M_H.M <- data1 %>% select(matches("ROI.H.M.Tumor|Protien"))
#Remove the 3 A columns that is wrongly present
data_M_H.M <- data_M_H.M %>% select(-matches("A.CD45"))
# write.xlsx(data_M_H.M, "Male_Hyperthermia.xlsx", showNA = FALSE)

#Male Radiation Hyperthermia
data_M_R.H.M <- data1 %>% select(matches("ROI.R.H.M.Tumor|Protien"))
# write.xlsx(data_M_R.H.M, "Male_Radiation_Hypethermia.xlsx", showNA = FALSE)

#Female Control
data_F_C.F <- data1 %>% select(matches("C.F.Tumor|Protien"))
# write.xlsx(data_F_C.F, "Female_Control.xlsx", showNA = FALSE)

#Female Radiation
data_F_R.F <- data1 %>% select(matches("R.F.Tumor|Protien"))
# write.xlsx(data_F_R.F, "Female_Radiation.xlsx", showNA = FALSE)

#Female Hyperthermia
data_F_H.F <- data1 %>% select(matches("ROI.H.F.Tumor|Protien"))
# write.xlsx(data_F_H.F, "Female_Hyperthermia.xlsx", showNA = FALSE)

#Female Radiation Hyperthermia
data_F_R.H.F <- data1 %>% select(matches("ROI.R.H.F.Tumor|Protien"))
# write.xlsx(data_F_R.H.F, "Female_Radiation_Hypethermia.xlsx", showNA = FALSE)

#### Heatmap of Protein Expressions ----

library(purrr)

# List of datasets
datasets_list <- list(data_M_C.M, data_M_R.M, data_M_H.M, data_M_R.H.M, data_F_C.F, data_F_R.F, data_F_H.F, data_F_R.H.F)

# Purrr: Reduce the list of datasets by left joining based on the "Protien" column
combined_data <- datasets_list %>%
  reduce(left_join, by = "Protien") %>% dplyr::arrange(Protien =="CD3e") #Move the CD3e to the end

# Set "Protien" column as row names
rownames(combined_data) <- combined_data$Protien
combined_data <- combined_data[, -which(colnames(combined_data) == "Protien")]
#Decide to Remove CD3e row
# combined_data <- combined_data[!(row.names(combined_data) %in% c("CD3e")),]

#SKIP -Create sample_info data template to modify from column names of commbined_data
# sample_info <- data.frame(row.names = colnames(combined_data))
# sample_info$Sex <- NA
# sample_info$Experiment <- NA
# sample_info$Sample <- row.names(sample_info)
# row.names(sample_info) <- NULL
# write.xlsx(sample_info, "sample_info.xlsx")
#SKIP-END

#Modify the written sample_info saved by above segment and read it back
sample_info <- read.xlsx("input/sample_info.xlsx")
#Manually add data in columns for plotting

# Set "Sample" column as row names
rownames(sample_info) <- sample_info$Sample
sample_info <- sample_info[, -which(colnames(sample_info) == "Sample")]


#Now I have two datasets to start with; combined_data and the sample_info tables

####PLOT-HEATMAP
library(ComplexHeatmap)
library(circlize)
library(grid)


create_heatmap <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), , drop = FALSE]
  
  # Define the desired order for Experiment
  experiment_order <- c("Control", "Hyperthermia", "Radiation", "Radiation+Hyperthermia")
  sex_order <- c("Male", "Female")
  
  # Convert Sex and Experiment to factors with specified levels
  sample_data$Sex <- factor(sample_data$Sex, levels = sex_order)
  sample_data$Experiment <- factor(sample_data$Experiment, levels = experiment_order)
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#707070", "#FF3079", "#049193", "#5B2897"), experiment_order)
  neuron_colors <- setNames(c("#007DEF", "#F08C00"), sex_order)
  
  # Create a diverging color palette for z-Score Normalized data
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  # max_abs <- max(abs(count_scores))
  max_abs <- max(5, max(abs(count_scores))) #To cap the range of colors at 5. (The range is till 9, more than 5 would have only 21 extreme values)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  ################START-For Grouping by Sex [INPUT_NEEDED]
  column_order <- order(sample_data$Sex, sample_data$Experiment) #Order columns first by Sex, then by Experiment
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
    Sex = sample_data$Sex,
    Experiment = sample_data$Experiment,
    col = list(Sex = neuron_colors, Experiment = slide_colors),
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
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                # row_title_side = c("left", "right"),
                row_names_gp = gpar(fontsize = 20, fontface = "bold"),
                top_annotation = ha_top,
                column_split = sample_data$Sex, #[INPUT_NEEDED] Change between sample_data$Sex or sample_data$Experiment
                # column_order = column_order, #Needed for matching legend order with column order
                column_gap = unit(2, "mm"), #Option
                border = TRUE, #Option
                show_heatmap_legend = FALSE  # Hide default heatmap legend
  )
  
  # Create custom legends
  sex_legend <- Legend(
    labels = sex_order,
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    legend_gp = gpar(fill = neuron_colors),
    column_gap = unit(5, "mm"), row_gap = unit(2, "mm"),
    title = "Sex",
    title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  )
  
  experiment_legend <- Legend(
    labels = experiment_order,
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    legend_gp = gpar(fill = slide_colors),
    row_gap = unit(2, "mm"),
    title = "Experiment",
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
    title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  )
  
  # Combine all legends into a single column
  combined_legend <- packLegend(
    sex_legend,
    experiment_legend,
    expression_legend,
    direction = "horizontal", #vertical or horizontal
    gap = unit(5, "mm")
  )
  
  num_rows <- nrow(count_scores)
  total_height <- unit(min(15, max(10, num_rows * 0.4)), "inch")
  
  # Draw the heatmap with only the combined legend on the left/bottom
  draw(ht, 
       annotation_legend_side = "bottom",
       annotation_legend_list = combined_legend,
       padding = unit(c(2, 20, 2, 10), "mm"),
       height = total_height)
}

plot_and_save_heatmap <- function(normalizedCountsData, sample_data_subset, pathway, output_file) {
  # Calculate height based on number of genes
  # height <- max(12, nrow(normalizedCountsData) * 0.2)  # Adjust the multiplier (0.2) as needed
  height <- 19
  width <- 10.5 #16
  
  pdf(output_file, width = width, height = height)  # Increased width to accommodate legends
  create_heatmap(normalizedCountsData, pathway, sample_data_subset)
  dev.off()
  print(create_heatmap(normalizedCountsData, pathway, sample_data_subset))
  print(paste("Heatmap for", pathway, "pathway saved to", output_file))
}


############################
#Create the Heatmap
############################
#Z score Normalize the combined_data
#Per Gene Scaling of filtered Normalized count ----
# Function to scale each row (gene) (Decide if to scale across all samples or after the filtering below)
scale_rows <- function(x) {
  (x - mean(x)) / sd(x)
}
# Apply scaling to each row
scaled_data <- t(apply(combined_data, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data[is.nan(scaled_data)] <- 0


# Plot and save the heatmap
plot_and_save_heatmap(
  scaled_data, #combined_data or scaled_data
  sample_info, 
  "Protein Expression", 
  "output/0_Protein_Expression.pdf"
)

##SKIP START(For Deciding Scale Range)----
#Get an estimate of score range in scaled data to cap the range
# Basic statistics
cat("Min value:", min(scaled_data), "\n")
cat("Max value:", max(scaled_data), "\n")
cat("Mean:", mean(scaled_data), "\n")
cat("Median:", median(scaled_data), "\n")
cat("Standard deviation:", sd(scaled_data), "\n")

# Quantiles
print(quantile(scaled_data, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)))

# Count of extreme values
cat("Number of values above 5:", sum(abs(scaled_data) > 5), "\n")
cat("Number of values above 3:", sum(abs(scaled_data) > 3), "\n")

hist(scaled_data, breaks = 50, main = "Distribution of Count Scores", xlab = "Score")
##SKIP-END(For Deciding Scale Range)----


####SKIP-START----
# Now second plot in incorporated in the code, but in case needed this is a backup
# For the second plot with grouping by Sex
# 
create_heatmap <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), , drop = FALSE]
  
  # Define the desired order for Experiment and Sex
  experiment_order <- c("Control", "Hyperthermia", "Radiation", "Radiation+Hyperthermia")
  sex_order <- c("Male", "Female")
  
  # Convert Sex and Experiment to factors with specified levels
  sample_data$Sex <- factor(sample_data$Sex, levels = sex_order)
  sample_data$Experiment <- factor(sample_data$Experiment, levels = experiment_order)
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#707070", "#FF3079", "#049193", "#5B2897"), experiment_order)
  neuron_colors <- setNames(c("#007DEF", "#F08C00"), sex_order)
  
  # Create a diverging color palette for z-Score Normalized data
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  max_abs <- max(5, max(abs(count_scores)))
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  # Order the data
  sample_data$OrderGroup <- paste(sample_data$Experiment, sample_data$Sex, sep="_")
  sample_data$OrderGroup <- factor(sample_data$OrderGroup, 
                                   levels = paste(rep(experiment_order, each=2), rep(sex_order, times=4), sep="_"))
  
  column_order <- order(sample_data$OrderGroup)
  
  # Reorder the data
  count_scores <- count_scores[, column_order]
  sample_data <- sample_data[column_order, , drop = FALSE]
  
  # Create top annotation with custom colors
  ha_top <- HeatmapAnnotation(
    Sex = sample_data$Sex,
    Experiment = sample_data$Experiment,
    col = list(Sex = neuron_colors, Experiment = slide_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    show_legend = FALSE  # Hide default legends for annotations
  )
  
  # Create main heatmap
  ht <- Heatmap(count_scores,
                name = paste(pathway, "Scores"),
                col = colorRamp2(breaks, colors),
                column_title = pathway,
                column_title_gp = gpar(fontsize = 20, fontface = "bold", col = "darkblue"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 12),
                top_annotation = ha_top,
                column_split = sample_data$Experiment,
                column_gap = unit(2, "mm"),
                border = TRUE,
                show_heatmap_legend = FALSE  # Hide default heatmap legend
  )
  
  # Create custom legends
  sex_legend <- Legend(
    labels = sex_order,
    legend_gp = gpar(fill = neuron_colors),
    title = "Sex"
  )
  
  experiment_legend <- Legend(
    labels = experiment_order,
    legend_gp = gpar(fill = slide_colors),
    title = "Experiment"
  )
  
  expression_legend <- Legend(
    col_fun = colorRamp2(breaks, colors),
    title = "Expression",
    at = c(-max_abs, 0, max_abs),
    labels = c("Low", "Medium", "High")
  )
  
  # Combine all legends into a single column
  combined_legend <- packLegend(
    sex_legend,
    experiment_legend,
    expression_legend,
    direction = "vertical",
    gap = unit(5, "mm")
  )
  
  num_rows <- nrow(count_scores)
  total_height <- unit(min(15, max(10, num_rows * 0.4)), "inch")
  
  # Draw the heatmap with only the combined legend on the left
  draw(ht, 
       annotation_legend_side = "left",
       annotation_legend_list = combined_legend,
       padding = unit(c(2, 20, 2, 10), "mm"),
       height = total_height)
}
####SKIP-END----