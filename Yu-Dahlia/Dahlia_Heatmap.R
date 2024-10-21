#The following code for plotting a heatmap from already normalized data of protein expression
#Rows have Gene/Protein names, columns have sample names
library(dplyr)
library(openxlsx)
library(tidyverse)
getwd()
setwd("/Users/i/Dropbox/Clinic3.0/Developer/RStudio/Visualizations/Yu-Dahlia")


########Real Data from Dahlia---- input: source csv output: gbm_metab_info, gbm_sample_info ----
# gbm_metab <- read.csv("input/GBM_Metabolite_SourceData.csv", header = TRUE, check.names = FALSE)
gbm_metab <- read.csv("input/GBM + Healthy Controls for ClustVis.csv", header = TRUE, check.names = FALSE)
progression_data <- read.csv("input/P vs NP vs HC.csv", header = TRUE, check.names = FALSE)
# Transpose the gbm_metab data
transposed_gbm_metab <- gbm_metab %>% t()



transposed_gbm_metab <- cbind(Samples = rownames(transposed_gbm_metab), transposed_gbm_metab)

# Reset the row names to NULL
rownames(transposed_gbm_metab) <- NULL


colnames(transposed_gbm_metab) <- as.character(unlist(transposed_gbm_metab[1, ]))
transposed_gbm_metab <- transposed_gbm_metab[-1, ]


colnames(transposed_gbm_metab)[1] <- "Samples"

transposed_gbm_metab <- as.data.frame(transposed_gbm_metab) %>% 
  mutate(GBM = as.numeric(GBM))

#Keep only first two columns
gbm_sample_info <- transposed_gbm_metab %>% select(1:2) %>% mutate(GBM =as.factor(GBM))

#Make the Sample column as the rownames and remove that column
rownames(gbm_sample_info) <- gbm_sample_info$Samples
# Remove the first column from the dataset
# gbm_sample_info <- gbm_sample_info[, -1] #Base R method does not work
# Alternative way Remove the first column from the dataset using select
gbm_sample_info <- gbm_sample_info %>% select(-Samples)

gbm_sample_info <- as.data.frame(gbm_sample_info)

# colnames(gbm_sample_info)[1] <- "GBM"



#Prepare gbm_metab_info
gbm_metab_info <- gbm_metab[-1, ]
rownames(gbm_metab_info) <- gbm_metab_info[, 1]

# Remove the first column from the dataset
gbm_metab_info <- gbm_metab_info[, -1]

########Sample Data---- input: sample data csv, output: data_metab_final, sample_metab ----
data_metab <- read.xlsx("input/Metabolomics_Sample_Dataset.xlsx", sheet="Sheet1", colNames= TRUE)
data_metab <- data_metab %>% select(-`0=.Progression;.1=no.progression`)
library(tidyr)
data_metab_new <- pivot_longer(data_metab, -Sample.Name, names_to = "Metabolite", values_to = "Value")
data_metab_final <- pivot_wider(data_metab_new, names_from = Sample.Name, values_from = Value)


data_metab_final <- as.data.frame(data_metab_final)
rownames(data_metab_final) <- data_metab_final$Metabolite
data_metab_final <- data_metab_final[, -which(colnames(data_metab_final) == "Metabolite")]

#SKIP -Create sample_info data template to modify from column names of combined_data
sample_metab <- data.frame(row.names = colnames(data_metab_final))
sample_metab$Sex <- NA
sample_metab$Experiment <- NA
# sample_metab$Progression <- data_metab$`0=.Progression;.1=no.progression`
sample_metab$Sample.Name <- row.names(sample_metab)
# row.names(sample_metab) <- NULL
write.xlsx(sample_metab, "input/sample_metab.xlsx")
#SKIP-END

#Modify the written sample_info saved by above segment and read it back
sample_metab <- read.xlsx("input/sample_metab.xlsx")
#Manually add data in columns for plotting

# Set "Sample" column as row names
rownames(sample_metab) <- sample_metab$Sample.Name
sample_metab <- sample_metab[, -which(colnames(sample_metab) == "Sample.Name")]

#Now I have two datasets to start with; data_metab_final and the sample_metab tables

####PLOT-HEATMAP
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
  slide_colors <- setNames(c("#707070", "#FF3079"), experiment_order)#, "#049193", "#5B2897"
  # neuron_colors <- setNames(c("#007DEF", "#F08C00"), sex_order)
  
  # Create a diverging color palette for z-Score Normalized data
  # colors <- colorRampPalette(c("#377eb8", "white", "#e41a1c"))(101)
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  # max_abs <- max(abs(count_scores))
  max_abs <- max(5, max(abs(count_scores))) #To cap the range of colors at 5. (The range is till 9, more than 5 would have only 21 extreme values)
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
                column_split = sample_data$GBM, #[INPUT_NEEDED] Change between sample_data$Sex or sample_data$Experiment
                column_order = column_order, #Needed for matching legend order with column order
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
    labels = experiment_order,
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    legend_gp = gpar(fill = slide_colors),
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
    title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  )
  
  # Combine all legends into a single column
  combined_legend <- packLegend(
    # sex_legend,
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
       padding = unit(c(2, 10, 2, 60), "mm"),#x, left, y, right
       height = total_height)
}
plot_and_save_heatmap_gbm <- function(normalizedCountsData, sample_data_subset, pathway, output_file) {
  # Calculate height based on number of genes
  # height <- max(12, nrow(normalizedCountsData) * 0.2)  # Adjust the multiplier (0.2) as needed
  height <- 19
  width <- 11 #16 or 10.5
  
  pdf(output_file, width = width, height = height)  # Increased width to accommodate legends
  create_heatmap_gbm(normalizedCountsData, pathway, sample_data_subset)
  dev.off()
  print(create_heatmap_gbm(normalizedCountsData, pathway, sample_data_subset))
  print(paste("Heatmap for", pathway, "pathway saved to", output_file))
}

#Z score Normalize the gbm_metab_info
#Per Gene Scaling of filtered Normalized count ----
# Function to scale each row (gene) (Decide if to scale across all samples or after the filtering below)
scale_rows <- function(x) {
  (x - mean(x)) / sd(x)
}
# Apply scaling to each row
scaled_data <- t(apply(gbm_metab_info, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data[is.nan(scaled_data)] <- 0


# Plot and save the heatmap
plot_and_save_heatmap_gbm(
  scaled_data, #data_metab_final or scaled_data
  gbm_sample_info, 
  "GBM Metabolites", 
  "25_GBM_Metabolite_Clustered_No Columns_No Scale Labels_Max Abs Range Full -Labels - Legend bottom -Column Split -Capped.pdf"
)


#Plotting the sample data----
create_heatmap <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), , drop = FALSE]
  
  # Define the desired order for Experiment
  experiment_order <- c("Progression", "No-Progression") #
  sex_order <- c("Male", "Female")
  
  # Convert Sex and Experiment to factors with specified levels
  sample_data$Sex <- factor(sample_data$Sex, levels = sex_order)
  sample_data$Experiment <- factor(sample_data$Experiment, levels = experiment_order)
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#707070", "#FF3079"), experiment_order)#, "#049193", "#5B2897"
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



#Create the Heatmap

#Z score Normalize the combined_data
#Per Gene Scaling of filtered Normalized count --
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
  data_metab_final, #data_metab_final or scaled_data
  sample_metab, 
  "Protein Expression", 
  "0_Protein_Expression.pdf"
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
dev.off()
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