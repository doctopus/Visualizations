library(dplyr)
library(openxlsx)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)


####Source Data----
progression_data_nov20 <- read.xlsx("input/Progressed vs Non-Progressed Metabolites (All 61).xlsx", check.names = FALSE)
##Wrangle Data
progression_nov20 <- progression_data_nov20 %>% rename(Samples = `Sample.Name`)
progression_nov20$Samples <- tools::toTitleCase(progression_nov20$Samples) #Sentence Case Samples
progression_nov20 <- progression_nov20 %>% mutate(Samples = ifelse(row_number() == 1, "PFS", Samples)) #First row name as PFS from "Progressed at 9 Months (0=yes, 1=no)"
progression_nov20 <- progression_nov20 %>% mutate(across(-c(Samples), ~as.numeric(as.character(.)))) # mutate_all(~as.numeric(as.character(.)))
#Make metabData
metabData_nov20 <- progression_nov20[-1, ]
rownames(metabData_nov20) <- metabData_nov20$Samples
metabData_nov20 <- metabData_nov20%>% select(-Samples)
#Make metabGroups
metabGroups_nov20 <- progression_nov20 %>% t() %>% as.data.frame()
colnames(metabGroups_nov20) <- as.character(unlist(metabGroups_nov20[1, ]))
metabGroups_nov20 <- metabGroups_nov20[-1, ]
metabGroups_nov20 <- metabGroups_nov20 %>% select(PFS) #Keep PFS Column Only
metabGroups_nov20 <- metabGroups_nov20 %>% mutate(PFS = as.numeric(PFS)) %>% mutate(PFS= as.factor(PFS))
summary(metabGroups_nov20)


#PLOT Functions
create_heatmap_pfs <- function(count_scores, pathway, sample_data) {
  # Ensure sample_data is in the same order as count_scores columns
  sample_data <- sample_data[colnames(count_scores), , drop = FALSE]
  
  # Define the desired order for Experiment
  experiment_order <- c(0, 1) #,2 (Progressed at 9 months (0=Yes, 1=No))
  # sex_order <- c("Male", "Female")
  
  # Convert Sex and Experiment to factors with specified levels
  # sample_data$Sex <- factor(sample_data$Sex, levels = sex_order)
  sample_data$PFS <- factor(sample_data$PFS, levels = experiment_order)
  
  # Create color vectors for each annotation
  slide_colors <- setNames(c("#F59F00", "#37B24D"), experiment_order)#, "#F4D03F", "#E67E22" "#049193", "#5B2897"
  # neuron_colors <- setNames(c("#007DEF", "#F08C00"), sex_order)
  
  # Create a diverging color palette for z-Score Normalized data
  # colors <- colorRampPalette(c("#377eb8", "white", "#e41a1c"))(101)
  colors <- colorRampPalette(c("blue", "white", "red"))(101)
  # max_abs <- max(abs(count_scores))
  # max_abs <- max(abs(count_scores[is.finite(count_scores)]), na.rm = TRUE)
  max_abs <- 3
  # max_abs <- max(5, max(abs(count_scores))) #To cap the range of colors at 5. (The range is till 9, more than 5 would have only 21 extreme values)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  # breaks <- seq(-3, 12, length.out = 101)
  
  ################START-For Grouping by Sex [INPUT_NEEDED]
  column_order <- order(sample_data$PFS)
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
    PFS = sample_data$PFS,
    # col = list(Sex = neuron_colors, Experiment = slide_colors),
    col = list(PFS = slide_colors),
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
                row_names_gp = gpar(fontsize = 20, fontface = "plain"),
                top_annotation = ha_top,
                clustering_distance_rows = "euclidean", #Default euclidean; other manhattan
                clustering_method_rows = "ward.D", #default is complete, ward is renamed to ward.D and there is ward.D2
                clustering_distance_columns = "euclidean", #Default euclidean; other manhattan
                clustering_method_columns = "ward.D",
                column_split = sample_data$PFS, #[INPUT_NEEDED] Change between sample_data$Sex or sample_data$Experiment
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
    # labels = experiment_order,
    labels = rev(c("Progressed at 9 months", "Not Progressed at 9 months")), #Reverse names (need to reverse legend_gp too)
    labels_gp = gpar(fontsize = 20, fontface='bold'),#Increase size of labels
    legend_gp = gpar(fill = rev(slide_colors)), #Reverse associated colors too (need to have labels rerversed too)
    row_gap = unit(2, "mm"),
    title = "PFS",
    title_gp = gpar(fontsize = 22, fontface='bold') #Increase size of legend label
  )
  
  expression_legend <- Legend(
    col_fun = colorRamp2(breaks, colors),
    at = c(-max_abs, 0, max_abs),
    # labels = c("Low", "Medium", "High"),
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
       # x = unit(1, "cm"), y = unit(1, "cm"), just = c("right", "top"),
       annotation_legend_list = combined_legend,
       padding = unit(c(2, 20, 2, 80), "mm"),#x, left, y, right
       height = total_height)
}
plot_and_save_heatmap_pfs <- function(normalizedCountsData, sample_data_subset, plot_title, base_filename) {
  # Calculate height based on number of genes
  # height <- max(12, nrow(normalizedCountsData) * 0.2)  # Adjust the multiplier (0.2) as needed
  height <- 20
  width <- 15 #16 or 10.5
  datetime <- format(Sys.time(), "%Y%m%d.%H%M%S")
  fileName <- paste(datetime, base_filename, sep="_")
  
  pdf(fileName, width = width, height = height)  # Increased width to accommodate legends
  create_heatmap_pfs(normalizedCountsData, plot_title, sample_data_subset)
  dev.off()
  print(create_heatmap_pfs(normalizedCountsData, plot_title, sample_data_subset))
  print(paste("Heatmap for", plot_title, "saved as", fileName))
}

#Z score Normalize the gbm_metab_info
#Per Gene Scaling of filtered Normalized count  Data
# Function to scale each row (gene) (Decide if to scale across all samples or after the filtering below)
scale_rows <- function(x) {
  (x - mean(x)) / sd(x)
}
# Apply scaling to each row
scaled_data <- t(apply(metabData_nov20, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data[is.nan(scaled_data)] <- 0

##START Alternative Scaling (Not-Used)----
#Log Scaling the whole data as an alterrnative approach
#Version1: Log2 Function if no negative values
log2_transform <- function(x) {
  log2(x + 1) # Add a small constant to avoid log(0); (1 is commonly used)
}

# #Version2: Log2 Function: Alternative versions if need to handle negative values:
# log2_transform <- function(x) {
#   sign(x) * log2(abs(x) + 1) # For negative values: sign(x) * log2(abs(x) + 1)
# }
# 
# #Version3: Log2 Functioon if need to handle zeros differently:
# log2_transform <- function(x) {
#   min_nonzero <- min(x[x > 0], na.rm = TRUE) # Replace zeros with minimum non-zero value divided by 2
#   x[x == 0] <- min_nonzero/2
#   log2(x)
# }

# Apply the desired log2 transformation to the entire dataset
log2_data <- as.matrix(log2_transform(metabData_nov20))

# Check for any infinite values or NAs if needed
sum(is.infinite(log2_data))
sum(is.na(log2_data))
##END Alternative Scaling (Not-Used)----

###START Plot only Top Contributing----
#Prepare Only the top contributing metabolite by filtering from the list generated from the PCA top loadings
#Source the combined group of top metabolites as combined_metabolite list. 
metabData_nov20_top <- metabData_nov20[rownames(metabData_nov20) %in% combined_top_metabolites, ]
#Metabolitewise scale the data using scale function
# Apply scaling to each row
scaled_data_top <- t(apply(metabData_nov20_top, 1, scale_rows))
# In case any rows have standard deviation of 0 (constant values),
# they will result in NaN. We can replace these with 0:
scaled_data_top[is.nan(scaled_data_top)] <- 0
###END Plot only Top Contributing----

#### Plot and save the heatmap
plot_and_save_heatmap_pfs(
  scaled_data, #log2_data or scaled_data or only the top metabolites (scaled_data_top)
  metabGroups_nov20, #Select the whole dataset of 60 metabolites or the top 20 contributing metabolites me
  "Metabolites by Progression", 
  "20250605_76_PFS Metabolites Nov20 zScaled data -rows cols euclidian ward -scale till 3 - size20x16.pdf"
)

