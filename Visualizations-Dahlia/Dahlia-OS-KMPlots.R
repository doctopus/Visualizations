library(tidyverse)
library(survival)
library(survminer)
library(lubridate)
library(openxlsx)
os_data <- read.xlsx("input/Characteristics of patients Glioblastoma & Microbiome Dataset_NR_5-9-21.xlsx", check.names = FALSE, detectDates = TRUE)
os_data_selected <- os_data %>% select(c("X1", "X2", "X13", "X14", "X15", "X16"))

colnames(os_data_selected) <- as.character(unlist(os_data_selected[1, ]))
os_data_selected <- os_data_selected[-1, ]

# First, identify and examine the problematic values
os_data_selected %>%
  filter(!is.na(DOD) & !str_detect(DOD, "^\\d{4}-\\d{2}-\\d{2}$")) %>%
  select(LABNUM, DOS, DOD, EXPIRED) %>%
  print()

# Clean the DOD column
os_data_selected_clean <- os_data_selected %>%
  mutate(
    # Create a cleaned DOD column
    DOD_clean = case_when(
      # Keep valid date format (YYYY-MM-DD)
      str_detect(DOD, "^\\d{4}-\\d{2}-\\d{2}$") ~ DOD,
      # Extract date from "2024-06-02-tx care" format
      str_detect(DOD, "^\\d{4}-\\d{2}-\\d{2}-") ~ str_extract(DOD, "^\\d{4}-\\d{2}-\\d{2}"),
      # Handle other non-date values
      DOD %in% c("<18", "transferred care") ~ NA_character_,
      # Keep original NA values as NA
      is.na(DOD) ~ NA_character_,
      # Any other unexpected format becomes NA
      TRUE ~ NA_character_
    ),
    # Also update EXPIRED status for non-death cases
    EXPIRED_clean = case_when(
      # If DOD indicates transfer of care, patient should be censored
      DOD %in% c("transferred care", "2024-06-02-tx care") ~ "N",
      # If DOD is "<18" (meaning died before 18 days?), keep as death
      DOD == "<18" ~ "Y",
      # Otherwise keep original EXPIRED value
      TRUE ~ EXPIRED
    )
  ) %>%
  # Replace the original columns with cleaned versions
  select(-DOD, -EXPIRED) %>%
  rename(DOD = DOD_clean, EXPIRED = EXPIRED_clean)

# Verify the cleaning
print("Summary of cleaned DOD column:")
os_data_selected_clean %>%
  summarise(
    total_rows = n(),
    valid_dates = sum(!is.na(DOD)),
    na_values = sum(is.na(DOD)),
    expired_y = sum(EXPIRED == "Y", na.rm = TRUE),
    expired_n = sum(EXPIRED == "N", na.rm = TRUE)
  )


# Find the latest date, removing NA values
max_death_date <- max(as.Date(os_data_selected_clean$DOD), na.rm = TRUE)
print(max_death_date)




# Calculate overall survival time and create survival object
os_data_survival <- os_data_selected_clean %>%
  # Convert date columns to Date type if they're not already
  mutate(
    DOS = as.Date(DOS),
    DOD = as.Date(DOD),
    DOB = as.Date(DOB)
  ) %>%
  # Calculate survival time in days from surgery to death or last follow-up
  mutate(
    # For patients who died (EXPIRED = 1), survival time is DOS to DOD
    # For patients who are alive (EXPIRED = 0), survival time is DOS to last known date
    survival_time = case_when(
      EXPIRED == "Y" & !is.na(DOD) ~ as.numeric(DOD - DOS),
      EXPIRED == "Y" & is.na(DOD) ~ NA_real_, #Death recorded but no date
      EXPIRED == "N" ~ as.numeric(max_death_date - DOS), # Alive/Censored: Use max death date for censoring or use a specific censoring date
      is.na(EXPIRED) ~ as.numeric(max_death_date -DOS), # Missing status - treat as censored
      TRUE ~ NA_real_
    ),
    # Create event indicator (1 = death, 0 = censored)
    event = case_when(
      EXPIRED == "Y" ~ 1,
      EXPIRED == "N" ~ 0,
      is.na(EXPIRED) ~ 0, # Missing status - treat as censored
      TRUE ~ NA_real_
    )
  ) %>%
  # Remove rows with missing or negative survival times
  filter(!is.na(survival_time) & survival_time >= 0) %>%
  # Remove duplicate rows (keep first occurrence of each duplicate)
  distinct() %>%
  # Keep only rows where LABNUM is in the rownames of metabGroups_nov20
  filter(LABNUM %in% rownames(metabGroups_nov20))

# Check the results
print(paste("Number of rows after deduplication and filtering:", nrow(os_data_survival)))
print(paste("Number of unique LABNUM values:", length(unique(os_data_survival$LABNUM))))

# Check the distribution of events
event_summary <- os_data_survival %>%
  count(event) %>%
  mutate(event_type = ifelse(event == 1, "Death", "Censored"))
print("Event distribution:")
print(event_summary)


# Find row names in metabGroups_nov20 that are not in os_data_survival
missing_labnums <- setdiff(rownames(metabGroups_nov20), os_data_survival$LABNUM)

print(paste("Number of missing LABNUM values:", length(missing_labnums)))
print("Missing LABNUM values:")
print(missing_labnums)

# Alternative way to see this information
print("\nComparison:")
print(paste("Rows in metabGroups_nov20:", nrow(metabGroups_nov20)))
print(paste("Rows in os_data_survival:", nrow(os_data_survival)))
print(paste("Missing from os_data_survival:", length(missing_labnums)))

# Show the missing LABNUMs in a more readable format
if(length(missing_labnums) > 0) {
  missing_df <- data.frame(Missing_LABNUM = missing_labnums)
  print(missing_df)
}






# Create survival object
surv_object <- Surv(time = os_data_survival$survival_time, 
                    event = os_data_survival$event)

# Fit Kaplan-Meier survival curve
# km_fit <- survfit(surv_object ~ 1, data = os_data_survival)
km_fit <- survfit(surv_object ~ GENDER, data = os_data_survival)

# Plot Kaplan-Meier curve using survminer
km_plot <- ggsurvplot(
  km_fit,
  data = os_data_survival,
  conf.int = FALSE,
  pval = TRUE,  # Show p-value for log-rank test comparing groups
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("blue", "red"),  # Different colors for M and F OR "blue",
  title = "Kaplan-Meier Survival Curve by Gender",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  legend.title = "Gender",
  legend.labs = c("Female", "Male")  # Customize labels (F will be first, M second)
)

# Display the plot
print(km_plot)

# Print summary statistics
print(km_fit)

# Perform log-rank test to compare survival between groups
logrank_test <- survdiff(surv_object ~ GENDER, data = os_data_survival)
print("Log-rank test results:")
print(logrank_test)

# Get median survival times for each group
median_survival <- surv_median(km_fit)
print("Median survival times by gender:")
print(median_survival)


##### Add Metabolite Data from metabData_nov20 to the os_data_survival----
metabolites_to_add <- c("Choline ", "Creatinine ", "Phenylalanine ", "Uric Acid ", 
                        "Pseudouridine ", "Phenylacetylglutamine ", "Tyrosine ", 
                        "Valine ", "Isoleucine ", "4-Hydroxyl-Phenyllactic Acid ")

# Check which metabolites are actually present in metabData_nov20
available_metabolites <- intersect(metabolites_to_add, rownames(metabData_nov20))
missing_metabolites <- setdiff(metabolites_to_add, rownames(metabData_nov20))

print(paste("Available metabolites:", length(available_metabolites)))
print("Available:", available_metabolites)
if(length(missing_metabolites) > 0) {
  print("Missing metabolites:")
  print(missing_metabolites)
}

# Extract and transpose the metabolite data
metab_subset <- metabData_nov20[available_metabolites, , drop = FALSE]

# Transpose so that LABNUM becomes rows and metabolites become columns
metab_transposed <- as.data.frame(t(metab_subset))

# Convert rownames to a column for merging
metab_transposed$LABNUM <- rownames(metab_transposed)

# Merge with os_data_survival
os_data_survival_expanded <- os_data_survival %>%
  left_join(metab_transposed, by = "LABNUM")

# Check the results
print(paste("Original os_data_survival columns:", ncol(os_data_survival)))
print(paste("Expanded os_data_survival columns:", ncol(os_data_survival_expanded)))
print(paste("Added metabolite columns:", length(available_metabolites)))

# Show the structure of the expanded data
print("First few rows of expanded data:")
print(head(os_data_survival_expanded))

# Check for any missing values in the metabolite columns
missing_summary <- os_data_survival_expanded %>%
  select(all_of(available_metabolites)) %>%
  summarise_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "Metabolite", values_to = "Missing_Count")

print("Missing values summary for metabolites:")
print(missing_summary)


#### Function to plot KM plot of Upper and Lower Quartile of Metabolites----
# Function to create quartile-based KM plots for a metabolite
# Alternative approach - create survival object within survfit
create_quartile_km_plot <- function(data, metabolite_name, title_prefix = "Kaplan-Meier Survival Curve") {
  
  # Check if metabolite exists in data
  if (!metabolite_name %in% colnames(data)) {
    stop(paste("Metabolite", metabolite_name, "not found in data"))
  }
  
  # Remove rows with missing metabolite values
  data_clean <- data %>%
    filter(!is.na(!!sym(metabolite_name)))
  
  if (nrow(data_clean) == 0) {
    stop(paste("No valid data for metabolite", metabolite_name))
  }
  
  # Calculate quartiles
  quartiles <- quantile(data_clean[[metabolite_name]], probs = c(0.25, 0.75), na.rm = TRUE)
  q1_value <- quartiles[1]
  q3_value <- quartiles[2]
  
  # Create quartile groups (top quartile vs bottom quartile)
  data_quartiles <- data_clean %>%
    mutate(
      quartile_group = case_when(
        !!sym(metabolite_name) <= q1_value ~ "Q1 (Lowest)",
        !!sym(metabolite_name) >= q3_value ~ "Q4 (Highest)",
        TRUE ~ "Q2-Q3 (Middle)"
      )
    ) %>%
    # Keep only Q1 and Q4 for comparison
    filter(quartile_group %in% c("Q1 (Lowest)", "Q4 (Highest)")) %>%
    # Convert to factor for proper ordering
    mutate(quartile_group = factor(quartile_group, levels = c("Q1 (Lowest)", "Q4 (Highest)")))
  
  # Check if we have data in both quartiles
  if (length(unique(data_quartiles$quartile_group)) < 2) {
    warning(paste("Insufficient data for quartile comparison for", metabolite_name))
    return(NULL)
  }
  
  # Fit Kaplan-Meier survival curve using formula directly
  km_fit <- survfit(Surv(survival_time, event) ~ quartile_group, data = data_quartiles)
  
  # Create the plot
  km_plot <- ggsurvplot(
    km_fit,
    data = data_quartiles,
    conf.int = FALSE,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_minimal(),
    palette = c("blue", "red"),
    title = paste(title_prefix, "by", metabolite_name, "Quartiles"),
    xlab = "Time (Days)",
    ylab = "Survival Probability",
    legend.title = paste(metabolite_name, "Level"),
    legend.labs = c("Q1 (Lowest 25%)", "Q4 (Highest 25%)")
  )
  
  # Print summary statistics
  cat("\n=== Summary for", metabolite_name, "===\n")
  cat("Q1 (25th percentile):", round(q1_value, 3), "\n")
  cat("Q3 (75th percentile):", round(q3_value, 3), "\n")
  cat("Sample sizes:\n")
  print(table(data_quartiles$quartile_group))
  cat("Events by quartile:\n")
  print(data_quartiles %>% count(quartile_group, event) %>% pivot_wider(names_from = event, values_from = n, values_fill = 0))
  
  return(list(plot = km_plot, data = data_quartiles, fit = km_fit))
}
# Or use the alternative approach
choline_result <- create_quartile_km_plot(os_data_survival_expanded, "Choline ")
if (!is.null(choline_result_alt)) {
  print(choline_result_alt$plot)
}




colnames(os_data_survival_expanded)
## USage of the Function
# Example 1: Create KM plot for Choline
choline_result <- create_quartile_km_plot(os_data_survival_expanded, "Choline ")
if (!is.null(choline_result)) {
  print(choline_result$plot)
}

# Example 2: Create KM plot for Creatinine
creatinine_result <- create_quartile_km_plot(os_data_survival_expanded, "4-Hydroxyl-Phenyllactic Acid ")
if (!is.null(creatinine_result)) {
  print(creatinine_result$plot)
}

# Example 3: Loop through all metabolites of interest
metabolites_to_analyze <- c("Choline ", "Creatinine", "Phenylalanine", "Uric Acid", 
                            "Pseudouridine", "Phenylacetylglutamine", "Tyrosine", 
                            "Valine", "Isoleucine", "4-Hydroxyl-Phenyllactic Acid")

# Store results for all metabolites
km_results <- list()

for (metabolite in metabolites_to_analyze) {
  if (metabolite %in% colnames(os_data_survival_expanded)) {
    cat("\n\nAnalyzing", metabolite, "...\n")
    result <- create_quartile_km_plot(os_data_survival_expanded, metabolite)
    if (!is.null(result)) {
      km_results[[metabolite]] <- result
      print(result$plot)
    }
  } else {
    cat("Metabolite", metabolite, "not found in data\n")
  }
}






### Create Upper and Lower Quartile OS metabgroups to use in heatmap as in PFS heatmap----
# Calculate quartiles of overall survival (survival_time)
survival_quartiles <- quantile(os_data_survival_expanded$survival_time, 
                               probs = c(0.25, 0.75), na.rm = TRUE)
q1_survival <- survival_quartiles[1]
q3_survival <- survival_quartiles[2]

print(paste("Q1 (25th percentile) survival time:", round(q1_survival, 1), "days"))
print(paste("Q3 (75th percentile) survival time:", round(q3_survival, 1), "days"))

# Create the quartile groups and filter to keep only Q1 and Q4
os_quartile_data <- os_data_survival_expanded %>%
  mutate(
    os_quartile_group = case_when(
      survival_time <= q1_survival ~ "Q1_Lowest",
      survival_time >= q3_survival ~ "Q4_Highest",
      TRUE ~ "Q2Q3_Middle"
    )
  ) %>%
  # Keep only Q1 and Q4 patients
  filter(os_quartile_group %in% c("Q1_Lowest", "Q4_Highest")) %>%
  # Create binary indicator: 1 for top quartile (Q4), 0 for bottom quartile (Q1)
  mutate(
    OS_Quartile = ifelse(os_quartile_group == "Q4_Highest", 1, 0)
  ) %>%
  # Select only LABNUM and OS_Quartile columns
  select(LABNUM, OS_Quartile)

# Create the final dataframe with LABNUM as row names
metabGroups_OS_Quartile <- os_quartile_data %>%
  column_to_rownames("LABNUM")

# Display summary
print("Summary of metabGroups_OS_Quartile:")
print(paste("Total patients:", nrow(metabGroups_OS_Quartile)))
print(paste("Bottom quartile (OS_Quartile = 0):", sum(metabGroups_OS_Quartile$OS_Quartile == 0)))
print(paste("Top quartile (OS_Quartile = 1):", sum(metabGroups_OS_Quartile$OS_Quartile == 1)))

# Show the structure
print("Structure of metabGroups_OS_Quartile:")
str(metabGroups_OS_Quartile)

# Show first few rows
print("First few rows:")
print(head(metabGroups_OS_Quartile))

# Show the row names (LABNUMs)
print("Row names (LABNUMs):")
print(head(rownames(metabGroups_OS_Quartile)))




## Create metabData_OS_Quartile to keep only the patients that are in the OS quartile data in metabGroups_OS_Quartile
# Get the LABNUMs from metabGroups_OS_Quartile rownames
os_quartile_labnums <- rownames(metabGroups_OS_Quartile)

# Check overlap between metabData_nov20 columns and metabGroups_OS_Quartile rownames
common_labnums <- intersect(colnames(metabData_nov20), os_quartile_labnums)
missing_in_metabData <- setdiff(os_quartile_labnums, colnames(metabData_nov20))
extra_in_metabData <- setdiff(colnames(metabData_nov20), os_quartile_labnums)

print(paste("Total LABNUMs in metabGroups_OS_Quartile:", length(os_quartile_labnums)))
print(paste("Total columns in metabData_nov20:", ncol(metabData_nov20)))
print(paste("Common LABNUMs:", length(common_labnums)))
print(paste("LABNUMs in OS_Quartile but not in metabData_nov20:", length(missing_in_metabData)))
print(paste("LABNUMs in metabData_nov20 but not in OS_Quartile:", length(extra_in_metabData)))

# Create metabData_OS_Quartile by selecting only the relevant columns
metabData_OS_Quartile <- metabData_nov20[, common_labnums, drop = FALSE]

# Display summary
print(paste("Dimensions of metabData_OS_Quartile:", nrow(metabData_OS_Quartile), "x", ncol(metabData_OS_Quartile)))

# Show the structure
print("Structure of metabData_OS_Quartile:")
str(metabData_OS_Quartile)

# Show first few rows and columns
print("First few rows and columns:")
print(metabData_OS_Quartile[1:5, 1:5])

# Verify that all column names are in metabGroups_OS_Quartile rownames
all_match <- all(colnames(metabData_OS_Quartile) %in% rownames(metabGroups_OS_Quartile))
print(paste("All column names match metabGroups_OS_Quartile rownames:", all_match))

## Now I have metabGroups_OS_Quartile and metabData_OS_Quartile to be used in PFS Heatmap script (Dahlia-PFS_Heatmap.R)


