# KM Plot for Mouse survival curves
library(survival)
library(survminer)
library(dplyr)
library(readxl) #To Read
library(openxlsx) #To Write

survival_data <- read.xlsx("Survival_Data.xlsx",
                           sheet = 1,          # can specify sheet number (1) or name ("Sheet1")
                           startRow = 1,       # start reading from first row
                           colNames = TRUE,    # first row contains column names
                           detectDates = TRUE) # automatically detect date columns

# View the first few rows of the data
head(survival_data)



survival_data <- read_excel("Survival_Data.xlsx", sheet = "Sheet 1")
# The column names should be: mouse_id, group, time, status
print(head(survival_data))

write.xlsx(survival_data, "Survival_Data.xlsx",
           overwrite = TRUE,      # overwrites existing file if it exists
           rowNames = FALSE,      # don't include row names
           colNames = TRUE,       # include column names
           asTable = FALSE)       # don't format as Excel table




# [SKIP] Create template dataset----
# Create data with:
# - Mouse ID
# - Group
# - Time (days until event/censoring)
# - Status (1 = event/death, 0 = censored)

set.seed(123) # for reproducibility

# Create data for each group
Tm1cxB6_data <- data.frame(
  mouse_id = 1:6,
  group = "Tm1cxB6",
  time = c(45, 52, 60, 60, 38, 60),  # some censored at 60 days
  status = c(1, 1, 0, 0, 1, 0)
)

LRIG1_KO_data <- data.frame(
  mouse_id = 7:14,
  group = "LRIG1-KO",
  time = c(30, 35, 40, 45, 50, 55, 60, 60),
  status = c(1, 1, 1, 1, 1, 0, 0, 0)
)

VISTA_KO_data <- data.frame(
  mouse_id = 15:19,
  group = "VISTA-KO",
  time = c(25, 30, 35, 60, 60),
  status = c(1, 1, 1, 0, 0)
)

# Combine all data
survival_data <- rbind(Tm1cxB6_data, LRIG1_KO_data, VISTA_KO_data)

# Create survival object----
surv_obj <- Surv(survival_data$time, survival_data$status)
fit <- survfit(surv_obj ~ group, data = survival_data)

# Create Kaplan-Meier plot
km_plot <- ggsurvplot(fit,
                      data = survival_data,
                      pval = TRUE,
                      conf.int = FALSE,
                      risk.table = TRUE,
                      risk.table.height = 0.25,
                      xlab = "Time (days)",
                      ylab = "Survival probability",
                      title = "Kaplan-Meier Survival Curve",
                      palette = c("#909090", "#2E9FDF", "#009669"),
                      legend.title = "Group",
                      legend.labs = c("Tm1cxB6", "LRIG1-KO", "VISTA-KO"),
                      xlim = c(0, max(survival_data$time)),
                      ylim = c(0, 1),
                      break.x.by = 15,     # Set x-axis breaks every 15 days
                      break.y.by = 0.2,    # Set y-axis breaks every 0.2
                      axes.offset = FALSE)

print(km_plot)
