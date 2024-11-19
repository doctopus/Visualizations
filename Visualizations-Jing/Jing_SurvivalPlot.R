library(ggplot2)
library(dplyr)
library(tidyr)
data_survival <- read.csv("input/cytokine/Survivial data table_10 pts on clinical trial (using survival data updated to Oct 9-2024).csv")

# Prepare the data
data_surv <- data.frame(
  `Subject #` = c(1,2,3,4,8,9,12,18,25,28),
  ttProg = c(15,2.3,4.3,NA,21.2,2.6,3.6,6.6,7.6,NA),
  ttDeath = c(30.5,14.1,9.7,2.2,29.4,3.7,5.1,NA,18.4,NA),
  ttLfp = c(26,13.5,9.7,2.2,14.9,3.7,4.6,26.3,18.4,15.4),
  Sex = c(1,1,1,1,1,1,1,0,1,1)  # 1=Male, 0=Female
)


data_surv <- data_survival %>%
  rename_with(~ gsub("\\.\\.|\\.+", "_", .x)) %>%  # Replace .. and . with _
  rename(Subject = Subject_) %>%                    # Fix Subject column name
  select(Subject, ttProg, ttDeath, ttLfp) %>%      # Select needed columns
  mutate(Sex = ifelse(Subject == 18, 0, 1)) %>%       # Add Sex column
  filter(!Subject %in% c(1,3,8))

# Define status attributes with stroke parameter
status_info <- data.frame(
  status = c("Death", "Last Followup", "Progression"),
  symbol = c(4, 1, 18),
  color = c("#495057", "darkgreen", "#E67700"),
  size = c(4, 5, 7),  # Reduced size for Death from 8 to 6
  stroke = c(2, 2, 1)
)

ggplot() +
  # Patient timeline
  geom_segment(data = data_surv,
               aes(x = 0, xend = ttLfp, 
                   y = reorder(factor(Subject), ttLfp),
                   yend = reorder(factor(Subject), ttLfp),
                   color = as.factor(Sex)),
               size = 2.5) +
  
  # Status markers
  geom_point(data = subset(data_surv, !is.na(ttProg)), 
             aes(x = ttProg, 
                 y = reorder(factor(Subject), ttLfp),
                 shape = "Progression"),
             color = "#E67700",
             size = status_info$size[3],
             stroke = status_info$stroke[3]) +

  
  geom_point(data = data_surv, 
             aes(x = ttLfp, 
                 y = reorder(factor(Subject), ttLfp),
                 shape = "Last Followup"),
             color = "darkgreen",
             size = status_info$size[2],
             stroke = status_info$stroke[2]) +
  
  geom_point(data = subset(data_surv, !is.na(ttDeath)), 
             aes(x = ttDeath, 
                 y = reorder(factor(Subject), ttLfp),
                 shape = "Death"),
             color = "#495057",
             size = status_info$size[1],
             stroke = status_info$stroke[1]) +
  
  # Scales
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  
  scale_shape_manual("Patient Status",
                     values = status_info$symbol,
                     breaks = c("Death", "Last Followup", "Progression")) +
  
  scale_color_manual("Sex",
                     values = c("0" = "#BE4BDB", "1" = "#4C6EF5"),
                     labels = c("Female", "Male")) +
  
  # Labels and theme
  labs(x = "Months since surgery",
       y = "Patient ID",
       title = "Patient Survival and Progression") +
  
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"), 
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(2, "lines")
  ) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2))
