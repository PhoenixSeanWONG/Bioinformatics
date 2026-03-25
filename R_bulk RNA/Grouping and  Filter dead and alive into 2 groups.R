## Grouping
```{r}
pfs_threshold <- 200  # days
os_threshold <- 400   # days

# Defining Grouping Functions
define_pfs_group <- function(pfs, pfs_threshold) {
  if (pfs >= pfs_threshold) {
    return("PFS_Long")
  } else {
    return("PFS_Short")
  }
}

define_os_group <- function(os, os_threshold) {
  if (os >= os_threshold) {
    return("OS_Long")
  } else {
    return("OS_Short")
  }
}
```

```{r}
combined_for_grouping$PFS_Group <- ifelse(
  combined_for_grouping$PFS.Event == 1, 
  mapply(define_pfs_group, 
         combined_for_grouping$PFS, 
         MoreArgs = list(pfs_threshold = pfs_threshold)), 
  NA
)

head(combined_for_grouping)
table(combined_for_grouping$PFS_Group, useNA = "ifany")
```

```{r}
combined_for_grouping$OS_Group <- ifelse(
  combined_for_grouping$OS.Event == 1, 
  mapply(define_os_group, 
         combined_for_grouping$OS, 
         MoreArgs = list(os_threshold = os_threshold)), 
  NA
)

head(combined_for_grouping)
table(combined_for_grouping$OS_Group, useNA = "ifany")
```

### Integration
```{r}
combined_gide <- combined_for_grouping %>%
  filter(PMID == "gide") %>%
  inner_join(ptmeta_subtype_purity, by = c("Patient_ID" = "NewPatientID"))
```

```{r}
combined_riaz <- combined_for_grouping %>%
  filter(PMID == "riaz") %>%
  inner_join(ptmeta_subtype_purity, by = c("Patient_ID" = "patient_id(pt)"))
```

```{r}
combined_for_grouping_withHugo <- combined_for_grouping_withHugo[combined_for_grouping_withHugo$PMID == "hugo", ]
combined_for_grouping_withHugo$OS_Group <- ifelse(combined_for_grouping_withHugo$OS >= os_threshold, "OS_Long", "OS_Short")
```

```{r}
library(stringr)
combined_for_grouping_withHugo <- combined_for_grouping_withHugo %>%
  mutate(Patient_ID = str_remove(Patient_ID, "Pt"))
```

```{r}
combined_hugo <- combined_for_grouping_withHugo %>%
  filter(PMID == "hugo") %>%
  inner_join(ptmeta_subtype_purity, by = c("Patient_ID" = "patient_id(pt)"))
```

```{r}
final_combined <- bind_rows(combined_gide, combined_riaz)
combined_hugo$PFS_Group <- NA
combined_hugo$`patient_id(pt)` <- combined_hugo$Patient_ID
final_combined_withHugo <- bind_rows(final_combined, combined_hugo)
```

```{r}
final_combined_withHugo_PFS <- final_combined_withHugo %>%
  filter(PFS.Event == 1)

final_combined_withHugo_PFS <- final_combined_withHugo_PFS %>%
  mutate(RECIST_Group = case_when(
    RECIST %in% c(0, 1) ~ "RECIST = 0,1",
    RECIST %in% c(2, 3) ~ "RECIST = 2,3",
    TRUE ~ NA_character_
  ))
```
```{r}
ggplot(final_combined_withHugo_PFS, aes(x = PFS, fill = RECIST_Group)) +
  geom_density(alpha = 0.5) + 
  labs(
    title = "Density Plot by RECIST Group",
    x = "PFS",
    y = "Density",
    fill = "RECIST Group"
  ) +
  theme_minimal()
```

```{r}
write.csv(final_combined_withHugo, file = "/Users/wangxin/Desktop/BINF90008/study_materials/final_combined_withHugo.csv", row.names = FALSE)
```

```{r}
RECIST_gide <- final_combined[final_combined$PMID.x == "riaz", ]
table(RECIST_gide$RECIST)
```

```{r}
# install.packages("dplyr")
library(dplyr)

# Create a crosstab showing the frequency of recist in each group
cross_table_Group_RECIST <- final_combined %>%
  group_by(PFS_Group,OS_Group, RECIST) %>%
  summarise(Frequency = n()) %>%
  ungroup()

# Show cross-tabulation
print(cross_table_Group_RECIST)
```

```{r}
write.csv(cross_table_Group_RECIST, "/Users/wangxin/Desktop/BINF90008/study_materials/cross_table_Group_RECIST.csv", row.names = FALSE)
```

## Modification of normalisation of life and death of patients in tables
### Hugo
```{r}
library(readxl)
library(writexl)
```

```{r}
file_path_hugo <- "/Users/wangxin/Desktop/BINF90008/study_materials/Hugo_mmc1.xlsx"
outcome_hugo <- read_excel(file_path_hugo)

print(colnames(outcome_hugo))

# Modify column name
colnames(outcome_hugo)[colnames(outcome_hugo) == "Vital Status"] <- "Dead/Alive(Dead = True)"

print(colnames(outcome_hugo))

print(unique(outcome_hugo$`Dead/Alive(Dead = True)`))

# Modify dead and alive to BOOL
outcome_hugo$`Dead/Alive(Dead = True)` <- ifelse(outcome_hugo$`Dead/Alive(Dead = True)` == "Alive", FALSE, 
                                         ifelse(outcome_hugo$`Dead/Alive(Dead = True)` == "Dead", TRUE, outcome_hugo$`Dead/Alive(Dead = True)`))

print(unique(outcome_hugo$`Dead/Alive(Dead = True)`))

write_xlsx(outcome_hugo, "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Hugo_mmc1.xlsx")
```

### Gide1
```{r}
file_path_gide1 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Gide_mmc2.xlsx"
outcome_gide1 <- read_excel(file_path_gide1)

print(colnames(outcome_gide1))

# Modify column contents
colnames(outcome_gide1)[colnames(outcome_gide1) == "Last Followup Status"] <- "Dead/Alive(Dead = True)"

print(colnames(outcome_gide1))

print(unique(outcome_gide1$`Dead/Alive(Dead = True)`))

# Modify dead and alive to BOOL
outcome_gide1$`Dead/Alive(Dead = True)` <- ifelse(outcome_gide1$`Dead/Alive(Dead = True)` == "Alive", FALSE, 
                                         ifelse(outcome_gide1$`Dead/Alive(Dead = True)` %in% c("Dead", "Dead, melanoma"), TRUE, outcome_gide1$`Dead/Alive(Dead = True)`))

print(unique(outcome_gide1$`Dead/Alive(Dead = True)`))

write_xlsx(outcome_gide1, "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Gide_mmc2.xlsx")
```

### Gide2
```{r}
file_path_gide2 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Gide_mmc3.xlsx"
outcome_gide2 <- read_excel(file_path_gide2)

print(colnames(outcome_gide2))

# Modify column contents
colnames(outcome_gide2)[colnames(outcome_gide2) == "Last Followup Status"] <- "Dead/Alive(Dead = True)"

print(colnames(outcome_gide2))

print(unique(outcome_gide2$`Dead/Alive(Dead = True)`))

# Modify dead and alive to BOOL
outcome_gide2$`Dead/Alive(Dead = True)` <- ifelse(outcome_gide2$`Dead/Alive(Dead = True)` == "Alive", FALSE, 
                                         ifelse(outcome_gide2$`Dead/Alive(Dead = True)` %in% c("Dead", "Dead, melanoma"), TRUE, outcome_gide2$`Dead/Alive(Dead = True)`))

print(unique(outcome_gide2$`Dead/Alive(Dead = True)`))

write_xlsx(outcome_gide2, "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Gide_mmc3.xlsx")
```

#### Modify Gide's Patient No.
```{r}
file_path_modified_gide1 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Gide_mmc2.xlsx"
modified_gide1 <- read_excel(file_path_modified_gide1)

print(colnames(modified_gide1))

# Modify data in Patient no. column by adding "PD1_" in front of it
modified_gide1$`Patient no.` <- paste0("PD1_", modified_gide1$`Patient no.`)

print(head(modified_gide1$`Patient no.`))

write_xlsx(modified_gide1, "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Patient_Gide_mmc2.xlsx")
```

```{r}
file_path_modified_gide2 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Gide_mmc3.xlsx"
modified_gide2 <- read_excel(file_path_modified_gide2)

print(colnames(modified_gide2))

# Modify data in Patient no. column by adding "PD1_" in front of it
modified_gide2$`Patient no.` <- paste0("ipiPD1_", modified_gide2$`Patient no.`)

print(head(modified_gide2$`Patient no.`))

write_xlsx(modified_gide2, "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Patient_Gide_mmc3.xlsx")
```

```{r}
file_path_modified_patient_gide1 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Patient_Gide_mmc2.xlsx"
modified_patient_gide1 <- read_excel(file_path_modified_patient_gide1)

print(colnames(modified_patient_gide1))

# Select Patient and Dead/Alive (Dead = True) columns
gide1_selected_columns <- modified_patient_gide1[, c("Patient no.", "Dead/Alive(Dead = True)")]
colnames(gide1_selected_columns)[colnames(gide1_selected_columns) == "Patient no."] <- "Patient"

# Add PMID column with "gide"
gide1_selected_columns$PMID <- "gide"

print(head(gide1_selected_columns))

write_xlsx(gide1_selected_columns, "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_gide2.xlsx")
```

```{r}
file_path_modified_patient_gide2 <- "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Patient_Gide_mmc3.xlsx"
modified_patient_gide2 <- read_excel(file_path_modified_patient_gide2)

print(colnames(modified_patient_gide2))

# Select Patient and Dead/Alive (Dead = True) columns
gide2_selected_columns <- modified_patient_gide2[, c("Patient no.", "Dead/Alive(Dead = True)")]
colnames(gide2_selected_columns)[colnames(gide2_selected_columns) == "Patient no."] <- "Patient"

# Add PMID column with "gide"
gide2_selected_columns$PMID <- "gide"

print(head(gide2_selected_columns))

write_xlsx(gide2_selected_columns, "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_gide3.xlsx")
```

```{r}
file_path_modified_hugo <- "/Users/wangxin/Desktop/BINF90008/study_materials/Modified_Hugo_mmc1.xlsx"
modified_hugo <- read_excel(file_path_modified_hugo)

print(colnames(modified_hugo))

# Select Patient and Dead/Alive (Dead = True) columns
hugo_selected_columns <- modified_hugo[, c("Patient ID", "Dead/Alive(Dead = True)")]
colnames(hugo_selected_columns)[colnames(hugo_selected_columns) == "Patient ID"] <- "Patient"

# Add PMID column with "gide"
hugo_selected_columns$PMID <- "hugo"

# Use regular expressions to preserve patient serial numbers beginning with ‘Pt’
hugo_selected_columns <- hugo_selected_columns %>%
  filter(grepl("^Pt\\d+$", Patient)) %>%
  mutate(Patient = str_remove(Patient, "Pt"))

print(head(hugo_selected_columns))

write_xlsx(hugo_selected_columns, "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_hugo.xlsx")
```

```{r}
file_path_modified_riaz <- "/Users/wangxin/Desktop/BINF90008/study_materials/Riaz_mmc2.xlsx"
modified_riaz <- read_excel(file_path_modified_riaz)

print(colnames(modified_riaz))

# Select Patient and Dead/Alive (Dead = True) columns
riaz_selected_columns <- modified_riaz[, c("Patient", "Dead/Alive(Dead = True)")]

# Add PMID column with "gide"
riaz_selected_columns$PMID <- "riaz"

print(head(riaz_selected_columns))

write_xlsx(riaz_selected_columns, "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_riaz.xlsx")
```

```{r}
file_paths_combined <- c(
  "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_gide2.xlsx",
  "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_gide3.xlsx",
  "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_hugo.xlsx",
  "/Users/wangxin/Desktop/BINF90008/study_materials/Selected_Patient_Dead_Alive_with_PMID_riaz.xlsx"
)

# Initialize an empty data frame
combined_data <- data.frame()

# Loop through each file and merge data
for (file_path in file_paths_combined) {
  data <- read_excel(file_path)
  combined_data <- rbind(combined_data, data)
}

print(head(combined_data))
print(dim(combined_data))
```

```{r}
output_path <- "/Users/wangxin/Desktop/BINF90008/study_materials/Patient_Vital_Status.xlsx"
write_xlsx(combined_data, output_path)
```

##  Filter dead and alive into 2 groups
```{r}
library(readxl)
library(dplyr)
library(writexl)
```

```{r}
file_path_dead_alive <- "/Users/wangxin/Desktop/BINF90008/study_materials/Patient_Vital_Status.xlsx"
data_dead_alive <- read_excel(file_path_dead_alive)

print(colnames(data_dead_alive))
print(colnames(final_combined_withHugo))

# Filter data_dead_alive to keep only rows in final_combined_data$Patient ID
filtered_data_dead_alive <- data_dead_alive[data_dead_alive$Patient %in% final_combined_withHugo$Patient_ID, ]

print(head(filtered_data_dead_alive))
print(dim(filtered_data_dead_alive))

output_path <- "/Users/wangxin/Desktop/BINF90008/study_materials/Filtered_Patient_Vital_Status.xlsx"
write_xlsx(filtered_data_dead_alive, output_path)
```

```{r}
# Separate living and dead patients
alive_patients <- filtered_data_dead_alive %>% filter(`Dead/Alive(Dead = True)` == FALSE)
dead_patients <- filtered_data_dead_alive %>% filter(`Dead/Alive(Dead = True)` == TRUE)

# Picked from final_combined_withHugo based on Patient
alive_group <- final_combined_withHugo %>% filter(Patient_ID %in% alive_patients$Patient)
dead_group <- final_combined_withHugo %>% filter(Patient_ID %in% dead_patients$Patient)

print(head(alive_group))
print(dim(alive_group))

print(head(dead_group))
print(dim(dead_group))

write_xlsx(alive_group, "/Users/wangxin/Desktop/BINF90008/study_materials/Alive_group.xlsx")
write_xlsx(dead_group, "/Users/wangxin/Desktop/BINF90008/study_materials/Dead_group.xlsx")
```

```{r}
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
```

```{r}
# Remove missing value of PFS
alive_group_pfs <- alive_group %>% drop_na(PFS)
dead_group_pfs <- dead_group %>% drop_na(PFS)

print(head(alive_group_pfs))
print(head(dead_group_pfs))
print(dim(alive_group_pfs))
print(dim(dead_group_pfs))
```

```{r}
group_PFS <- bind_rows(alive_group %>% drop_na(PFS), 
                                dead_group %>% drop_na(PFS))
```

```{r}
# Mapping of PFS density in alive_group
library(ggplot2)
ggplot(alive_group_pfs, aes(x = PFS)) + 
  geom_density(fill = "blue", alpha = 0.7) +
  labs(title = "PFS Density Plot", x = "PFS (days)", y = "Density")
```

```{r}
# Mapping of OS density in alive group
ggplot(alive_group, aes(x = OS)) + 
  geom_density(fill = "green", alpha = 0.7) +
  labs(title = "OS Density Plot", x = "OS (days)", y = "Density")
```

```{r}
# Mapping of PFS density in dead_group
library(ggplot2)
ggplot(dead_group_pfs, aes(x = PFS)) + 
  geom_density(fill = "blue", alpha = 0.7) +
  labs(title = "PFS Density Plot", x = "PFS (days)", y = "Density")
```

```{r}
# Mapping of OS density in dead_group
ggplot(dead_group, aes(x = OS)) + 
  geom_density(fill = "green", alpha = 0.7) +
  labs(title = "OS Density Plot", x = "OS (days)", y = "Density")
```

```{r}
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
```

```{r}
alive_group_pfsRECIST <- alive_group_pfs %>% filter(RECIST %in% c(0, 1, 2, 3))
alive_group_osRECIST <- alive_group %>% filter(RECIST %in% c(0, 1, 2, 3))

# Define function that plots density
plot_density <- function(data, variable, group_name, recist_value) {
  ggplot(data, aes_string(x = variable)) +
    geom_density(fill = "blue", alpha = 0.5) +
    ggtitle(paste(group_name, "- RECIST", recist_value, "-", variable, "Density Plot")) +
    xlab(variable) +
    ylab("Density")
}

# Plot and save density maps for PFS and OS
variables <- c("PFS", "OS")
recist_values <- c(0, 1, 2, 3)

# Initialize lists to store drawing objects
plots_pfs <- list()
plots_os <- list()

for (recist_value in recist_values) {
  data_pfs <- alive_group_pfs %>% filter(RECIST == recist_value)
  data_os <- alive_group %>% filter(RECIST == recist_value)
  
  # Map density of PFS and store
  plot_pfs <- plot_density(data_pfs, "PFS", "Alive Group", recist_value)
  plots_pfs[[length(plots_pfs) + 1]] <- plot_pfs
  
  # Map density of OS and store
  plot_os <- plot_density(data_os, "OS", "Alive Group", recist_value)
  plots_os[[length(plots_os) + 1]] <- plot_os
}

# Using patchwork to merge PFS diagrams
combined_alive_pfs_plot <- wrap_plots(plots_pfs, ncol = 2)
ggsave("/Users/wangxin/Desktop/BINF90008/study_materials/Combined_PFS_density_plot.png", combined_alive_pfs_plot, width = 10, height = 8)

# Using patchwork to merge OS diagrams
combined_alive_os_plot <- wrap_plots(plots_os, ncol = 2)
ggsave("/Users/wangxin/Desktop/BINF90008/study_materials/Combined_OS_density_plot.png", combined_alive_os_plot, width = 10, height = 8)
```

```{r}
dead_group_pfsRECIST <- dead_group_pfs %>% filter(RECIST %in% c(0, 1, 2, 3))
dead_group_osRECIST <- dead_group %>% filter(RECIST %in% c(0, 1, 2, 3))

# Define function that plots density
plot_density <- function(data, variable, group_name, recist_value) {
  ggplot(data, aes_string(x = variable)) +
    geom_density(fill = "green", alpha = 0.5) +
    ggtitle(paste(group_name, "- RECIST", recist_value, "-", variable, "Density Plot")) +
    xlab(variable) +
    ylab("Density")
}

# Plot and save density maps for PFS and OS
variables <- c("PFS", "OS")
recist_values <- c(0, 1, 2, 3)

# Initialize the list to store drawing objects
plots_pfs <- list()
plots_os <- list()

for (recist_value in recist_values) {
  data_pfs <- dead_group_pfs %>% filter(RECIST == recist_value)
  data_os <- dead_group %>% filter(RECIST == recist_value)
  
  # Map density of PFS and store
  plot_pfs <- plot_density(data_pfs, "PFS", "Alive Group", recist_value)
  plots_pfs[[length(plots_pfs) + 1]] <- plot_pfs
  
  # Map density of OS and store
  plot_os <- plot_density(data_os, "OS", "Alive Group", recist_value)
  plots_os[[length(plots_os) + 1]] <- plot_os
}

# Using patchwork to merge PFS diagrams
combined_dead_pfs_plot <- wrap_plots(plots_pfs, ncol = 2)
ggsave("/Users/wangxin/Desktop/BINF90008/study_materials/Combined_PFS_dead_density_plot.png", combined_dead_pfs_plot, width = 10, height = 8)

# Using patchwork to merge OS diagrams
combined_dead_os_plot <- wrap_plots(plots_os, ncol = 2)
ggsave("/Users/wangxin/Desktop/BINF90008/study_materials/Combined_OS_dead_density_plot.png", combined_dead_os_plot, width = 10, height = 8)
```

```{r}
# RECIST items 0 and 1 as one group, 2 and 3 as another group
alive_group_pfsRECIST_2 <- alive_group_pfs %>%
  mutate(RECIST_group = ifelse(RECIST %in% c(0, 1), "RECIST 0,1", "RECIST 2,3"))

dead_group_pfsRECIST_2 <- dead_group_pfs %>%
  mutate(RECIST_group = ifelse(RECIST %in% c(0, 1), "RECIST 0,1", "RECIST 2,3"))

alive_group_osRECIST_2 <- alive_group %>%
  mutate(RECIST_group = ifelse(RECIST %in% c(0, 1), "RECIST 0,1", "RECIST 2,3"))

dead_group_osRECIST_2 <- dead_group %>%
  mutate(RECIST_group = ifelse(RECIST %in% c(0, 1), "RECIST 0,1", "RECIST 2,3"))

# Define function that plots density
plot_density <- function(data, variable, group_name) {
  ggplot(data, aes_string(x = variable, fill = "RECIST_group", color = "RECIST_group")) +
    geom_density(alpha = 0.5) +
    ggtitle(paste(group_name, "-", variable, "Density Plot")) +
    xlab(variable) +
    ylab("Density")
}

# Plot and save PFS density maps
pfs_group_names <- list(Alive = alive_group_pfsRECIST_2, Dead = dead_group_pfsRECIST_2)
for (group_name in names(pfs_group_names)) {
  data <- pfs_group_names[[group_name]]
  density_plot <- plot_density(data, "PFS", group_name)
  ggsave(paste0("/Users/wangxin/Desktop/BINF90008/study_materials/", group_name, "_PFS_density_plot.png"), density_plot)
}

# Plot and save OS density maps
os_group_names <- list(Alive = alive_group_osRECIST_2, Dead = dead_group_osRECIST_2)
for (group_name in names(os_group_names)) {
  data <- os_group_names[[group_name]]
  density_plot <- plot_density(data, "OS", group_name)
  ggsave(paste0("/Users/wangxin/Desktop/BINF90008/study_materials/", group_name, "_OS_density_plot.png"), density_plot)
}
```
