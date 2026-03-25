## Cox analysis and Kaplan plots
```{r}
library(survival)
library(survminer)
library(coxphf)

# OS time: OS
# OS statement: vital_status (0 = occurrence，1 = eventless)
# Other variables: Group, subtypes, tumour_purity, disease_status, anat_loc, Treatment, Response, AvgSpotLen, RECIST
```

### Cox regression analysis of PFS & OS in relation to each variable
```{r}
final_combined_withHugo_no_na_PFS <- final_combined_withHugo[complete.cases(final_combined_withHugo[, c("PFS", "PFS_Group", "PMID.x")]), ]
```

```{r}
final_combined_withHugo_no_na_OS <- final_combined_withHugo[complete.cases(final_combined_withHugo[, c("OS", "OS_Group", "PMID.x")]), ]
```


#### Group
```{r}
cox_model_PFS_Group <- coxphf(Surv(PFS, vital_status) ~ PFS_Group, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_PFS_Group)
```

```{r}
cox_model_OS_Group <- coxphf(Surv(OS, vital_status) ~ OS_Group, data = final_combined_withHugo_no_na_OS)
summary(cox_model_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Groups
surv_fit_PFS_Group <- survfit(Surv(PFS, vital_status) ~ PFS_Group, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "PFS_Group")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Groups
surv_fit_OS_Group <- survfit(Surv(OS, vital_status) ~ OS_Group, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "OS_Group")      # title
```

#### subtypes
```{r}
cox_model_subtypes_PFS_Group <- coxph(Surv(PFS, vital_status) ~ subtypes, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_subtypes_PFS_Group)
```

```{r}
cox_model_subtypes_OS_Group <- coxph(Surv(OS, vital_status) ~ subtypes, data = final_combined_withHugo_no_na_OS)
summary(cox_model_subtypes_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by subtypes
surv_fit_subtypes_PFS_Group <- survfit(Surv(PFS, vital_status) ~ subtypes, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_subtypes_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by subtypes PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "subtypes")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by subtypes
surv_fit_subtypes_OS_Group <- survfit(Surv(OS, vital_status) ~ subtypes, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_subtypes_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by subtypes OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "subtypes")      # title
```

#### tumour_purity
```{r}
final_combined_withHugo_no_na_PFS$tumour_purity <- as.numeric(final_combined_withHugo_no_na_PFS$tumour_purity)
```

```{r}
cox_model_tumour_purity_PFS_Group <- coxph(Surv(PFS, vital_status) ~ tumour_purity, data = final_combined_withHugo_no_na_PFS)  # tumour_purity as continuous variable
summary(cox_model_tumour_purity_PFS_Group)
```

```{r}
final_combined_withHugo_no_na_OS$tumour_purity <- as.numeric(final_combined_withHugo_no_na_OS$tumour_purity)
```

```{r}
cox_model_tumour_purity_OS_Group <- coxph(Surv(OS, vital_status) ~ tumour_purity, data = final_combined_withHugo_no_na_OS)  # tumour_purity as continuous variable
summary(cox_model_tumour_purity_OS_Group)
```

#### disease_status
```{r}
cox_model_disease_status_PFS_Group <- coxph(Surv(PFS, vital_status) ~ disease_status, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_disease_status_PFS_Group)
```

```{r}
cox_model_disease_status_OS_Group <- coxph(Surv(OS, vital_status) ~ disease_status, data = final_combined_withHugo_no_na_OS)
summary(cox_model_disease_status_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by disease_status
surv_fit_disease_status_PFS_Group <- survfit(Surv(PFS, vital_status) ~ disease_status, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_disease_status_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by disease_status PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "disease_status")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by disease_status
surv_fit_disease_status_OS_Group <- survfit(Surv(OS, vital_status) ~ disease_status, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_disease_status_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by disease_status OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "disease_status")      # title
```

#### anat_loc
```{r}
cox_model_anat_loc_PFS_Group <- coxph(Surv(PFS, vital_status) ~ anat_loc, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_anat_loc_PFS_Group)
```

```{r}
cox_model_anat_loc_OS_Group <- coxph(Surv(OS, vital_status) ~ anat_loc, data = final_combined_withHugo_no_na_OS)
summary(cox_model_anat_loc_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by anat_loc
surv_fit_anat_loc_PFS_Group <- survfit(Surv(PFS, vital_status) ~ anat_loc, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_anat_loc_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by anat_loc PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "anat_loc")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by anat_loc
surv_fit_anat_loc_OS_Group <- survfit(Surv(OS, vital_status) ~ anat_loc, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_anat_loc_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by anat_loc OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "anat_loc")      # title
```

#### Treatment
```{r}
cox_model_Treatment_PFS_Group <- coxph(Surv(PFS, vital_status) ~ Treatment, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_Treatment_PFS_Group)
```

```{r}
cox_model_Treatment_OS_Group <- coxph(Surv(OS, vital_status) ~ Treatment, data = final_combined_withHugo_no_na_OS)
summary(cox_model_Treatment_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Treatment
surv_fit_Treatment_PFS_Group <- survfit(Surv(PFS, vital_status) ~ Treatment, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_Treatment_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by Treatment PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "Treatment")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Treatment
surv_fit_Treatment_OS_Group <- survfit(Surv(OS, vital_status) ~ Treatment, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_Treatment_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by Treatment OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "Treatment")      # title
```

#### Response
```{r}
cox_model_Response_PFS_Group <- coxph(Surv(PFS, vital_status) ~ Response, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_Response_PFS_Group)
```

```{r}
cox_model_Response_OS_Group <- coxph(Surv(OS, vital_status) ~ Response, data = final_combined_withHugo_no_na_OS)
summary(cox_model_Response_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Treatment
surv_fit_Response_PFS_Group <- survfit(Surv(PFS, vital_status) ~ Response, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_Response_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by Response PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "Response")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by Treatment
surv_fit_Response_OS_Group <- survfit(Surv(OS, vital_status) ~ Response, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_Response_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by Response OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "Response")      # title
```

#### AvgSpotLen
```{r}
cox_model_AvgSpotLen_PFS_Group <- coxph(Surv(PFS, vital_status) ~ AvgSpotLen, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_AvgSpotLen_PFS_Group)
```

```{r}
cox_model_AvgSpotLen_OS_Group <- coxph(Surv(OS, vital_status) ~ AvgSpotLen, data = final_combined_withHugo_no_na_OS)
summary(cox_model_AvgSpotLen_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by AvgSpotLen
surv_fit_AvgSpotLen_PFS_Group <- survfit(Surv(PFS, vital_status) ~ AvgSpotLen, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_AvgSpotLen_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by AvgSpotLen PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "AvgSpotLen")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by AvgSpotLen
surv_fit_AvgSpotLen_OS_Group <- survfit(Surv(OS, vital_status) ~ AvgSpotLen, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_AvgSpotLen_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by AvgSpotLen OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "AvgSpotLen")      # title
```

#### RECIST
```{r}
cox_model_RECIST_PFS_Group <- coxph(Surv(PFS, vital_status) ~ RECIST, data = final_combined_withHugo_no_na_PFS)
summary(cox_model_RECIST_PFS_Group)
```

```{r}
cox_model_RECIST_OS_Group <- coxph(Surv(OS, vital_status) ~ RECIST, data = final_combined_withHugo_no_na_OS)
summary(cox_model_RECIST_OS_Group)
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by RECIST
surv_fit_RECIST_PFS_Group <- survfit(Surv(PFS, vital_status) ~ RECIST, data = final_combined_withHugo_no_na_PFS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_RECIST_PFS_Group, 
           data = final_combined_withHugo_no_na_PFS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by RECIST PFS_Group",
           xlab = "Time (PFS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "RECIST")      # title
```

```{r}
# Surv and survfit generate Kaplan-Meier survival curves by RECIST
surv_fit_RECIST_OS_Group <- survfit(Surv(OS, vital_status) ~ RECIST, data = final_combined_withHugo_no_na_OS)

# Mapping Kaplan-Meier survival curve
ggsurvplot(surv_fit_RECIST_OS_Group, 
           data = final_combined_withHugo_no_na_OS, 
           risk.table = TRUE,           # Show risk table
           pval = TRUE,                 # Display Log-rank test p-value
           conf.int = TRUE,             # Show confidence intervals
           title = "Kaplan-Meier Survival Curves by RECIST OS_Group",
           xlab = "Time (OS)",          # x-axis labels
           ylab = "Survival Probability", # y-axis labels
           legend.title = "RECIST")      # title
```
