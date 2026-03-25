## Open files and standardize Patient_ID
### ptmeta
```{r}
library(readr)
# Read CSV file:ptmeta_subtype_purity
ptmeta_subtype_purity <- read_csv("/Users/wangxin/Desktop/BINF90008/jess_code_v1/ptmeta_subtype_purity.csv")

table(ptmeta_subtype_purity$PMID)
```   

```{r}
# Relevant column name as SampleID"
# Add a new column to remove the last four digits
library(dplyr)
ptmeta_subtype_purity <- ptmeta_subtype_purity %>%
  mutate(NewPatientID = substr(`patient_id(pt)`, 1, nchar(`patient_id(pt)`) - 4))

head(ptmeta_subtype_purity)
```

### Gide
```{r}
library(readxl)
library(readr)
library(data.table)
```

```{r}
clinical_data_Gide1 <- read_delim("/Users/wangxin/Desktop/BINF90008/study_materials/Gide2019_PD1_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab_Melanoma.clinical", delim = "\t")

head(clinical_data_Gide1)
```

```{r}
clinical_data_Gide2 <- read_delim("/Users/wangxin/Desktop/BINF90008/study_materials/Gide2019_PD1+CTLA4_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.clinical", delim = "\t")

head(clinical_data_Gide2)
```

#### Gide Patient_ID standrization
```{r}
# Put Treatment prefix on ID
add_treatment_prefix <- function(df, treatment_prefix) {
  df$...1 <- paste0(treatment_prefix, df$...1)
  return(df)
}

# Add a prefix to two data frames
clinical_data_Gide1 <- add_treatment_prefix(clinical_data_Gide1, "PD1_")
clinical_data_Gide2 <- add_treatment_prefix(clinical_data_Gide2, "ipiPD1_")

clinical_data_Gide1$PMID <- "gide"
clinical_data_Gide2$PMID <- "gide"

head(clinical_data_Gide1)
head(clinical_data_Gide2)
```

### Hugo
```{r}
clinical_data_Hugo <- read_delim("/Users/wangxin/Desktop/BINF90008/study_materials/Hugo2016_PD1_Melanoma_RNASeq/ICB.Hugo2016_Pembrolizumab_Melanoma.clinical", delim = "\t")

head(clinical_data_Hugo)
```

#### Hugo Patient_ID standrization
```{r}
clinical_data_Hugo$PMID <- "hugo"
```

### Riaz
```{r}
# Read Excel file
clinical_data_Riaz_N <- read_delim("/Users/wangxin/Desktop/BINF90008/study_materials/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Naive/ICB.Riaz2017_Nivolumab_Melanoma_Naive.clinical", delim = "\t")

head(clinical_data_Riaz_N)
```

```{r}
clinical_data_Riaz_P <- read_delim("/Users/wangxin/Desktop/BINF90008/study_materials/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Prog/ICB.Riaz2017_Nivolumab_Melanoma_Prog.clinical", delim = "\t")

head(clinical_data_Riaz_P)
```

#### Riaz Patient_ID standrization
```{r}
clinical_data_Riaz_N$PMID <- "riaz"
clinical_data_Riaz_P$PMID <- "riaz"
```


## Selection of PFS, OS and RECIST programs
```{r}
# Extract and rename ID columns
extract_from_clinical_files <- function(df) {
  df_subset <- df[, c("...1", "OS", "PFS", "RECIST","PMID","PFS.Event","OS.Event")]
  colnames(df_subset)[colnames(df_subset) == "...1"] <- "Patient_ID"
  return(df_subset)
}

# Extract data
data_Gide1 <- extract_from_clinical_files(clinical_data_Gide1)
data_Gide2 <- extract_from_clinical_files(clinical_data_Gide2)
data_Riaz_N <- extract_from_clinical_files(clinical_data_Riaz_N)
data_Riaz_P <- extract_from_clinical_files(clinical_data_Riaz_P)

# Merge all data frames
combined_for_grouping<- rbind(data_Gide1, data_Gide2, data_Riaz_N, data_Riaz_P)

head(combined_for_grouping)
```

```{r}
# Extract and rename the columns of clinical_data_Hugo
clinical_data_Hugo_subset <- clinical_data_Hugo[, c("OS", "OS.Event", "RECIST", "PMID")]
colnames(clinical_data_Hugo_subset) <- c("Patient_ID", "OS", "RECIST", "PMID")

# Add PFS column, populate NA
clinical_data_Hugo_subset$PFS <- NA
clinical_data_Hugo_subset$PFS_Group <- NA
clinical_data_Hugo_subset$PFS.Event <- NA
clinical_data_Hugo_subset$OS.Event <- 1

# Reorder columns
clinical_data_Hugo_subset <- clinical_data_Hugo_subset[, c("Patient_ID", "OS", "PFS", "RECIST", "PMID", "PFS.Event", "OS.Event")]
```

```{r}
# Merging data frames
combined_for_grouping_withHugo <- rbind(combined_for_grouping, clinical_data_Hugo_subset)

head(combined_for_grouping_withHugo)
```

```{r}
combined_for_grouping_withHugo <- subset(combined_for_grouping_withHugo, !is.na(OS))
any(is.na(combined_for_grouping_withHugo$OS))
```
