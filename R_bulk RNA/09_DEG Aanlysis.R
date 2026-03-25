## DEG Aanlysis
### PFS Group and OS Group
```{r}
library(dplyr)
library(tidyr)
```

#### PFS Group
```{r}
# Convert ruv_iii_prps_normal_PFS to long format
ruv_iii_prps_normal_PFS_long <- as.data.frame(ruv_iii_prps_normal_PFS) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Run", values_to = "Expression")
```
```{r}
# Merge grouped information ( use Run column)
ruv_iii_prps_normal_PFS_DEG <- left_join(
  ruv_iii_prps_normal_PFS_long, 
  final_combined_withHugo[, c("Run", "anat_loc", "Treatment", "Response", "RECIST", "AvgSpotLen", "disease_status", "subtypes", "PFS_Group", "OS_Group")], 
  by = "Run"
)

head(ruv_iii_prps_normal_PFS_DEG)
```

#### OS Group
```{r}
# Convert ruv_iii_prps_normal_PFS to long format
ruv_iii_prps_normal_OS_long <- as.data.frame(ruv_iii_prps_normal_OS) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Run", values_to = "Expression")
```
```{r}
# Merge grouped information ( use Run column)
ruv_iii_prps_normal_OS_DEG <- left_join(
  ruv_iii_prps_normal_OS_long, 
  final_combined_withHugo[, c("Run", "anat_loc", "Treatment", "Response", "RECIST", "AvgSpotLen", "disease_status", "subtypes", "PFS_Group", "OS_Group")], 
  by = "Run"
)

head(ruv_iii_prps_normal_OS_DEG)
```


#### PFS DEG Analysis
```{r}
# Filter out rows with PFS_Group of NA
filtered_ruv_iii_prps_normal_PFS_DEG <- ruv_iii_prps_normal_PFS_DEG %>% 
  filter(!is.na(PFS_Group)) %>% 
  mutate(PFS_Group = as.factor(PFS_Group))
```

```{r}
# Prepare expression matrix
expression_matrix_ruv_iii_prps_normal_PFS_DEG <- filtered_ruv_iii_prps_normal_PFS_DEG %>% 
  select(Gene, Run, Expression) %>% 
  pivot_wider(names_from = Run, values_from = Expression) %>% 
  column_to_rownames(var = "Gene")

# Replace negative values with 0
expression_matrix_ruv_iii_prps_normal_PFS_DEG[expression_matrix_ruv_iii_prps_normal_PFS_DEG < 0] <- 0

# Prepare grouping information
expression_matrix_group_PFS <- filtered_ruv_iii_prps_normal_PFS_DEG %>%
  select(Run, PFS_Group) %>%  # Extract sample name and grouping information
  distinct() %>%  # Make sure there is only one subgroup for each sample
  arrange(match(Run, colnames(expression_matrix_ruv_iii_prps_normal_PFS_DEG))) %>%  # Ensure consistent order
  pull(PFS_Group)  # Extract as vector
```

```{r}
library(edgeR)
```
```{r}
# Create a DGEList object
dge_filtered_ruv_iii_prps_normal_PFS_DEG <- DGEList(counts = expression_matrix_ruv_iii_prps_normal_PFS_DEG, group = expression_matrix_group_PFS)

# Create a design matrix
design_filtered_ruv_iii_prps_normal_PFS_DEG <- model.matrix(~expression_matrix_group_PFS)

# Estimate dispersion
dge_dispersion_PFS_DEG <- estimateDisp(dge_filtered_ruv_iii_prps_normal_PFS_DEG, design_filtered_ruv_iii_prps_normal_PFS_DEG)
```

```{r}
# Fit model with GLM
fit_filtered_ruv_iii_prps_normal_PFS_DEG <- glmFit(dge_dispersion_PFS_DEG, design_filtered_ruv_iii_prps_normal_PFS_DEG)

# Perform differential expression test
lrt_filtered_ruv_iii_prps_normal_PFS_DEG <- glmLRT(fit_filtered_ruv_iii_prps_normal_PFS_DEG, coef = 2)  # coef=2 is to compare the second grouping (Short vs Long)

# Extract significantly differentially expressed genes
top_tags_PFS_DEG <- topTags(lrt_filtered_ruv_iii_prps_normal_PFS_DEG, n = Inf)  # Extract all genes
DEG_PFS_results <- top_tags_PFS_DEG$table
```

```{r}
DEG_PFS_results_significant <- DEG_PFS_results[DEG_PFS_results$FDR < 0.05, ]

DEG_PFS_results_final <- DEG_PFS_results_significant[, c("logFC", "PValue", "FDR")]
```

```{r}
significant_genes_DEG_PFS_results <- DEG_PFS_results %>% filter(FDR < 0.05)

nrow(significant_genes_DEG_PFS_results)
head(significant_genes_DEG_PFS_results)
significant_genes_DEG_PFS_results_names <- rownames(significant_genes_DEG_PFS_results)
```

```{r}
library(EnhancedVolcano)
```
```{r}
# Ensure `DEG_PFS_results` contains the columns: logFC (log2 Fold Change) and PValue (p-value)
EnhancedVolcano(
  DEG_PFS_results,
  lab = rownames(DEG_PFS_results),       # Gene names as labels
  x = 'logFC',                           # log2 Fold Change
  y = 'PValue',                          # P-value
  xlab = bquote(~Log[2]~ 'Fold Change'), # X-axis label
  ylab = bquote(~-Log[10]~italic(P)),    # Y-axis label
  title = 'PFS Differential Expression Analysis', # Plot title
  subtitle = 'EnhancedVolcano',
  pCutoff = 0.05,                        # P-value cutoff
  FCcutoff = 1,                          # Fold Change cutoff
  pointSize = 3,                         # Size of points
  labSize = 3,                           # Size of labels
  ylim = c(0, 8),                        # Limit y-axis range
  xlim = c(-2.5, 2.5),                   # Limit x-axis range
  legendPosition = 'right',              # Legend position
  legendLabSize = 10,                    # Legend label size
  legendIconSize = 3,                    # Legend icon size
  col = c('grey', 'green', 'blue', 'red'), # Colors for different significance levels
  colAlpha = 0.6,                        # Transparency of points
  gridlines.major = TRUE,                # Display major gridlines
  gridlines.minor = FALSE                # Hide minor gridlines
)
```


#### OS DEG Analysis
```{r}
# Filter out rows with PFS_Group of NA
filtered_ruv_iii_prps_normal_OS_DEG <- ruv_iii_prps_normal_OS_DEG %>% 
  filter(!is.na(OS_Group)) %>% 
  mutate(OS_Group = as.factor(OS_Group))
```

```{r}
# Prepare expression matrix
expression_matrix_ruv_iii_prps_normal_OS_DEG <- filtered_ruv_iii_prps_normal_OS_DEG %>% 
  select(Gene, Run, Expression) %>% 
  pivot_wider(names_from = Run, values_from = Expression) %>% 
  column_to_rownames(var = "Gene")

# Replace negative values with 0
expression_matrix_ruv_iii_prps_normal_OS_DEG[expression_matrix_ruv_iii_prps_normal_OS_DEG < 0] <- 0

# Prepare grouping information
expression_matrix_group_OS <- filtered_ruv_iii_prps_normal_OS_DEG %>%
  select(Run, OS_Group) %>%  # Extract sample name and grouping information
  distinct() %>%  # Make sure there is only one subgroup for each sample
  arrange(match(Run, colnames(expression_matrix_ruv_iii_prps_normal_OS_DEG))) %>%  # Ensure consistent order
  pull(OS_Group)  # Extract as vector
```

```{r}
library(edgeR)
```
```{r}
# Create a DGEList object
dge_filtered_ruv_iii_prps_normal_OS_DEG <- DGEList(counts = expression_matrix_ruv_iii_prps_normal_OS_DEG, group = expression_matrix_group_OS)

# Create a design matrix
design_filtered_ruv_iii_prps_normal_OS_DEG <- model.matrix(~expression_matrix_group_OS, data = filtered_ruv_iii_prps_normal_OS_DEG)

# Estimate dispersion
dge_dispersion_OS_DEG <- estimateDisp(dge_filtered_ruv_iii_prps_normal_OS_DEG, design_filtered_ruv_iii_prps_normal_OS_DEG)
```

```{r}
# Fit model with GLM
fit_filtered_ruv_iii_prps_normal_OS_DEG <- glmFit(dge_dispersion_OS_DEG, design_filtered_ruv_iii_prps_normal_OS_DEG)

# Perform differential expression test
lrt_filtered_ruv_iii_prps_normal_OS_DEG <- glmLRT(fit_filtered_ruv_iii_prps_normal_OS_DEG, coef = 2)  # coef=2 is to compare the second grouping (Short vs Long)

# Extract significantly differentially expressed genes
top_tags_OS_DEG <- topTags(lrt_filtered_ruv_iii_prps_normal_OS_DEG, n = Inf)  # Extract all genes
DEG_OS_results <- top_tags_OS_DEG$table
```

```{r}
significant_genes_DEG_OS_results <- DEG_OS_results %>% filter(FDR < 0.05)

nrow(significant_genes_DEG_OS_results)
head(significant_genes_DEG_OS_results)
significant_genes_DEG_OS_results_names <- rownames(significant_genes_DEG_OS_results)
```

```{r}
# Ensure `DEG_PFS_results` contains the columns: logFC (log2 Fold Change) and PValue (p-value)
EnhancedVolcano(
  DEG_OS_results,
  lab = rownames(DEG_OS_results),       # Gene names as labels
  x = 'logFC',                           # log2 Fold Change
  y = 'PValue',                          # P-value
  xlab = bquote(~Log[2]~ 'Fold Change'), # X-axis label
  ylab = bquote(~-Log[10]~italic(P)),    # Y-axis label
  title = 'OS Differential Expression Analysis', # Plot title
  subtitle = 'EnhancedVolcano',
  pCutoff = 0.05,                        # P-value cutoff
  FCcutoff = 1,                          # Fold Change cutoff
  pointSize = 3,                         # Size of points
  labSize = 3,                           # Size of labels
  ylim = c(0, 8),                        # Limit y-axis range
  xlim = c(-2.5, 2.5),                   # Limit x-axis range
  legendPosition = 'right',              # Legend position
  legendLabSize = 10,                    # Legend label size
  legendIconSize = 3,                    # Legend icon size
  col = c('grey', 'green', 'blue', 'red'), # Colors for different significance levels
  colAlpha = 0.6,                        # Transparency of points
  gridlines.major = TRUE,                # Display major gridlines
  gridlines.minor = FALSE                # Hide minor gridlines
)
```
