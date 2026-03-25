## RUV-III
### RLE before RUV-III
```{r}
library(ggplot2)
library(reshape2)
```

```{r}
# Data ‘rawCounts_expressing’ is a frame of data with rows for genes and columns for samples
# In order to plot RLE, we need to calculate median for each gene and then calculate difference between each sample and this median
# ‘RawCounts_express’ is already log2 converted data.

# 
log2_rawCounts_expressing <- log2(rawCounts_expressing + 1)  # Add 1 to avoid taking a log of 0

# Calculate median for each gene
gene_medians_before <- apply(log2_rawCounts_expressing, 1, stats::median)

# Calculate expression difference relative to median gene (i.e., relative log expression value)
rle_values_before <- sweep(log2_rawCounts_expressing, 1, gene_medians_before)

# Converting data to long format for plotting with ggplot2
rle_values_long_before <- melt(rle_values_before)
```

```{r}
# Merge grouping information in pdata_2 directly based on row names
rle_values_long_before$PFS_Group <- pdata_2[rle_values_long_before$Var2, "PFS_Group"]
rle_values_long_before$OS_Group <- pdata_2[rle_values_long_before$Var2, "OS_Group"]
```

#### RLE for PFS
```{r}
rle_values_long_PFS_before <- subset(rle_values_long_before, !is.na(PFS_Group))
```

```{r}
ggplot(rle_values_long_PFS_before, aes(x = Var2, y = value, color = PFS_Group)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  labs(title = "RLE Plot with PFS Grouping before RUV-III", x = "Samples", y = "Relative Log Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange", size = 0.8)+
  coord_cartesian(ylim = c(-5, 5)) +
  scale_color_manual(values = c("PFS_Long" = "darkgreen", "PFS_Short" = "darkred"))
```

#### RLE for OS
```{r}
rle_values_long_OS_before <- subset(rle_values_long_before, !is.na(OS_Group))
```

```{r}
ggplot(rle_values_long_OS_before, aes(x = Var2, y = value, color = OS_Group)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
  labs(title = "RLE Plot with OS Grouping before RUV-III", x = "Samples", y = "Relative Log Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange", size = 0.8) +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_color_manual(values = c("OS_Long" = "purple", "OS_Short" = "navy"))
```

### RUV-III precedure
```{r}
library(ruv)  # RUV package, which provides functions such as RUV1, residop, etc.
library(Matrix)  # Matrix manipulation functions are provided
library(stats)  # Includes basic statistics-related functions such as svd, solve, etc.
```

```{r}
RUV_III_PRPS <- function(
  Y,
  M,
  ctl,
  k = NULL, 
  eta = NULL, 
  include.intercept = TRUE,
  average = FALSE, 
  fullalpha = NULL, 
  return.info = FALSE, 
  inputcheck = TRUE) {

  print("Starting RUV_III_PRPS function...")
  
  if (is.data.frame(Y)) {
    print("Converting Y to matrix if necessary...")
    Y <- data.matrix(Y)
  }
  
  m <- nrow(Y)
  n <- ncol(Y)
  print(paste("Dimensions of Y:", dim(Y)))
  
  M <- ruv::replicate.matrix(M)
  print(paste("Dimensions of M:", dim(M)))

  # Modify here by replacing tological with as.logical
  ctl <- as.logical(ctl)

  if (inputcheck) {
    if (m > n) {
      warning("m is greater than n! This may indicate that you need to transpose your data matrix.")
    }
    if (sum(is.na(Y)) > 0) {
      warning("Y contains missing values. This is not supported.")
    }
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0) {
      warning("Y contains infinities. This is not supported.")
    }
  }
  
  Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
  print("Applying RUV1...")

  mu <- colMeans(Y)
  mu_mat <- rep(1, m) %*% t(mu)
  Y_stand <- Y - mu_mat
  print(paste("Y_stand created successfully:", dim(Y_stand)))
  
  if (ncol(M) > m) {
    print("Condition met: ncol(M) >= m, setting newY to Y.")
    newY <- Y
  } else if (is.null(k)) {
    print("Calculating newY using default k...")
    ycyctinv <- solve(Y[, ctl] %*% t(Y[, ctl]))
    newY <- (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv)) %*% Y
    fullalpha <- NULL
  } else if (k == 0) {
    print("k is 0, setting newY to Y.")
    newY <- Y
    fullalpha <- NULL
  } else {
    print("Calculating Y0...")
    Y0 <- tryCatch({
      ruv::residop(Y, M)
    }, error = function(e) {
      print("Error in generating Y0")
      print(e)
      return(NULL)
    })
    if (is.null(Y0)) return(NULL)
    
    print(paste("Dimensions of Y0:", dim(Y0)))

    fullalpha <- tryCatch({
      t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M), sum(ctl)), drop = FALSE]) %*% Y
    }, error = function(e) {
      print("Error in generating fullalpha")
      print(e)
      return(NULL)
    })
    if (is.null(fullalpha)) return(NULL)
    
    print(paste("Dimensions of fullalpha:", dim(fullalpha)))

    alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]

    print(paste("Dimensions of ac:", dim(ac)))
    
    W <- tryCatch({
      Y_stand[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    }, error = function(e) {
      print("Error in calculating W")
      print(e)
      return(NULL)
    })
    if (is.null(W)) return(NULL)
    
    print("W calculated successfully.")
    
    newY <- Y - W %*% alpha
  }

  if (average) {
    print("Averaging newY...")
    newY <- ((1/apply(M, 2, sum)) * t(M)) %*% newY
  }

  if (!return.info) {
    return(newY)
  } else {
    return(list(newY = newY, M = M, fullalpha = fullalpha, W = W))
  }
}
```

```{r}
# Handle prps_combined column names to ensure that there are no errors in merge
colnames(prps_LP$ps_batch) <- gsub("ps_batch\\.", "", colnames(prps_LP$ps_batch))

colnames(prps_SP$ps_batch) <- gsub("ps_batch\\.", "", colnames(prps_SP$ps_batch))

colnames(prps_LO$ps_batch) <- gsub("ps_batch\\.", "", colnames(prps_LO$ps_batch))

colnames(prps_SO$ps_batch) <- gsub("ps_batch\\.", "", colnames(prps_SO$ps_batch))
```


#### PFS RUV-III precedure
```{r}
# Preparation of matrices for RUV-III
ruv_iii_input_PFS <- t(log2(cbind(
  rawCounts_expressing,
  prps_LP$ps_batch,
  prps_SP$ps_batch
) + 1)) # log2 transformed data
```

```{r}
# Filter out lines with non-sample identifiers
real_sample_names_PFS <- grep("^SRR|^ERR", row.names(ruv_iii_input_PFS), value = TRUE)

# Generate replication matrices
ruv_iii_replicate_matrix_PFS <- ruv::replicate.matrix(real_sample_names_PFS)

# Check generated matrices are normal
print(ruv_iii_replicate_matrix_PFS)
```

```{r}
BiocManager::install("RUVSeq") 
```
```{r}
library(RUVSeq)
```

```{r}
# Identify the rows that contain the pattern "gide" or "riaz"
rows_to_remove_PFS <- grepl("gide|riaz", row.names(ruv_iii_input_PFS))
```

```{r}
# Remove these rows from your matrix
ruv_iii_input_filtered_PFS <- ruv_iii_input_PFS[!rows_to_remove_PFS, ]
```

```{r}
# Now you should have correct dimensions: 228 samples
print(dim(ruv_iii_input_filtered_PFS))  # Should print 228 x 19423
```

```{r}
# Select negative control genes
ruv_iii_NCGs_PFS <- colnames(ruv_iii_input_PFS) %in% all_set_ncg[[6]]
```

```{r}
# RUV-III standardisation
ruv_iii_normal_PFS  <- RUV_III_PRPS(
  Y = ruv_iii_input_filtered_PFS,
  M = ruv_iii_replicate_matrix_PFS,
  ctl = ruv_iii_NCGs_PFS,
  k = 10,
  eta = NULL,
  return.info = TRUE)
```

```{r}
# Extraction of normalised expression matrix
ruv_iii_prps_normal_PFS <- t(ruv_iii_normal_PFS$newY[1:ncol(set_mk_se), ])
```

#### OS RUV-III precedure
```{r}
# Preparation of matrices for RUV-III
ruv_iii_input_OS <- t(log2(cbind(
  rawCounts_expressing,
  prps_LO$ps_batch,
  prps_SO$ps_batch
) + 1)) # log2 transformed data
```

```{r}
# Filter out lines with non-sample identifiers
real_sample_names_OS <- grep("^SRR|^ERR", row.names(ruv_iii_input_OS), value = TRUE)

# Generate replication matrices
ruv_iii_replicate_matrix_OS <- ruv::replicate.matrix(real_sample_names_OS)

# Check generated matrices are normal
print(ruv_iii_replicate_matrix_OS)
```

```{r}
# Identify the rows that contain the pattern "gide" or "riaz"
rows_to_remove_OS <- grepl("gide|riaz|hugo", row.names(ruv_iii_input_OS))
```

```{r}
# Remove these rows from your matrix
ruv_iii_input_filtered_OS <- ruv_iii_input_OS[!rows_to_remove_OS, ]
```

```{r}
# Now you should have correct dimensions: 228 samples
print(dim(ruv_iii_input_filtered_OS))  # Should print 228 x 19423
```

```{r}
# Select negative control genes
ruv_iii_NCGs_OS <- colnames(ruv_iii_input_OS) %in% all_set_ncg[[6]]
```

```{r}
# RUV-III standardisation
ruv_iii_normal_OS  <- RUV_III_PRPS(
  Y = ruv_iii_input_filtered_OS,
  M = ruv_iii_replicate_matrix_OS,
  ctl = ruv_iii_NCGs_OS,
  k = 10,
  eta = NULL,
  return.info = TRUE)
```

```{r}
# Extraction of normalised expression matrix
ruv_iii_prps_normal_OS <- t(ruv_iii_normal_OS$newY[1:ncol(set_mk_se), ])
```

### RLE after RUV-III
```{r}
library(ggplot2)
library(reshape2)
```

#### PFS of RLE after RUV-III
```{r}
# ‘ruv_iii_prps_normal’ is a data frame with rows for genes and columns for samples

# Calculate the median for each gene
gene_medians_PFS_after <- apply(ruv_iii_prps_normal_PFS, 1, median)

# Calculate expression difference relative to median gene (i.e., relative log expression value)
rle_values_PFS_after <- sweep(ruv_iii_prps_normal_PFS, 1, gene_medians_PFS_after)

# Converting data to long format for plotting with ggplot2
rle_values_long_PFS_after <- melt(rle_values_PFS_after)
```

```{r}
# Merge grouping information in pdata_2 directly based on row names
rle_values_long_PFS_after$PFS_Group <- pdata_2[rle_values_long_PFS_after$Var2, "PFS_Group"]
```

```{r}
# Remove NA row in PFS_Group
rle_values_long_PFS_after <- subset(rle_values_long_PFS_after, !is.na(PFS_Group))
```

```{r}
ggplot(rle_values_long_PFS_after, aes(x = Var2, y = value, color = PFS_Group)) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) + 
  labs(title = "RLE Plot with PFS NO Grouping after RUV-III", x = "Samples", y = "Relative Log Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange", size = 0.8) + 
  coord_cartesian(ylim = c(-5, 5)) +
  scale_color_manual(values = c("PFS_Long" = "darkblue", "PFS_Short" = "#C71585"))
```

#### OS of RLE after RUV-III
```{r}
# ‘ruv_iii_prps_normal’ is a data frame with rows for genes and columns for samples

# Calculate the median for each gene
gene_medians_OS_after <- apply(ruv_iii_prps_normal_OS, 1, median)

# Calculate expression difference relative to median gene (i.e., relative log expression value)
rle_values_OS_after <- sweep(ruv_iii_prps_normal_OS, 1, gene_medians_OS_after)

# Converting data to long format for plotting with ggplot2
rle_values_long_OS_after <- melt(rle_values_OS_after)
```

```{r}
# Merge grouping information in pdata_2 directly based on row names
rle_values_long_OS_after$OS_Group <- pdata_2[rle_values_long_OS_after$Var2, "OS_Group"]
```

```{r}
# Remove NA row in PFS_Group
rle_values_long_OS_after <- subset(rle_values_long_OS_after, !is.na(OS_Group))
```

```{r}
ggplot(rle_values_long_OS_after, aes(x = Var2, y = value, color = OS_Group)) +
  geom_boxplot(fill = "lightblue", alpha = 0.6, outlier.shape = NA) + 
  labs(title = "RLE Plot with OS Grouping after RUV-III", x = "Samples", y = "Relative Log Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange", size = 0.8) + 
  coord_cartesian(ylim = c(-5, 5)) +
  scale_color_manual(values = c("OS_Long" = "#8B008B", "OS_Short" = "#007373"))
```

### Whitney-Mann-U test
#### PFS
```{r}
wmu_test_PFS <- wilcox.test(rle_values_long_PFS_before$value, rle_values_long_PFS_after$value, alternative = "two.sided", paired = FALSE)

print(wmu_test_PFS)
```

#### OS
```{r}
wmu_test_OS <- wilcox.test(rle_values_long_OS_before$value, rle_values_long_OS_after$value, alternative = "two.sided", paired = FALSE)

print(wmu_test_OS)
```

### Median absolute deviation, MAD
#### PFS
```{r}
median_dev_PFS_before <- median(abs(rle_values_long_PFS_before$value))
median_dev_PFS_after <- median(abs(rle_values_long_PFS_after$value))

cat("PFS - Absolute median deviation before RLE:", median_dev_PFS_before, "\n")
cat("PFS - Absolute median deviation after RLE:", median_dev_PFS_after, "\n")
cat("PFS - Median reduction:", median_dev_PFS_before - median_dev_PFS_after, "\n\n")
```

#### OS
```{r}
median_dev_OS_before <- median(abs(rle_values_long_OS_before$value))
median_dev_OS_after <- median(abs(rle_values_long_OS_after$value))

cat("OS - Absolute median deviation before RLE:", median_dev_PFS_before, "\n")
cat("OS - Absolute median deviation after RLE:", median_dev_PFS_after, "\n")
cat("OS - Median reduction:", median_dev_OS_before - median_dev_OS_after, "\n\n")
```

### Check before DEG Analysis
```{r}
# Extract colData to data.frame
check_set_mk_se <- as.data.frame(colData(set_mk_se))

head(check_set_mk_se)

colnames(check_set_mk_se)

pfs_long_samples <- check_set_mk_se[check_set_mk_se$PFS_Group == "PFS_Long", ]

print(pfs_long_samples)

table(sample_metadata$PFS_Group, sample_metadata$PMID)
```

