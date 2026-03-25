## PRPS
```{r}
.create_pseudo_samples_for_ls_purity_batch <- function(
    expr_data,   # Gene Expression Data Matrix
    sample_info, # Sample annotation information
    library_size, # Library Size Column Name
    biology,     # Biological labelling of listings
    batch,       # Column names for batch information
    purity,      # Column name for purity
    include_ls = FALSE,  # Whether to consider batch information to generate pseudo-samples
    include_purity = FALSE, # Whether to consider purity information to generate pseudo-samples
    min_samples_per_batch_ps = 3, # Minimum number of samples required to create a batch of pseudo-samples
    min_samples_for_purity_per_biology = 12,  # Minimum number of samples required to create a purity pseudo-sample for each biological label
    min_samples_for_purity_ps = 3,   # Minimum number of samples required to create a purity pseudo-sample
    min_samples_for_library_size_per_batch = 10,  # Minimum number of samples required to create library size pseudo-samples for each batch
    min_samples_for_library_size_ps = 3  # Minimum number of samples to create PS for library
){
  
  ### Check
  if(include_purity & min_samples_for_purity_ps > min_samples_for_purity_per_biology){
    stop('error: min_samples_for_purity_ps can not be smaller than min_samples_for_purity_per_biology')
  } else if(include_purity & min_samples_for_purity_per_biology < 2*min_samples_for_purity_ps){
    stop('error: min_samples_for_purity_per_biology should be at least two times larger than min_samples_for_purity_ps')
  } else if(include_ls & min_samples_for_library_size_ps > min_samples_for_library_size_per_batch) {
    stop('error: min_samples_for_library_size_per_batch can not be smaller than min_samples_for_library_size_ps')
  } else if(include_ls & min_samples_for_library_size_per_batch < 2*min_samples_for_library_size_ps ){
    stop('error: min_samples_for_library_size_per_batch should be at least two times larger than min_samples_for_library_size_ps')
  }

  ### Biology
  row.names(sample_info) <- colnames(expr_data)
  sample_info$biology <- apply(
    sample_info[ , biology, drop = FALSE],
    1,
    paste,
    collapse = "-"
  )

  ### Biology - Batch
  sample_info$biology_batch <- apply(
    sample_info[, c(biology, batch)],
    1,
    paste,
    collapse = "_"
  )

  ### removing batch effects
  # create PS per biology/batch
  selected_biology_ps_batch <- unlist(lapply(
    unique(sample_info$biology), 
    function(x){
      index <- sample_info$biology == x
      if(sum( table(sample_info$biology_batch[index] ) >= min_samples_per_batch_ps) > 1 ){
        x
      }
    }))
  if(length(selected_biology_ps_batch) > 0){
    message('Generated PRPS for batch effects successfully')
  }else{
    message('error: there are not enough samples to create 
            pseudo-samples for batch effects removal, you may want to lower min_samples_per_batch_ps')
  }
  sample_info_ps_batch <- sample_info[sample_info$biology %in% selected_biology_ps_batch , ]
  expr_data_ps_batch <- expr_data[, row.names(sample_info_ps_batch)]

  ### sort samples
  selected_batches <- names(which(table(sample_info_ps_batch$biology_batch) >= min_samples_per_batch_ps))
  ps_batch <- sapply(
    selected_batches,
    function(x) {
      index <- sample_info_ps_batch$biology_batch == x
      Matrix::rowMeans(expr_data_ps_batch[, index])
    })

  if(include_ls){
    selected_batches_ls <- names(
      which(table(sample_info$biology_batch) >= min_samples_for_library_size_per_batch)
    )
    if(length(selected_batches_ls) > 0){
      message('Generated PRPS for library size successfully')
      sample_info <- sample_info[
        with(sample_info,
             order(sample_info[, 'biology_batch'],
                   sample_info[, library_size])), ]
      expr_data <- expr_data[, row.names(sample_info)]
      ps_ls <- lapply(
        selected_batches_ls, 
        function(x){
          index <- sample_info$biology_batch == x
          ls_data <- expr_data[ , index]
          low_ls <- Matrix::rowMeans(ls_data[ , 1:min_samples_for_library_size_ps])
          high_ls <- Matrix::rowMeans(ls_data[ , c(ncol(ls_data)-(min_samples_for_library_size_ps - 1)):ncol(ls_data) ])
          all_ls <- cbind(low_ls, high_ls)
          colnames(all_ls) <- rep(paste(x, 'LS', sep = '_'), 2)
          all_ls
        })
      ps_ls <- do.call(cbind, ps_ls)
      
    }else{
      message('error: there are not enough samples to create pseudo-samples for 
              removal of library size effects, you may want to lower min_samples_for_library_size_per_batch')
    }
  }else if (! include_ls){
    print('PRPS is not generated for library size effects')
    ps_ls = list()
  }

  if(include_purity ){
    selected_biology_purity <- names(
      which(table(sample_info$biology) >= min_samples_for_purity_per_biology)
    ) 
    if(length(selected_biology_purity) > 0){
      message('Generated PRPS for purity effects successfully.')
      sample_info <- sample_info[
        with(sample_info,
             order(sample_info[, 'biology'],
                   sample_info[, purity])), ]
      expr_data <- expr_data[, row.names(sample_info)]
      ps_purity <- lapply(
        selected_biology_purity,
        function(x) {
          index <- sample_info$biology == x
          purity_data <- expr_data[, index]
          low_pur <- Matrix::rowMeans(purity_data[, 1:min_samples_for_purity_ps])
          high_pur <- Matrix::rowMeans(purity_data[, c(ncol(purity_data) - (min_samples_for_purity_ps - 1)):ncol(purity_data)])
          all_purity <- cbind(low_pur, high_pur)
          colnames(all_purity) <- rep(paste(x, 'purity', sep = '_'), 2)
          all_purity
        })
      ps_purity <- do.call(cbind, ps_purity)
    }else{
      message('error: there are not enough samples to make pseudo-samples 
              for purity variation, you may want to lower min_samples_for_purity_per_biology')
    }
  } else if (!include_purity){
    print('PRPS is not generated for purity effects')
    ps_purity = list()
  }
  return(list(ps_batch = ps_batch, ps_ls = ps_ls, ps_purity = ps_purity))
}
```

```{r}
samples_with_PFS <- !is.na(pdata_2$PFS_Group)

rawCounts_expressing_PFS <- rawCounts_expressing[ , samples_with_PFS]

pdata_2_PFS <- pdata_2[samples_with_PFS, ]
```

```{r}
# Select relevant columns in pdata_2 and merge them with set_sampleAnnot
set_sampleAnnot_2 <- merge(
  set_sampleAnnot, 
  pdata_2[, c("PMID", "Biosample", "PFS", "OS", "PFS_Group", "OS_Group")], 
  by = c("PMID", "Biosample"), 
  all.x = TRUE
)

head(set_sampleAnnot_2)
```

```{r}
row.names(set_sampleAnnot_2) <- set_sampleAnnot_2$Run
```

```{r}
# Split into PFS_Long and PFS_Short.
pdata_2_LP <- pdata_2[pdata_2$PFS_Group == "PFS_Long" & !is.na(pdata_2$PFS_Group), ]

pdata_2_SP <- pdata_2[pdata_2$PFS_Group == "PFS_Short" & !is.na(pdata_2$PFS_Group), ]
```

```{r}
# Filter samples of PFS_Long and PFS_Short in set_sampleAnnot_2
set_sampleAnnot_2_LP <- droplevels(set_sampleAnnot_2[set_sampleAnnot_2$PFS_Group == "PFS_Long" & !is.na(set_sampleAnnot_2$PFS_Group), ])
set_sampleAnnot_2_SP <- droplevels(set_sampleAnnot_2[set_sampleAnnot_2$PFS_Group == "PFS_Short" & !is.na(set_sampleAnnot_2$PFS_Group), ])
```

```{r}
# Filter columns for rawCounts_expressing using row names from pdata_2_LP
sample_names_LP <- rownames(pdata_2_LP)
rawCounts_expressing_LP <- rawCounts_expressing[ , colnames(rawCounts_expressing) %in% sample_names_LP]

sample_names_SP <- rownames(pdata_2_SP)
rawCounts_expressing_SP <- rawCounts_expressing[ , colnames(rawCounts_expressing) %in% sample_names_SP]
```

```{r}
# Generate pseudo-samples
prps_LP <- .create_pseudo_samples_for_ls_purity_batch(
  expr_data = rawCounts_expressing_LP,
  sample_info = set_sampleAnnot_2_LP,
  library_size = 'libSize',
  batch = 'PMID', 
  biology = 'PFS_Group',
  purity = 'tumour_purity',
  include_ls = TRUE,
  include_purity = TRUE,
  min_samples_per_batch_ps = 3,
  min_samples_for_purity_per_biology = 10,
  min_samples_for_purity_ps = 3,
  min_samples_for_library_size_per_batch = 10,
  min_samples_for_library_size_ps = 3
)
```

```{r}
# Generate pseudo-samples
prps_SP <- .create_pseudo_samples_for_ls_purity_batch(
  expr_data = rawCounts_expressing_SP,
  sample_info = set_sampleAnnot_2_SP,
  library_size = 'libSize',
  batch = 'PMID', 
  biology = 'PFS_Group',
  purity = 'tumour_purity',
  include_ls = TRUE,
  include_purity = TRUE,
  min_samples_per_batch_ps = 3,
  min_samples_for_purity_per_biology = 10,
  min_samples_for_purity_ps = 3,
  min_samples_for_library_size_per_batch = 10,
  min_samples_for_library_size_ps = 3
)
```
```{r}
prps_all <- .create_pseudo_samples_for_ls_purity_batch(
    expr_data = rawCounts_expressing,
    sample_info = set_sampleAnnot_2,
    library_size = 'libSize',
    batch = 'PMID', 
    biology = 'PFS_Group',
    purity = 'tumour_purity',
    include_ls = TRUE,
    include_purity = TRUE,
    min_samples_per_batch_ps = 3,
    min_samples_for_purity_per_biology = 10,
    min_samples_for_purity_ps = 3,
    min_samples_for_library_size_per_batch = 10,
    min_samples_for_library_size_ps = 3
)
```
```{r}
# Split into PFS_Long and PFS_Short.
pdata_2_LO <- pdata_2[pdata_2$OS_Group == "OS_Long" & !is.na(pdata_2$OS_Group), ]

pdata_2_SO <- pdata_2[pdata_2$OS_Group == "OS_Short" & !is.na(pdata_2$OS_Group), ]
```

```{r}
# Filter samples of PFS_Long and PFS_Short in set_sampleAnnot_2
set_sampleAnnot_2_LO <- droplevels(set_sampleAnnot_2[set_sampleAnnot_2$OS_Group == "OS_Long" & !is.na(set_sampleAnnot_2$OS_Group), ])
set_sampleAnnot_2_SO <- droplevels(set_sampleAnnot_2[set_sampleAnnot_2$OS_Group == "OS_Short" & !is.na(set_sampleAnnot_2$OS_Group), ])
```

```{r}
# Filter columns for rawCounts_expressing using row names from pdata_2_LP
sample_names_LO <- rownames(pdata_2_LO)
rawCounts_expressing_LO <- rawCounts_expressing[ , colnames(rawCounts_expressing) %in% sample_names_LO]

sample_names_SO <- rownames(pdata_2_SO)
rawCounts_expressing_SO <- rawCounts_expressing[ , colnames(rawCounts_expressing) %in% sample_names_SO]
```

```{r}
# Generate pseudo-samples
prps_LO <- .create_pseudo_samples_for_ls_purity_batch(
  expr_data = rawCounts_expressing_LO,
  sample_info = set_sampleAnnot_2_LO,
  library_size = 'libSize',
  batch = 'PMID', 
  biology = 'OS_Group',
  purity = 'tumour_purity',
  include_ls = TRUE,
  include_purity = TRUE,
  min_samples_per_batch_ps = 3,
  min_samples_for_purity_per_biology = 10,
  min_samples_for_purity_ps = 3,
  min_samples_for_library_size_per_batch = 10,
  min_samples_for_library_size_ps = 3
)
```

```{r}
# Generate pseudo-samples
prps_SO <- .create_pseudo_samples_for_ls_purity_batch(
  expr_data = rawCounts_expressing_SO,
  sample_info = set_sampleAnnot_2_SO,
  library_size = 'libSize',
  batch = 'PMID', 
  biology = 'OS_Group',
  purity = 'tumour_purity',
  include_ls = TRUE,
  include_purity = TRUE,
  min_samples_per_batch_ps = 3,
  min_samples_for_purity_per_biology = 10,
  min_samples_for_purity_ps = 3,
  min_samples_for_library_size_per_batch = 10,
  min_samples_for_library_size_ps = 3
)
```
