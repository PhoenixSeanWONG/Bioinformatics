## NCGs Choosing
### Installation and loading of required packages
```{r}
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DelayedArray")
```
```{r}
library(readr)
library(tibble)
library(SummarizedExperiment)
library(cowplot)
```

### data import
```{r}
pdata_csv <- file.choose()    # ptmeta_subtype_purity.csv, patient metadata
counts_csv <- file.choose()  # all_counts.tsv
```
```{r}
pdata <- read_csv(pdata_csv)
counts <- read_tsv(counts_csv)
```

### data filtered
```{r}
filtered_pdata <- pdata[pdata$PMID != "30013197", ]
```

```{r}
valid_runs <- filtered_pdata$Run
filtered_counts <- counts[, c("Run", valid_runs)]
```

```{r}
pdata <- filtered_pdata
counts <- filtered_counts
```

```{r}
# Extract Run, PFS and OS columns from final_combined_withHugo to create new dataframe
before_final_combined <- final_combined_withHugo[, c("Run", "PFS", "OS","PMID.x", "PFS_Group", "OS_Group")]
```

```{r}
# Modify the column name by renaming PMID.x to PMID
colnames(before_final_combined)[colnames(before_final_combined) == "PMID.x"] <- "PMID"
```

```{r}
# Merge dataframes, using Run columns as merge criteria
# Perform a left join, preserving all rows in pdata
# Assume both dataframes have ‘Run’ and ‘PMID’ columns.
pdata_2 <- merge(pdata, before_final_combined, by = c("Run", "PMID"), all.x = TRUE)

head(pdata_2)
```

### Data conversion and SummarizedExperiment Object creation
```{r}
pdata_2 <- column_to_rownames(pdata_2, "Run")
counts <- column_to_rownames(counts, "Run")
t_counts <- as.data.frame(t(counts))
```

```{r}
# This command will generate SE object
trial_se <- SummarizedExperiment(assays=list(counts=t(t_counts)), colData = pdata)
saveRDS(trial_se, file='proj_se_latest.rds')
```

```{r}
set_mk_se <- readRDS('proj_se_latest.rds')
set_data_se <- colnames(set_mk_se)
```

### Filter low-expressed genes
```{r}
mlnma_tissue_check <- SummarizedExperiment::colData(set_mk_se)$PMID != '30013197'
keep_set_se <- apply(
  SummarizedExperiment::assay(set_mk_se[, mlnma_tissue_check], 'counts'),
  1,
  function(x) length(x[x > 15]) >= round(.1 * sum(mlnma_tissue_check), digits = 0)
)
keep_set_se <- names(keep_set_se[keep_set_se == TRUE])
```

### Calculate library size
```{r}
set_mk_se$libSize <- log2(
  Matrix::colSums(SummarizedExperiment::assay(set_mk_se, 'counts'))
)
```

### Extract raw count data and sample information
```{r}
set_geneAnnot_ncg <- SummarizedExperiment::rowData(set_mk_se)
set_sampleAnnot <- as.data.frame(SummarizedExperiment::colData(set_mk_se))
rawCounts_expressing <- SummarizedExperiment::assay(set_mk_se, 'counts')
```

## Correlation analysis
### Definition of correlation analysis function
```{r}
.correlation_gene <- function(
    express_data, 
    is_log, 
    variable, 
    calc_method, 
    n_cores, 
    grouping)
{
  if(is_log) express_data <- express_data
  else express_data <- log2(express_data + 1)
  
  # Add noise to break ties
  set.seed(316)
  express_data <- express_data + runif(length(express_data), 0, 1e-10)
  
  rho <- parallel::mclapply(
    1:nrow(express_data),
    function(x){
      round(cor.test(
        x = express_data[x, ], 
        y = variable, 
        method = calc_method)[[4]], 6)},
    mc.cores = n_cores
  )
  p_value <- parallel::mclapply(
    1:nrow(express_data),
    function(x){
      cor.test(
        x = express_data[x, ], 
        y = variable,
        method = calc_method)[[3]]},
    mc.cores = n_cores)
  
  correlation_results <- data.frame(
    genes = row.names(express_data),
    rho = unlist(rho), 
    p_value = unlist(p_value), 
    adj_p_value = p.adjust(unlist(p_value), 'BH')
  )
  colnames(correlation_results) <- paste(
    grouping, 
    colnames(correlation_results), 
    sep = '_'
  )
  return(correlation_results)
}
```

### Calculated correlation with library size
```{r}
correlation_ls_ncg <- .correlation_gene(
  express_data  = rawCounts_expressing,
  variable  = set_sampleAnnot$libSize,
  is_log = FALSE,
  calc_method = 'spearman',
  n_cores = 1,
  grouping = 'library_size'
)
```

```{r}
set_geneAnnot_ncg$corr_ls <- correlation_ls_ncg$library_size_rho
```

### Calculated correlation with tumor purity
```{r}
correlation_tp_ncg <- .correlation_gene(
  express_data  = rawCounts_expressing,
  variable  = set_sampleAnnot$tumour_purity,
  is_log = FALSE,
  calc_method = 'spearman',
  n_cores = 1,
  grouping = 'tumor_purity'
)
```

```{r}
set_geneAnnot_ncg$corr_tp <- correlation_tp_ncg$tumor_purity_rho
```

### Calculated correlation between library size and tumor purity
```{r}
corr_ls_tp <- cor(
  set_sampleAnnot$libSize, 
  set_sampleAnnot$tumour_purity, 
  method = "spearman"
)
```

```{r}
set_geneAnnot_ncg$corr_ls_tp <- corr_ls_tp
```

## NCG Selection
### Criteria for selection of NCG
```{r}
Genes_threholds <- c(
  1000,
  2000,
  4000,
  6000,
  8000,
  10000,
  12000
)
tp_correlation <- .8
ls_correlation <- .5
```

### Selection of NCG based on correlation
#### Tumor purity
```{r}
tp_set_ncg <- lapply(
  Genes_threholds,
  function(x) {
    set_ncg <- set_geneAnnot_ncg$corr_tp > tp_correlation
    return(set_ncg)
  })
```

#### Library size
```{r}
ls_set_ncg <- lapply(
  Genes_threholds,
  function(x) {
    set_ncg <- set_geneAnnot_ncg$corr_ls > ls_correlation
    return(set_ncg)
  })
```

### 
```{r}
all_counts_data <- read.table("/Users/wangxin/Desktop/BINF90008/jess_code_v1/all_counts.tsv", header = TRUE, row.names = 1, sep = "\t")

set_geneAnnot_ncg$gene_name <- rownames(all_counts_data)

all_set_ncg <- lapply(
  c(1:length(Genes_threholds)),
  function(x) {
    set_ncg <- c(
      set_geneAnnot_ncg$gene_name[ls_set_ncg[[x]]],
      set_geneAnnot_ncg$gene_name[tp_set_ncg[[x]]]
    )
    return(set_ncg)
  })
```

### Define PCA analysis function
```{r}
.pca <- function(data, is_log) {
  if (is_log)
    data <- data
  else
    data <- log2(data + 1)
  svd <- base::svd(scale(
    x = t(data),
    center = TRUE,
    scale = FALSE
  ))
  percent <- svd$d ^ 2 / sum(svd$d ^ 2) * 100
  percent <-
    sapply(seq_along(percent),
           function(i) {
             round(percent[i], 1)
           })
  return(list(
    svd_value = svd,
    variation = percent))
}
```

```{r}
perform_pca_analysis <- function(data, is_log = FALSE) {
  .pca(
    data = data,
    is_log = is_log
  )
}
```

### Define PCA scatter plot generation function
```{r}
.scatter_density_pca <- function(
    pcs, 
    pcs_var, 
    group_name, 
    group, 
    color, 
    strokeSize, 
    pointSize, 
    strokeColor,
    alpha,
    title
){
  pair.pcs <- utils::combn(ncol(pcs), 2)
  pList <- list()
  
  # Determining whether a group is continuous or discrete Determining whether a group is continuous or discrete
  if (is.numeric(group)) {
    color_scale <- scale_color_gradient2(low = "green", mid = "white", high = "yellow")
  } else {
    color_scale <- scale_fill_manual(values = color, name = group_name)
  }
  
  for(i in 1:ncol(pair.pcs)){
    x <- pair.pcs[1,i]
    y <- pair.pcs[2,i]
    p <- ggplot(mapping = aes(
        x = pcs[,x], 
        y = pcs[,y], 
        fill = group)) +  
      xlab(paste0('PC', x, ' (', pcs_var[x], '%)')) +
      ylab(paste0('PC', y, ' (', pcs_var[y], '%)')) +
      geom_point(
        pch = 21,  
        color = strokeColor, 
        stroke = strokeSize, 
        size = pointSize,
        alpha = alpha) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      ggtitle(title) +
      theme(
        legend.position = "right",
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1.1),
        legend.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
      guides(fill = guide_legend(override.aes = list(size = 4))) +
      color_scale
    
    if(i == 1){
      le <- ggpubr::get_legend(p)
    }
    
    xdens <- cowplot::axis_canvas(p, axis = "x") +
      geom_density(
        mapping = aes(
          x = pcs[,x], 
          fill = group),  
        alpha = 0.5,  # Set transparency to ensure fill effect
        size = 0.2
      ) +
      theme(legend.position = "none") +
      color_scale
    
    ydens <- cowplot::axis_canvas(
      p, 
      axis = "y", 
      coord_flip = TRUE) +
      geom_density(
        mapping = aes(
          x = pcs[,y],
          fill = group),  
        alpha = 0.5,
        size = 0.2) +
      theme(legend.position = "none") +
      color_scale +
      coord_flip()
    
    p1 <- insert_xaxis_grob(
      p,
      xdens,
      grid::unit(.2, "null"),
      position = "top"
    )
    p2 <- insert_yaxis_grob(
      p1,
      ydens,
      grid::unit(.2, "null"),
      position = "right"
    )
    pList[[i]] <- ggdraw(p2)
  }
  
  pList[[i+1]] <- le
  return(pList)
}
```

```{r}
generate_pca_plot <- function(
    pca_result, 
    group_name, 
    group, 
    color, 
    title
  ) {
  
  # Determine whether a group is continuous or discrete
if (is.numeric(group)) {
  color_scale <- scale_color_gradient2(low = "green", mid = "white", high = "yellow")
} else {
  color_scale <- scale_fill_manual(values = color, name = group_name)
}
  
  p <- .scatter_density_pca(
    pcs = pca_result$svd_value$u[, 1:3],
    pcs_var = pca_result$variation,
    group_name = group_name,
    group = group,
    color = color,
    strokeSize = .2,
    pointSize = 2,
    strokeColor = 'gray30',
    alpha = .6,
    title = title
  )
  
   # Apply color scale to each plot in list
  pp <- lapply(p, function(plot) {
    plot + color_scale
  })
  
  return(pp)
}
```

### Perform PCA analysis for library size and tumor purity
```{r}
ls_tp_pca <- perform_pca_analysis(
  data = rawCounts_expressing,
  is_log = FALSE
)
```

### Define color palettes
```{r}
library(ggplot2)
```

```{r}
library_size_colors <- scale_color_manual(values = c("red",  "blue", "purple", "orange"))
purity_colors <- scale_color_gradient2(low = "green", mid = "white", high = "yellow")
subtype_colors <- scale_color_manual(values = c("red", "green", "blue", "purple", "orange"))
```

### Generate PCA scatter plots
#### Library size
```{r}
#pca_plot_ls <- generate_pca_plot(
#  pca_result = ls_tp_pca,
#  group_name = 'Library Size',
#  group = set_sampleAnnot$libSize,
#  color = NULL,
#  title = 'Library Size PCA'
#)
```
```{r}
lb_colours <- c(
  'red',
  'cyan3', 
  'darkorange') 
names(lb_colours) <- c(
  "gide",
  "hugo",
  "riaz"
)
```

```{r}
pca_plot_ls <- generate_pca_plot(
  pca_result = ls_tp_pca,
  group_name = 'Library Size',
  group = set_sampleAnnot$PMID,
  color = lb_colours,
  title = 'Library Size PCA'
)
```
```{r}
pca_plot_ls
```

#### Tumor purity
```{r}
pca_plot_tp <- generate_pca_plot(
  pca_result = ls_tp_pca,
  group_name = 'Tumor Purity',
  group = set_sampleAnnot$tumour_purity,
  color = NULL,
  title = 'Tumor Purity PCA'
)
```
```{r}
pca_plot_tp
```

#### Subtype
```{r}
subtype_colours <- c(
  'red',
  'darkgreen',
  'navy',
  'cyan3', 
  'darkorange',
  'purple') 
names(subtype_colours) <- c(
  "1",
  "2",
  "3",
  "0", 
  "M",
  "NA")
```
```{r}
pca_plot_st <- generate_pca_plot(
  pca_result = ls_tp_pca,
  group_name = 'Subtype',
  group = set_sampleAnnot$subtypes,
  color = subtype_colours,
  title = 'Subtype PCA'
)
```
```{r}
pca_plot_st
```


### Visualize the PCA results side by side
```{r}
visualize_pca_results <- function(
    pca_plot_ls, 
    pca_plot_tp,
    pca_plot_st
) {
  gridExtra::grid.arrange(
    pca_plot_ls[[1]],
    pca_plot_tp[[1]],
    pca_plot_st[[1]],
    pca_plot_ls[[2]],
    pca_plot_tp[[2]],
    pca_plot_st[[2]],
    pca_plot_ls[[3]],
    pca_plot_tp[[3]],
    pca_plot_st[[3]],
    ncol = 3
  )
}
```
```{r}
visualize_pca_results(
  pca_plot_ls,
  pca_plot_tp,
  pca_plot_st)
```

## Assessment of PCA analysis on NCGs
### NCG function
```{r}
ncg_pca <- perform_pca_analysis(
  data = rawCounts_expressing[all_set_ncg[[6]], ],
  is_log = FALSE
)
```

### Generate PCA Scatter Plots
#### Library Size
```{r}
pca_plot_ls_ncg <- generate_pca_plot(
  pca_result = ncg_pca,
  group_name = 'Library Size',
  group = set_sampleAnnot$PMID,
  color = lb_colours,
  title = 'NCG PCA by Library Size'
)
```
```{r}
pca_plot_ls_ncg
```

#### Tumor Purity
```{r}
pca_plot_tp_ncg <- generate_pca_plot(
  pca_result = ncg_pca,
  group_name = 'Tumor Purity',
  group = set_sampleAnnot$tumour_purity,
  color = NULL,
  title = 'NCG PCA by Tumor Purity'
)
```
```{r}
pca_plot_tp_ncg
```

#### Subtype
```{r}
pca_plot_st_ncg <- generate_pca_plot(
  pca_result = ncg_pca,
  group_name = 'Subtype',
  group = set_sampleAnnot$subtypes,
  color = subtype_colours,
  title = 'NCG PCA by Subtype'
)
```
```{r}
pca_plot_st_ncg
```

### Visualize the PCA Results Side by Side
```{r}
visualize_pca_results <- function(
    pca_plot_ls_ncg, 
    pca_plot_tp_ncg, 
    pca_plot_st_ncg) {
  gridExtra::grid.arrange(
    pca_plot_ls_ncg[[1]],
    pca_plot_tp_ncg[[1]],
    pca_plot_st_ncg[[1]],
    pca_plot_ls_ncg[[2]],
    pca_plot_tp_ncg[[2]],
    pca_plot_st_ncg[[2]],
    pca_plot_ls_ncg[[3]],
    pca_plot_tp_ncg[[3]],
    pca_plot_st_ncg[[3]],
    ncol = 3
  )
}
```
```{r}
visualize_pca_results(
  pca_plot_ls_ncg, 
  pca_plot_tp_ncg, 
  pca_plot_st_ncg)
```
