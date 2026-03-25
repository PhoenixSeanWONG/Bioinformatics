### Pathway studays
#### PFS
```{r}
library(clusterProfiler)
library(org.Hs.eg.db) # Use the appropriate annotation package for your species (e.g., mouse, zebrafish)
library(ReactomePA)  # For Reactome pathway enrichment
```
```{r}
# Convert gene names to Entrez IDs (if needed)
entrez_DEG_PFS <- bitr(
  significant_genes_DEG_PFS_results_names, 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db)

# Run KEGG pathway enrichment analysis
kegg_result_DEG_PFS <- enrichKEGG(
  gene = entrez_DEG_PFS$ENTREZID,
  organism = 'hsa',           # 'hsa' is for humans; change to 'mmu' for mouse, etc.
  pvalueCutoff = 0.1
)

head(kegg_result_DEG_PFS)
```

```{r}
# Visualize enriched pathways
dotplot(kegg_result_DEG_PFS, showCategory = 20)   # Dotplot for KEGG pathways
```

```{r}
# Run Reactome pathway enrichment analysis
reactome_result_DEG_PFS <- enrichPathway(
  gene = entrez_DEG_PFS$ENTREZID,
  organism = 'human',         # Replace 'human' with the correct species name if needed
  pvalueCutoff = 0.05
)

head(reactome_result_DEG_PFS)
```

```{r}
dotplot(reactome_result_DEG_PFS, showCategory = 20, font.size = 8) # Dotplot for Reactome pathways
```

```{r}
# Perform GO enrichment analysis for Biological Process (BP)
go_bp_DEG_PFS <- enrichGO(
  gene = entrez_DEG_PFS$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",                # Options: BP (Biological Process), CC (Cellular Component), MF (Molecular Function)
  pvalueCutoff = 0.05
)
```

```{r}
# Visualize top enriched GO terms
dotplot(go_bp_DEG_PFS, showCategory = 20, font.size = 8) # Dotplot for enriched GO terms
```

```{r}
library(pathview)
```

```{r}
DEG_PFS_gene_data <- DEG_PFS_results$logFC

names(DEG_PFS_gene_data) <- rownames(DEG_PFS_results) 

DEG_PFS_gene_data <- DEG_PFS_gene_data[entrez_DEG_PFS$SYMBOL] 

names(DEG_PFS_gene_data) <- entrez_DEG_PFS$ENTREZID
```

```{r}
# Select pathway of interest
pathway_id_DEG_PFS <- kegg_result_DEG_PFS@result$ID[1]

pathview(
  gene.data = DEG_PFS_gene_data,
  pathway.id = pathway_id_DEG_PFS,
  species = "hsa",
  out.suffix = "PFS_DEG",
  kegg.native = TRUE,
  low = list(gene = "navy"),
  mid = list(gene = "yellow"),
  high = list(gene = "#FF00FF")
)
```

```{r}
target_genes <- c("MS4A1", "CR2", "CD19")

expression_data_PFS_DEG_target <- expression_matrix_ruv_iii_prps_normal_PFS_DEG[target_genes, ]

print(expression_data_PFS_DEG_target)
```

```{r}
expression_data_PFS_DEG_long <- as.data.frame(expression_data_PFS_DEG_target)

expression_data_PFS_DEG_long$Gene <- rownames(expression_data_PFS_DEG_long)

expression_data_PFS_DEG_long <- tidyr::pivot_longer(expression_data_PFS_DEG_long, cols = -Gene, names_to = "Sample", values_to = "Expression")

expression_data_PFS_DEG_long <- merge(expression_data_PFS_DEG_long, filtered_ruv_iii_prps_normal_PFS_DEG[, c("Run", "PFS_Group")], by.x = "Sample", by.y = "Run")
```

```{r}
library(ggplot2)
```
```{r}
ggplot(expression_data_PFS_DEG_long, aes(x = PFS_Group, y = Expression, fill = PFS_Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Gene, scales = "free") +
  theme_minimal() +
  labs(title = "Expression Levels of MS4A1, CR2, CD19", x = "PFS Group", y = "Expression")
```
