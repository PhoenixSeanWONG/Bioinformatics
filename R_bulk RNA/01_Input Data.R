# Construct SummarizedExperiment object
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DelayedArray")
library(readr)
library(tibble)
library(SummarizedExperiment)

pdata_file <- file.choose()    # ptmeta_subtype_purity.csv, patient metadata
counts_file <- file.choose()  # all_counts.tsv
fpkm_file <- file.choose() # fpkm_deseq2.csv
uq_file <- file.choose() # fpkmUQ_lastest_definition.csv

pdata <- read_csv(pdata_file)

counts <- read_tsv(counts_file)
fpkm <- read_csv(fpkm_file)
uq <- read_csv(uq_file)

pdata <- column_to_rownames(pdata, "Run")
counts <- column_to_rownames(counts, "Run")
t_counts <- as.data.frame(t(counts))
fpkm <- column_to_rownames(fpkm, "Run")
t_fpkm <- as.data.frame(t(fpkm))

uq <- column_to_rownames(uq, "Run")
t_uq <- as.data.frame(t(uq))


# this command will generate the SE object stored in a file [proj_se_latestUQ.rds]
trial_se <- SummarizedExperiment(assays=list(counts=t(t_counts), fpkm=t(t_fpkm), fpkm_uq=t(t_uq)), colData = pdata)
saveRDS(trial_se, file='proj_se_latestUQ.rds')

write_rds(trial_se, 
          "/Users/proj_se_latestUQ.rds")
