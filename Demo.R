library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)

# load processed data matrices for each assay
rna <- Read10X("data/GSE126074_AdBrainCortex_rna/", gene.column = 1)
atac <- Read10X("data/GSE126074_AdBrainCortex_atac/", gene.column = 1)
fragments <- "data/fragments.sort.bed.gz"

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(seurat[["ATAC"]]) <- annotations





DefaultAssay(seurat) <- "ATAC"
seurat <- TSSEnrichment(seurat)
seurat <- NucleosomeSignal(seurat)
seurat$blacklist_fraction <- FractionCountsInRegion(
  object = seurat,
  assay = 'ATAC',
  regions = blacklist_mm10
)



Idents(seurat) <- "all"  # group all cells together, rather than by replicate
VlnPlot(
  seurat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)



seurat <- subset(
  x = seurat,
  subset = blacklist_fraction < 0.04 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)
seurat





DefaultAssay(seurat) <- "RNA"

seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- NormalizeData(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 30)
seurat <- RunUMAP(seurat, dims = 1:30, reduction.name = "umap.rna")
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5, algorithm = 3)



p1 <- DimPlot(seurat, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")

p1



DefaultAssay(seurat) <- 'ATAC'

seurat <- FindTopFeatures(seurat, min.cutoff = 10)
seurat <- RunTFIDF(seurat)
seurat <- RunSVD(seurat)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')
p2 <- DimPlot(seurat, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")



p1+p2






DefaultAssay(seurat) <- "ATAC"
CoveragePlot(seurat, region = "chr11-73192851-73193767",features = "Ctns")

rownames(seurat)[21000:22000]

nrow(seurat)
# Load necessary libraries
library(tidyverse)

# Example peaks matrix
# Replace this with your actual data
peaks_matrix <- matrix(sFL1ATACSeurat@assays$ATAC$counts)

# Convert matrix to data frame for manipulation
peaks_df <- as.data.frame(peaks_matrix)

# Add row names as a separate column
peaks_df <- peaks_df %>%
  rownames_to_column(var = "peak")

# Separate 'peak' column into chromosome, start, and end
peaks_df <- peaks_df %>%
  separate(peak, into = c("chr", "coords"), sep = ":") %>%
  separate(coords, into = c("start", "end"), sep = "-")

# Convert start and end to numeric
peaks_df <- peaks_df %>%
  mutate(across(c(start, end), as.numeric))

# Pivot longer to have one row per cell barcode per peak
bed_df <- peaks_df %>%
  pivot_longer(cols = starts_with("Cell"),
               names_to = "cell_barcode",
               values_to = "read_count") %>%
  filter(read_count > 0)  # Keep only rows with non-zero read counts

# Select and reorder columns for BED format
bed_df <- bed_df %>%
  select(chr, start, end, cell_barcode, read_count)

# Write to a BED file
write.table(bed_df, "peaks_matrix.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


remotes::install_github("10XGenomics/loupeR")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LoupeR")
