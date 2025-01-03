source("Libraries.R")

source("Functions.R")


source("Demo.R")


#Read the file paths

files <- list.files("data/matrix_files/",full.names = TRUE)

filePaths <- files

fileIndexes <- c(1,3,5,7)

fileNames <- c(
  
  "1" = "SFL1",
  
  "3" = "SFL2",
  
  "5" = "SKM1",
  
  "7" = "SKM2"
  
)



#Build seurat object SFL1

data  <- Read10X_h5(filename = filePaths[1])

fragments <- "data/SFL1/outs/atac_fragments.tsv.gz"

ATACSeurat <- CreateSeuratObject(counts = data[["Peaks"]], assay = "ATAC", project = fileNames[as.character(1)])

seuratSFL1 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(1)])

seuratSFL1[['ATAC']] <- CreateChromatinAssay(
  
  counts = ATACSeurat@assays$ATAC$counts,
  
  sep = c(":", "-"),
  
  genome = "mm10",
  
  fragments = fragments
  
)





#Build seurat object SFL2

data  <- Read10X_h5(filename = filePaths[3])

fragments <- "data/SFL2/outs/atac_fragments.tsv.gz"

ATACSeurat <- CreateSeuratObject(counts = data[["Peaks"]], assay = "ATAC", project = fileNames[as.character(3)])

seuratSFL2 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(3)])

seuratSFL2[['ATAC']] <- CreateChromatinAssay(
  
  counts = ATACSeurat@assays$ATAC$counts,
  
  sep = c(":", "-"),
  
  genome = "mm10",
  
  fragments = fragments
  
)







#Build seurat object SKM1

data  <- Read10X_h5(filename = filePaths[5])

fragments <- "data/SKM1/outs/atac_fragments.tsv.gz"

ATACSeurat <- CreateSeuratObject(counts = data[["Peaks"]], assay = "ATAC", project = fileNames[as.character(5)])

seuratSKM1 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(5)])

seuratSKM1[['ATAC']] <- CreateChromatinAssay(
  
  counts = ATACSeurat@assays$ATAC$counts,
  
  sep = c(":", "-"),
  
  genome = "mm10",
  
  fragments = fragments
  
)





#Build seurat object SKM2

data  <- Read10X_h5(filename = filePaths[7])

fragments <- "data/SKM2/outs/atac_fragments.tsv.gz"

ATACSeurat <- CreateSeuratObject(counts = data[["Peaks"]], assay = "ATAC", project = fileNames[as.character(7)])

seuratSKM2 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(7)])

seuratSKM2[['ATAC']] <- CreateChromatinAssay(
  
  counts = ATACSeurat@assays$ATAC$counts,
  
  sep = c(":", "-"),
  
  genome = "mm10",
  
  fragments = fragments
  
)







# extract gene annotations from EnsDb

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)



# change to UCSC style since the data was mapped to mm10

seqlevels(annotations) <- paste0('chr', seqlevels(annotations))

genome(annotations) <- "mm10"





#Add gene annotations to all individual seurat objects

seuratSFL1 <- AddAnnotations(seuratSFL1,annotations)

seuratSFL2 <- AddAnnotations(seuratSFL2,annotations)

seuratSKM1 <- AddAnnotations(seuratSKM1,annotations)

seuratSKM2 <- AddAnnotations(seuratSKM2,annotations)









seuratSFL1 <- FindDoublets(seuratSFL1)

seuratSFL1$CellType <- seuratSFL1$DF.classifications_0.25_0.3_565

seuratSFL2 <- FindDoublets(seuratSFL2)

seuratSFL2$CellType <- seuratSFL2$DF.classifications_0.25_0.1_428

seuratSKM1 <- FindDoublets(seuratSKM1)

seuratSKM1$CellType <- seuratSKM1$DF.classifications_0.25_0.13_430

seuratSKM2 <- FindDoublets(seuratSKM2)

seuratSKM2$CellType <- seuratSKM2$DF.classifications_0.25_0.18_481





SFL1Plot <- DimPlot(seuratSFL1,reduction = "umap",group.by = "CellType") + ggtitle("SFL1")

SFL2Plot <- DimPlot(seuratSFL2,reduction = "umap",group.by = "CellType") + ggtitle("SFL2")

SKM1Plot <- DimPlot(seuratSKM1,reduction = "umap",group.by = "CellType") + ggtitle("SKM1")

SKM2Plot <- DimPlot(seuratSKM2,reduction = "umap",group.by = "CellType") + ggtitle("SKM2")



SFL1Plot+SFL2Plot+SKM1Plot+SKM2Plot











seurat <- merge(seuratSFL1,y=c(seuratSFL2,seuratSKM1,seuratSKM2), project = "MultiOmics")







seurattest <- DoSignac(seurat = seurat)







seurattest@meta.data <- seurattest@meta.data %>%
  
  dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))



p1 <- DimPlot(seurattest, reduction = 'umap.rna', label = TRUE, split.by = "SampleType") + ggtitle("RNA UMAP")

p2 <- DimPlot(seurattest, reduction = 'umap.atac', label = TRUE, split.by = "SampleType") + ggtitle("ATAC UMAP")



p1+p2





SaveSeuratRds(seurattest, file = "data/BethSeurat.RDS")



seurattest <- readRDS("data/BethSeurat.RDS")



install.packages("clustree")



DefaultAssay(seurattest) <- "RNA"

seurattest <- FindClusters(seurattest, resolution = 0.3, algorithm = 3,cluster.name = "cca_clusters_0.3")

seurattest <- FindClusters(seurattest, resolution = 0.7, algorithm = 3,cluster.name = "cca_clusters_0.7")



# Create a data frame with cluster assignments for each resolution

cluster_data <- data.frame(
  
  cell_id = colnames(seurattest),
  
  resolution_0.3 = seurattest$cca_clusters_0.3,  # Cluster assignments at resolution 0.3
  
  resolution_0.5 = seurattest$cca_clusters,  # Cluster assignments at resolution 0.5
  
  resolution_0.7 = seurattest$cca_clusters_0.7   # Cluster assignments at resolution 0.7
  
)





# Convert to long format

cluster_data_long <- cluster_data %>%
  
  gather(key = "resolution", value = "cluster", -cell_id) %>%
  
  arrange(resolution, cell_id)  # Make sure the data is sorted



# Use clustree to visualize cluster stability across resolutions

clustree(cluster_data, prefix = "resolution_")



data(nba_clusts)

clustree(nba_clusts, prefix = "K")



head(seurattest@meta.data)



allMarkers <- FindAllMarkers(seurattest)



seurattest <- JoinLayers(seurattest)



Idents(seurattest) <- "seurat_clusters"

Idents(seurattest) <- "SampleType"





Idents(seurattest)



DefaultAssay(seurattest) <- "RNA"



markers <- FindAllMarkers(
  
  object = seurattest, 
  
  test.use = "wilcox",
  
  verbose = TRUE
  
)







cluster0seurat <- subset(seurattest,idents="0")







cluster0Markers <- FindMarkers(
  
  object = cluster0seurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType"
  
)





cluster0seurat@assays$ATAC@annotation[cluster0seurat@assays$ATAC@annotation$gene_name == "Pelo",]



Idents(cluster0seurat) <- "SampleType"



cluster0seuratWT <- subset(cluster0seurat, idents = "WildType")



cluster0seuratKO <- subset(cluster0seurat, idents = "KnockOut")





DefaultAssay(cluster0seuratWT) <- "ATAC"



CoveragePlot(cluster0seuratWT,region = "chr16-57326609-57327532",features = "Cmss1")



DefaultAssay(cluster0seuratKO) <- "ATAC"



CoveragePlot(cluster0seuratKO,region = "chr16-57326609-57327532",features = "Cmss1")



DefaultAssay(cluster0seurat) <- "ATAC"



CoveragePlot(cluster0seurat,region = "chr13-115088530-115088961",features = "Pelo")







rownames(cluster0seuratWT)[24000:25000]













seuratSFL1 <- DoSignac(seuratSFL1,annotations)

seuratSFL2 <- DoSignac(seuratSFL2,annotations)

seuratSKM1 <- DoSignac(seuratSKM1,annotations)

seuratSKM2 <- DoSignac(seuratSKM2,annotations)



p1 <- DimPlot(seuratSFL1, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SFL1")

p2 <- DimPlot(seuratSFL1, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SFL1")



SFL1Plots <- p1 + p2



SFL1Plots





p1 <- DimPlot(seuratSFL2, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SFL2")

p2 <- DimPlot(seuratSFL2, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SFL2")



SFL2Plots <- p1 + p2



SFL2Plots







p1 <- DimPlot(seuratSKM1, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SKM1")

p2 <- DimPlot(seuratSKM1, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SKM1")



SKM1Plots <- p1 + p2



SKM1Plots





p1 <- DimPlot(seuratSKM2, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SKM2")

p2 <- DimPlot(seuratSKM2, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SKM2")



SKM2Plots <- p1 + p2



SKM2Plots





ggsave("data/SFL1.jpeg",SFL1Plots)

ggsave("data/SFL2.jpeg",SFL2Plots)

ggsave("data/SKM1.jpeg",SKM1Plots)

ggsave("data/SKM2.jpeg",SKM2Plots)



seuratSFL1@meta.data



# add the gene information to the object

Annotation(seurat[['ATAC']]) <- annotations



ATACregions <- rownames(ATACSeurat)

colnames(ATACSeurat)



rownames(RNASeurat)

colnames(RNASeurat)



bethSeurat <- CreateProjectSeurat(filePaths = filePaths,fileIndexes = fileIndexes,fileNames = fileNames)



bethSeurat@meta.data <- bethSeurat@meta.data %>%
  
  dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))



bethSeurat@meta.data



DefaultAssay(bethSeurat) <- 'ATAC'



bethSeurat <- FindTopFeatures(bethSeurat, min.cutoff = 10)

bethSeurat <- RunTFIDF(bethSeurat)

bethSeurat <- RunSVD(bethSeurat)

bethSeurat <- RunUMAP(bethSeurat, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')

p2 <- DimPlot(bethSeurat, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")



bethSeurat <- DoStandardSingleCellAnalysis(bethSeurat)





DimPlot(bethSeurat, reduction = "uamp.unintegrated",split.by = "SampleType", label = TRUE)





bethSeurat <- IntegrateSeuratObjects(bethSeurat)





bethSeurat <- PerformAfterIntegrationStandardAnalysis(bethSeurat)





plot1 <- DimPlot(bethSeurat, reduction = "uamp.unintegrated",split.by = "SampleType",label = TRUE)

plot2 <- DimPlot(bethSeurat, reduction = "umap.cca", label = TRUE, repel = TRUE, split.by = "SampleType")

plot1+plot2



plot2



rownames(seuratSFL1)



CoveragePlot(seuratSFL1, region = "chr1-37504883-37505816",features = "Mgat4a")





head(seurat@meta.data)

head(seuratSFL1@meta.data)



rownames(seurat)

rownames(seuratSFL1)