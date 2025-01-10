source("Libraries.R")

source("Functions.R")



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

seuratSFL1$CellType <- seuratSFL1$DF.classifications_0.25_0.3_563

seuratSFL2 <- FindDoublets(seuratSFL2)

seuratSFL2$CellType <- seuratSFL2$DF.classifications_0.25_0.1_428

seuratSKM1 <- FindDoublets(seuratSKM1)

seuratSKM1$CellType <- seuratSKM1$DF.classifications_0.25_0.13_430

seuratSKM2 <- FindDoublets(seuratSKM2)

seuratSKM2$CellType <- seuratSKM2$DF.classifications_0.25_0.22_481





SFL1Plot <- DimPlot(seuratSFL1,reduction = "umap",group.by = "CellType") + ggtitle("SFL1")

SFL2Plot <- DimPlot(seuratSFL2,reduction = "umap",group.by = "CellType") + ggtitle("SFL2")

SKM1Plot <- DimPlot(seuratSKM1,reduction = "umap",group.by = "CellType") + ggtitle("SKM1")

SKM2Plot <- DimPlot(seuratSKM2,reduction = "umap",group.by = "CellType") + ggtitle("SKM2")



SFL1Plot+SFL2Plot+SKM1Plot+SKM2Plot











seurat <- merge(seuratSFL1,y=c(seuratSFL2,seuratSKM1,seuratSKM2), project = "MultiOmicsMacrophage")







seurattest <- DoSignac(seurat = seurat)







seurattest@meta.data <- seurattest@meta.data %>%
  
  dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))



p1 <- DimPlot(seurattest, reduction = 'umap.rna', label = TRUE, split.by = "SampleType") + ggtitle("RNA UMAP")

p2 <- DimPlot(seurattest, reduction = 'umap.atac', label = TRUE, split.by = "SampleType") + ggtitle("ATAC UMAP")



p1+p2





SaveSeuratRds(seurattest, file = "data/BethSeurat.RDS")



seurattest <- readRDS("data/BethSeurat.RDS")


seurattest <- VisualizeClustree(0,1,0.1,seurattest)


seuratJoined <- JoinLayers(seurattest)

counts <- seuratJoined@assays$RNA$counts

new_clusters = testClusters(counts, 
                            cluster_ids = as.character(seuratJoined$cca_clusters_0_1_0.2),
                            alpha = 0.05, #FWER control, can be relaxed if needed
                            num_features = 2500,#default number
                            num_PCs = 30, #default number
                            parallel = FALSE, #can set to FALSE
                            cores = 1)
table(new_clusters[[1]],seuratJoined$cca_clusters_0_1_0.2)
class(new_clusters)
new_clusters[[2]]#
seuratJoined$Seurat <- seuratJoined$cca_clusters_0_1_0.2
seuratJoined$Corrected_Seurat <- new_clusters[[1]]
ggarrange(#DimPlot(scnorm_nocc_Pfncntss,group.by='scSHC')+NoLegend(),
  DimPlot(seuratJoined,group.by='Seurat', label=T)+NoLegend(),
  DimPlot(seuratJoined,group.by='Corrected_Seurat', label=T)+NoLegend(),
  nrow=1)
d1 = DimPlot(seuratJoined,group.by='Seurat', label=T) #+
d2 = DimPlot(seuratJoined,group.by='Corrected_Seurat', label=T) #+
#NoLegend()
d1 + d2
ggsave("compare_clusters_Seurat.png", width=10, height=8, units="in", dpi=588)
ggsave("compare_clusters_Seurat.pdf", width=10, height=8)


seuratJoined$SampleType

DimPlot(seuratJoined,group.by='Corrected_Seurat', label=T, split.by = "orig.ident")

SFL1New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL1" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
SFL2New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL2" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
SKM1New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM1" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
SKM2New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM2" & seuratJoined@meta.data$Corrected_Seurat == "new1",])


SFL1New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL1" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
SFL2New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL2" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
SKM1New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM1" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
SKM2New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM2" & seuratJoined@meta.data$Corrected_Seurat == "new2",])


cellCountDf <- data.frame(
  SFL1 = c(SFL1New1Count,SFL1New2Count),
  SFL2 = c(SFL2New1Count,SFL2New2Count),
  SKM1 = c(SKM1New1Count,SKM1New2Count),
  SKM2 = c(SKM2New1Count,SKM2New2Count)
)

rownames(cellCountDf) <- c("new1","new2")

normalizedCellCountDf <- sweep(cellCountDf, 2, colSums(cellCountDf), FUN = "/")

normalizedCellCountDf$ClusterType <- rownames(normalizedCellCountDf)
normalizedCellCountDf[] <- lapply(normalizedCellCountDf, as.character)
normalizedCellCountDfLong <- pivot_longer(normalizedCellCountDf, cols = -"ClusterType")
colnames(normalizedCellCountDfLong) <- c("ClusterType","SampleType","CellCountProportion")

Plot1 <- ggplot(normalizedCellCountDfLong[grep("1",normalizedCellCountDfLong$SampleType),], aes(x = SampleType, y = CellCountProportion, fill = factor(ClusterType))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster Type", y = "Cell Count Proportion", fill = "Sample Type") +
  theme_minimal()

Plot2 <- ggplot(normalizedCellCountDfLong[grep("2",normalizedCellCountDfLong$SampleType),], aes(x = SampleType, y = CellCountProportion, fill = factor(ClusterType))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster Type", y = "Cell Count Proportion", fill = "Sample Type") +
  theme_minimal()

Plot1/Plot2


Idents(seuratJoined) <- "Corrected_Seurat"




allMarkers <- FindAllMarkers(seuratJoined)











clusterNew1seurat <- subset(seuratJoined,idents="new1")
clusterNew2seurat <- subset(seuratJoined,idents="new2")






clusterNew1Markers <- FindMarkers(
  
  object = clusterNew1seurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  logfc.threshold = 1,
  
  test.use = "wilcox"
  
)


clusterNew2Markers <- FindMarkers(
  
  object = clusterNew2seurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  logfc.threshold = 1,
  
  test.use = "wilcox"
  
)

geneName <- "Erich6"

geneRegions <- GetGenomeRanges(seuratJoined,geneName)

GenerateCoveragePlots(seuratJoined, geneRegions,geneName)

Idents(seuratJoined) <- "SampleType"
ranges.show <- StringToGRanges(geneRegions[15:20])
CoveragePlot(seuratJoined,assay = "ATAC",
             region = geneRegions[15:20],
             features = "Cmss1", 
             expression.assay = "RNA",
             region.highlight = ranges.show, 
             split.by = "Corrected_Seurat")

clusterNew1seurat@assays$ATAC@annotation[clusterNew1seurat@assays$ATAC@annotation$gene_name == "Siah2",]

Idents(clusterNew1seurat) <- "SampleType"

clusterNew1seuratWT <- subset(clusterNew1seurat, idents = "WildType")

clusterNew1seuratKO <- subset(clusterNew1seurat, idents = "KnockOut")





DefaultAssay(clusterNew1seuratWT) <- "ATAC"



p1 <- CoveragePlot(clusterNew1seuratWT,assay = "ATAC",region = c("chr3-58688799-58689631","chr3-58690316-58691181","chr3-58691816-58692742"),features = "Siah2", expression.assay = "RNA")



DefaultAssay(clusterNew1seuratKO) <- "ATAC"



p2 <- CoveragePlot(clusterNew1seuratKO,assay = "ATAC",region = c("chr3-58688799-58689631","chr3-58690316-58691181","chr3-58691816-58692742"),features = "Siah2", expression.assay = "RNA")



p1/p2





DefaultAssay(cluster0seurat) <- "ATAC"



CoveragePlot(cluster0seurat,region = "chr13-115088530-115088961",features = "Pelo")







rownames(clusterNew1seuratWT)[24000:25000]













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