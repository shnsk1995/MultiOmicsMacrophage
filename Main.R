source("Libraries.R")
source("Functions.R")

####################################################################Seurat Objects##########################################################################

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



#Process ATAC data first to create a unified peak set across all the samples

#Extract genome ranges from peaks of each sample
peaksSFL1 <- read.table(
  file = "data/SFL1/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
grSFL1 <- makeGRangesFromDataFrame(peaksSFL1)


peaksSFL2 <- read.table(
  file = "data/SFL2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
grSFL2 <- makeGRangesFromDataFrame(peaksSFL2)


peaksSKM1 <- read.table(
  file = "data/SKM1/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
grSKM1 <- makeGRangesFromDataFrame(peaksSKM1)


peaksSKM2 <- read.table(
  file = "data/SKM2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
grSKM2 <- makeGRangesFromDataFrame(peaksSKM2)



#Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(grSFL1,grSFL2,grSKM1,grSKM2))



# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]


# load metadata for each sample
mdSFL1 <- read.table(
  file = "data/SFL1/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

mdSFL2 <- read.table(
  file = "data/SFL2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

mdSKM1 <- read.table(
  file = "data/SKM1/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

mdSKM2 <- read.table(
  file = "data/SKM2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row



# perform an initial filtering of low count cells
mdSFL1 <- mdSFL1[mdSFL1$atac_raw_reads > 500, ]
mdSFL2 <- mdSFL2[mdSFL2$atac_raw_reads > 500, ]
mdSKM1 <- mdSKM1[mdSKM1$atac_raw_reads > 500, ]
mdSKM2 <- mdSKM2[mdSKM2$atac_raw_reads > 500, ]


# create fragment objects
fragsSFL1 <- CreateFragmentObject(
  path = "data/SFL1/outs/atac_fragments.tsv.gz",
  cells = rownames(mdSFL1)
)

fragsSFL2 <- CreateFragmentObject(
  path = "data/SFL2/outs/atac_fragments.tsv.gz",
  cells = rownames(mdSFL2)
)

fragsSKM1 <- CreateFragmentObject(
  path = "data/SKM1/outs/atac_fragments.tsv.gz",
  cells = rownames(mdSKM1)
)

fragsSKM2 <- CreateFragmentObject(
  path = "data/SKM2/outs/atac_fragments.tsv.gz",
  cells = rownames(mdSKM2)
)



#Create Counts objects
countsSFL1 <- FeatureMatrix(
  fragments = fragsSFL1,
  features = combined.peaks,
  cells = rownames(mdSFL1)
)

countsSFL2 <- FeatureMatrix(
  fragments = fragsSFL2,
  features = combined.peaks,
  cells = rownames(mdSFL2)
)

countsSKM1 <- FeatureMatrix(
  fragments = fragsSKM1,
  features = combined.peaks,
  cells = rownames(mdSKM1)
)

countsSKM2 <- FeatureMatrix(
  fragments = fragsSKM2,
  features = combined.peaks,
  cells = rownames(mdSKM2)
)



#Finally create assays and respective seurat objects
assaySFL1 <- CreateChromatinAssay(countsSFL1, fragments = fragsSFL1)
atacseuratSFL1 <- CreateSeuratObject(assaySFL1, assay = "ATAC", meta.data=mdSFL1)

assaySFL2 <- CreateChromatinAssay(countsSFL2, fragments = fragsSFL2)
atacseuratSFL2 <- CreateSeuratObject(assaySFL2, assay = "ATAC", meta.data=mdSFL2)

assaySKM1 <- CreateChromatinAssay(countsSKM1, fragments = fragsSKM1)
atacseuratSKM1 <- CreateSeuratObject(assaySKM1, assay = "ATAC", meta.data=mdSKM1)

assaySKM2 <- CreateChromatinAssay(countsSKM2, fragments = fragsSKM2)
atacseuratSKM2 <- CreateSeuratObject(assaySKM2, assay = "ATAC", meta.data=mdSKM2)



#Now construct a seurat object for each sample with both RNA and ATAC data
#Keep only those cells that are common to RNA and ATAC data in each sample

#Build seurat object SFL1
data  <- Read10X_h5(filename = filePaths[1])
rnaseuratSFL1 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(1)])
commonCells <- intersect(colnames(rnaseuratSFL1),colnames(atacseuratSFL1)) 
rnaseuratSFL1 <- subset(rnaseuratSFL1, cells = commonCells)
atacseuratSFL1 <- subset(atacseuratSFL1, cells = commonCells)
seuratSFL1 <- rnaseuratSFL1
seuratSFL1[['ATAC']] <- CreateChromatinAssay(
  counts = atacseuratSFL1@assays$ATAC$counts,
  genome = "mm10",
  fragments = fragsSFL1
)





#Build seurat object SFL2
data  <- Read10X_h5(filename = filePaths[3])
rnaseuratSFL2 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(3)])
commonCells <- intersect(colnames(rnaseuratSFL2),colnames(atacseuratSFL2))
rnaseuratSFL2 <- subset(rnaseuratSFL2, cells = commonCells)
atacseuratSFL2 <- subset(atacseuratSFL2, cells = commonCells)
seuratSFL2 <- rnaseuratSFL2
seuratSFL2[['ATAC']] <- CreateChromatinAssay(
  counts = atacseuratSFL2@assays$ATAC$counts,
  genome = "mm10",
  fragments = fragsSFL2
)
SaveSeuratRds(seuratSFL2,"data/SFL2/outs/seuratSFL2.RDS")






#Build seurat object SKM1
data  <- Read10X_h5(filename = filePaths[5])
rnaseuratSKM1 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(5)])
commonCells <- intersect(colnames(rnaseuratSKM1),colnames(atacseuratSKM1))
rnaseuratSKM1 <- subset(rnaseuratSKM1, cells = commonCells)
atacseuratSKM1 <- subset(atacseuratSKM1, cells = commonCells)
seuratSKM1 <- rnaseuratSKM1
seuratSKM1[['ATAC']] <- CreateChromatinAssay(
  counts = atacseuratSKM1@assays$ATAC$counts,
  genome = "mm10",
  fragments = fragsSKM1
)
SaveSeuratRds(seuratSKM1,"data/SKM1/outs/seuratSKM1.RDS")





#Build seurat object SKM2
data  <- Read10X_h5(filename = filePaths[7])
rnaseuratSKM2 <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(7)])
commonCells <- intersect(colnames(rnaseuratSKM2),colnames(atacseuratSKM2))
rnaseuratSKM2 <- subset(rnaseuratSKM2, cells = commonCells)
atacseuratSKM2 <- subset(atacseuratSKM2, cells = commonCells)
seuratSKM2 <- rnaseuratSKM2
seuratSKM2[['ATAC']] <- CreateChromatinAssay(
  counts = atacseuratSKM2@assays$ATAC$counts,
  genome = "mm10",
  fragments = fragsSKM2
)
SaveSeuratRds(seuratSKM2,"data/SKM2/outs/seuratSKM2.RDS")


####################################################################Annotations#########################################################################


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



####################################################################Doublet Finder#########################################################################

#Call FindDoublets function to perform quality control and detect doublets
#Ensure to select the correct classifications column when assigning to cell type
seuratSFL1 <- FindDoublets(seuratSFL1)
seuratSFL1$CellType <- seuratSFL1$DF.classifications_0.25_0.28_557

seuratSFL2 <- FindDoublets(seuratSFL2)
seuratSFL2$CellType <- seuratSFL2$DF.classifications_0.25_0.1_415

seuratSKM1 <- FindDoublets(seuratSKM1)
seuratSKM1$CellType <- seuratSKM1$DF.classifications_0.25_0.27_418

seuratSKM2 <- FindDoublets(seuratSKM2)
seuratSKM2$CellType <- seuratSKM2$DF.classifications_0.25_0.21_473

#Save all Seurat objects for backup
SaveSeuratRds(seuratSFL1,"data/SFL1/outs/seuratDoubletsSFL1.RDS")
SaveSeuratRds(seuratSFL2,"data/SFL2/outs/seuratDoubletsSFL2.RDS")
SaveSeuratRds(seuratSKM1,"data/SKM1/outs/seuratDoubletsSKM1.RDS")
SaveSeuratRds(seuratSKM2,"data/SKM2/outs/seuratDoubletsSKM2.RDS")

#Visualize the doublets
SFL1Plot <- DimPlot(seuratSFL1,reduction = "umap",group.by = "CellType") + ggtitle("SFL1")
SFL2Plot <- DimPlot(seuratSFL2,reduction = "umap",group.by = "CellType") + ggtitle("SFL2")
SKM1Plot <- DimPlot(seuratSKM1,reduction = "umap",group.by = "CellType") + ggtitle("SKM1")
SKM2Plot <- DimPlot(seuratSKM2,reduction = "umap",group.by = "CellType") + ggtitle("SKM2")

SFL1Plot+SFL2Plot+SKM1Plot+SKM2Plot



####################################################################Merge Seurat Objects#########################################################################


# Rename cells for each sample
seuratSFL1 <- RenameCells(seuratSFL1, add.cell.id = "SFL1")
seuratSFL2 <- RenameCells(seuratSFL2, add.cell.id = "SFL2")
seuratSKM1 <- RenameCells(seuratSKM1, add.cell.id = "SKM1")
seuratSKM2 <- RenameCells(seuratSKM2, add.cell.id = "SKM2")

#Merge and save the Seurat Object
seurat <- merge(seuratSFL1,y=c(seuratSFL2,seuratSKM1,seuratSKM2), project = "MultiOmicsMacrophage")
SaveSeuratRds(seurat,"data/mergedseurat.RDS")


####################################################################Signac#########################################################################



#Using Signac package do quality control of ATAC Assays and then preprocess both RNA (Integration) and ATAC Assays to compute dimension reduction
#Refer DoSignac function in Functions.R for more info
seurat <- DoSignac(seurat = seurat)

#Assign sample type to the samples
seurat@meta.data <- seurat@meta.data %>%
  dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))


#Join layers in RNA assay for further processing
DefaultAssay(seurat) <- "RNA"
seurat <- JoinLayers(seurat)


#Save Seurat Object
SaveSeuratRds(seurat,"data/MergedSignacSeurat.RDS")


####################################################################Visualize Umaps#########################################################################

#RNA UMAP
p1 <- DimPlot(seurat, reduction = 'umap.rna', label = TRUE, group.by = "orig.ident") + ggtitle("RNA UMAP")
p2 <- DimPlot(seurat, reduction = 'umap.rna', label = TRUE, split.by = "orig.ident") + ggtitle("RNA UMAP")
p1+p2

#ATAC UMAP
p1 <- DimPlot(seurat, reduction = 'umap.atac', label = TRUE, group.by = "orig.ident") + ggtitle("ATAC UMAP")
p2 <- DimPlot(seurat, reduction = 'umap.atac', label = TRUE, split.by = "orig.ident") + ggtitle("ATAC UMAP")
p1+p2


####################################################################Link peaks to genes#########################################################################


DefaultAssay(seurat) <- "ATAC"

# first compute the GC content for each peak
seurat <- RegionStats(seurat, genome = BSgenome.Mmusculus.UCSC.mm10)


# link peaks to genes
seurat <- LinkPeaks(
  object = seurat,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)

#Save the seurat object
SaveSeuratRds(seurat, file = "data/GeneLinkedSeurat.RDS")



####################################################################Clustering resolution testing using Clustree#########################################################################

DefaultAssay(seurat) <- "RNA"

#Call VisualizeClustree function to cluster the cells at different resolutions passed to the function
seurat <- VisualizeClustree(0,1,0.1,seurat)


#Use the resolution that is determined to be the best after visualizing the cluster tree


####################################################################Testing Clustering Significance#########################################################################

#Extract counts data from RNA assay to test the clustering significance
counts <- seurat@assays$RNA$counts

#Perform the significance test
new_clusters = testClusters(counts, 
                            cluster_ids = as.character(seurat$cca_clusters_0_1_0.2), # Choose the correct cluster column as determined in Clustree step
                            alpha = 0.05, #FWER control, can be relaxed if needed
                            num_features = 2500,#default number
                            num_PCs = 30, #default number
                            parallel = FALSE, #can set to FALSE
                            cores = 1)

#Compare the new clusters with the existing clusters
table(new_clusters[[1]],seurat$cca_clusters_0_1_0.2)

#Add old clusters as Seurat and new clusters as Corrected Seurat to the Seurat Object
seurat$Seurat <- seurat$cca_clusters_0_1_0.2
seurat$Corrected_Seurat <- new_clusters[[1]]


#Visualize the old and new clusters
ggarrange(DimPlot(seurat,group.by='Seurat', label=T)+NoLegend(),
  DimPlot(seurat,group.by='Corrected_Seurat', label=T)+NoLegend(),
  nrow=1)


d1 = DimPlot(seurat,group.by='Seurat', label=T) #+
d2 = DimPlot(seurat,group.by='Corrected_Seurat', label=T) #+
d1 + d2


ggsave("compare_clusters_Seurat.png", width=10, height=8, units="in", dpi=588)
ggsave("compare_clusters_Seurat.pdf", width=10, height=8)


#Visualize UMAPs by imposing the clusters
rna <- DimPlot(seurat,group.by='Corrected_Seurat', label=T, split.by = "orig.ident", reduction = "umap.rna")
atac <- DimPlot(seurat,group.by='Corrected_Seurat', label=T, split.by = "orig.ident", reduction = "umap.atac")
rna + atac


####################################################################Find Cluster Markers#########################################################################

#Call WriteClusterMarkerGenes function to find markers for each cluster and save them in a text file
WriteClusterMarkerGenes(seurat, unique(seurat$Corrected_Seurat))

####################################################################Annotate Cell Types#########################################################################

#Call DoEnrichRAnnotation function to use the cluster markers and assign cell types to the clusters
#Use Mouse Cell Atlas and Tabula Muris # Call the function twice with one of these databases each time  #Has to be changed inside the function
seurat <- DoEnrichRAnnotation(seurat, unique(seurat$Corrected_Seurat))

####################################################################Visualize Cell Type Annotations#########################################################################

#Plot cell type annotations based on significance according to different metrics
Idents(seurat) <- "er_adj_pval_annotation"
adjPValAnn <- DimPlot(seurat, label=T, reduction = "umap.rna") + ggtitle("er_adj_pval_annotation")

Idents(seurat) <- "er_odds_ratio_annotation"
oddsratioAnn <- DimPlot(seurat, label=T, reduction = "umap.rna") + ggtitle("er_odds_ratio_annotation")

Idents(seurat) <- "er_combined_score_annotation"
combinedscoreAnn <- DimPlot(seurat, label=T, reduction = "umap.rna") + ggtitle("er_combined_score_annotation")


adjPValAnn/oddsratioAnn/combinedscoreAnn


#Save the plots
ggsave(filename = "data/Annotations.jpeg", plot = adjPValAnn/oddsratioAnn/combinedscoreAnn, width = 15, height = 12)
ggsave(filename = "data/OddsRatioAnnotations.jpeg", plot = oddsratioAnn)
ggsave(filename = "data/CombinedScoreAnnotations.jpeg", plot = combinedscoreAnn)

#Save the seurat object
SaveSeuratRds(seurat, file = "data/EnrichRAnnotatedJoinedSeurat.RDS")


###############################################################################DE & DA Analysis ##############################################################

#Do differential gene expression analysis between KO and WT samples
RNADGE <- FindMarkers(
  
  object = seurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "RNA"
  
)

RNADGE$Gene <- rownames(RNADGE)



#Do differential accessability analysis between KO and WT samples
ATACDA <- FindMarkers(
  
  object = seurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "ATAC"
  
)

ATACDA$Region <- rownames(ATACDA)



#Map the differentially accessable regions to the nearest gene
DefaultAssay(seurat) <- "ATAC"
regionGeneMap <- ClosestFeature(object = seurat, regions = rownames(ATACDA), annotation = NULL)

regionGeneMap$Gene_PValue <- NA
regionGeneMap$Gene_LFC <- NA
regionGeneMap$Gene_AdjPVal <- NA
regionGeneMap$Region_PValue <- NA
regionGeneMap$Region_LFC <- NA
regionGeneMap$Region_AdjPVal <- NA


#Add the gene differential expression metrics to region gene map data frame if they exists in RNADGE results
for (gene in unique(regionGeneMap$gene_name)) {
  
  if(length(RNADGE$p_val[RNADGE$Gene == gene]) != 0){
    
    regionGeneMap$Gene_PValue[regionGeneMap$gene_name == gene] <- RNADGE$p_val[RNADGE$Gene == gene]
    regionGeneMap$Gene_LFC[regionGeneMap$gene_name == gene] <- RNADGE$avg_log2FC[RNADGE$Gene == gene]
    regionGeneMap$Gene_AdjPVal[regionGeneMap$gene_name == gene] <- RNADGE$p_val_adj[RNADGE$Gene == gene]
    
  }
  
  
}



#Add the region differential accessability metrics to the region gene map data frame
for (region in unique(regionGeneMap$query_region)) {
  
  if(length(ATACDA$p_val[ATACDA$Region == region]) != 0){
    
    regionGeneMap$Region_PValue[regionGeneMap$query_region == region] <- ATACDA$p_val[ATACDA$Region == region]
    regionGeneMap$Region_LFC[regionGeneMap$query_region == region] <- ATACDA$avg_log2FC[ATACDA$Region == region]
    regionGeneMap$Region_AdjPVal[regionGeneMap$query_region == region] <- ATACDA$p_val_adj[ATACDA$Region == region]
    
  }
  
  
}


#Save the region gene map data
write.csv(regionGeneMap, file = "data/RegionGeneMap.csv", row.names = TRUE)



#Filter the rows if there are no differential expression or accessability metrics
regionGeneMap <- regionGeneMap[!is.na(regionGeneMap$Gene_PValue),]
regionGeneMap <- regionGeneMap[!is.na(regionGeneMap$Region_PValue),]


#Filter the rows for significant region and gene maps
regionGeneMap <- regionGeneMap[
  regionGeneMap$Gene_AdjPVal < 0.05 &
    regionGeneMap$Region_AdjPVal < 0.05
  # abs(regionGeneMap$Gene_LFC) > 1 &
  # abs(regionGeneMap$Region_LFC) > 1
  ,]


#Filter out the region gene maps if the direction of change is different for regions and maps
regionGeneMap <- regionGeneMap[
  (regionGeneMap$Gene_LFC < 0 &
     regionGeneMap$Region_LFC < 0) |
    (regionGeneMap$Gene_LFC > 0 &
       regionGeneMap$Region_LFC > 0),]



#Filter for distance less than 2Kb
regionGeneMap <- regionGeneMap[
  (regionGeneMap$distance < 2000 &
     regionGeneMap$distance > -2000),]




#Filter for type exon, UTR and cds
regionGeneMap <- regionGeneMap[
  (regionGeneMap$type == "exon" |
     regionGeneMap$type == "cds" |
     regionGeneMap$type == "utr"),]


#Save the filtered list of genes to use for pathway enrichment analysis
write(paste0(unique(regionGeneMap$gene_name),"\n"), file = paste0("data/GenesOfInterest.txt"))


#Generate Coverage plots for the significant genes as these are showing changes in both RNA and ATAC data in the same direction
GenerateCoveragePlots(EnrichRAnnotatedJoinedSeurat, unique(regionGeneMap$gene_name))


#Do pathway enrichment using the significant genes
#No code for doing this. Do manually using Enrich R website and download the results table to the correct location
Pathwayresults <- fread(file = "data/WikiPathways_2024_Mouse_table2.txt")

#Filter for significant pathways
Pathwayresults <- Pathwayresults[Pathwayresults$`P-value` < 0.01,]

#Visualize the bar plot of the pathways
ggplot(Pathwayresults, aes(x = reorder(`Term`, `Combined Score`), y = `Combined Score`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Pathway", y = "Combined Score", title = "Pathway Enrichment") +
  theme_minimal()


#Visualize the dot plot of the pathways
ggplot(Pathwayresults, aes(x=`Overlap`, y=reorder(`Term`, -log10(`Adjusted P-value`)), 
                           size=-log10(`Adjusted P-value`), color=`Odds Ratio`)) +
  geom_point() +
  scale_color_gradient(low="blue", high="red") +
  labs(x="Gene Overlap", y="Pathway", title="Pathway Enrichment Dot Plot") +
  theme_minimal()


####################################################################THE END#########################################################################