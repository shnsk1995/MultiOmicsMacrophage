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
SaveSeuratRds(seuratSFL1,"data/SFL1/outs/seuratSFL1.RDS")


 


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





seuratSFL1
seuratSFL2
seuratSKM1
seuratSKM2




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

seuratSFL1$CellType <- seuratSFL1$DF.classifications_0.25_0.28_557

seuratSFL2 <- FindDoublets(seuratSFL2)

seuratSFL2$CellType <- seuratSFL2$DF.classifications_0.25_0.1_415

seuratSKM1 <- FindDoublets(seuratSKM1)

seuratSKM1$CellType <- seuratSKM1$DF.classifications_0.25_0.27_418

seuratSKM2 <- FindDoublets(seuratSKM2)

seuratSKM2$CellType <- seuratSKM2$DF.classifications_0.25_0.21_473




SaveSeuratRds(seuratSFL1,"data/SFL1/outs/seuratDoubletsSFL1.RDS")
SaveSeuratRds(seuratSFL2,"data/SFL2/outs/seuratDoubletsSFL2.RDS")
SaveSeuratRds(seuratSKM1,"data/SKM1/outs/seuratDoubletsSKM1.RDS")
SaveSeuratRds(seuratSKM2,"data/SKM2/outs/seuratDoubletsSKM2.RDS")


SFL1Plot <- DimPlot(seuratSFL1,reduction = "umap",group.by = "CellType") + ggtitle("SFL1")
SFL2Plot <- DimPlot(seuratSFL2,reduction = "umap",group.by = "CellType") + ggtitle("SFL2")
SKM1Plot <- DimPlot(seuratSKM1,reduction = "umap",group.by = "CellType") + ggtitle("SKM1")
SKM2Plot <- DimPlot(seuratSKM2,reduction = "umap",group.by = "CellType") + ggtitle("SKM2")

SFL1Plot+SFL2Plot+SKM1Plot+SKM2Plot




# Rename cells for each sample
seuratSFL1 <- RenameCells(seuratSFL1, add.cell.id = "SFL1")
seuratSFL2 <- RenameCells(seuratSFL2, add.cell.id = "SFL2")
seuratSKM1 <- RenameCells(seuratSKM1, add.cell.id = "SKM1")
seuratSKM2 <- RenameCells(seuratSKM2, add.cell.id = "SKM2")






seurat <- merge(seuratSFL1,y=c(seuratSFL2,seuratSKM1,seuratSKM2), project = "MultiOmicsMacrophage")
SaveSeuratRds(seurat,"data/mergedseurat.RDS")






seurat <- DoSignac(seurat = seurat)
SaveSeuratRds(seurat,"data/mergedsignacseurat.RDS")

seurattest <- seurat




seurattest@meta.data <- seurattest@meta.data %>%
  
  dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))



p1 <- DimPlot(seurattest, reduction = 'umap.rna', label = TRUE, group.by = "orig.ident") + ggtitle("RNA UMAP")

p2 <- DimPlot(seurattest, reduction = 'umap.rna', label = TRUE, split.by = "orig.ident") + ggtitle("RNA UMAP")



p1+p2


p1 <- DimPlot(seurattest, reduction = 'umap.atac', label = TRUE, group.by = "orig.ident") + ggtitle("ATAC UMAP")

p2 <- DimPlot(seurattest, reduction = 'umap.atac', label = TRUE, split.by = "orig.ident") + ggtitle("ATAC UMAP")



p1+p2



DefaultAssay(seurattest) <- "RNA"

seuratJoined <- JoinLayers(seurattest)





#Linking peaks to genes
DefaultAssay(seuratJoined) <- "ATAC"

# first compute the GC content for each peak
seuratJoined <- RegionStats(seuratJoined, genome = BSgenome.Mmusculus.UCSC.mm10)



# link peaks to genes
seuratJoined <- LinkPeaks(
  object = seuratJoined,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)







SaveSeuratRds(seuratJoined, file = "data/GeneLinkedSeurat.RDS")


DefaultAssay(seuratJoined) <- "ATAC"

rownames(seuratJoined)

seuratJoined <- readRDS("data/GeneLinkedSeurat.RDS")


#Clustering resolution significance testing

DefaultAssay(seuratJoined) <- "RNA"
seuratJoined <- VisualizeClustree(0,1,0.1,seuratJoined)


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

rna <- DimPlot(seuratJoined,group.by='Corrected_Seurat', label=T, split.by = "orig.ident", reduction = "umap.rna")
atac <- DimPlot(seuratJoined,group.by='Corrected_Seurat', label=T, split.by = "orig.ident", reduction = "umap.atac")

rna + atac







# SFL1New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL1" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
# SFL2New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL2" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
# SKM1New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM1" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
# SKM2New1Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM2" & seuratJoined@meta.data$Corrected_Seurat == "new1",])
# 
# 
# SFL1New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL1" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
# SFL2New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SFL2" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
# SKM1New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM1" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
# SKM2New2Count <- nrow(seuratJoined@meta.data[seuratJoined@meta.data$orig.ident == "SKM2" & seuratJoined@meta.data$Corrected_Seurat == "new2",])
# 
# 
# cellCountDf <- data.frame(
#   SFL1 = c(SFL1New1Count,SFL1New2Count),
#   SFL2 = c(SFL2New1Count,SFL2New2Count),
#   SKM1 = c(SKM1New1Count,SKM1New2Count),
#   SKM2 = c(SKM2New1Count,SKM2New2Count)
# )
# 
# rownames(cellCountDf) <- c("new1","new2")
# 
# normalizedCellCountDf <- sweep(cellCountDf, 2, colSums(cellCountDf), FUN = "/")
# 
# normalizedCellCountDf$ClusterType <- rownames(normalizedCellCountDf)
# normalizedCellCountDf[] <- lapply(normalizedCellCountDf, as.character)
# normalizedCellCountDfLong <- pivot_longer(normalizedCellCountDf, cols = -"ClusterType")
# colnames(normalizedCellCountDfLong) <- c("ClusterType","SampleType","CellCountProportion")
# 
# Plot1 <- ggplot(normalizedCellCountDfLong[grep("1",normalizedCellCountDfLong$SampleType),], aes(x = SampleType, y = CellCountProportion, fill = factor(ClusterType))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Cluster Type", y = "Cell Count Proportion", fill = "Sample Type") +
#   theme_minimal()
# 
# Plot2 <- ggplot(normalizedCellCountDfLong[grep("2",normalizedCellCountDfLong$SampleType),], aes(x = SampleType, y = CellCountProportion, fill = factor(ClusterType))) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Cluster Type", y = "Cell Count Proportion", fill = "Sample Type") +
#   theme_minimal()
# 
# Plot1/Plot2
#####################Rough work from here###############################

# Idents(seuratJoined) <- "Corrected_Seurat"
# 
# Idents(seuratJoined) <- "SampleType"
# 
# Idents(SKMseuratJoined) <- "Corrected_Seurat"
# 
# Idents(seuratJoined) <- "Corrected_Seurat"
# New1seuratJoined <- subset(seuratJoined, idents = "new1")
# Idents(New1seuratJoined) <- "Corrected_Seurat"
# 
# allRNAMarkers <- FindAllMarkers(New1seuratJoined, assay = "RNA")
# allATACMarkers <- FindAllMarkers(seuratJoined, assay = "ATAC")
# colnas <- colnames(allATACMarkers)
# colnas[7] <- "peak"
# colnames(allATACMarkers) <- colnas
# 
# links <- Links(seuratJoined[["ATAC"]])
# 
# # Convert links to a data frame for merging
# linksDF <- as.data.frame(links)
# 
# # Merge DA results with gene links
# allATACMarkers <- allATACMarkers %>%
#   left_join(linksDF, by = "peak", relationship = "many-to-many")
# 
# unique(seuratJoined$Corrected_Seurat)
# 
# WriteClusterSignificantGenes(seuratJoined,unique(seuratJoined$Corrected_Seurat))

WriteClusterMarkerGenes(seuratJoined, unique(seuratJoined$Corrected_Seurat))

EnrichRAnnotatedJoinedSeurat <- DoEnrichRAnnotation(seuratJoined, unique(seuratJoined$Corrected_Seurat))

Idents(EnrichRAnnotatedJoinedSeurat) <- "er_adj_pval_annotation"
adjPValAnn <- DimPlot(EnrichRAnnotatedJoinedSeurat, label=T, reduction = "umap.rna") + ggtitle("er_adj_pval_annotation")

Idents(EnrichRAnnotatedJoinedSeurat) <- "er_odds_ratio_annotation"
oddsratioAnn <- DimPlot(EnrichRAnnotatedJoinedSeurat, label=T, reduction = "umap.rna") + ggtitle("er_odds_ratio_annotation")


Idents(EnrichRAnnotatedJoinedSeurat) <- "er_combined_score_annotation"
combinedscoreAnn <- DimPlot(EnrichRAnnotatedJoinedSeurat, label=T, reduction = "umap.rna") + ggtitle("er_combined_score_annotation")


adjPValAnn/oddsratioAnn/combinedscoreAnn



ggsave(filename = "data/Annotations.jpeg", plot = adjPValAnn/oddsratioAnn/combinedscoreAnn, width = 15, height = 12)
ggsave(filename = "data/OddsRatioAnnotations.jpeg", plot = oddsratioAnn)
ggsave(filename = "data/CombinedScoreAnnotations.jpeg", plot = combinedscoreAnn)

SaveSeuratRds(EnrichRAnnotatedJoinedSeurat, file = "data/EnrichRAnnotatedJoinedSeurat.RDS")


###############################################################################DE & DA Analysis ##############################################################

#Do differential gene expression analysis between KO and WT samples


RNADGE <- FindMarkers(
  
  object = EnrichRAnnotatedJoinedSeurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "RNA"
  
)

RNADGE$Gene <- rownames(RNADGE)


#Do differential accessability analysis between KO and WT samples


ATACDA <- FindMarkers(
  
  object = EnrichRAnnotatedJoinedSeurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "ATAC"
  
)

ATACDA$Region <- rownames(ATACDA)



#Map the differentially accessable regions to the nearest gene

DefaultAssay(EnrichRAnnotatedJoinedSeurat) <- "ATAC"
regionGeneMap <- ClosestFeature(object = EnrichRAnnotatedJoinedSeurat, regions = rownames(ATACDA), annotation = NULL)

regionGeneMap$Gene_PValue <- NA
regionGeneMap$Gene_LFC <- NA
regionGeneMap$Gene_AdjPVal <- NA
regionGeneMap$Region_PValue <- NA
regionGeneMap$Region_LFC <- NA
regionGeneMap$Region_AdjPVal <- NA



for (gene in unique(regionGeneMap$gene_name)) {
  
  if(length(RNADGE$p_val[RNADGE$Gene == gene]) != 0){
    
    regionGeneMap$Gene_PValue[regionGeneMap$gene_name == gene] <- RNADGE$p_val[RNADGE$Gene == gene]
    regionGeneMap$Gene_LFC[regionGeneMap$gene_name == gene] <- RNADGE$avg_log2FC[RNADGE$Gene == gene]
    regionGeneMap$Gene_AdjPVal[regionGeneMap$gene_name == gene] <- RNADGE$p_val_adj[RNADGE$Gene == gene]
    
  }
  
  
}






for (region in unique(regionGeneMap$query_region)) {
  
  if(length(ATACDA$p_val[ATACDA$Region == region]) != 0){
    
    regionGeneMap$Region_PValue[regionGeneMap$query_region == region] <- ATACDA$p_val[ATACDA$Region == region]
    regionGeneMap$Region_LFC[regionGeneMap$query_region == region] <- ATACDA$avg_log2FC[ATACDA$Region == region]
    regionGeneMap$Region_AdjPVal[regionGeneMap$query_region == region] <- ATACDA$p_val_adj[ATACDA$Region == region]
    
  }
  
  
}















write.csv(regionGeneMap, file = "data/RegionGeneMap.csv", row.names = TRUE)

regionGeneMap <- read.csv(file = "data/RegionGeneMap.csv", row.names = 1)

regionGeneMap <- regionGeneMap[!is.na(regionGeneMap$Gene_PValue),]
regionGeneMap <- regionGeneMap[!is.na(regionGeneMap$Region_PValue),]

regionGeneMap <- regionGeneMap[
  regionGeneMap$Gene_AdjPVal < 0.05 &
    regionGeneMap$Region_AdjPVal < 0.05
    # abs(regionGeneMap$Gene_LFC) > 1 &
    # abs(regionGeneMap$Region_LFC) > 1
    ,]


regionGeneMap <- regionGeneMap[
  (regionGeneMap$Gene_LFC < 0 &
     regionGeneMap$Region_LFC < 0) |
    (regionGeneMap$Gene_LFC > 0 &
       regionGeneMap$Region_LFC > 0),]



#Filter for distance less than 2Kb
regionGeneMap <- regionGeneMap[
  (regionGeneMap$distance < 2000 &
     regionGeneMap$distance > -2000),]



unique(regionGeneMap$type)

#Filter for type exon and cds
regionGeneMap <- regionGeneMap[
  (regionGeneMap$type == "exon" |
     regionGeneMap$type == "cds" |
     regionGeneMap$type == "utr"),]



write(paste0(unique(regionGeneMap$gene_name),"\n"), file = paste0("data/GenesOfInterest.txt"))


GenerateCoveragePlots(EnrichRAnnotatedJoinedSeurat, unique(regionGeneMap$gene_name))


Pathwayresults <- fread(file = "data/WikiPathways_2024_Mouse_table2.txt")

Pathwayresults <- Pathwayresults[Pathwayresults$`P-value` < 0.01,]



enrichplot::dotplot(Pathwayresults)



ggplot(Pathwayresults, aes(x = reorder(`Term`, `Combined Score`), y = `Combined Score`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Pathway", y = "Combined Score", title = "Pathway Enrichment") +
  theme_minimal()





ggplot(Pathwayresults, aes(x=`Overlap`, y=reorder(`Term`, -log10(`Adjusted P-value`)), 
                   size=-log10(`Adjusted P-value`), color=`Odds Ratio`)) +
  geom_point() +
  scale_color_gradient(low="blue", high="red") +
  labs(x="Gene Overlap", y="Pathway", title="Pathway Enrichment Dot Plot") +
  theme_minimal()




DefaultAssay(EnrichRAnnotatedJoinedSeurat)
atactest <- ClosestFeature(object = EnrichRAnnotatedJoinedSeurat, regions = rownames(ATACDA), annotation = NULL)

interestingGenes <- rownames(RNADGE)[rownames(RNADGE) %in% atactest$gene_name]
GenerateCoveragePlots(EnrichRAnnotatedJoinedSeurat, interestingGenes)

g <- "Filip1l"
atactest$query_region[atactest$gene_name == g]
ATACDA$Region <- rownames(ATACDA)
ATACDA$avg_log2FC[ATACDA$Region== atactest$query_region[atactest$gene_name == g]]
RNADGE$Genes <- rownames(RNADGE)

gLFCpLFC <- data.frame(matrix(ncol=3))
colnames(gLFCpLFC) <- c("Gene","DE_LFC","DA_LFC")

for (gene in interestingGenes){
  
  row <- data.frame(matrix(ncol=3, nrow = 1))
  colnames(row) <- c("Gene","DE_LFC","DA_LFC")
  
  row[1,1] <- gene
  row[1,2] <- RNADGE$avg_log2FC[RNADGE$Genes==gene]
  if(length(ATACDA$avg_log2FC[ATACDA$Region== atactest$query_region[atactest$gene_name == gene]])==0) next
  if(length(ATACDA$avg_log2FC[ATACDA$Region== atactest$query_region[atactest$gene_name == gene]]) > 1){
    row[1,3] <- ATACDA$avg_log2FC[ATACDA$Region== atactest$query_region[atactest$gene_name == gene]][1]
  }else{
    row[1,3] <- ATACDA$avg_log2FC[ATACDA$Region== atactest$query_region[atactest$gene_name == gene]]
  }
  
  
  gLFCpLFC <- rbind(gLFCpLFC,row)
  
}

gLFCpLFC <- na.omit(gLFCpLFC)

ggplot(gLFCpLFC, aes(x = DE_LFC, y = DA_LFC))+
  geom_smooth(method = "lm")


write.csv(RNADGE, file = "data/RNADGE.csv",row.names = TRUE, col.names = TRUE)
write.csv(ATACDA, file = "data/ATACDA.csv",row.names = TRUE, col.names = TRUE)
write.csv(atactest, file = "data/MappedATACRegions.csv",row.names = TRUE, col.names = TRUE)

cor(gLFCpLFC$DE_LFC, gLFCpLFC$DA_LFC, method = "pearson")

Idents(seuratJoined) <- "Corrected_Seurat"
New2seuratJoined <- subset(seuratJoined, idents = "new2")
Idents(New2seuratJoined) <- "Corrected_Seurat"

RNADGE <- FindMarkers(
  
  object = EnrichRAnnotatedJoinedSeurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "RNA"
  
)

sigRNADGE <- RNADGE[RNADGE$p_val_adj < 0.05,]

write(paste0(rownames(sigRNADGE),"\n"), file = "data/New2DGESIGGENES.txt")

ATACDA <- FindMarkers(
  
  object = EnrichRAnnotatedJoinedSeurat,
  
  ident.1 = "KnockOut",
  
  ident.2 = "WildType",
  
  group.by = "SampleType",
  
  test.use = "wilcox",
  
  assay = "ATAC"
  
)

RNADGE <- RNADGE[RNADGE$p_val_adj < 0.05,]
RNADGE$gene <- rownames(RNADGE)

commonSignificantDGEDA <- FindCommonDGEDA(RNADGE, ATACDA)


GenerateCoveragePlots(seuratJoined, RNADGE, 91)


clusterMarkersNew1 <- allRNAMarkers$gene[allRNAMarkers$cluster == "new1" & allRNAMarkers$p_val_adj < 0.05 & abs(allRNAMarkers$avg_log2FC) > 1]
genes <- paste0(clusterNew1Markers, collapse = ",")

# geneName <- "Erich6"
# 
# geneRegions <- GetGenomeRanges(seuratJoined,geneName)
# 
# GenerateCoveragePlots(seuratJoined, geneRegions,geneName)
# 
# Idents(seuratJoined) <- "SampleType"
# ranges.show <- StringToGRanges(geneRegions[15:20])
# CoveragePlot(seuratJoined,assay = "ATAC",
#              region = geneRegions[15:20],
#              features = "Cmss1", 
#              expression.assay = "RNA",
#              region.highlight = ranges.show, 
#              split.by = "Corrected_Seurat")
# 
# clusterNew1seurat@assays$ATAC@annotation[clusterNew1seurat@assays$ATAC@annotation$gene_name == "Siah2",]
# 
# Idents(clusterNew1seurat) <- "SampleType"
# 
# clusterNew1seuratWT <- subset(clusterNew1seurat, idents = "WildType")
# 
# clusterNew1seuratKO <- subset(clusterNew1seurat, idents = "KnockOut")
# 
# 
# 
# 
# 
# DefaultAssay(clusterNew1seuratWT) <- "ATAC"
# 
# 
# 
# p1 <- CoveragePlot(clusterNew1seuratWT,assay = "ATAC",region = c("chr3-58688799-58689631","chr3-58690316-58691181","chr3-58691816-58692742"),features = "Siah2", expression.assay = "RNA")
# 
# 
# 
# DefaultAssay(clusterNew1seuratKO) <- "ATAC"
# 
# 
# 
# p2 <- CoveragePlot(clusterNew1seuratKO,assay = "ATAC",region = c("chr3-58688799-58689631","chr3-58690316-58691181","chr3-58691816-58692742"),features = "Siah2", expression.assay = "RNA")
# 
# 
# 
# p1/p2
# 
# 
# 
# 
# 
# DefaultAssay(cluster0seurat) <- "ATAC"
# 
# 
# 
# CoveragePlot(cluster0seurat,region = "chr13-115088530-115088961",features = "Pelo")
# 
# 
# 
# 
# 
# 
# 
# rownames(clusterNew1seuratWT)[24000:25000]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# seuratSFL1 <- DoSignac(seuratSFL1,annotations)
# 
# seuratSFL2 <- DoSignac(seuratSFL2,annotations)
# 
# seuratSKM1 <- DoSignac(seuratSKM1,annotations)
# 
# seuratSKM2 <- DoSignac(seuratSKM2,annotations)
# 
# 
# 
# p1 <- DimPlot(seuratSFL1, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SFL1")
# 
# p2 <- DimPlot(seuratSFL1, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SFL1")
# 
# 
# 
# SFL1Plots <- p1 + p2
# 
# 
# 
# SFL1Plots
# 
# 
# 
# 
# 
# p1 <- DimPlot(seuratSFL2, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SFL2")
# 
# p2 <- DimPlot(seuratSFL2, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SFL2")
# 
# 
# 
# SFL2Plots <- p1 + p2
# 
# 
# 
# SFL2Plots
# 
# 
# 
# 
# 
# 
# 
# p1 <- DimPlot(seuratSKM1, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SKM1")
# 
# p2 <- DimPlot(seuratSKM1, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SKM1")
# 
# 
# 
# SKM1Plots <- p1 + p2
# 
# 
# 
# SKM1Plots
# 
# 
# 
# 
# 
# p1 <- DimPlot(seuratSKM2, reduction = 'umap.rna', label = TRUE) + NoLegend() + ggtitle("RNA UMAP SKM2")
# 
# p2 <- DimPlot(seuratSKM2, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP SKM2")
# 
# 
# 
# SKM2Plots <- p1 + p2
# 
# 
# 
# SKM2Plots
# 
# 
# 
# 
# 
# ggsave("data/SFL1.jpeg",SFL1Plots)
# 
# ggsave("data/SFL2.jpeg",SFL2Plots)
# 
# ggsave("data/SKM1.jpeg",SKM1Plots)
# 
# ggsave("data/SKM2.jpeg",SKM2Plots)
# 
# 
# 
# seuratSFL1@meta.data
# 
# 
# 
# # add the gene information to the object
# 
# Annotation(seurat[['ATAC']]) <- annotations
# 
# 
# 
# ATACregions <- rownames(ATACSeurat)
# 
# colnames(ATACSeurat)
# 
# 
# 
# rownames(RNASeurat)
# 
# colnames(RNASeurat)
# 
# 
# 
# bethSeurat <- CreateProjectSeurat(filePaths = filePaths,fileIndexes = fileIndexes,fileNames = fileNames)
# 
# 
# 
# bethSeurat@meta.data <- bethSeurat@meta.data %>%
#   
#   dplyr::mutate(SampleType = ifelse(orig.ident %in% c("SFL1","SFL2"),"WildType","KnockOut"))
# 
# 
# 
# bethSeurat@meta.data
# 
# 
# 
# DefaultAssay(bethSeurat) <- 'ATAC'
# 
# 
# 
# bethSeurat <- FindTopFeatures(bethSeurat, min.cutoff = 10)
# 
# bethSeurat <- RunTFIDF(bethSeurat)
# 
# bethSeurat <- RunSVD(bethSeurat)
# 
# bethSeurat <- RunUMAP(bethSeurat, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')
# 
# p2 <- DimPlot(bethSeurat, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")
# 
# 
# 
# bethSeurat <- DoStandardSingleCellAnalysis(bethSeurat)
# 
# 
# 
# 
# 
# DimPlot(bethSeurat, reduction = "uamp.unintegrated",split.by = "SampleType", label = TRUE)
# 
# 
# 
# 
# 
# bethSeurat <- IntegrateSeuratObjects(bethSeurat)
# 
# 
# 
# 
# 
# bethSeurat <- PerformAfterIntegrationStandardAnalysis(bethSeurat)
# 
# 
# 
# 
# 
# plot1 <- DimPlot(bethSeurat, reduction = "uamp.unintegrated",split.by = "SampleType",label = TRUE)
# 
# plot2 <- DimPlot(bethSeurat, reduction = "umap.cca", label = TRUE, repel = TRUE, split.by = "SampleType")
# 
# plot1+plot2
# 
# 
# 
# plot2
# 
# 
# 
# rownames(seuratSFL1)
# 
# 
# 
# CoveragePlot(seuratSFL1, region = "chr1-37504883-37505816",features = "Mgat4a")
# 
# 
# 
# 
# 
# head(seurat@meta.data)
# 
# head(seuratSFL1@meta.data)
# 
# 
# 
# rownames(seurat)
# 
# rownames(seuratSFL1)