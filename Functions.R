CreateProjectSeurat <- function(filePaths,fileIndexes,fileNames){
  
  
  
  seuratObjects <- list()
  
  
  
  for (fileIndex in fileIndexes) {
    
    
    
    data  <- Read10X_h5(filename = filePaths[fileIndex])
    
    
    
    #RNASeurat <- CreateSeuratObject(counts = data[["Gene Expression"]], assay = "RNA", project = fileNames[as.character(fileIndex)], min.cells = 3, min.features = 200)
    
    ATACSeurat <- CreateSeuratObject(counts = data[["Peaks"]], assay = "ATAC", project = fileNames[as.character(fileIndex)], min.cells = 3, min.features = 200)
    
    
    
    #seurat <- merge(RNASeurat, y = ATACSeurat, add.cell.ids = c(paste0(fileNames[as.character(fileIndex)],"_RNA"),paste0(fileNames[as.character(fileIndex)],"_ATAC")), project = fileNames[as.character(fileIndex)])
    
    seurat <- ATACSeurat
    
    
    
    if(fileIndex == "1"){
      
      bethSeurat <- seurat
      
    }else{
      
      bethSeurat <- merge(bethSeurat,y = seurat)
      
    }
    
    
    
  }
  
  
  
  return(bethSeurat)
  
  
  
}







DoStandardSingleCellAnalysis <- function(seuratObject){
  
  
  
  
  
  
  
  
  
  seuratObject <- PercentageFeatureSet(seuratObject, pattern = "^mt-",col.name = "percent.mt", assay = "RNA")
  
  
  
  
  
  
  
  
  
  
  
  VlnPlot(seuratObject,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
  
  
  
  
  
  
  
  
  
  plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  plot2 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  plot1 + plot2
  
  
  
  
  
  seuratObject <- subset(seuratObject,subset = nFeature_RNA > 200 &
                           
                           nFeature_RNA < 4000 &
                           
                           nCount_RNA < 20000 &
                           
                           percent.mt < 10)
  
  
  
  
  
  VlnPlot(seuratObject,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
  
  
  
  
  
  
  
  
  
  
  
  seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  
  
  
  
  
  
  
  seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
  
  
  
  top10 <- head(VariableFeatures(seuratObject),10)
  
  
  
  plot1 <- VariableFeaturePlot(seuratObject)
  
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  plot1 + plot2
  
  
  
  
  
  seuratObject <- ScaleData(seuratObject)
  
  
  
  
  
  
  
  
  
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  
  
  
  VizDimLoadings(seuratObject, dims = 1:2, reduction = "pca")
  
  
  
  
  
  DimPlot(seuratObject, reduction = "pca") + NoLegend()
  
  
  
  
  
  DimHeatmap(seuratObject, dims = 1, cells = 500, balanced = TRUE)
  
  
  
  
  
  ElbowPlot(seuratObject)
  
  
  
  
  
  
  
  seuratObject <- FindNeighbors(seuratObject,dims = 1:20,reduction = "pca")
  
  seuratObject <- FindClusters(seuratObject, resolution = 1.0, cluster.name = "unintegrated_clusters")
  
  
  
  
  
  seuratObject <- RunUMAP(seuratObject, dims = 1:20, reduction = "pca", reduction.name = "uamp.unintegrated")
  
  
  
  
  
  return(seuratObject)
  
  
  
  
  
  
  
  
  
}







IntegrateSeuratObject <- function(seuratObject){
  
  
  
  integratedSeuratObject <- IntegrateLayers(seuratObject,
                                            
                                            method = CCAIntegration,
                                            
                                            orig.reduction = "pca",
                                            
                                            assay = "RNA",
                                            
                                            new.reduction = "integrated.cca")
  
  
  
  return(integratedSeuratObject)
  
  
  
}









PerformAfterIntegrationStandardAnalysis <- function(seuratObject){
  
  
  
  
  
  
  
  
  
  
  
  
  
  # seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  
  
  
  # VizDimLoadings(seuratObject, dims = 1:2, reduction = "pca")
  
  #
  
  #
  
  # DimPlot(seuratObject, reduction = "pca") + NoLegend()
  
  #
  
  #
  
  # DimHeatmap(seuratObject, dims = 1, cells = 500, balanced = TRUE)
  
  #
  
  #
  
  # ElbowPlot(seuratObject)
  
  
  
  
  
  
  
  seuratObject <- FindNeighbors(seuratObject,dims = 1:20,reduction="integrated.cca")
  
  seuratObject <- FindClusters(seuratObject, resolution = 1.0,cluster.name = "cca_clusters")
  
  
  
  
  
  seuratObject <- RunUMAP(seuratObject, dims = 1:20, reduction = "integrated.cca", reduction.name = "umap.cca")
  
  
  
  
  
  return(seuratObject)
  
  
  
}





AddAnnotations <- function(seurat,annotations){
  
  
  
  
  
  # add the gene information to the object
  
  Annotation(seurat[["ATAC"]]) <- annotations
  
  
  
  return(seurat)
  
  
  
}





DoSignac <- function(seurat){
  
  
  
  
  
  
  
  
  
  
  
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
      
      CellType == "Singlet" &
      
      nCount_ATAC > 500
    
  )
  
  seurat
  
  
  
  
  
  
  
  
  
  
  
  DefaultAssay(seurat) <- "RNA"
  
  
  
  seurat <- FindVariableFeatures(seurat, nfeatures = 2000)
  
  seurat <- NormalizeData(seurat)
  
  seurat <- ScaleData(seurat)
  
  seurat <- RunPCA(seurat, npcs = 20, reduction.name = "pca.2")
  
  seurat <- IntegrateLayers(seurat,
                            
                            method = CCAIntegration,
                            
                            orig.reduction = "pca.2",
                            
                            assay = "RNA",
                            
                            new.reduction = "integrated.cca")
  
  seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "umap.rna", reduction = "integrated.cca")
  
  seurat <- FindNeighbors(seurat, dims = 1:20, reduction = "integrated.cca")
  
  seurat <- FindClusters(seurat, resolution = 0.5, algorithm = 3,cluster.name = "cca_clusters")
  
  
  
  
  
  
  
  
  
  
  
  DefaultAssay(seurat) <- 'ATAC'
  
  
  
  seurat <- FindTopFeatures(seurat, min.cutoff = 10)
  
  seurat <- RunTFIDF(seurat)
  
  seurat <- RunSVD(seurat)
  
  seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')
  
  
  
  
  
  return(seurat)
  
  
  
  
  
}





FindDoublets <- function(seurat){
  
  
  
  #Set default Assay as RNA
  
  DefaultAssay(seurat) <- "RNA"
  
  
  
  #QC
  
  seurat$mitoPercent <- PercentageFeatureSet(seurat,pattern = "^mt-")
  
  qcPlot <- VlnPlot(seurat,features = c("nCount_RNA","nFeature_RNA","mitoPercent"),ncol = 3)
  
  plot(qcPlot)
  
  cat("Select the lower threshold of nCount_RNA: \n")
  
  lowerThresholdRNA <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of nCount_RNA: ",as.numeric(lowerThresholdRNA)))
  
  cat("Select the upper threshold of nCount_RNA: \n")
  
  upperThresholdRNA <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of nCount_RNA: ",as.numeric(upperThresholdRNA)))
  
  cat("Select the lower threshold of nFeature_RNA: \n")
  
  lowerThresholdFeatures <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of nFeature_RNA: ",as.numeric(lowerThresholdFeatures)))
  
  cat("Select the upper threshold of nFeature_RNA: \n")
  
  upperThresholdFeatures <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of nFeature_RNA: ",as.numeric(upperThresholdFeatures)))
  
  cat("Select the upper threshold of mitoPercent: \n")
  
  mitoThreshold <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of mitoPercent: ",as.numeric(mitoThreshold)))
  
  seurat <- subset(seurat, subset = nCount_RNA > lowerThresholdRNA &
                     
                     nCount_RNA < upperThresholdRNA &
                     
                     nFeature_RNA > lowerThresholdFeatures &
                     
                     nFeature_RNA < upperThresholdFeatures &
                     
                     mitoPercent < mitoThreshold)
  
  
  
  
  
  #Remove ambient RNA
  
  #seurat <- RemoveAmbientRNA(seurat)
  
  
  
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  
  seurat <- NormalizeData(seurat)
  
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  
  seurat <- ScaleData(seurat)
  
  seurat <- RunPCA(seurat)
  
  seurat <- FindNeighbors(object = seurat)
  
  seurat <- FindClusters(object = seurat)
  
  elbow <- ElbowPlot(seurat)
  
  plot(elbow)
  
  cat("Select the number of dimensions: \n")
  
  dimen <- scan(nmax = 1)
  
  print(paste0("Selected number of dimensions: ",as.numeric(dimen)))
  
  seurat <- RunUMAP(seurat, dims = 1:as.numeric(dimen))
  
  
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  
  sweep.res.seurat <- paramSweep(seurat, PCs = 1:as.numeric(dimen), sct = FALSE)
  
  sweep.stats_seurat <- summarizeSweep(sweep.res.seurat, GT = FALSE)
  
  bcmvn_seurat <- find.pK(sweep.stats_seurat)
  
  pK <- bcmvn_seurat$pK[bcmvn_seurat$BCmetric == max(bcmvn_seurat$BCmetric)]
  
  pK <- as.numeric(as.character(pK[[1]]))
  
  
  
  
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  
  annotations <- seurat@meta.data$seurat_clusters
  
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seurat@meta.data$ClusteringResults
  
  nExp_poi <- round(0.075*nrow(seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  
  seurat <- doubletFinder(seurat, PCs = 1:as.numeric(dimen), pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  seurat <- doubletFinder(seurat, PCs = 1:as.numeric(dimen), pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  
  
  return(seurat)
  
  
  
}





RemoveAmbientRNA <- function(seurat){
  
  
  
  rawSeurat <- NormalizeData(seurat)
  
  rawSeurat <- FindVariableFeatures(rawSeurat, selection.method = "vst", nfeatures = 2000)
  
  rawSeurat <- ScaleData(rawSeurat)
  
  rawSeurat <- RunPCA(rawSeurat, npcs = 20)
  
  rawSeurat <- FindNeighbors(object = rawSeurat, dims = 1:20)
  
  rawSeurat <- FindClusters(object = rawSeurat, resolution = 0.5)
  
  
  
  seurat$soup_group <- rawSeurat@meta.data[['seurat_clusters']]
  
  
  
  path <- paste0("data/",Project(seurat), "/outs/raw_feature_bc_matrix/")
  
  print(path)
  
  raw <- Read10X(data.dir = path)
  
  raw <- raw$`Gene Expression`
  
  sc <- SoupChannel(raw,seurat@assays$RNA@layers$counts, metaData = seurat@assays$RNA@cells)
  
  sc <- setClusters(sc,seurat$soup_group)
  
  sc <- autoEstCont(sc, doPlot=FALSE)
  
  out <- adjustCounts(sc, roundToInt = TRUE)
  
  
  
  seurat[["original.counts"]] <- CreateAssayObject(counts = seurat@assays$RNA@layers$counts)
  
  
  
  seurat@assays$RNA@layers$counts <- out
  
  
  
  return(seurat)
  
}