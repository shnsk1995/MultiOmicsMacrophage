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

IntegrateSeuratObjectRNAv5 <- function(seurat){
  
  DefaultAssay(seurat) <- "RNA"
  
  seurat <- IntegrateLayers(seurat,
                            
                            method = CCAIntegration,
                            
                            orig.reduction = "pca.rna",
                            
                            assay = "RNA",
                            
                            new.reduction = "integrated.cca")
  
  
  seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "umap.integrated_rna", reduction = "integrated.cca")
  
  
  
  return(seurat)
  
  
  
}


IntegrateSeuratObjectRNAv4 <- function(seurat){
  
  DefaultAssay(seurat) <- "RNA"
  
  seurat.list <- SplitObject(seurat,split.by = "orig.ident")
  
  for (i in 1:length(seurat.list)) {
    
    seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], nfeatures = 10000)
    seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]])
    seurat.list[[i]] <- ScaleData(object = seurat.list[[i]])
    
  }
  
  features <- SelectIntegrationFeatures(object.list = seurat.list)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                    anchor.features = features)
  
  seurat.integrated <- IntegrateData(anchorset = anchors, new.assay.name = "Integrated_RNA")
  
  seurat.integrated <- ScaleData(object = seurat.integrated)
  
  seurat.integrated <- RunPCA(object = seurat.integrated)
  
  seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:30, reduction.name = "umap.integrated_rna")
  
  return(seurat.integrated)
  
  
}

IntegrateSeuratObjectATACSeurat <- function(seurat){
  
  browser()
  
  DefaultAssay(seurat) <- "ATAC"
  # 
  # seurat <- FindTopFeatures(seurat, min.cutoff = 10)
  # 
  # seurat <- RunTFIDF(seurat)
  # 
  # seurat <- RunSVD(seurat)
  # 
  # seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:30, reduction.name = 'umap.atac')
  
  seurat.list <- SplitObject(seurat,split.by = "orig.ident")
  
  for (i in 1:length(seurat.list)) {
    
    DefaultAssay(seurat.list[[i]]) <- "ATAC"
    seurat.list[[i]] <- FindTopFeatures(seurat.list[[i]], min.cutoff = 10)
    seurat.list[[i]] <- RunTFIDF(seurat.list[[i]])
  }
  
  
  features <- SelectIntegrationFeatures(seurat.list, nfeatures = 2000)
  
  anchors <- FindIntegrationAnchors(
    object.list = seurat.list,
    anchor.features = features,
    assay = c('ATAC','ATAC','ATAC','ATAC')
  )
  
  seurat.integrated <- IntegrateData(anchorset = anchors, new.assay.name = "Integrated_ATAC")
  
  # seurat <- RunUMAP(seurat,dims=2:30,reduction.name ='umap.integrated_atac', reduction = 'harmony')
  # 
  # 
  # 
  # 
  # # 
  # # # peaks.use <- Reduce(intersect, list(VariableFeatures(individualATACObjects[[1]]),
  # # #                                     VariableFeatures(individualATACObjects[[2]]),
  # # #                                     VariableFeatures(individualATACObjects[[3]]),
  # # #                                     VariableFeatures(individualATACObjects[[4]])))
  # # 
  # 
  # # 
  # # browser()
  # # 
  # 
  # # 
  # # any(is.na(rownames(GetAssayData(individualATACObjects[[4]], assay = "ATAC", slot = "counts"))))
  # # 
  # # # integrate LSI embeddings
  # # integrated <- IntegrateEmbeddings(
  # #   anchorset = integration.anchors,
  # #   reductions = seurat[["lsi"]],
  # #   new.reduction.name = "integrated_lsi",
  # #   dims.to.integrate = 1:30
  # # )
  # # 
  # # # create a new UMAP using the integrated embeddings
  # # integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 1:30)
  
  
  
  return(seurat.integrated)
  
  
  
}

IntegrateSeuratObjectATACHarmony <- function(seurat){
  
  DefaultAssay(seurat) <- "ATAC"
  
  seurat <- RunHarmony(object=seurat,group.by.vars='orig.ident',reduction.use='lsi',assay.use='ATAC', project.dim=F)
  seurat <- RunUMAP(seurat,dims=1:30,reduction.name ='umap.integrated_atac', reduction = 'harmony')
  
  

  # for (i in 1:length(individualATACObjects)) {
  # 
  #   DefaultAssay(individualATACObjects[[i]]) <- "ATAC"
  #   individualATACObjects[[i]] <- FindTopFeatures(individualATACObjects[[i]], min.cutoff = 50)
  #   individualATACObjects[[i]] <- RunTFIDF(individualATACObjects[[i]])
  # 
  # }
  # 
  # # peaks.use <- Reduce(intersect, list(VariableFeatures(individualATACObjects[[1]]),
  # #                                     VariableFeatures(individualATACObjects[[2]]),
  # #                                     VariableFeatures(individualATACObjects[[3]]),
  # #                                     VariableFeatures(individualATACObjects[[4]])))
  # 
  # peaks.use <- Reduce(intersect, list(VariableFeatures(individualATACObjects[[1]]),
  #                                     VariableFeatures(individualATACObjects[[2]]),
  #                                     VariableFeatures(individualATACObjects[[3]]),
  #                                     VariableFeatures(individualATACObjects[[4]])))
  # peaks.use <- sample(peaks.use, 20000)
  # 
  # browser()
  # 
  # integration.anchors <- FindIntegrationAnchors(
  #   object.list = individualATACObjects,
  #   anchor.features = peaks.use,
  #   assay = c('ATAC','ATAC','ATAC','ATAC'),
  #   scale = TRUE
  # )
  # 
  # any(is.na(rownames(GetAssayData(individualATACObjects[[4]], assay = "ATAC", slot = "counts"))))
  # 
  # # integrate LSI embeddings
  # integrated <- IntegrateEmbeddings(
  #   anchorset = integration.anchors,
  #   reductions = seurat[["lsi"]],
  #   new.reduction.name = "integrated_lsi",
  #   dims.to.integrate = 1:30
  # )
  # 
  # # create a new UMAP using the integrated embeddings
  # integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 1:30)
  
  
  
  return(seurat)
  
  
  
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

CalculateATACQCMetrics <- function(seurat){
  
  DefaultAssay(seurat) <- "ATAC"
  
  seurat <- TSSEnrichment(seurat)
  
  seurat <- NucleosomeSignal(seurat)
  
  seurat$blacklist_fraction <- FractionCountsInRegion(
    
    object = seurat,
    
    assay = 'ATAC',
    
    regions = blacklist_mm10
    
  )
  
  return(seurat)
  
}

DoQCForATACAssay <- function(seurat){
  
  DefaultAssay(seurat) <- "ATAC"
  
  Idents(seurat) <- "orig.ident"
  
  
  qcPlot <- VlnPlot(
    
    seurat,
    
    features = c("nCount_ATAC", "TSS.enrichment",
                 
                 "nucleosome_signal", "blacklist_fraction"),
    
    pt.size = 0.1,
    
    ncol = 4,
    
    group.by = NULL
    
  )
  
  
  ggsave(filename = paste0("data/",seurat@project.name,"/outs/ATACQCBeforeFiltering.jpeg"), plot = qcPlot)
  
  plot(qcPlot)
  
  cat("Select the lower threshold of nCount_ATAC: \n")
  
  lowerThresholdATAC <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of nCount_ATAC: ",as.numeric(lowerThresholdATAC)))
  
  cat("Select the upper threshold of nCount_ATAC: \n")
  
  upperThresholdATAC <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of nCount_ATAC: ",as.numeric(upperThresholdATAC)))
  
  cat("Select the lower threshold of TSSErichment: \n")
  
  lowerThresholdTSSErichment <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of TSSErichment: ",as.numeric(lowerThresholdTSSErichment)))
  
  cat("Select the upper threshold of TSSErichment: \n")
  
  upperThresholdTSSErichment <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of TSSErichment: ",as.numeric(upperThresholdTSSErichment)))
  
  cat("Select the lower threshold of NucleosomeSignal: \n")
  
  lowerThresholdNucleosomeSignal <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of NucleosomeSignal: ",as.numeric(lowerThresholdNucleosomeSignal)))
  
  cat("Select the upper threshold of NucleosomeSignal: \n")
  
  upperThresholdNucleosomeSignal <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of NucleosomeSignal: ",as.numeric(upperThresholdNucleosomeSignal)))
  
  cat("Select the lower threshold of BlacklistFraction: \n")
  
  lowerThresholdBlacklistFraction <- scan(nmax = 1)
  
  print(paste0("Selected lower threshold of BlacklistFraction: ",as.numeric(lowerThresholdBlacklistFraction)))
  
  cat("Select the upper threshold of BlacklistFraction: \n")
  
  upperThresholdBlacklistFraction <- scan(nmax = 1)
  
  print(paste0("Selected upper threshold of BlacklistFraction: ",as.numeric(upperThresholdBlacklistFraction)))
  
  
  seurat <- subset(
    
    x = seurat,
    
    subset = blacklist_fraction < upperThresholdBlacklistFraction &
      
      blacklist_fraction > lowerThresholdBlacklistFraction &
      
      TSS.enrichment < upperThresholdTSSErichment &
      
      TSS.enrichment > lowerThresholdTSSErichment &
      
      CellType == "Singlet" &
      
      nCount_ATAC < upperThresholdATAC &
      
      nCount_ATAC > lowerThresholdATAC &
      
      nucleosome_signal < upperThresholdNucleosomeSignal &
      
      nucleosome_signal > lowerThresholdNucleosomeSignal
    
  )
  
  
  
  qcPlot <- VlnPlot(
    
    seurat,
    
    features = c("nCount_ATAC", "TSS.enrichment",
                 
                 "nucleosome_signal", "blacklist_fraction"),
    
    pt.size = 0.1,
    
    ncol = 4,
    
    group.by = NULL
    
  )
  
  
  ggsave(filename = paste0("data/",seurat@project.name,"/outs/ATACQCAfterFiltering.jpeg"), plot = qcPlot)
  
  
  return(seurat)
  
}

Preprocess <- function(seurat){
  
  DefaultAssay(seurat) <- "RNA"
  
  
  
  seurat <- FindVariableFeatures(seurat, nfeatures = 10000)
  
  seurat <- NormalizeData(seurat)
  
  seurat <- ScaleData(seurat)
  
  seurat <- RunPCA(seurat, npcs = 20, reduction.name = "pca.rna")
  
  seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "umap.rna", reduction = "pca.rna")
  
  #seurat <- RunUMAP(seurat, dims = 1:20, reduction.name = "umap.rna", reduction = "integrated.cca")
  
  #seurat <- FindNeighbors(seurat, dims = 1:20, reduction = "integrated.cca")
  
  
  DefaultAssay(seurat) <- "ATAC"
  
  seurat <- FindTopFeatures(seurat, min.cutoff = 10)
  
  seurat <- RunTFIDF(seurat)
  
  seurat <- RunSVD(seurat)
  
  seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 1:30, reduction.name = 'umap.atac')
  
  
  
  return(seurat)
  
  
  
}

FindDoublets <- function(seurat){
  
  
  
  #Set default Assay as RNA
  
  DefaultAssay(seurat) <- "RNA"
  
  
  
  #QC
  
  seurat$mitoPercent <- PercentageFeatureSet(seurat,pattern = "^mt-")
  
  qcPlot <- VlnPlot(seurat,features = c("nCount_RNA","nFeature_RNA","mitoPercent"),ncol = 3)
  
  plot(qcPlot)
  
  ggsave(filename = paste0("data/",seurat@project.name,"/outs/QCMetricsBeforeFiltering.jpeg"), plot = qcPlot)
  
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
  
  
  
  qcPlot <- VlnPlot(seurat,features = c("nCount_RNA","nFeature_RNA","mitoPercent"),ncol = 3)
  
  
  
  ggsave(filename = paste0("data/",seurat@project.name,"/outs/QCMetricsAfterFiltering.jpeg"), plot = qcPlot)
  
  
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
  
  ggsave(filename = paste0("data/",seurat@project.name,"/outs/Elbow.jpeg"), plot = elbow)
  
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

VisualizeClustree <- function(lowerResolutionThreshold, higherResolutionThreshold, thresholdIncrement, seuratObject){
  
  
  DefaultAssay(seuratObject) <- "RNA"
  
  seuratObject <- FindNeighbors(seuratObject, reduction="integrated.cca")
  
  for (i in seq(lowerResolutionThreshold,higherResolutionThreshold,thresholdIncrement)) {
    
    seuratObject <- FindClusters(seuratObject, resolution = i, algorithm = 3,cluster.name = paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",i))
    
  }
  
  
  cluster_data <- data.frame(
    
    cell_id = colnames(seuratObject),
    
    resolution_1 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold)]], 
    
    resolution_2 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + thresholdIncrement)]], 
    
    resolution_3 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (2 * thresholdIncrement))]],
    
    resolution_4 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (3 * thresholdIncrement))]], 
    
    resolution_5 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (4 * thresholdIncrement))]], 
    
    resolution_6 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (5 * thresholdIncrement))]],
    
    resolution_7 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (6 * thresholdIncrement))]], 
    
    resolution_8 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (7 * thresholdIncrement))]], 
    
    resolution_9 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (8 * thresholdIncrement))]],
    
    resolution_10 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (9 * thresholdIncrement))]],
    
    resolution_11 = seuratObject[[paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_",lowerResolutionThreshold + (10 * thresholdIncrement))]]
)
  
  
  # # # Convert to long format
  # # 
  # cluster_data_long <- cluster_data %>%
  # 
  #   gather(key = "resolution", value = "cluster", -cell_id) %>%
  # 
  #   arrange(resolution, cell_id)  # Make sure the data is sorted
  
  
  
  # Use clustree to visualize cluster stability across resolutions
  
  clust <- clustree(cluster_data, prefix = paste0("cca_clusters_",lowerResolutionThreshold,"_",higherResolutionThreshold,"_"))
  
  
  jpeg(filename = paste0("Clustree_",lowerResolutionThreshold,"_",higherResolutionThreshold,".jpeg"), quality = 100,height = 2000, width = 2000)
  
  print(clust)
  
  dev.off()
  
  
  return(seuratObject)
  
}

GetGenomeRanges <- function(seuratObject, genes){
  
  
  
  if(length(genes) == 1){
    
    geneRanges <- seuratObject@assays$ATAC@annotation[seuratObject@assays$ATAC@annotation$gene_name == genes,]
    
    DefaultAssay(seuratObject) <- "ATAC"
    
    allRegions <- rownames(seuratObject)
    
    allRegions <- strsplit(allRegions, "[--]")
    
    allRegionsDF <- do.call(rbind, lapply(allRegions, function(x) {
      data.frame(Chromosome = x[1], Start = as.numeric(x[2]), End = as.numeric(x[3]))
    }))
    
    allGeneRanges <- GRanges(
      seqnames = allRegionsDF$Chromosome,
      ranges = IRanges(start = allRegionsDF$Start, end = allRegionsDF$End)
    )
    
    overlaps <- findOverlaps(geneRanges, allGeneRanges)
    
    overlappingRanges <- allGeneRanges[subjectHits(overlaps)]
    
    overlappingRanges
    
    genomicRanges <- paste(seqnames(overlappingRanges), start(overlappingRanges), end(overlappingRanges), sep = "-")
    
    return(genomicRanges)
  }
  
}

GenerateCoveragePlots <- function(seuratObject, genesOfInterest){
  
  Idents(seuratObject) <- "SampleType"
  
  
  if(!dir.exists("data/CoveragePlots")){
           dir.create("data/CoveragePlots")
  }
  
  for (gene in genesOfInterest) {
    
      
      tryCatch({
      c <- CoveragePlot(seuratObject,assay = "ATAC",
                   region = gene,
                   features = gene,
                   expression.assay = "RNA",
                   extend.upstream = 5000,
                   extend.downstream = 5000,
                   split.by = "orig.ident"
                   )
      
      
        ggsave(filename = paste0("data/CoveragePlots/",gene,".jpeg"), plot = c, height = 15, width = 25, units = "cm")
        # jpeg(filename = paste0("data/CoveragePlots/",gene,".jpeg"),quality = 100)
        # 
        # print(c)
        # 
        # dev.off()
      }, error = function(e) {
        # Handle errors by printing the message and skipping the gene
        message("Skipping ", gene, ": ", e$message)
      })
    
  }
  
  # if(length(geneRegions) < 6){
  #   
  #   ranges.show <- StringToGRanges(geneRegions)
  #   
  #   c <- CoveragePlot(seuratObject,assay = "ATAC",
  #                region = geneRegions,
  #                features = geneName, 
  #                expression.assay = "RNA",
  #                region.highlight = ranges.show, 
  #                split.by = "Corrected_Seurat")
  #   
  #   if(!dir.exists("data/CoveragePlots")){
  #     dir.create("data/CoveragePlots")
  #   }
  #   
  #   if(!dir.exists(paste0("data/CoveragePlots/",geneName))){
  #     dir.create(paste0("data/CoveragePlots/",geneName))
  #   }
  #   
  #   jpeg(filename = paste0("data/CoveragePlots/",geneName,"/",geneName,".jpeg"),quality = 100, height = 800, width = 1000)
  #   
  #   print(c)
  #   
  #   dev.off()
  #   
  # }else{
  #   
  #   iterations <- (length(geneRegions)/6) + 1
  #   
  #   j <- 1
  #   k <- 6
  #   for (i in 1:iterations) {
  #     
  #     if(k>length(geneRegions)) k <- length(geneRegions)
  #     
  #     ranges.show <- StringToGRanges(geneRegions[j:k])
  #     
  #     c <- CoveragePlot(seuratObject,assay = "ATAC",
  #                       region = geneRegions[j:k],
  #                       features = geneName, 
  #                       expression.assay = "RNA",
  #                       region.highlight = ranges.show, 
  #                       split.by = "Corrected_Seurat")
  #     
  #     if(!dir.exists("data/CoveragePlots")){
  #       dir.create("data/CoveragePlots")
  #     }
  #     
  #     if(!dir.exists(paste0("data/CoveragePlots/",geneName))){
  #       dir.create(paste0("data/CoveragePlots/",geneName))
  #     }
  #     
  #     jpeg(filename = paste0("data/CoveragePlots/",geneName,"/",geneName,"_",j,"_",k,".jpeg"),quality = 100, height = 800, width = 1000)
  #     
  #     print(c)
  #     
  #     dev.off()
  #     
  #     j <- j+6
  #     k <- k+6
  #     
  #     
  #   }
  #   
  # }
  
  
  
}

FindCommonDGEDA <- function(RNADGE, ATACDA){
  
  combinedResults <- merge(
    x = RNADGE,
    y = ATACDA,
    by.x = "row.names", # Gene names in DGE
    by.y = "gene",      # Gene names in DA
    suffixes = c("_RNA", "_ATAC")
  )
  
  
  consistentResults <- combinedResults[
    (combinedResults$logFC_RNA > 0 & combinedResults$logFC_ATAC > 0) |
      (combinedResults$logFC_RNA < 0 & combinedResults$logFC_ATAC < 0),
  ]
  
  
  #Filter the rows based on adj p value significance of both RNA and ATAC data
  consistentResults <- consistentResults[
    consistentResults$adjpvalue_RNA < 0.05 &
      consistentResults$adjpvalue_ATAC < 0.05,
  ]
  
  
  consistentResults <- consistentResults[order(-abs(consistentResults$logFC_RNA)), ]
  
  return(consistentResults)
  
}

WriteClusterSignificantGenes <- function(seuratObject, ids){
  
  Idents(seuratObject) <- "Corrected_Seurat"
  
  for (id in ids) {
    
    
    subSeurat <- subset(seuratObject, idents = id)
    
    RNADGE <- FindMarkers(
      
      object = subSeurat,
      
      ident.1 = "KnockOut",
      
      ident.2 = "WildType",
      
      group.by = "SampleType",
      
      logfc.threshold = 1,
      
      test.use = "wilcox",
      
      assay = "RNA"
      
    )
    
    sigRNADGE <- RNADGE[RNADGE$p_val_adj < 0.05,]
    
    write(paste0(rownames(sigRNADGE),"\n"), file = paste0("data/",id,"DGESIGGENES.txt"))
    
    
  }
  
  
  
}

WriteClusterMarkerGenes <- function(seuratObject, ids){
  
  Idents(seuratObject) <- "Corrected_Seurat"
  
  allMarkers <- FindAllMarkers(seuratObject, assay = "RNA", logfc.threshold = 1)
  
  write.csv(allMarkers, file = 'data/AllMarkers.csv', row.names = TRUE)
    
  sigMarkers <- allMarkers[allMarkers$p_val_adj < 0.05,]
    
  for (id in ids) {
    
    write(paste0(rownames(sigMarkers[sigMarkers$cluster==id,]),"\n"), file = paste0("data/",id,"MARKSIGGENES.txt"))
    
  }
    
    
}

DoEnrichRAnnotation <- function(seuratObject, ids){
  
  Idents(seuratObject) <- "Corrected_Seurat"
  
  seuratObject@meta.data$er_adj_pval_annotation <- NA
  
  seuratObject@meta.data$er_odds_ratio_annotation <- NA
  
  seuratObject@meta.data$er_combined_score_annotation <- NA
  
  for (id in ids) {
    
    genes <- read.table(paste0("data/",id,"MARKSIGGENES.txt"))
    enrichmentResults <- enrichr(genes = genes$V1 , databases = "Tabula_Muris")
    res <- enrichmentResults$Tabula_Muris
    adjpvalTerm <- res$Term[res$Adjusted.P.value == min(res$Adjusted.P.value)]
    if(length(adjpvalTerm)>1) adjpvalTerm <- adjpvalTerm[1]
    oddsratioTerm <- res$Term[res$Odds.Ratio == max(res$Odds.Ratio)]
    if(length(oddsratioTerm)>1) oddsratioTerm <- oddsratioTerm[1]
    combinedscoreTerm <- res$Term[res$Combined.Score == max(res$Combined.Score)]
    if(length(combinedscoreTerm)>1) combinedscoreTerm <- combinedscoreTerm[1]
    seuratObject@meta.data[seuratObject@meta.data$Corrected_Seurat == id, "er_adj_pval_annotation"] <- adjpvalTerm
    seuratObject@meta.data[seuratObject@meta.data$Corrected_Seurat == id, "er_odds_ratio_annotation"] <- oddsratioTerm
    seuratObject@meta.data[seuratObject@meta.data$Corrected_Seurat == id, "er_combined_score_annotation"] <- combinedscoreTerm
    
    
  }
  
  return(seuratObject)
  
}

# Create a new Seurat object with only one assay
SubsetSeurat <- function(seurat, assayName) {
  
  newObj <- CreateSeuratObject(counts = seurat[[assayName]], assay = assayName)
  return(newObj)
}

CapitalizeFirst <- function(genes) {
  sapply(genes, function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
}


#Function to perform cell cycle scoring and plot marker genes in each cycle
DoCellCycleScoring <- function(seurat){
  
  DefaultAssay(seurat) <- "RNA"
  
  
  # Manually curated mouse orthologs of Seurat's cell cycle genes
  s.genes.mouse <- cc.genes$s.genes
  
  g2m.genes.mouse <- cc.genes$g2m.genes
  
  # Convert the cell cycle genes
  s.genes.mouse <- CapitalizeFirst(s.genes.mouse)
  g2m.genes.mouse <- CapitalizeFirst(g2m.genes.mouse)
  
  
  seurat <- CellCycleScoring(seurat, 
                                 s.features = s.genes.mouse, 
                                 g2m.features = g2m.genes.mouse, 
                                 set.ident = TRUE)
  
  Idents(seurat) <- "Phase"
  
  
  G1Markers <- FindMarkers(seurat, ident.1 = "G1", only.pos = TRUE)
  SMarkers <- FindMarkers(seurat, ident.1 = "S", only.pos = TRUE)
  G2MMarkers <- FindMarkers(seurat, ident.1 = "G2M", only.pos = TRUE)
  
  
  Markers <- c(head(rownames(G1Markers),3),head(rownames(G2MMarkers),2),head(rownames(SMarkers),1))
  
  
  ridgePlot <- RidgePlot(seurat, features = Markers)
  
  ggsave(filename = "data/TopCellPhaseMarkers.jpeg", plot = ridgePlot)
  
  
  return(seurat)
  
}