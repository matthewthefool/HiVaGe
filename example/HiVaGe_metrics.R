### Set-up libraries for metrics
library(scRNAseq)
library(DropletUtils)
library(scater)
library(scDblFinder)
library(celldex)
library(SingleR)
library(Rtsne)
library(scry)
library(scran)
library(igraph)
library(clusterSim)
library(fossil)
library(RBGL)
library(tibble)
library(ROGUE)

### Set-up libraries for HVGs
library(reticulate)
library(magrittr)
library(Seurat)


## Message printing function
replaceCat <- function(x, width = 50)
{
  cat("\r",rep(" ", times = width - length(x)), "\r")
  cat(x)
}


### Dependency with mean expression testing
getOverlapWithHighLowExpressed <- function(hvgs, sce_object, amount_of_genes_to_check = length(hvgs)) {
  # Get pseudo-bulk data from sce_objectExperiment
  pseudo_bulk = matrix(apply(counts(sce_object), 1, sum))
  rownames(pseudo_bulk) = rownames(counts(sce_object))
  # Log-normalized counts of pseudo-bulk data
  pseudo_bulk = normalizeCounts(pseudo_bulk)
  # Order based on expression
  pseudo_bulk = pseudo_bulk[order(pseudo_bulk, decreasing = TRUE),]
  
  # Get the set amount of highly/lowly expressed genes to test overlap
  hegs = names(pseudo_bulk[1:amount_of_genes_to_check])
  legs = names(pseudo_bulk[(length(pseudo_bulk)-amount_of_genes_to_check):length(pseudo_bulk)])
  
  # Calculate and return overlap
  overlap_h = table(unique(c(hvgs)) %in% hegs)
  overlap_l = table(unique(c(hvgs)) %in% legs)
  
  expr_var_cor = cor(log(apply(counts(sce_object), 1, mean)), log(apply(counts(sce_object), 1, var)), method = "pearson")
  
  return(list(highly_overlap = unname(overlap_h["TRUE"]/sum(overlap_h)), lowly_overlap = unname(overlap_l["TRUE"]/sum(overlap_l)), correlation_meanVariance = expr_var_cor))
}


### Purity testing
calcPurity = function(clusters, classes) {
  # Just the purity math formula
  sum(apply(table(classes, clusters), 1, max))/length(clusters)
}

getAvrPurity = function(hvgs, sce_object, labels, tSNE_count = 5) {
  set.seed(42)
  replaceCat(paste("Chosen t-SNE count:", tSNE_count, ".\n", sep = ""))
  # Retrieve relevant data
  sce_logcounts = logcounts(sce_object)
  # Vector for storing purity for each clustering
  purities = c(1:tSNE_count)
  
  # Getting expression data for genes of HVGs
  sce_logcounts_genes = as.data.frame(t(sce_logcounts[hvgs,]))
  # Adding cell-type as separate column
  sce_logcounts_genes$CellType = labels
  # Removing duplicates (ignoring the CellType column)
  sce_logcounts_genes = sce_logcounts_genes[!duplicated(subset(sce_logcounts_genes, select=-c(CellType))), ]
  
  for (i in 1:tSNE_count) {
    replaceCat(paste("Running t-SNE #", i, "...", sep = ""))
    # Running t-SNE (with the library from the paper) on expression data
    sce_tSNE = Rtsne(subset(sce_logcounts_genes, select=-c(CellType)), dims = 2) #, perplexity=50, max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000)
    
    # K-means clustering. I chose 28 clusters, since that's the amount of predicted
    # cell types. If that's not how that should be done, pls change
    sce_cluster <- kmeans(sce_tSNE$Y, center=28, nstart=20)
    
    # Calculating purity for (i+1)th clustering
    purities[i] = calcPurity(sce_logcounts_genes$CellType, sce_cluster$cluster)
  }
  
  # Final purity score for this method (seurat) with this database (MairPBMCData)
  purity_score = mean(purities)
  return(list(purities = purities, purity_score = purity_score))
}


### Calculate ROGUE scores for each cluster for each sample
getROGUEScore <- function(hvgs, sce_object, clustering, sampling, platform, assay.type = "counts", min.cell.n = 10) {
  # Get all the cluster names
  clusters = unique(clustering)
  # Get all the sample names
  samples = unique(sampling)
  # Prepare results matrix
  res = matrix(nrow = length(clusters), ncol = length(samples))
  rownames(res) = clusters
  colnames(res) = samples
  
  replaceCat(paste("Cluster count:", length(clusters), ", sample count:", length(samples), sep = ""))
  
  for (i in 1:length(clusters)) {
    for (j in 1:length(samples)) {
      replaceCat(paste("Working with cluster: ", clusters[i], ", Sample: ", samples[j], sep = ""))
      
      # Get relevant data (cells from cluster i-th and sample j-th)
      if (assay.type == "counts"){
        tmp = counts(sce_object[,clustering == clusters[i] & sampling == samples[j]])
      } else if (assay.type == "logcounts"){
        tmp = logcounts(sce_object[,clustering == clusters[i] & sampling == samples[j]])
      } else {
        stop("Innappropriate assay type")
      }
      
      # Filter out data will too little cells
      if (dim(tmp)[2] >= min.cell.n) {
        # Running the S-E model
        tmp.res <- SE_fun(tmp)
        # Running ROGUE score calculation
        res[i,j] = CalculateRogue(tmp.res, features = hvgs, platform = platform)
      } else {
        # Too little cells!
        res[i,j] = NA
      }
    }
  }
  # Return result matrix
  return(res)
}


### Calculate correlation between HVGs
HiVaGe <- function(sce_object, flavour, percentile) {
  # Should already be preinstalled rPython, scLVM, scVEGs (file "scVEGs.R" should be in same directory as "HiVaGe" function), SIEVE
  num_HVGs <- as.integer(nrow(assays(sce_object)$counts) * percentile)
  sce_object <- logNormCounts(sce_object)
  library(scran)
  library(DropletUtils)
  library(scRNAseq)
  library(magrittr)
  library(Seurat)
  library(scater)
  library(scDblFinder)
  valid_flavours <- c("BASiCS", "M3Drop","M3Drop_Basic", "M3Drop_Brennecke", "ROGUE", "ROGUE_n", "Seurat_vst", "Seurat_sct", "Seurat_disp", "scVEGs", "SCHS", "scmap", "SIEVE_Scmap", "SIEVE_Scran", "SIEVE_ROGUE", "SIEVE_M3Drop", "SIEVE_Seurat_vst", "SIEVE_Seurat_disp", "scLVM_log", "scLVM_logvar")
  if (!(flavour %in% valid_flavours)) {
    stop(paste("Invalid flavour. Allowed values are:", paste(valid_flavours, collapse = ", ")))
  }
  
  # Perform actions based on the specified flavour
  if (flavour == "BASiCS") {
    # BASiCS
    # Could use spike-ins
    library(BASiCS)  # Load the BASiCS package if not already loaded
    BASiCS <- newBASiCS_Data(counts(sce_object), Tech = FALSE, SpikeInfo = NULL, BatchInfo = sce_object$Sample_Name)
    Chain <- BASiCS_MCMC(Data = BASiCS, N = 1000, Thin = 10, Burn = 500, WithSpikes = FALSE, Regression = TRUE, PrintProgress = FALSE)
    HVGs_BASiCS <- BASiCS_DetectVG(Chain, Task = c("HVG"), PercentileThreshold = percentile, VarThreshold = NULL, ProbThreshold = 0.5, EpsilonThreshold = NULL, EFDR = 0.1, Plot = FALSE, MinESS = 100)
    HVGs_BASiCS <- unlist(as.data.frame(HVGs_BASiCS)$GeneName)
    return(HVGs_BASiCS)
    
  } else if (flavour == "M3Drop") {
    # M3Drop 
    # Raw UMI counts
    # Variance-based Feature Selection. Ranks genes by residual dispersion from mean-dispersion power-law relationship.
    # Uses linear regression on log-transformed mean expression and fitted dispersions. Ranks genes by the residuals of this fit, negative = high variance, positive = low variance.
    library(M3Drop)
    HVGs_M3Drop <- sce_object %>%
      counts() %>%
      NBumiConvertData() %>%
      NBumiFitModel() %>%
      NBumiFeatureSelectionHighVar() %>%
      as.data.frame() %>%
      rownames() %>%
      unlist() %>%
      .[1:num_HVGs]
    return(HVGs_M3Drop)
    
  } else if (flavour == "M3Drop_Basic") {
    # M3Drop 
    # Raw UMI counts
    library(M3Drop)
    HVGs_M3Drop <- sce_object %>%
      counts() %>%
      NBumiConvertData() %>%
      NBumiFitBasicModel() %>%
      NBumiFeatureSelectionHighVar() %>%
      as.data.frame() %>%
      rownames() %>%
      unlist() %>%
      .[1:num_HVGs]
    
    return(HVGs_M3Drop)
    
  } else if (flavour == "M3Drop_Brennecke") {
    # Spike-ins is possible
    # Normalized or raw (not log-transformed) expression values, columns = samples, rows = genes.
    library(M3Drop)
    HVGs_M3Drop_Brennecke <- sce_object %>%
      counts() %>%
      BrenneckeGetVariableGenes() %>%
      rownames()
    return(HVGs_M3Drop_Brennecke)
    
  } else if (flavour == "ROGUE") {
    # ROGUE
    # Raw counts
    library(ROGUE)
    count_matrix <- assay(sce_object, "counts") %>%
      as.matrix()
    rownames(count_matrix) <- rownames(sce_object)
    colnames(count_matrix) <- colnames(sce_object)
    ds_ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
    HVGs_ROGUE <- unlist(ds_ROGUE$Gene)[1:num_HVGs]
    return(HVGs_ROGUE)
    
  } else if (flavour == "ROGUE_n") {
    # ROGUE_n
    # Normalized counts
    library(ROGUE)
    count_matrix <- assay(sce_object, "logcounts") %>%
      as.matrix()
    rownames(count_matrix) <- rownames(sce_object)
    colnames(count_matrix) <- colnames(sce_object)
    ds_ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
    HVGs_ROGUE_n <- unlist(ds_ROGUE$Gene)[1:num_HVGs]
    return(HVGs_ROGUE_n)
    
  } else if (flavour == "Seurat_vst") {
    # Seurat_vst
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = num_HVGs, verbose = TRUE)
    HVGs_seurat_vst <- ds_seurat@assays$RNA@var.features
    return(HVGs_seurat_vst)
    
  } else if (flavour == "Seurat_sct") {
    # Seurat_sct
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      SCTransform(assay="RNA",variable.features.n = num_HVGs)
    HVGs_seurat_sct <- ds_seurat@assays$SCT@var.features
    return(HVGs_seurat_sct)
    
  } else if (flavour == "Seurat_disp") {
    # Seurat_disp
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      FindVariableFeatures(selection.method = "disp", nfeatures = num_HVGs, verbose = TRUE)
    HVGs_seurat_disp <- ds_seurat@assays$RNA@var.features
    return(HVGs_seurat_disp)
    
  } else if (flavour == "scVEGs") {
    source("scVEGs.R")
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    pVal <- 0.01
    pFlag <- 1
    species <- 'hs'
    cellSum <- apply(ds_counts_df, 2, sum)
    scaleNum <- mean(cellSum)/cellSum
    scaleFactor <- t(kronecker( matrix(1, 1, dim(ds_counts_df)[1]), scaleNum))
    normData <- scaleFactor * ds_counts_df
    outputName <- 'HVGs_scVEGs'
    HVGs_scVEGs <- scVEGs(normData, pVal, pFlag, species, outputName) %>%
      rownames()
    return(HVGs_scVEGs)
    
  } else if (flavour == "SCHS") {
    library("singleCellHaystack")
    set.seed(42)
    HVGs_SCHS <- reducedDim(sce_object, "PCA") %>%
      haystack(expression = logcounts(sce_object)) %>%
      show_result_haystack(n = num_HVGs) %>%
      rownames()
    return(HVGs_SCHS)
    
  } else if (flavour == "scmap") {
    library("scmap")
    rowData(sce_object)$feature_symbol <- rowData(sce_object)$Symbol
    sce_object <- selectFeatures(sce_object, n_features = num_HVGs,
                                 suppress_plot = TRUE) 
    HVGs_scmap <- rowData(sce_object)[rowData(sce_object)$scmap_features == TRUE, ]
    return(HVGs_scmap)
    
  } else if (flavour == "SIEVE_Scmap") {
    # If your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap". If your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".
    # SCHS, Seurat_sct doesn't work
    library(SIEVE)
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Scmap <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Scmap",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Scmap)
    
  } else if (flavour == "SIEVE_Scran") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Scran <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Scran",n=num_HVGs) %>%
      SIEVE(n=50)
    return (HVGs_SIEVE_Scran)
    
  } else if (flavour == "SIEVE_ROGUE") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_ROGUE <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="ROGUE",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_ROGUE)
    
  } else if (flavour == "SIEVE_M3Drop") {
    library(SIEVE)
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_M3Drop <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="M3Drop",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_M3Drop)
    
  } else if (flavour == "SIEVE_Seurat_vst") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Seurat_vst <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Seurat_vst",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Seurat_vst)
    
  } else if (flavour == "SIEVE_Seurat_disp") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Seurat_disp <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Seurat_disp",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Seurat_disp)
    
  } else if (flavour == "scLVM_log") {
    # could use spike-ins and size Factors. As threshhold used default
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 
    
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_log <- as.data.frame(is_hetLog[is_hetLog == TRUE]) %>%
      rownames()
    return(HVGs_scLVM_log)
    
  } else if (flavour == "scLVM_logvar") {
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoiseLogFit = fitTechnicalNoise(norm_counts, fit_type = 'logvar', use_ERCC = FALSE, plot=FALSE) 
    
    is_hetLog = getVariableGenes(norm_counts, techNoiseLogFit$fit, plot=FALSE)
    HVGs_scLVM_logvar <- as.data.frame(is_hetLog[is_hetLog == TRUE]) %>%
      rownames()
    return(HVGs_scLVM_logvar)
  }
}

correlatedHVGs <- function(hvgs, sce_object) {
  cor_genes <- correlatePairs(sce_object)
  cor_genes <- cor_genes[(cor_genes$gene1 %in% hvgs) & (cor_genes$gene2 %in% hvgs),]
  # Get significant correlated pairs
  sig_cor_genes <- cor_genes$FDR <= 0.05
  #summary(sig_cor_genes)
  # Build undirected graph of correlated genes
  g_cor_genes <- ftM2graphNEL(cbind(cor_genes$gene1, cor_genes$gene2)[sig_cor_genes,],
                              W=NULL, V=NULL, edgemode="undirected")
  # Get clusters of correlated genes
  cl_cor_genes <- highlyConnSG(g_cor_genes)$clusters
  # Order them by size
  cl_cor_genes <- cl_cor_genes[order(lengths(cl_cor_genes), decreasing=TRUE)]
  
  # Prepare result list
  res = list()
  res$gene_correlation_df = cor_genes
  res$correlated_gene_clusters = cl_cor_genes
  # Return results
  return(res)
}


### Get HVGs based on method
