# Set-up libraries
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


# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

# Create mapping vector of rownames to Symbol
symbol_map = rowData(sce)$Symbol
names(symbol_map) = rownames(sce)

symbol_map

# Fix rownames (duplicates and use of dashes)
rownames(sce) <- rowData(sce)$Symbol
new_rownames = gsub("-", "_", rownames(sce))
rownames(sce) = make.names(new_rownames, unique=TRUE)

standardizeGeneSymbols = function(map, genes) {
  # Voodoo magic that makes all gene names the same across all datasets.. hopefully
  standard_gene_symbols = unname(map[genes])
  standard_gene_symbols[is.na(standard_gene_symbols)] = genes[is.na(standard_gene_symbols)]
  standard_gene_symbols = gsub("-", "_", standard_gene_symbols)  
  return(standard_gene_symbols)
}


### Data preparation (copied from profiling)
# Filtering for empty droplets
# Set random seed for reproducibility
set.seed(42)

# Identify empty droplets using 'emptyDrops' function
e.out <- emptyDrops(counts(sce))

# Summary of empty droplet filtering results
summary(e.out$FDR <= 0.001)

# Create a table showing the counts threshold and FDR filtering results
table(colSums(counts(sce)) > 100, e.out$FDR <= 0.001, useNA = "ifany")

# Apply the filtering based on FDR threshold to 'sce'
sce <- sce[, which(e.out$FDR <= 0.001)]
sce

# Summary of columns and rows with all-zero counts
summary(colMeans(counts(sce) == 0))
summary(rowMeans(counts(sce) == 0))

# Identify cells with all-zero counts for all genes
allzero <- rowMeans(counts(sce) == 0) == 1

# Create a table showing the number of cells with all-zero counts
table(allzero)

# Remove cells with all-zero counts for all genes from 'sce'
sce <- sce[which(!allzero),]

# Log-normalization
# Perform log-normalization on 'sce' using logNormCounts function
sce <- logNormCounts(sce)

# Dimensionality reduction with PCA
# Run PCA on 'sce'
sce <- runPCA(sce)

# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(sce, d = ncol(reducedDim(sce)))

# Add doublet scores to the 'sce' object
sce$DoubletScore <- dbl.dens

# Summary of doublet scores
summary(dbl.dens)

# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")

# Remove doublets from 'sce' based on the identified doublet calls
sce <- sce[, dbl.calls == "singlet"]

# Perform PCA again
sce <- runPCA(sce)

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)



### Metrics
# Retrieve HVGs (M3Drop and seurat)
data <- read.csv("M3Drop_matrix.csv", header = FALSE)
M3Drop_matrix <- as.matrix(data)
M3Drop_matrix <- M3Drop_matrix[2:31, ]
# Standardize names for HVGs
M3Drop_matrix = standardizeGeneSymbols(symbol_map, M3Drop_matrix[,1])

data <- read.csv("seurat_matrix.csv", header = FALSE)
seurat_matrix <- as.matrix(data)
seurat_matrix <- seurat_matrix[2:31, ]
# Standardize names for HVGs
seurat_matrix = standardizeGeneSymbols(symbol_map, seurat_matrix[,1])

### Purity testing
## Define functions
getPurity = function(clusters, classes) {
  # Just the purity math formula
  sum(apply(table(classes, clusters), 1, max))/length(clusters)
}

getPurityByMethod = function(hvgs, singleCell, tSNE_count = 5) {
  print(paste("Chosen t-SNE count:", tSNE_count, sep = ""))
  # Retrieve relevant data
  sce_logcounts = logcounts(singleCell)
  # Vector for storing purity for each clustering
  purities = c(1:tSNE_count)
  
  # Getting expression data for genes of HVGs
  sce_logcounts_genes = as.data.frame(t(sce_logcounts[hvgs,]))
  # Adding cell-type as separate column
  sce_logcounts_genes$CellType = singleCell$CellTypes
  # Removing duplicates (ignoring the CellType column)
  sce_logcounts_genes = sce_logcounts_genes[!duplicated(subset(sce_logcounts_genes, select=-c(CellType))), ]
  
  for (i in 1:tSNE_count) {
    print(paste("Running t-SNE #", i, "...", sep = ""))
    # Running t-SNE (with the library from the paper) on expression data
    sce_tSNE = Rtsne(subset(sce_logcounts_genes, select=-c(CellType)), dims = 2) #, perplexity=50, max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000)
    
    # K-means clustering. I chose 28 clusters, since that's the amount of predicted
    # cell types. If that's not how that should be done, pls change
    sce_cluster <- kmeans(sce_tSNE$Y, center=28, nstart=20)
    
    # Calculating purity for (i+1)th clustering
    purities[i] = getPurity(sce_logcounts_genes$CellType, sce_cluster$cluster)
  }
  
  # Final purity score for this method (seurat) with this database (MairPBMCData)
  purity_score = mean(purities)
  return(list(purities = purities, purity_score = purity_score))
}

## Example use
# Run purity scoring function on HVGs
M3Drop_purity_score = getPurityByMethod(M3Drop_matrix, sce)
seurat_purity_score = getPurityByMethod(seurat_matrix, sce)

# First list is 5 purity scores of 5 t-SNE runs, second list is the average value
# and thus it's the final score
M3Drop_purity_score
seurat_purity_score


### Dependency with mean expression testing
## Define function
getOverlapWithHighlyExpressed <- function(hvgs, singlecell, amount_of_genes_to_check = 30) {
  # Get pseudo-bulk data from SingleCellExperiment
  pseudo_bulk = matrix(apply(counts(singlecell), 1, sum))
  rownames(pseudo_bulk) = rownames(counts(singlecell))
  # Log-normalized counts of pseudo-bulk data
  pseudo_bulk = normalizeCounts(pseudo_bulk)
  # Order based on expression
  pseudo_bulk = pseudo_bulk[order(pseudo_bulk, decreasing = TRUE),]
  
  # Get the set amount of highly expressed genes to test overlap
  hegs = names(pseudo_bulk[1:amount_of_genes_to_check])

  # Calculate and return overlap
  overlap = table(unique(c(hvgs)) %in% hegs)
  return(unname(overlap["TRUE"]/sum(overlap)))
}

## Example use
# Run the dependency test
M3Drop_meanDependence = getOverlapWithHighlyExpressed(M3Drop_matrix, sce)
seurat_meanDependence = getOverlapWithHighlyExpressed(seurat_matrix, sce)

M3Drop_meanDependence
seurat_meanDependence


### Calculate CH and DB indexes
# Get data from SingleCellExperiment for the HVGs
filtered <- sce[M3Drop_matrix, ]
# Run GLMPCA for dimensionality reduction
filtered <- GLMPCA(filtered, L=10, minibatch="stochastic")
# Build shared nearest-neighbor graph (SNNGraph)
g <- buildSNNGraph(filtered, k=10, use.dimred = 'GLMPCA')
# Perform Louvain clustering on the SNNGraph
clust <- cluster_louvain(g)
# Add cluster data as a factor to "filtered" variable
filtered$Louvain <- factor(membership(clust))

# Calculate Calinski-Harabasz index
M3Drop_CH_index = round(index.G1(t(counts(filtered)), as.integer(filtered$Louvain)), digits=2)
# Calculate Davies-Bouldin index
M3Drop_DB_index = round(index.DB(t(counts(filtered)), as.integer(filtered$Louvain))$DB, digits=2)

M3Drop_CH_index
M3Drop_DB_index
