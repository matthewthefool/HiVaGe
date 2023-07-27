
# Set-up libraries
library(scRNAseq)
library(scater)
library(scDblFinder)
library(scry)
library(SingleR)
library(celldex)
library(DropletUtils)

library(Rtsne)

# Define functions
getPurity = function(clusters, classes) {
  sum(apply(table(classes, clusters), 1, max))/length(clusters)
}


### Data preparation
# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

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

# GLMPCA clustering
# Perform GLMPCA clustering on 'sce' with L=10 and minibatch="stochastic"
sce <- GLMPCA(sce, L = 10, minibatch = "stochastic")

# tSNE visualization
# Run tSNE on 'sce' using the GLMPCA dimensions
sce <- runTSNE(sce, dimred = "GLMPCA")

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Assign row names to 'sce' using the gene symbols
rownames(sce) <- rowData(sce)$Symbol

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)



### Purity testing
# Retrieve relevant data
sce_logcounts = logcounts(sce)

# Fix rownames (duplicates and use of underscores)
new_rownames = gsub("-", "_", rownames(sce_logcounts))
new_rownames = make.names(new_rownames, unique=TRUE)
rownames(sce_logcounts) = gsub("_", "-", new_rownames)

# Retrieve HVGs (seurat)
data <- read.csv("seurat_matrix.csv", header = FALSE)
seurat_matrix <- as.matrix(data)
seurat_matrix <- seurat_matrix[2:31, ]

# According to the paper, using PCs 2 to 6 is best ¯\_(ツ)_/¯
# and i don't know what else it could mean, but these components

# Vector for storing purity for each clustering
purities = c(1:5)

for (i in 1:5) {
  # Getting expression data for genes of (i+1)th PC 
  sce_logcounts_genes = as.data.frame(t(sce_logcounts[seurat_matrix[,i+1],]))
  # Adding cell-type as separate column
  sce_logcounts_genes$CellType = sce$CellTypes
  # Removing duplicates (ignoring the CellType column)
  sce_logcounts_genes = sce_logcounts_genes[!duplicated(subset(sce_logcounts_genes, select=-c(CellType))), ]
  
  # Running t-SNE (with the library from the paper) on expression data
  sce_tSNE = Rtsne(subset(sce_logcounts_genes, select=-c(CellType)), perplexity=50, max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000)
  # K-means clustering. I chose 28 clusters, since that's the amount of predicted
  # cell types. If that's not how that should be done, pls change
  sce_cluster <- kmeans(sce_tSNE$Y, center=28, nstart=20)
  
  # Calculating purity for (i+1)th clustering
  purities[i] = getPurity(sce_logcounts_genes$CellType, sce_cluster$cluster)
}

# Final purity score for this method (seurat) with this database (MairPBMCData)
purity_score = mean(purities)
purity_score
