library(scRNAseq)
Richard_ds <- RichardTCellData(location = TRUE)

# Making a plot that shows that there's clear distinction in cells that have big UMIs and low UMIs. Inflection is 100 UMIs so we could use default threshold. Big difference in genes expression we consider normal, as we have PBMC cells, most of which are specific to some function and don't express big amount of genes.
bcrank <- barcodeRanks(counts(Richard_ds))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)




allzero <- rowMeans(counts(Richard_ds) == 0) == 1
Richard_ds <- Richard_ds[which(!allzero), ]

Richard_ds <- logNormCounts(Richard_ds)
Richard_ds <- runPCA(Richard_ds)
library(scDblFinder)
dbl.dens <- computeDoubletDensity(Richard_ds, d=ncol(reducedDim(Richard_ds)))
Richard_ds$DoubletScore <- dbl.dens
dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
Richard_ds <- Richard_ds[, dbl.calls == "singlet"]





library(ensembldb)
library(EnsDb.Mmusculus.v75)

ens.mm.v75 <- EnsDb.Mmusculus.v75


#annotation
rowData(Richard_ds)$location <- mapIds(ens.mm.v75,
                                    keys=rownames(rowData(Richard_ds)), 
                                    column="SEQNAME", keytype="GENEID")
stats <- perCellQCMetrics(Richard_ds, subsets=list(Mito=which(rowData(Richard_ds)$location=="MT")), use.altexps=TRUE)
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
Richard_ds <- Richard_ds[,!high.mito]

# Discarding cells using perCellQCMetrics
df <- perCellQCMetrics(Richard_ds)
reasons <- perCellQCFilters(stats, sub.fields = TRUE)
Richard_ds <- Richard_ds[,!reasons$discard]


