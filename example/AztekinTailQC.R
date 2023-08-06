library(scRNAseq)
library(ensembldb)
library(scater)
library(xlaevis.db)
library(org.Xl.eg.db)

#loading dataset
sce.aztekin <- AztekinTailData()

unfiltered <- sce.aztekin
#qc

ah <- AnnotationHub()

annotation = ah[['AH111581']]


# qcstats <- perCellQCMetrics(sce.aztekin)
# filtered <- quickPerCellQC(qcstats)


# sce.aztekin <- sce.aztekin[, !filtered$discard]
annotated_dataset <- select(annotation, keys = rownames(sce.aztekin), keytype = "SYMBOL", column = 'GENENAME')

is.mito <- grepl("mitoch", annotated_dataset$GENENAME, ignore.case = TRUE)

stats <- perCellQCMetrics(sce.aztekin, subsets=list(Mito=is.mito))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.aztekin <- sce.aztekin[,!high.mito]


#Remove transcripts with low total counts
total_counts_threshold <- 10  # Set your desired total counts threshold here, originally was 10

#Calculate total counts for each transcript
transcript_total_counts <- rowSums(counts(sce.aztekin))

# Subset transcripts with total counts above the threshold
sce <- sce.aztekin[transcript_total_counts >= total_counts_threshold, ]


# 
# colData(unfiltered) <- cbind(colData(unfiltered), stats)
# unfiltered$discard <- high.mito
# 
# gridExtra::grid.arrange(
#   plotColData(unfiltered, y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count"),
#   plotColData(unfiltered, y="detected", colour_by="discard") +
#     scale_y_log10() + ggtitle("Detected features"),
#   plotColData(unfiltered, y="subsets_Mito_percent",
#               colour_by="discard") + ggtitle("Mito percent"),
#   ncol=2
# )



