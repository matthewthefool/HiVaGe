library(scRNAseq)
library(ensembldb)
library(scater)

#loading dataset
sce.lawlor =  LawlorPancreasData()



#annotation
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(edb, keys=rownames(sce.lawlor), keytype="GENEID", 
               columns=c("SYMBOL", "SEQNAME"))
rowData(sce.lawlor) <- anno[match(rownames(sce.lawlor), anno[,1]),-1]

#unfiltered data
unfiltered <- sce.lawlor

#qc and filtering
stats <- perCellQCMetrics(sce.lawlor, 
                          subsets=list(Mito=which(rowData(sce.lawlor)$SEQNAME=="MT")))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch=sce.lawlor$`islet unos id`)
sce.lawlor <- sce.lawlor[,!qc$discard]

#Remove transcripts with low total counts
total_counts_threshold <- 10  # Set your desired total counts threshold here, originally was 10

#Calculate total counts for each transcript
transcript_total_counts <- rowSums(counts(sce.lawlor))

# Subset transcripts with total counts above the threshold
sce <- sce.lawlor[transcript_total_counts >= total_counts_threshold, ]


# colData(unfiltered) <- cbind(colData(unfiltered), stats)
# unfiltered$discard <- qc$discard
# 
# gridExtra::grid.arrange(
#   plotColData(unfiltered, x="islet unos id", y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count") +
#     theme(axis.text.x = element_text(angle = 90)),
#   plotColData(unfiltered, x="islet unos id", y="detected", 
#               colour_by="discard") + scale_y_log10() + ggtitle("Detected features") +
#     theme(axis.text.x = element_text(angle = 90)), 
#   plotColData(unfiltered, x="islet unos id", y="subsets_Mito_percent",
#               colour_by="discard") + ggtitle("Mito percent") +
#     theme(axis.text.x = element_text(angle = 90)),
#   ncol=2
# )
# 
# 
# plotColData(unfiltered, x="sum", y="subsets_Mito_percent",
#             colour_by="discard") + scale_x_log10()
