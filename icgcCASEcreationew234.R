
setwd("/home/vant/tfmdir/signats/citvalidation/icgcCA/")
load("newtransicgcCA.RData")
library(DGEobj.utils)

icgcnew <- read.csv("icgc_ca_raw0.csv", stringsAsFactors = FALSE, row.names = 1)

icgcnewtpm <- icgcnew[rownames(genGClenCAED),]

icgcnewtpm <- as.matrix(icgcnewtpm)
tpmnorm4 <- convertCounts(countsMatrix = icgcnewtpm, unit = "TPM", geneLength = genGClenCAED[,1])

all(rownames(genGClenCAED)==rownames(icgcnewtpm))

coldatas <- DataFrame(donorED)

all(rownames(tpmnorm4)==rownames(setpmXCA))
seicgcnew234 <- SummarizedExperiment(assays = list(counts=icgcnewtpm), colData = coldatas)
rowData(seicgcnew234) <- rowData(setpmXCA)

assays(seicgcnew234)$tpmnorm <- tpmnorm4
assays(seicgcnew234)$logtpm <- log2(tpmnorm4 + 1)

save(seicgcnew234, file = "icgcnewtpm.RData")

write.csv(icgcnewtpm,"icgcnewtpm.csv")
write.csv(tpmnorm4,"tpmicgcnew.csv")

