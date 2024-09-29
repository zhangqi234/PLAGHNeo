library(Seurat)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(viridis)
library(copykat)
options(stringsAsFactors=FALSE)
load("/data5/zhangq/ThrombusST/data/scRNA/301_RCC_Thrombus_FirstAnnotation.rda")
setwd("/data5/zhangq/ThrombusST/data/scRNA/Tumor")
D301obj <- subset(RCC_SCT, subset = bcluster %in% c("Epi", "T") & nFeature_RNA >= 500)

listCountMtx <- list(); listNormCells <- list()
for(samplei in unique(as.character(D301obj$PatientID)))
{
    tmp <- subset(D301obj, subset = PatientID == samplei)
    listCountMtx[[samplei]] <- tmp@assays$RNA@counts
    listNormCells[[samplei]] <- rownames(tmp@meta.data)[which(tmp@meta.data$bcluster == "T")]
    print(samplei)
}

for(samplei in unique(as.character(D301obj$PatientID)))
{
    #samplei <- "LMD"
    exp.rawdata <- listCountMtx[[samplei]]
    norm_cell <- listNormCells[[samplei]]
    copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", cell.line="no",  ngene.chr=5, win.size=25, 
                        norm.cell.names=norm_cell, KS.cut=0.1,  sam.name=samplei, distance="euclidean",  n.cores=4)
    save(copykat.test, file=paste0(samplei, "_Filter_copykat.rda"))
    print(samplei)
}
