#!/usr/bin/env Rscript
library(tximport)
library(readr)
library(DESeq2)
library(knitr)
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)
samples <- read.csv("/scratch/SampleDataFiles/Samples.csv", header=TRUE)
head(samples)
files <- file.path("quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ Menthol + Vibrio)
dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
    if(result != 'Intercept'){
        res <- results(dds, alpha=.05, name=result)
        dfRes <- as.data.frame(res)
        dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
        dfRes$Factor <- result
        dfAll <- rbind(dfAll, dfRes)
    }
}
head(dfAll)
write.csv(dfAll,file="dfAll.csv")
dfall<-read.csv("dfAll.csv",header=T)
colnames(dfall)<-c("ko","log2FoldChange","padj","Factor")
dfall<-subset(dfall,dfall$padj<.05)


#Getting the path and description from the path.txtand ko files
path <- read.table("/scratch/SampleDataFiles/Annotation/path.txt", sep="\t")
colnames(path) <- c("ko","pathway")
descript<-read.table("/scratch/SampleDataFiles/Annotation/ko",sep="\t")
colnames(descript)<- c("pathway","description")
#merging of the pathways
pathko <- merge(dfall,path)
pathdesc <- merge (pathko,descript)
#writing in to the csv file
write.csv(pathdesc,file="deAnnotated.csv")

