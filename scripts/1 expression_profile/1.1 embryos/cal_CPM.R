suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

rm(list = ls())  
options(stringsAsFactors = F)

#################################################### import data
filelist<-paste0("rawData/counts/",dir("rawData/counts/"))
file1.count<-read.table(filelist[1],header = T,sep = "\t")
count<-matrix(ncol = length(filelist)+1,nrow = nrow(file1.count))
count<-as.data.frame(count)
count[,1:2]<-file1.count[,c(1,3)]
for(i in 2:length(filelist)){
  file.count<-read.table(filelist[i],header = T,sep = "\t")
  count[,i+1]<-file.count[,3]
}
colnames(count) <- c("Geneid",filelist)
colnames(count) <- gsub(pattern="rawData/counts/",replacement = "",colnames(count))
colnames(count) <- gsub(pattern=".repeats.count.simple",replacement = "",colnames(count))
colnames(count) <- gsub(pattern=".repeats.Aligned.sortedByCoord.out",replacement = "",colnames(count))
rownames(count) <- count$Geneid
count <- count[,-1]

#################################################### PCA
library(xlsx)
pheno <- read.xlsx("rawData/pheno.xlsx",sheetIndex=1)
counts <- count[,colnames(count) %in% pheno$SRR]
rownames(pheno) <- pheno$SRR
pheno <- pheno[colnames(counts),]
colnames(counts) <- pheno$name

colData <- data.frame(stage=factor(pheno$stage,
                                   levels = c("oocytes","pronuclei","zygote","2-cell","4-cell","8-cell",
                                              "morula","blastocyst ","hESC.passage#0","hESC.passage#10")),
                      sequence=as.factor(pheno$sequencing),
                      project=as.factor(pheno$project),
                      row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ stage) 
dds <- DESeq(dds)
cpm <- fpm(dds,robust = FALSE)
cpm.raw <- as.data.frame(cpm)

chosen.cpm <- as.data.frame(t(cpm.raw))
chosen.cpm$name <- rownames(chosen.cpm)
chosen <- merge(pheno[,c("name","stage")],chosen.cpm,by="name")

#################################################### get mean of replicate samples in the same stage.
chosen.mean <-  aggregate(chosen[,3:ncol(chosen)],list(chosen$stage),mean)
rownames(chosen.mean) <- chosen.mean$Group.1
chosen.mean <- chosen.mean[,-1]
chosen.mean <- as.data.frame(chosen.mean)
chosen.mean.t <- t(chosen.mean)
data <- chosen.mean.t[,c(8,9,10,1,2,3,7,4,5,6)]
write.table(data, file = "repeats_CPM_embryogenesis.txt",sep="\t")