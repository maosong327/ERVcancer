suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))

rm(list = ls())  
options(stringsAsFactors = F)

####################################################### load counts 
filelist <- paste0("rawData/counts/",dir("rawData/counts/"))
file1.count<-read.table(filelist[1],header = T,sep = "\t")
raw.count<-matrix(ncol = length(filelist)+1,nrow = nrow(file1.count))
raw.count<-as.data.frame(raw.count)
raw.count[,1:(ncol(file1.count)-1)]<-file1.count[,c(1,3)]
for(i in 2:length(filelist)){
  file.count<-read.table(filelist[i],header =T,sep = "\t")
  raw.count[,i+1]<-file.count[,3]
}
colnames(raw.count) <- c("Geneid",filelist)
colnames(raw.count) <- gsub(pattern="rawData/counts/",replacement = "",colnames(raw.count))
colnames(raw.count) <- gsub(pattern=".count.brief.txt",replacement = "",colnames(raw.count))
colnames(raw.count) <- gsub(pattern=".counts.simple",replacement = "",colnames(raw.count))
colnames(raw.count) <- gsub(pattern=".counts.brief.txt",replacement = "",colnames(raw.count))
rownames(raw.count) <- raw.count$Geneid
raw.count <- raw.count[,-1]

####################################################### input phenotype information for each tissue sample  
raw.pheno <- read.table(file="rawData/bam_files_submitter_id.txt",header = TRUE,sep = "\t")
raw.ffpe <- read.table(file = "rawData/bam_files_is_ffpe.txt",header = TRUE,sep = "\t")
raw.type <- read.table(file="rawData/all_samples_and_its_type.txt",header = TRUE)
raw.pheno.merge1 <- merge(raw.pheno,raw.ffpe,by="bam_file_name")
raw.pheno.merge2 <- merge(raw.pheno.merge1,raw.type,by="bam_file_name")
raw.pheno.merge2 <- raw.pheno.merge2[,c("file_id.x","bam_file_name","submitter_id","group","is_ffpe","type")]
raw.pheno.merge2$type.group <- paste(raw.pheno.merge2$type,raw.pheno.merge2$group,sep=".")

####################################################### delete samples that are not normal solid or primary tumor tissues.
filter.pheno <- raw.pheno.merge2[raw.pheno.merge2$group %in% c("solid.tissue.normal","primary.tumor"),]

####################################################### delete FFPE(Formalin Fixed Paraffin-Embedded) samples.
filter.pheno <- filter.pheno[filter.pheno$is_ffpe =="FALSE",]
rownames(filter.pheno) <- filter.pheno$file_id.x

####################################################### the final count matrix of samples filtered.
filter.count <- raw.count[,colnames(raw.count) %in% filter.pheno$file_id.x]
final.pheno <- filter.pheno[colnames(filter.count),]
table(rownames(final.pheno)==colnames(filter.count))  
colnames(filter.count) <- final.pheno$submitter_id

####################################################### start differential expression analysis 
colData <- data.frame(group=final.pheno$group,row.names = final.pheno$submitter_id,type=final.pheno$type,type.group=final.pheno$type.group)
colData$patient <- substr(rownames(colData),1,12)
colData$type.group <- factor(colData$type.group)

# Not all the counts in the expression matrix is interger because of the "--fraction" used in the featureCounts quantification step. 
# So the round function was used.
dds <- DESeqDataSetFromMatrix(countData = round(filter.count), colData = colData,design = ~ type.group ) 
dds <- DESeq(dds)

all_types <- as.data.frame(unique(raw.type$type))
colnames(all_types) <- "type"
all_types$normal <- paste(all_types$type,".primary.tumor",sep="")
all_types$tumor <- paste(all_types$type,".solid.tissue.normal",sep="")
all_types$order <- rownames(all_types)

####################################################### write a loop for differential expressed repeats.
for (i in 1:nrow(all_types))
{
res <- results(dds,contrast = c("type.group",all_types[i,2],all_types[i,3]))
res <- as.data.frame(res)
res$repeatMasker.repName <- rownames(res)
res$cancer.type <- all_types[i,1]
write.table(res,file=paste(all_types[i,1],"DEG_res.txt",sep="."),sep = "\t")
}

#merge all the differential expression result.
#all.res <- rbind(res,res4,res6,res8,res9,res11,res12)
#write.table(all.res[,c(8,7,1,2,5,6)],file="all_cancer_types_res.txt",sep="\t",row.names=F)

####################################################### caculate CPM(Counts Per Million) for each tissue sample in each patient
res.cpm <- fpm(dds,robust = FALSE)
res.cpm <- as.data.frame(res.cpm)
write.table(res.cpm,file="cpm.txt",sep="\t")