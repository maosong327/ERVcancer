#1、import tissues sample CPM data
raw <- read.table("TCGA_tissue_samples_ERV_CPM.txt",header=TRUE,sep="\t")
#2、only tumor tissues were included
cancer <- raw[raw$cancer.type %in% chosen & raw$group %in% "primary.tumor",c(-1,-3)]
rownames(cancer) <- cancer$submitter_id
cancer <- cancer[,-1]
library(Hmisc)
#3、correlation analysis using pearson method based on log2(CPM+1)
cor.res <- rcorr(as.matrix(log2(cancer+1),method="pearson"))
#4、output correlation analysis result
write.table(as.data.frame(cor.res$P),file=paste(chosen,"P_value.txt",sep="_"),sep="\t")
write.table(as.data.frame(cor.res$r),file=paste(chosen,"r_value.txt",sep="_"),sep="\t")
#5、output raw data, i.e. expression value (log2(CPM+1)) for correlation analysis
raw <- read.table("TCGA_tissue_samples_ERV_CPM.txt",header=TRUE,sep="\t")
raw$log2_cpm_plus_1 <- log2(raw$CPM+1)
write.table(raw[,-5],file="TCGA_tumor_tissue_ERV_log2_CPM_plus_1.txt",sep="\t",row.names=F)