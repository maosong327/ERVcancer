#collect genes expression matrix from the UCSC Xena platform (https://xena.ucsc.edu/)
# reference:
# Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting
# cancer genomics data via the Xena platform. Nat Biotechnol (2020).
# https://doi.org/10.1038/s41587-020-0546-8
#1、input genes expression matrix.
genes.raw <- read.table("rawData/tcga_RSEM_Hugo_norm_count",header=TRUE,sep="\t")
genes.raw <- as.data.frame(genes.raw)
rownames(genes.raw) <- genes.raw$sample
genes.raw <- genes.raw[,-1]
genes.raw <- as.data.frame(genes.raw)
genes.raw.log2 <- log2(genes.raw+1)
genes.raw.log2 <- as.data.frame(genes.raw.log2)
#transpose
genes.raw.log2.t <- t(genes.raw.log2)
genes.raw.log2.t <- as.data.frame(genes.raw.log2.t)
genes.raw.log2.t$sample <- rownames(genes.raw.log2.t)
#2、input human endogenous retroviral sequences (HERVs) expression matrix
erv.raw <- read.table("rawData/TCGA_tissue_samples_ERV_CPM.txt",header=TRUE,sep="\t")
#include only tumor tissue samples.
cancer <- erv.raw[ erv.raw$group %in% "primary.tumor",-c(1,3)]
# for patients with multiple cancer tissue sample, their expression were merged by mean
cancer$sample <- substr(cancer$submitter_id,1,15)
cancer$sample <- gsub("-",".",cancer$sample)
cancer.mean <- aggregate(cancer[,2:(ncol(cancer)-1)],list(cancer$sample),mean)
colnames(cancer.mean)[1] <- "sample"
rownames(cancer.mean) <- cancer.mean$sample
cancer.mean <- cancer.mean[,-1]
cancer.mean.log2 <- log2(cancer.mean+1)
cancer.mean.log2$sample <- rownames(cancer.mean.log2)
cancer.mean.log2 <- as.data.frame(cancer.mean.log2)
#3、merge genes and HERVs expression matrix
erv.genes <- merge(genes.raw.log2.t, cancer.mean.log2,by="sample")
rownames(erv.genes) <- erv.genes$sample
erv.genes <- erv.genes[,-1]
# when include all ~60,000 genes and ERV subfamilies,
# the R function “rcorr” cannot calculate so many genes or repeats
# so only include genes annotated by both UCSC genes annotation file and Ensemble
# genes annotation file
# there are ~24,000 genes annotated by both annotation files.
# input genes list
chosen.gene <- read.table("rawData/anno/share.txt",header=F,sep="\t")
chosen.gene$V1 <- gsub(";","",chosen.gene$V1)
chosen.erv.genes <- erv.genes[,colnames(erv.genes) %in% c(chosen.gene$V1,colnames(erv.raw))]
chosen.erv.genes <- as.data.frame(chosen.erv.genes)
chosen.erv.genes$sample <- rownames(chosen.erv.genes)
#4、correlation analysis and output result
library(Hmisc)
cor.res <- rcorr(as.matrix(chosen.erv.genes,method="pearson"))
res.p.raw <- as.data.frame(cor.res$P)
res.r.raw <- as.data.frame(cor.res$r)
#only include result of genes and HERVs,
#remove genes and genes, HERVs and HERVs
res.p <- res.p.raw[rownames(res.p.raw) %in% colnames(erv.raw),colnames(res.p.raw) %in% chosen.gene$V1]
write.table(as.data.frame(res.p),file=paste(chosen,"P_value.txt",sep="_"),sep="\t")
# only include result of genes and HERVs,
# remove genes and genes, HERVs and HERVs
res.r <- res.r.raw[rownames(res.r.raw) %in% colnames(erv.raw),colnames(res.p.raw) %in% chosen.gene$V1]
write.table(as.data.frame(res.r),file=paste(chosen,"r_value.txt",sep="_"),sep="\t")
# 5、output raw data, i.e. expression value for correlation analysis
anno <- erv.raw[,c(1,2)]
anno$sample <- substr(anno$submitter_id,1,15)
anno$sample <- gsub("-",".",anno$sample)
data <- merge(anno[,c("sample","cancer.type")],chosen.erv.genes,by="sample")
write.table(as.data.frame(data[,c(2,1,3:ncol(data))]),row.names=F,file="log2_cpm_plus_1.txt",sep="\t")