# data frame with rows as cell lines and columns as repeat names
a <- read.csv("all_repeats_cpm_in_all_ccle_cancer_cell_lines.csv")
rownames(a) <- a$repName 
a <- a[,-1] 
a1 <- t(a) 
a1 <- as.data.frame(a1)
a1$cell <- rownames(a1)

anno <- read.table("cells_annotation.txt",header=TRUE,sep="\t")
a1.anno <- merge(anno,a1,by="cell")

write.table(a1.anno, file="repeats_CPM_in_cancer_cells.txt",sep="\t")
