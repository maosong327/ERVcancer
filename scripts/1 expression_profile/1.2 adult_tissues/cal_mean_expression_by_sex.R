library(MultiAssayExperiment)
library(data.table)

# download GTEx (Genotype-Tissue Expression) data from previous paper:
# Bogu, G.K., Reverter, F., Marti-Renom, M.A., Snyder, M.P., Guigó, R., 2019. Atlas of  # transcriptionally active transposable elements in human adult tissues. bioRxiv, 
# 714212. https://doi.org/10.1101/714212. 
#import GTEx data
exp <- readRDS("111222_tes_8051_samples.matrix.main.ann.backup.rds")

#caculate HERVs (human endogenous retroviral sequences) expression for female or #male respectively
# for example, for female
#for a specific ERV subfamily, there is many copies in the genome, 
# so their expression value were summarized, 
# the median of all copies was kept, which belonging to the same HERV subfamily 
female.exp <- exp[exp$Gender_name %in% "Female",]
female.exp.median <- aggregate(female.exp[,3:ncol(female.exp)], list(female.exp $repName), median)

#caculate HERV subfamily expression for a type of tissue, 
# for example, there are many sigmoid colon tissues for female
# so their expression value were summarized,
# the mean of all replicate tissue samples was kept

female.exp.mean <- aggregate(female.exp[,2:ncol(female.exp)], list(female.exp$Tissue), mean) 

