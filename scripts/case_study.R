library(tidyverse)

# pwd
# /home/maosong/project/xylei/EEPIA/1.single_ERV_analysis/5.correlation/cancer_tissues/1.TCGA.update
# erv_all <- read_tsv("1.correlated_ERVs/TCGA_tumor_tissue_ERV_log2_CPM_plus_1.txt")
# "cancer.type"   "tissue_sample_id"      "group" "repName"       "log2_cpm_plus_1"
# "UCEC"  "TCGA-EY-A1GI-01A-11R-A13S-07"  "primary tumor" "MSTA-int"      8.6976916150430

# 2.correlated_genes/2.log2_cpm_plus_1_uniq_merge.txt
# "cancer.type"   "patient_sample_id"     "repName_or_geneName"   "log2_cpm_plus_1"
# "BLCA"  "TCGA.2F.A9KO.01"       "A1BG"  2.68198890998603

# load data -----------------
all <- read_tsv("2.correlated_genes/2.log2_cpm_plus_1_uniq_merge.txt")
all <- all %>% 
  mutate(patient_sample_id = str_replace_all(patient_sample_id, "\\.", "-")) %>%
  mutate(repName_or_geneName = str_replace_all(repName_or_geneName, "\\.", "-"))
all_sysbols = all %>%
  distinct(repName_or_geneName) %>%
  pull(repName_or_geneName)

znf_genes <- read_tsv("2.correlated_genes/all_znf.txt") %>%
  distinct(assigned_gene, .keep_all = T) %>%
  pull(assigned_gene)
#467 znfs

znfs <- intersect(znf_genes, all_sysbols)
#391 znfs
data.frame(znfs)
write.csv(data.frame(znfs), "znf_analysis_391.csv")

ervs <- read_tsv("/home/maosong/project/xylei/EEPIA/1.single_ERV_analysis/1.general/ERV_subfamilies_names.txt",
                col_names = F) %>%
  pull(X1)
#580 ervs
hervhs <- c("HERVH-int", "LTR7", "LTR7A", "LTR7B", "LTR7C", "LTR7Y")

# call correlation -----------------
choose <- all %>%
  filter(repName_or_geneName %in% c(znfs, ervs))
ervs_znfs <- choose %>%
  pivot_wider(names_from=repName_or_geneName, values_from=log2_cpm_plus_1) %>%
  select(all_of(ervs), everything(), -cancer.type, -patient_sample_id)
ervs_znfs_cor <- cor(ervs_znfs)
ervs_znfs_cor <- ervs_znfs_cor[1:length(ervs),(length(ervs)+1):(length(ervs)+length(znfs))]
write.csv(ervs_znfs_cor, "ervs_znfs_cor.csv")

znfs_cor_ge03 <- apply(ervs_znfs_cor, 2, function(x) sum(x>=0.3))
znfs_cor_le03 <- apply(ervs_znfs_cor, 2, function(x) -sum(x <= -0.3))
top_cor_znf_plus <- head(names(znfs_cor_ge03)[order(znfs_cor_ge03, decreasing = T)], 10)
top_cor_znf_minus <- tail(names(znfs_cor_le03)[order(znfs_cor_le03, decreasing = T)], 10)
znf_cor_stat = cbind(as.data.frame(znfs_cor_ge03), as.data.frame(znfs_cor_le03))
pdf("Pancancer_KZFPs_cor_HERVs_barplot.pdf")
znf_cor_stat %>%
  as_tibble(rownames = "gene") %>%
  filter(rownames(znf_cor_stat) %in% c(top_cor_znf_plus, top_cor_znf_minus)) %>%
  pivot_longer(-gene, names_to = "type", values_to = "number") %>%
  mutate(gene = factor(gene, levels = rev(c(top_cor_znf_plus, top_cor_znf_minus)))) %>%
  # = fct_reorder(gene, number)
  ggplot(aes(x = gene, fill = type, y = number)) +
  geom_bar(stat = "identity") +
  # scale_y_continuous(labels = abs, limits = c(-200,200)) +
  labs(x = "", y = "Number of HERVs targeting a KZFP", fill = "",
       title = "Pan-cancer analysis of KZFP-related HERVs") +
  scale_fill_hue(labels=c("Positive correlation","Negative correlation")) +
  coord_flip() +
  theme_classic()
dev.off() 

ervs_cor_ge03 <- apply(ervs_znfs_cor, 1, function(x) sum(x>=0.3))
ervs_cor_le03 <- apply(ervs_znfs_cor, 1, function(x) -sum(x <= -0.3))
top_cor_ervs_plus <- head(names(ervs_cor_ge03)[order(ervs_cor_ge03, decreasing = T)], 10)
top_cor_ervs_minus <- tail(names(ervs_cor_le03)[order(ervs_cor_le03, decreasing = T)], 10)
ervs_cor_stat = cbind(as.data.frame(ervs_cor_ge03), as.data.frame(ervs_cor_le03))
pdf("Pancancer_HERVs_cor_KZFPs_barplot.pdf")
ervs_cor_stat %>%
  as_tibble(rownames = "gene") %>%
  filter(rownames(ervs_cor_stat) %in% c(top_cor_ervs_plus, top_cor_ervs_minus)) %>%
  pivot_longer(-gene, names_to = "type", values_to = "number") %>%
  mutate(gene = factor(gene, levels = rev(c(top_cor_ervs_plus, top_cor_ervs_minus)))) %>%
  # = fct_reorder(gene, number)
  ggplot(aes(x = gene, fill = type, y = number)) +
  geom_bar(stat = "identity") +
  # scale_y_continuous(labels = abs, limits = c(-200,200)) +
  labs(x = "", y = "Number of KZFPs targeting a subfamily", fill = "",
       title = "Pan-cancer analysis of HERV-related KZFPs") +
  scale_fill_hue(labels=c("Positive correlation","Negative correlation")) +
  coord_flip() +
  theme_classic()
dev.off() 

pdf("Pancancer_HERVHs_cor_KZFPs_barplot.pdf")
ervs_cor_stat[hervhs,] %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(-gene, names_to = "type", values_to = "number") %>%
  mutate(gene = factor(gene, levels = 
                         c("LTR7Y", "HERVH-int", "LTR7", "LTR7A", "LTR7B", "LTR7C"))) %>%
  ggplot(aes(x = gene, fill = type, y = number)) +
  geom_bar(stat = "identity") +
  # scale_y_continuous(labels = abs, limits = c(-200,200)) +
  labs(x = "", y = "Number of KZFPs targeting a subfamily", fill = "",
       title = "Pan-cancer analysis of HERV-related KZFPs") +
  scale_fill_hue(labels=c("Positive correlation","Negative correlation")) +
  coord_flip() +
  theme_classic()
dev.off() 

HERVH_int_cor <- ervs_znfs_cor["HERVH-int",]
topn = 20
pdf(paste0("Pancancer_HERVH_int_cor_KZFPs_barplot_top", topn ,".pdf"))
data.frame(HERVH_int_cor) %>%
  as_tibble(rownames = "gene") %>% 
  arrange(HERVH_int_cor) %>%
  filter(HERVH_int_cor <= -0.3) %>%
  head(topn) %>%
  ggplot(aes(x = fct_reorder(gene, HERVH_int_cor), y = HERVH_int_cor)) +
  geom_point(size = 4, shape = 21,
             fill = "blue", color = "black") +
  geom_segment(aes(xend = after_stat(x), yend = -0.3)) +
  labs(x = "", y = "Correlation between KZFP and HERVH-int expression", 
       title = "Pan-cancer analysis of KZFPs targeting HERVH-int") +
  coord_flip() +
  theme_classic()
dev.off()

LTR7Y_cor <- ervs_znfs_cor["LTR7Y",]
pdf(paste0("Pancancer_LTR7Y_cor_KZFPs_barplot_top", topn ,".pdf"))
data.frame(LTR7Y_cor) %>%
  as_tibble(rownames = "gene") %>% 
  arrange(LTR7Y_cor) %>%
  filter(LTR7Y_cor <= -0.3) %>%
  head(topn) %>%
  ggplot(aes(x = fct_reorder(gene, LTR7Y_cor), y = LTR7Y_cor)) +
  geom_point(size = 4, shape = 21,
             fill = "red", color = "black") +
  geom_segment(aes(xend = after_stat(x), yend = -0.3)) +
  labs(x = "", y = "Correlation between KZFP and LTR7Y expression", 
       title = "Pan-cancer analysis of KZFPs targeting LTR7Y") +
  coord_flip() +
  theme_classic()
dev.off()

LTR7_cor <- ervs_znfs_cor["LTR7",]
pdf(paste0("Pancancer_LTR7_cor_KZFPs_barplot_top", topn ,".pdf"))
data.frame(LTR7_cor) %>%
  as_tibble(rownames = "gene") %>% 
  arrange(LTR7_cor) %>%
  filter(LTR7_cor <= -0.3) %>%
  head(topn) %>%
  ggplot(aes(x = fct_reorder(gene, LTR7_cor), y = LTR7_cor)) +
  geom_point(size = 4, shape = 21,
             fill = "green", color = "black") +
  geom_segment(aes(xend = after_stat(x), yend = -0.3)) +
  labs(x = "", y = "Correlation between KZFP and LTR7_cor expression", 
       title = "Pan-cancer analysis of KZFPs targeting LTR7_cor") +
  coord_flip() +
  theme_classic()
dev.off()

#
library(eulerr)
# vd <- euler(c(A = 0.3, B = 0.3, C = 1.1,
#               "A&B" = 0, "A&C" = 0.2, "B&C" = 0.1,
#               "A&B&C" = 0.1))
a = names(LTR7Y_cor)[LTR7Y_cor < -0.3]
b = names(HERVH_int_cor)[HERVH_int_cor < -0.3]
c = names(LTR7_cor)[LTR7_cor < -0.3]

vd <- euler(c("LTR7Y" = length(setdiff(setdiff(a, b), c)), 
              "HERVH-int" = length(setdiff(b, a)),
              "LTR7" = length(setdiff(c, a)),
              "LTR7Y&HERVH-int" = length(intersect(a,b)),
              "LTR7Y&LTR7" = length(intersect(a,c)))
            )

pdf("Pancancer_LTR7Y_HERVH_int_intersect_veen.pdf")
plot(vd,
     fills = list(fill = c("red", "blue", "green"), alpha = 0.6),
     labels = list(col = "white", font = 4), 
     edges = FALSE,
     quantities = TRUE)
dev.off()

# group by cancer type ------------------------
znfs_intersect <- intersect(a, b)
choose2 <- choose %>%
  filter(repName_or_geneName %in% c(znfs_intersect, "HERVH-int", "LTR7Y")) %>%
  pivot_wider(names_from = "repName_or_geneName", 
              values_from = "log2_cpm_plus_1") %>%
  # nest_by(cancer.type)
  group_by(cancer.type) %>%
  split(group_indices(.))

choose %>%
  filter(repName_or_geneName %in% c(znfs_intersect, "HERVH-int", "LTR7Y")) %>%
  pivot_wider(names_from = "repName_or_geneName", 
              values_from = "log2_cpm_plus_1") %>%
  select(-`cancer.type`, -patient_sample_id, -`HERVH-int`) %>%
  cor()  %>%
  as.data.frame() %>%
  slice(-1) %>%
  pull(LTR7Y)

choose %>%
  filter(repName_or_geneName %in% c(znfs_intersect, "HERVH-int", "LTR7Y")) %>%
  pivot_wider(names_from = "repName_or_geneName", 
              values_from = "log2_cpm_plus_1") %>%
  select(-`cancer.type`, -patient_sample_id, -LTR7Y) %>%
  cor()  %>%
  as.data.frame() %>%
  slice(-1) %>%
  pull(`HERVH-int`)

choose2[[1]] %>%
  ungroup %>%
  select(-`cancer.type`, -patient_sample_id, -`HERVH-int`) %>%
  cor() %>%
  as.data.frame() 

cal_cor_LTR7Y <- function(tb){
  tb %>%
    ungroup %>%
    select(-`cancer.type`, -patient_sample_id, -`HERVH-int`) %>%
    cor() %>%
    as.data.frame() %>%
    slice(-1) %>%
    pull(LTR7Y)
}

LTR7_cor_df <- map_dfc(choose2, cal_cor_LTR7Y) %>%
  set_names(unique(choose$cancer.type)) %>%
  mutate(gene = znfs_intersect) %>%
  column_to_rownames(var = "gene")

pdf("LTR7_cor_by_cancer_heatmap.pdf")
pheatmap::pheatmap(LTR7_cor_df,
                   scale = "none", 
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize_number = 10,
                   colorRampPalette(colors = c("red", "yellow","white"))(100)
                   )
dev.off()

cal_cor_HERVH_int <- function(tb){
  tb %>%
    ungroup %>%
    select(-`cancer.type`, -patient_sample_id, -`LTR7Y`) %>%
    cor() %>%
    as.data.frame() %>%
    slice(-1) %>%
    pull(`HERVH-int`)
}

HERVH_int_cor_df <- map_dfc(choose2, cal_cor_HERVH_int) %>%
  set_names(unique(choose$cancer.type)) %>%
  mutate(gene = znfs_intersect) %>%
  column_to_rownames(var = "gene")

pdf("HERVH_int_cor_by_cancer_heatmap.pdf")
pheatmap::pheatmap(HERVH_int_cor_df,
                   scale = "none", 
                   display_numbers = T,
                   number_format = "%.2f",
                   fontsize_number = 10,
                   colorRampPalette(colors = c("red", "yellow","white"))(100)
)
dev.off()

# data for cytoscape ---------------------------
ervs_znfs_cor[, znfs_intersect] %>%
  as_tibble(rownames = "HERVs") %>%
  pivot_longer(cols = -HERVs, names_to = "KZFP", values_to = "cor") %>%
  filter(abs(cor) >= 0.3) %>%
  select(2,1,3) %>%
  set_names(c("source", "target", "weight")) %>%
  write_csv("ervs_znfs_network.csv")

bind_rows(tibble(Id = ervs, Lable = ervs, Type = "HERV"),
                     tibble(Id = znfs_intersect, Lable = znfs_intersect, Type = "KZFP")
) %>%
  write_csv("dot_type.csv")

ervs_znfs_cor[, znfs_intersect] %>%
  as_tibble(rownames = "HERVs") %>%
  pivot_longer(cols = -HERVs, names_to = "KZFP", values_to = "cor") %>%
  filter(cor <= -0.3) %>%
  select(2,1,3) %>%
  set_names(c("source", "target", "weight")) %>%
  write_csv("ervs_znfs_network_minus.csv")