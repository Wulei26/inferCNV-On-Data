#!/public/home/wulei/miniconda3/envs/seurat/bin/Rscript
library(infercnv)
library(Seurat)
library(qs)
library(dplyr)
library(stringr)
crct <- qread("/public2022/wulei/research/counts/tumor/CRC.qs")
crct <- crct[, crct@meta.data$dataset == 8]

epit <- crct[, crct@meta.data$major == "Epithelial"]
tumor_sub <- subset(epit, downsample = 1500)
anno1 <- data.frame(Idents(tumor_sub))
colnames(anno1) <- "subtype"


reft <- crct[, crct@meta.data$major %in% c("T/NK","Myeloid")]
Idents(reft) <- "major"
ref_sub <- subset(reft, downsample = 5000)

anno2 <- data.frame(Idents(ref_sub))
colnames(anno2) <- "subtype"
anno <- rbind(anno1, anno2)
write.table(anno, file = "/public2022/wulei/research/InferCNV/CRC/reftnk/anno.txt",
                 sep = "\t", row.names = TRUE, col.names = FALSE)

scemerge <- merge(x = tumor_sub, y = ref_sub)
counts <- GetAssayData(scemerge, slot = 'counts')

gene_order <- read.table("/public2022/wulei/research/InferCNV/scripts/human_genes_pos_10X.txt", sep = "\t", header = FALSE)
counts <- counts[rownames(counts) %in% gene_order$V1, ]

anno <- "/public2022/wulei/research/InferCNV/CRC/reftnk/anno.txt"
gene_order <- "/public2022/wulei/research/InferCNV/scripts/human_genes_pos_10X.txt"

options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                        annotations_file = anno,
                                        delim="\t",
                                        gene_order_file = gene_order,
                                        min_max_counts_per_cell = c(100, +Inf),
                                        ref_group_names = c("T/NK", "Myeloid"))
setwd("/public2022/wulei/research/InferCNV/CRC/reftnk")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="/public2022/wulei/research/InferCNV/CRC/reftnk",
                             cluster_by_groups=TRUE, 
                             #analysis_mode="subclusters", #默认是"samples"
                             denoise=TRUE,
                             num_threads=48
                             )
