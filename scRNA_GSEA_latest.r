# 为每个细胞类型运行GSEA
# 为每个细胞类型运行GSEA
library(Seurat)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(fgsea)
load("F:/workspace/CDH11/seurat_stro_rename.rda")
load("F:/workspace/CDH11/seurat_T_rename.rda")
# 确认细胞类型注释
Idents(seurat_stro_rename) <- "celltype"
table(Idents(seurat_stro_rename))

cell_types <- unique(Idents(seurat_stro_rename))
gsea_by_celltype <- list()
# 从MSigDB获取基因集
# 可以选择不同的基因集类别
msigdbr_species <- "Homo sapiens"  # 根据物种调整

# 获取Hallmark基因集
hallmark_sets <- msigdbr(species = msigdbr_species, 
                         category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# 转换为GSEA需要的列表格式
hallmark_list <- split(hallmark_sets$gene_symbol, 
                       hallmark_sets$gs_name)

# 首先查看所有可用的子集合
available_subcollections <- msigdbr_collections()
print(available_subcollections)

# 或者查看C2类别下的子集合
c2_subcollections <- msigdbr_collections() %>%
  filter(collection == "C2")
print(c2_subcollections)
# 获取KEGG通路
kegg_sets <- msigdbr(species = msigdbr_species,
                     category = "C2", 
                     subcategory = "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)
kegg_list <- split(kegg_sets$gene_symbol, kegg_sets$gs_name)
for(celltype in cell_types) {
  # 提取该细胞类型与其他所有细胞的差异表达
  de_celltype <- FindMarkers(seurat_stro_rename,
                           ident.1 = celltype,
                           ident.2 = NULL,
                           min.pct = 0.1)
  
  # 准备排序基因列表
  gene_rank_celltype <- de_celltype$avg_log2FC
  names(gene_rank_celltype) <- rownames(de_celltype)
  gene_rank_celltype <- sort(gene_rank_celltype, decreasing = TRUE)
  
  # 运行GSEA
  gsea_celltype <- fgsea(pathways = hallmark_list,
                        stats = gene_rank_celltype,
                        minSize = 15,
                        maxSize = 500)
  
  gsea_by_celltype[[celltype]] <- gsea_celltype
}

# 比较不同细胞类型的富集结果
compare_gsea_results <- function(gsea_list, pathway) {
  results <- data.frame()
  for(celltype in names(gsea_list)) {
    pathway_result <- gsea_list[[celltype]][pathway == pathway, ]
    if(nrow(pathway_result) > 0) {
      results <- rbind(results, 
                      data.frame(celltype = celltype,
                                NES = pathway_result$NES,
                                padj = pathway_result$padj))
    }
  }
  return(results)
}

# 示例：比较炎症反应通路
inflammatory_comparison <- compare_gsea_results(gsea_by_celltype, 
                                               "HALLMARK_INFLAMMATORY_RESPONSE")
# 提取所有细胞类型的GSEA结果并合并
# prepare_comparison_data <-function(gsea_list, top_n = 10) {
#   all_results <- data.frame()

#   for(celltype in names(gsea_list)[c(1,2,4,5)]) {
#     gsea_df <- gsea_list[[celltype]]
#     if(nrow(gsea_df) > 0) {
#       # 选择每个细胞类型中NES最高和最低的top_n个通路
#       top_pos <- gsea_df[gsea_df$NES > 0, ][1:min(top_n, sum(gsea_df$NES > 0)), ]
#       top_neg <- gsea_df[gsea_df$NES < 0, ][1:min(top_n, sum(gsea_df$NES < 0)), ]

#       selected <- rbind(top_pos, top_neg)
#       selected$celltype <- celltype
#       all_results <- rbind(all_results, selected)
#     }
#   }

#   # 移除重复和NA值
#   all_results <- na.omit(all_results)
#   return(all_results)
# }
prepare_comparison_data <- function(gsea_list, top_n = 10) {
  all_results <- data.frame()
  
  for(celltype in names(gsea_list)) {
    gsea_df <- gsea_list[[celltype]]
    if(nrow(gsea_df) > 0) {
      gsea_df_sort=gsea_df[order(gsea_df$padj),]
      # 选择每个细胞类型中NES最高和最低的top_n个通路
      top_pos <- gsea_df_sort[gsea_df_sort$NES > 0, ][1:min(top_n, sum(gsea_df_sort$NES > 0)), ]
    #  top_neg <- gsea_df_sort[gsea_df_sort$NES < 0, ][1:min(top_n, sum(gsea_df_sort$NES < 0)), ]
      
      selected <- top_pos
      selected$celltype <- celltype
      all_results <- rbind(all_results, selected)
    }
  }
  
  # 移除重复和NA值
  all_results <- na.omit(all_results)
  return(all_results)
}

# 准备比较数据
comparison_data <- prepare_comparison_data(gsea_by_celltype, top_n = 5)
# 气泡图：显示NES和显著性
library(ggplot2)

bubble_plot <- ggplot(comparison_data, 
                     aes(x = celltype, y = pathway, 
                         size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "-log10(FDR)", 
                       range = c(2, 8),
                       breaks = c(1, 2, 3, 4)) +
  scale_color_gradient2(name = "NES",
                       low = "blue", 
                       mid = "white", 
                       high = "red",
                       midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = "grey80")) +
  labs(x = "Cell Type", y = "Pathway",
       title = "Cell Type Specific Pathway Enrichment")

pdf("F:/workspace/CDH11/scRNA_enrich_diff_top5.pdf",8,8)
print(bubble_plot)
dev.off()
save.image(file="F:/workspace/CDH11/scRNA_GSEA.rda")
load("F:/workspace/CDH11/scRNA_GSEA.rda")
