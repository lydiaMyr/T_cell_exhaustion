load("F:/workspace/CDH11/CDH11_group_gsea.rda")
library(GSVA)
library(ggplot2)
library(patchwork)
library(ggpubr)


# 1. 根据CDH11表达将样本分为高低两组


tpm_matrix <- read.table("F:/workspace/CDH11/TCGA-COAD.star_tpm.tsv",header=T,sep="\t",check.names = F)
print("原始数据维度:")
print(dim(tpm_matrix))
print("前几个Ensembl_ID:")
print(head(tpm_matrix$Ensembl_ID))

# 连接到Ensembl数据库
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# 提取Ensembl ID（去掉版本号）
ensembl_ids <- gsub("\\.[0-9]+$", "", tpm_matrix$Ensembl_ID)

# 获取基因生物类型信息
gene_biotypes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# 查看基因类型分布
cat("基因类型分布:\n")
gene_type_table <- table(gene_biotypes$gene_biotype)
print(gene_type_table)

# 筛选蛋白编码基因
protein_coding_ensembl <- gene_biotypes %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(ensembl_gene_id)

print(paste("找到", length(protein_coding_ensembl), "个蛋白编码基因"))

# 在原始count_matrix中筛选
tpm_matrix_pc <- tpm_matrix %>%
  mutate(ensembl_id_no_version = gsub("\\.[0-9]+$", "", Ensembl_ID)) %>%
  filter(ensembl_id_no_version %in% protein_coding_ensembl)
 # select(-ensembl_id_no_version)
tpm_matrix_pc=tpm_matrix_pc[,-516]
print(paste("筛选后矩阵维度:", dim(tpm_matrix_pc)))

library(biomaRt)

# 连接到Ensembl数据库
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# 提取Ensembl ID（去掉版本号）
ensembl_ids <- gsub("\\.[0-9]+$", "", tpm_matrix_pc$Ensembl_ID)

# 获取基因符号映射
gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# 查看映射结果
head(gene_annotation)


# 合并转换结果到原始数据
tpm_matrix_with_symbol <- tpm_matrix_pc

# 添加去版本化的Ensembl ID用于匹配
tpm_matrix_with_symbol$ensembl_id_no_version <- gsub("\\.[0-9]+$", "", tpm_matrix_with_symbol$Ensembl_ID)

# 方法一：使用biomaRt结果
tpm_matrix_with_symbol_final <- merge(tpm_matrix_with_symbol, gene_annotation, 
                           by.x = "ensembl_id_no_version", 
                           by.y = "ensembl_gene_id", all.x = TRUE)


# 查看合并后的数据
head(tpm_matrix_with_symbol_final)


# 处理重复的基因符号
library(dplyr)
tpm_matrix_with_symbol_final <- tpm_matrix_with_symbol_final %>%
  group_by(hgnc_symbol) %>%  # 或 SYMBOL，根据使用的包
  mutate(dup_count = n()) %>%
  ungroup()

# 移除空基因符号的行
tpm_matrix_with_symbol_final <- tpm_matrix_with_symbol_final[!is.na(tpm_matrix_with_symbol_final$hgnc_symbol) & 
                                         tpm_matrix_with_symbol_final$hgnc_symbol != "", ]

# 处理重复基因符号：取表达量平均值或最大值
tpm_matrix_aggregated <- tpm_matrix_with_symbol_final %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# 设置行名为基因符号
rownames(tpm_matrix_aggregated) <- tpm_matrix_aggregated$hgnc_symbol

# 移除不再需要的列
# tpm_matrix_final <- tpm_matrix_aggregated %>%
#   select(-hgnc_symbol, -dup_count)
tpm_matrix_final=tpm_matrix_aggregated[-which(colnames(tpm_matrix_aggregated)%in%c("hgnc_symbol","dup_count"))]

# 查看最终数据
head(tpm_matrix_final[, 1:5])  # 显示前5列
row.names(tpm_matrix_final)=tpm_matrix_aggregated$hgnc_symbol


tpm_tumor <- tpm_matrix_final[, tumor_samples]
row.names(tpm_tumor)=row.names(tpm_matrix_final)


# T细胞耗竭基因
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TIM3", 
                      "TOX", "ENTPD1", "CD39", "BATF", "TCF7", "EOMES", 
                      "TBX21", "CXCL13", "IFNG", "GZMB", "PRF1")

# 纤维化相关基因
fibrosis_genes <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "ACTA2", "TGFB1",
                   "TGFB2", "TGFB3", "MMP2", "MMP9", "TIMP1", "TIMP2",
                   "POSTN", "FAP", "PDGFRA", "VIM", "S100A4", "LOXL2")

# 目标基因：CDH11 + 耗竭基因 + 纤维化基因
target_genes <- c("CDH11", exhaustion_genes, fibrosis_genes)

# 过滤在TPM矩阵中存在的基因
available_genes <- target_genes[target_genes %in% rownames(tpm_matrix_final)]
print(paste("可分析的基因数量:", length(available_genes)))

# 提取目标基因的表达矩阵
expr_matrix <- tpm_matrix_final[available_genes, ]
row.names(expr_matrix) = available_genes
# 转置矩阵（样本在行，基因在列）
expr_df <- as.data.frame(t(expr_matrix))

# 计算相关性矩阵
cor_matrix <- cor(expr_df, method = "spearman")  # 使用Spearman相关性
#row.names(cor_matrix)=
# 提取CDH11与其他基因的相关性
cdh11_cor <- cor_matrix["CDH11", ]
cdh11_cor <- cdh11_cor[names(cdh11_cor) != "CDH11"]  # 移除CDH11与自身的相关性

print("CDH11与其他基因的相关性:")
print(sort(cdh11_cor, decreasing = TRUE))

library(pheatmap)
library(viridis)

# Remove CDH11 from the correlation matrix
genes_without_cdh11 <- available_genes[!available_genes %in% c("CDH11","TCF7","ENTPD1")]
# genes_without_cdh11 <- available_genes[available_genes != "TCF7"]
# genes_without_cdh11 <- available_genes[available_genes != "ENTPD1"]
plot_cor_matrix_noCDH11 <- cor_matrix[genes_without_cdh11, genes_without_cdh11]

# Prepare annotation information
gene_types_noCDH11 <- data.frame(
  GeneType = ifelse(rownames(plot_cor_matrix_noCDH11) %in% exhaustion_genes, "T_cell_exhaustion", "Fibrosis"),
  row.names = rownames(plot_cor_matrix_noCDH11)
)

# Define colors
annotation_colors <- list(
  GeneType = c("T_cell_exhaustion" = "#2e9df7", "Fibrosis" = "#ffc907")
)

# Draw heatmap with English labels
pdf("F:/workspace/CDH11/result_new/CDH11_exhausted_Fib_cor.pdf",8,8)
pheatmap(plot_cor_matrix_noCDH11,
         main = "Correlation Heatmap: T Cell Exhaustion vs Fibrosis Genes",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 7,
         fontsize_col = 7,
         annotation_row = gene_types_noCDH11,
         annotation_colors = annotation_colors,
         display_numbers = FALSE,
         border_color = "grey60",
         cellwidth = 12,
         cellheight = 12)
dev.off()
# Alternative with viridis color scheme
pheatmap(plot_cor_matrix_noCDH11,
         main = "Gene Correlation Matrix\n(T Cell Exhaustion vs Fibrosis Genes)",
         color = viridis::viridis(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_row = gene_types_noCDH11,
         annotation_colors = annotation_colors,
         annotation_legend = TRUE,
         legend = TRUE,
         border_color = NA)
library(GSVA)


cdh11_expression <- as.numeric(tpm_tumor["CDH11", ]) # 假设expr_matrix是您的表达矩阵
cdh11_median <- median(cdh11_expression, na.rm = TRUE)
cdh11_group <- ifelse(cdh11_expression > cdh11_median, "CDH11_High", "CDH11_Low")

# 2. 定义T细胞耗竭和纤维化相关的基因集
# T细胞耗竭基因集
t_exhaustion_genes <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", 
                        "HAVCR2", "TIGIT", "BTLA", "CD160", "VSIR")

# 纤维化相关基因集
fibrosis_genes <- c("ACTA2", "COL1A1", "COL1A2", "COL3A1", "FN1", "TGFB1", "TGFBR1", "TGFBR2", "MMP2", "MMP9", "TIMP1", "CTGF")

# 3. 计算基因集评分（使用ssGSEA）
gene_sets <- list(
  T_cell_exhaustion = t_exhaustion_genes,
  Fibrosis = fibrosis_genes
)

# 计算ssGSEA分数
gsva_scores <- gsva(as.matrix(log(tpm_tumor+1)), gene_sets, method = "ssgsea", verbose = FALSE)
cor.test(gsva_scores[1,],gsva_scores[2,])

# 提取两个基因集的GSVA得分
score1 <- gsva_scores[1, ]
score2 <- gsva_scores[2, ]

# 获取基因集名称
geneset1_name <- rownames(gsva_scores)[1]
geneset2_name <- rownames(gsva_scores)[2]

# 计算相关性
cor_test <- cor.test(score1, score2, method = "spearman")
cor_coef <- round(cor_test$estimate, 3)
cor_pvalue <- format(cor_test$p.value, scientific = TRUE, digits = 3)

# 绘制相关性点图
library(ggplot2)
library(ggpubr)
pdf("F:/workspace/CDH11/result_new/CDH11_exhausted_fibro_cor.pdf",4,4)
ggplot(data = data.frame(Score1 = score1, Score2 = score2), 
       aes(x = Score1, y = Score2)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = paste("Correlation between", geneset1_name, "and", geneset2_name),
       x = paste(geneset1_name, "GSVA Score"),
       y = paste(geneset2_name, "GSVA Score"),
       subtitle = paste("Spearman r =", cor_coef, ", p =", cor_pvalue)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()
# 转换为数据框便于绘图
score_df <- data.frame(
  Sample = colnames(tpm_tumor),
  CDH11_Group = cdh11_group,
  T_exhaustion_Score = as.numeric(gsva_scores["T_cell_exhaustion", ]),
  Fibrosis_Score = as.numeric(gsva_scores["Fibrosis", ])
)


# 执行GSEA分析
library(clusterProfiler)
library(org.Hs.eg.db)
expr_matrix=tpm_tumor
# 准备排名基因列表（基于与CDH11的相关性）
cor_results <- apply(expr_matrix, 1, function(x) cor(x, cdh11_expression, method = "spearman"))
gene_rank <- sort(cor_results, decreasing = TRUE)

# 准备基因集
gene_set_list <- list(
  T_cell_exhaustion = t_exhaustion_genes,
  Fibrosis = fibrosis_genes
)

# 执行GSEA
gsea_results <- GSEA(gene_rank, 
                     TERM2GENE = data.frame(term = rep(names(gene_set_list), 
                                                     lengths(gene_set_list)),
                                          gene = unlist(gene_set_list)),
                     pvalueCutoff = 1, # 先不设阈值查看所有结果
                     verbose = FALSE)


# 绘制T细胞耗竭评分比较
p1 <- ggplot(score_df, aes(x = CDH11_Group, y = T_exhaustion_Score, fill = CDH11_Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_fill_manual(values = c("CDH11_High" = "#E74C3C", "CDH11_Low" = "#3498DB")) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     label.x = 1.5, label.y = max(score_df$T_exhaustion_Score) * 0.95) +
  theme_classic() +
  labs(x = "", y = "T cell Exhaustion Score", 
       title = "T Cell Exhaustion Signature Enrichment") +
  theme(legend.position = "none")

# 绘制纤维化评分比较
p2 <- ggplot(score_df, aes(x = CDH11_Group, y = Fibrosis_Score, fill = CDH11_Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  scale_fill_manual(values = c("CDH11_High" = "#E74C3C", "CDH11_Low" = "#3498DB")) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     label.x = 1.5, label.y = max(score_df$Fibrosis_Score) * 0.95) +
  theme_classic() +
  labs(x = "", y = "Fibrosis Score", 
       title = "Fibrosis Signature Enrichment") +
  theme(legend.position = "none")

# 组合图形
combined_plot <- p1 + p2 + 
  plot_annotation(title = "CDH11 High Group Shows Enhanced T Cell Exhaustion and Fibrosis Signatures",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

print(combined_plot)
ggsave("CDH11_Enrichment_Analysis.pdf", combined_plot, width = 12, height = 5)


# 执行GSEA分析
library(clusterProfiler)
library(org.Hs.eg.db)

# 准备排名基因列表（基于与CDH11的相关性）
cor_results <- apply(expr_matrix, 1, function(x) cor(x, cdh11_expression, method = "spearman"))
gene_rank <- sort(cor_results, decreasing = TRUE)

# 准备基因集
gene_set_list <- list(
  T_cell_exhaustion = t_exhaustion_genes,
  Fibrosis = fibrosis_genes
)
# 按功能分组，便于后续分析
t_exhaustion_gene_sets <- list(
  # 主要免疫检查点
  immune_checkpoints = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", 
                        "HAVCR2", "TIGIT", "BTLA", "CD160", "VSIR"),
  
  # 转录调控网络
  transcription_factors = c("TOX", "TOX2", "NFATC1", "NFATC2", "PRDM1", 
                           "BATF", "IRF4", "TCF7", "EOMES", "TBX21"),
  
  # 代谢和信号
  metabolic_signaling = c("ENTPD1", "NT5E", "IDO1", "IDO2", "HK2", "LDHA",
                         "SLC2A1", "ENO1", "PKM", "ACLY", "DGKA"),
  
  # 细胞因子环境
  cytokine_environment = c("IL10", "IL10RA", "IL10RB", "TGFB1", "TGFBR1",
                          "TGFBR2", "IL2RA", "IL7R", "IL2RB"),
  
  # 抑制性受体
  inhibitory_receptors = c("CD244", "KLRG1", "KLRB1", "LAYN", "CD38",
                          "CD101", "CD96", "SIGLEC7", "SIGLEC9", "LAIR1"),
  
  # 效应功能（通常下调）
  effector_functions = c("GZMB", "IFNG", "TNF", "IL2", "PRF1")
)

# 按功能分组纤维化相关marker
fibrosis_gene_sets <- list(
  # 细胞外基质成分 - 核心纤维化标志物
  extracellular_matrix = c("COL1A1", "COL1A2", "COL3A1", "FN1", "LAMB1", 
                          "LAMC1", "ELN", "TNC", "SPARC", "BGN"),
  
  # 基质重塑酶 - ECM降解和修饰
  matrix_remodeling_enzymes = c("MMP2", "MMP9", "MMP14", "TIMP1", "TIMP2",
                               "LOX", "LOXL1", "LOXL2", "PLOD2", "HSPG2"),
  
  # 促纤维化细胞因子和生长因子
  profibrotic_cytokines = c("TGFB1", "TGFB2", "TGFB3", "CTGF", "PDGFA", 
                           "PDGFB", "FGF2", "VEGFA", "IL6", "IL11", "IL13"),
  
  # 肌成纤维细胞标志物 - 主要效应细胞
  myofibroblast_markers = c("ACTA2", "TAGLN", "MYH11", "MYLK", "DES", 
                           "VIM", "FAP", "PDGFRB", "THY1", "CNN1"),
  
  # 纤维化转录调控因子
  transcription_regulators = c("TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2",
                              "SMAD2", "SMAD3", "SMAD4", "STAT3", "YAP1", "TAZ"),
  
  # 整合素和粘附分子
  integrins_adhesion = c("ITGA1", "ITGA2", "ITGA5", "ITGAV", "ITGB1", 
                        "ITGB3", "ITGB5", "CDH2", "CDH11"),
  
  # 炎症和免疫调节 - 纤维化驱动因素
  inflammation_immune = c("IL1B", "IL17A", "TNF", "CCL2", "CCL5", 
                         "CXCL8", "CXCL12", "NFKB1", "NFKB2"),
  
  # 内质网应激和凋亡相关
  stress_apoptosis = c("HSPA5", "XBP1", "ATF4", "ATF6", "DDIT3", 
                      "BCL2", "BAX", "CASP3", "CASP8")
)

gene_set_list_all=c(t_exhaustion_gene_sets,fibrosis_gene_sets)
library(GSVA)
gsea_results_exhausted <- GSEA(gene_rank, 
                     TERM2GENE = data.frame(term = rep(names(t_exhaustion_gene_sets), 
                                                     lengths(t_exhaustion_gene_sets)),
                                          gene = unlist(t_exhaustion_gene_sets)),
                     pvalueCutoff = 1, # 先不设阈值查看所有结果
                     verbose = FALSE)

gsea_results_fibro <- GSEA(gene_rank, 
                     TERM2GENE = data.frame(term = rep(names(fibrosis_gene_sets), 
                                                     lengths(fibrosis_gene_sets)),
                                          gene = unlist(fibrosis_gene_sets)),
                     pvalueCutoff = 1, # 先不设阈值查看所有结果
                     verbose = FALSE)     
save(gsea_results_exhausted,file="F:/workspace/CDH11/gsea_results_exhausted.rda")  
save(gsea_results_fibro,file="F:/workspace/CDH11/gsea_results_fibro.rda")  


load("gsea_results_exhausted.rda")
load("gsea_results_fibro.rda")
library(enrichplot)   
gsea_plot <- gseaplot2(gsea_results_exhausted,
                       geneSetID = 1:length(gsea_results_exhausted$Description),  # 所有通路
                       title = "GSEA Enrichment Plot - Exhausted Pathways",
                       color = scales::hue_pal()(length(gsea_results_exhausted$Description)), 
                       base_size = 12)
gsea_plot2 <- gseaplot2(gsea_results_fibro,
                       geneSetID = 1:length(gsea_results_fibro$Description),  # 所有通路
                       title = "GSEA Enrichment Plot - Exhausted Pathways",
                       color = scales::hue_pal()(length(gsea_results_fibro$Description)), 
                       base_size = 12)
pdf("GSEA_plot.pdf",5,4)
print(gsea_plot)
print(gsea_plot2)
dev.off()

# # 山脊图展示所有通路
# ridge_plot <- ridgeplot(gsea_results_exhausted,
#                        showCategory = length(gsea_results_exhausted$Description),
#                        fill = "p.adjust",  # 根据校正p值填充颜色
#                        core_enrichment = TRUE,
#                        label_format = 30)  # 标签长度

# print(ridge_plot)
# 执行GSEA
gsea_results <- GSEA(gene_rank, 
                     TERM2GENE = data.frame(term = rep(names(gene_set_list), 
                                                     lengths(gene_set_list)),
                                          gene = unlist(gene_set_list)),
                     pvalueCutoff = 1, # 先不设阈值查看所有结果
                     verbose = FALSE)

gsea_results_2 <- GSEA(gene_rank, 
                     TERM2GENE = data.frame(term = rep(names(t_exhaustion_gene_sets), 
                                                     lengths(t_exhaustion_gene_sets)),
                                          gene = unlist(t_exhaustion_gene_sets)),
                     pvalueCutoff = 1, # 先不设阈值查看所有结果
                     verbose = FALSE)

# 方法一：使用enrichplot包的内置函数（推荐）
library(enrichplot)
# 尝试修复的gseaplot2

# 检查并安装最新版本的enrichplot
if (!require("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot", update = TRUE)
}
library(enrichplot)

# 手动绘制GSEA图（完全自定义）
# 安装并加载必要的包
if (!require("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}
library(enrichplot)
library(ggplot2)

save(gsea_results,file="F:/workspace/CDH11/result_new/gsea_results.rda")
# 绘制GSEA富集图

gsea_plot <- gseaplot2(gsea_results,
                       geneSetID = 1,  # 绘制第一个富集通路
                       title = gsea_results$Description[1],
                       color = "red",
                       base_size = 12)
gsea_plot2 <- gseaplot2(gsea_results,
                       geneSetID = 2,  # 绘制第一个富集通路
                       title = gsea_results$Description[2],
                       color = "red",
                       base_size = 12)
pdf("GSEA_plot.pdf",4,4)
print(gsea_plot)
print(gsea_plot2)
dev.off()
ggsave("GSEA_Plot_Single.pdf", gsea_plot, width = 10, height = 8)

