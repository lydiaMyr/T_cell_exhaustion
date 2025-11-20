# 安装并加载必要的包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
BiocManager::install("limma")

library(TCGAbiolinks)
library(DESeq2)
library(limma)
library(ggplot2)
library(dplyr)
library(pheatmap)

# 1. 查询和下载TCGA-COAD数据
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 下载数据
GDCdownload(query)

# 准备表达矩阵
data <- GDCprepare(query)

# 获取原始counts矩阵
counts_matrix <- assay(data)



counts_matrix <- read.table("F:/workspace/CDH11/TCGA-COAD.star_counts.tsv",header=T,sep="\t",check.names = F)
tpm_matrix <- read.table("F:/workspace/CDH11/TCGA-COAD.star_tpm.tsv",header=T,sep="\t",check.names = F)


library(biomaRt)
library(dplyr)

# 查看原始数据
print("原始数据维度:")
print(dim(counts_matrix))
print("前几个Ensembl_ID:")
print(head(counts_matrix$Ensembl_ID))

# 连接到Ensembl数据库
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# 提取Ensembl ID（去掉版本号）
ensembl_ids <- gsub("\\.[0-9]+$", "", counts_matrix$Ensembl_ID)

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
count_matrix_pc <- counts_matrix %>%
  mutate(ensembl_id_no_version = gsub("\\.[0-9]+$", "", Ensembl_ID)) %>%
  filter(ensembl_id_no_version %in% protein_coding_ensembl)
count_matrix_pc=count_matrix_pc[,-516]
print(paste("筛选后矩阵维度:", dim(count_matrix_pc)))

# library(biomaRt)

# # 连接到Ensembl数据库
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# # 提取Ensembl ID（去掉版本号）
# ensembl_ids <- gsub("\\.[0-9]+$", "", count_matrix_pc$Ensembl_ID)

# # 获取基因符号映射
# gene_annotation <- getBM(
#   attributes = c("ensembl_gene_id", "hgnc_symbol"),
#   filters = "ensembl_gene_id",
#   values = ensembl_ids,
#   mart = ensembl
# )

# 查看映射结果
head(gene_annotation)


# 合并转换结果到原始数据
count_matrix_with_symbol <- count_matrix_pc

# 添加去版本化的Ensembl ID用于匹配
count_matrix_with_symbol$ensembl_id_no_version <- gsub("\\.[0-9]+$", "", count_matrix_with_symbol$Ensembl_ID)

# 方法一：使用biomaRt结果
count_matrix_final <- merge(count_matrix_with_symbol, gene_annotation, 
                           by.x = "ensembl_id_no_version", 
                           by.y = "ensembl_gene_id", all.x = TRUE)


# 查看合并后的数据
head(count_matrix_final)


# 处理重复的基因符号
library(dplyr)
count_matrix_final <- count_matrix_final %>%
  group_by(hgnc_symbol) %>%  # 或 SYMBOL，根据使用的包
  mutate(dup_count = n()) %>%
  ungroup()

# 移除空基因符号的行
count_matrix_final <- count_matrix_final[!is.na(count_matrix_final$hgnc_symbol) & 
                                         count_matrix_final$hgnc_symbol != "", ]

# 处理重复基因符号：取表达量平均值或最大值
count_matrix_aggregated <- count_matrix_final %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# 设置行名为基因符号
rownames(count_matrix_aggregated) <- count_matrix_aggregated$hgnc_symbol

library(dplyr)
# 移除不再需要的列
# count_matrix_final <- count_matrix_aggregated %>%
#   select(-hgnc_symbol, -dup_count)
count_matrix_final=count_matrix_aggregated[-which(colnames(count_matrix_aggregated)%in%c("hgnc_symbol","dup_count"))]
# 查看最终数据
head(count_matrix_final[, 1:5])  # 显示前5列
row.names(count_matrix_final)=count_matrix_aggregated$hgnc_symbol





# 加载必要的包
library(DESeq2)
library(dplyr)
library(tibble)

# 查看数据
print("TPM矩阵维度:")
print(dim(tpm_matrix))
print("Count矩阵维度:")
print(dim(count_matrix_final))

# CDH11的Ensembl ID (根据您的数据调整)
cdh11_ensembl <- "ENSG00000134575"  # CDH11的Ensembl ID

# 在TPM矩阵中提取肿瘤样本（假设肿瘤样本名包含"01A"）
tumor_samples <- grep("01A$", colnames(tpm_matrix), value = TRUE)
print(paste("找到", length(tumor_samples), "个肿瘤样本"))

# 提取肿瘤样本的TPM矩阵
tpm_tumor <- tpm_matrix[, tumor_samples]

cdh11_tpm <- as.numeric(tpm_tumor[grep(cdh11_ensembl,tpm_matrix$Ensembl_ID), ])

# 根据CDH11表达中位数分组
cdh11_median <- median(cdh11_tpm, na.rm = TRUE)
cdh11_group <- ifelse(cdh11_tpm > cdh11_median, "CDH11_High", "CDH11_Low")

# 创建样本信息数据框
sample_info <- data.frame(
  sample = tumor_samples,
  group = cdh11_group,
  cdh11_expression = cdh11_tpm,
  row.names = tumor_samples
)

print("样本分组信息:")
print(table(sample_info$group))
print(head(sample_info))


# 确保count矩阵只包含肿瘤样本
count_tumor <- count_matrix_final[, tumor_samples]

# 检查矩阵是否包含负值或非整数值（DESeq2需要count数据）
print("Count数据摘要:")
print(summary(as.numeric(as.matrix(count_tumor))))

# 如果有负值或小数，可能需要处理（取整等）
if(any(count_tumor < 0)) {
  warning("Count矩阵包含负值,这可能不是原始的count数据")
}

# 确保count数据是数值矩阵
count_tumor <- apply(count_tumor, 2, as.numeric)
rownames(count_tumor) <- rownames(count_matrix_final)

print(paste("Count肿瘤矩阵维度:", dim(count_tumor)))


# 确保样本顺序一致
common_samples <- intersect(colnames(count_tumor), rownames(sample_info))
count_tumor <- count_tumor[, common_samples]
sample_info <- sample_info[common_samples, ]

print(paste("共同样本数:", length(common_samples)))

# 创建DESeq2数据对象
dds <- DESeqDataSetFromMatrix(
  countData = round(count_tumor),  # 确保是整数
  colData = sample_info,
  design = ~ group
)

print("DESeq2数据对象创建成功:")
print(dds)

# 运行DESeq2分析
dds <- DESeq(dds)

# 获取结果
results <- results(dds, contrast = c("group", "CDH11_High", "CDH11_Low"))
results_df <- as.data.frame(results) %>%
  rownames_to_column("gene") %>%
  arrange(padj)  # 按调整p值排序

print("差异表达分析完成!")
print(paste("总基因数:", nrow(results_df)))
print(paste("显著差异基因数 (padj < 0.05):", sum(results_df$padj < 0.05, na.rm = TRUE)))
print(paste("显著上调基因数:", sum(results_df$padj < 0.05 & results_df$log2FoldChange > 0, na.rm = TRUE)))
print(paste("显著下调基因数:", sum(results_df$padj < 0.05 & results_df$log2FoldChange < 0, na.rm = TRUE)))

# 安装并加载可视化包
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
library(ggplot2)
library(EnhancedVolcano)


# 首先提取top5上调和下调的基因
# 获取log2FoldChange最大的5个上调基因
top_up <- results_df %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  head(5) %>%
  pull(gene)

# 获取log2FoldChange最小的5个下调基因
top_down <- results_df %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  head(5) %>%
  pull(gene)

top_genes <- c(top_up, top_down)



# 创建火山图
volcano_plot <- EnhancedVolcano(results_df,
                lab = results_df$gene,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top_genes,
                title = 'CDH11 High vs Low - Differential Expression',
                subtitle = 'Top 5 up and down genes by log2FoldChange',
                pCutoff = 0.05,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 0.6,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 20)

print(volcano_plot)
# 火山图
ggsave("F:/workspace/CDH11/CDH11_DE_Volcano.pdf", volcano_plot, width = 4, height = 6)

# MA图
plotMA(results, main = "MA Plot: CDH11 High vs Low", ylim = c(-5, 5))

# 提取显著差异基因
significant_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.58) %>%
  arrange(desc(log2FoldChange))

print("Top 10上调基因:")
print(head(significant_genes[significant_genes$log2FoldChange > 0, ], 10))

print("Top 10下调基因:")
print(head(significant_genes[significant_genes$log2FoldChange < 0, ], 10))


#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# 获取上下调基因（按log2FoldChange排序）
up_genes <- results_df %>%
  filter(log2FoldChange > 0.5) %>%filter(padj < 0.05)%>% 
  arrange(desc(log2FoldChange)) %>%
  head(500)  # 取前200个上调基因用于富集分析

down_genes <- results_df %>%
  filter(log2FoldChange < (-0.3)) %>%filter(padj < 0.05)%>% 
  arrange(log2FoldChange) %>%
  head(500)  # 取前200个下调基因用于富集分析

print(paste("上调基因数量:", nrow(up_genes)))
print(paste("下调基因数量:", nrow(down_genes)))


# 将基因名转换为Entrez ID
up_entrez <- bitr(up_genes$gene, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

down_entrez <- bitr(down_genes$gene, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

print(paste("成功转换的上调基因:", nrow(up_entrez)))
print(paste("成功转换的下调基因:", nrow(down_entrez)))

# 上调基因KEGG分析
kegg_up <- enrichKEGG(gene         = up_entrez$ENTREZID,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1)

# 下调基因KEGG分析
kegg_down <- enrichKEGG(gene         = down_entrez$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1)
if(nrow(kegg_up) > 0) {
  dotplot(kegg_up, showCategory = 15, title = "上调基因 - KEGG通路")
}
if(nrow(kegg_down) > 0) {
  dotplot(kegg_down, showCategory = 15, title = "下调基因 - KEGG通路")
}

#GO富集
# 生物过程
ego_up_bp <- enrichGO(gene          = up_entrez$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE)

# 细胞组分
ego_up_cc <- enrichGO(gene          = up_entrez$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE)

# 分子功能
ego_up_mf <- enrichGO(gene          = up_entrez$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE)


# 生物过程
ego_down_bp <- enrichGO(gene          = down_entrez$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1,
                       readable      = TRUE)

# 细胞组分
ego_down_cc <- enrichGO(gene          = down_entrez$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1,
                       readable      = TRUE)

# 分子功能
ego_down_mf <- enrichGO(gene          = down_entrez$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1,
                       readable      = TRUE)


pdf("F:/workspace/CDH11/result_new/CDH11_DEG_enrich.pdf",5.5,6)
if(nrow(ego_up_bp) > 0) {
  dotplot(ego_up_bp, showCategory = 15, title = "up-regulated genes") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 下调基因GO气泡图
if(nrow(ego_down_bp) > 0) {
  dotplot(ego_down_bp, showCategory = 15, title = "down-regulated genes") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
dev.off()
save.image("F:/workspace/CDH11/CDH11_DEG.rda",compress = T,ascii = F)
load("F:/workspace/CDH11/CDH11_DEG.rda")
