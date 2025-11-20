ME5=read.csv("F:/workspace/CDH11/result_new/WGCNA_geneInfo_ME5_module.csv")
ME9=read.csv("F:/workspace/CDH11/result_new/WGCNA_geneInfo_ME9_module.csv")

#hub gene 筛选并展示
# 加载包
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# 1. 提取Hub基因
target_GS <- "GS_Stromal"
hub_genes <- ME5 %>%
  filter(abs(MM) > 0.4 & (!!sym(target_GS)) > 0.4) %>%
  arrange(desc(abs(!!sym(target_GS))))

# 2. 打印Hub基因列表
print("Hub Genes:")
print(hub_genes)

# 3. 绘制散点图
p1 <- ggplot(ME5, aes(x = MM, y = !!sym(target_GS))) +
  geom_point(color = "grey80") +
  geom_point(data = hub_genes, color = "red") +
  geom_text_repel(data = hub_genes, aes(label = Gene), color = "red", size = 3, max.overlaps = 15) +
  labs(title = "Module Membership vs. Gene Significance", x = "MM", y = target_GS) +
  theme_minimal()
pdf("F:/workspace/CDH11/result_new/ME5_strommal.pdf",4,4)
print(p1)
dev.off()

#ME5 CDH11
# 1. 提取Hub基因
target_GS <- "GS_CDH11"
hub_genes <- ME5 %>%
  filter(abs(MM) > 0.4 & (!!sym(target_GS)) > 0.4) %>%
  arrange(desc(abs(!!sym(target_GS))))

# 2. 打印Hub基因列表
print("Hub Genes:")
print(hub_genes)

# 3. 绘制散点图
p1 <- ggplot(ME5, aes(x = MM, y = !!sym(target_GS))) +
  geom_point(color = "grey80") +
  geom_point(data = hub_genes, color = "red") +
  geom_text_repel(data = hub_genes, aes(label = Gene), color = "red", size = 3, max.overlaps = 15) +
  labs(title = "Module Membership vs. Gene Significance", x = "MM", y = target_GS) +
  theme_minimal()
pdf("F:/workspace/CDH11/result_new/ME5_CDH11.pdf",4,4)
print(p1)
dev.off()

# 4. 绘制热图 (需要 datExpr 数据)
# hub_expr <- datExpr[, hub_genes$Gene]
# pheatmap(t(hub_expr), scale = "row", main = "Hub Genes Expression Heatmap")

# 5. 绘制Top Hub基因条形图
top_hub <- head(hub_genes, 15)
p2 <- ggplot(top_hub, aes(x = reorder(Gene, !!sym(target_GS)), y = !!sym(target_GS))) +
  geom_col(aes(fill = MM)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.85) +
  coord_flip() +
  labs(title = paste("Top Hub Genes for", target_GS), x = "Gene") +
  theme_minimal()

print(p2)


# 1. 提取Hub基因
target_GS <- "GS_Exhausted"
hub_genes <- ME9 %>%
  filter(abs(MM) > 0.4 & abs(!!sym(target_GS)) > 0.3) %>%
  arrange(desc(abs(!!sym(target_GS))))

# 2. 打印Hub基因列表
print("Hub Genes:")
print(hub_genes)

# 3. 绘制散点图
p1 <- ggplot(ME9, aes(x = MM, y = !!sym(target_GS))) +
  geom_point(color = "grey80") +
  geom_point(data = hub_genes, color = "red") +
  geom_text_repel(data = hub_genes, aes(label = Gene), color = "red", size = 3, max.overlaps = 15) +
  labs(title = "Module Membership vs. Gene Significance", x = "MM", y = target_GS) +
  theme_minimal()

pdf("F:/workspace/CDH11/result_new/ME9_exhausted.pdf",4,4)
print(p1)
dev.off()


target_GS <- "GS_CDH11"
hub_genes <- ME9 %>%
  filter(abs(MM) > 0.4 & abs(!!sym(target_GS)) > 0.4) %>%
  arrange(desc(abs(!!sym(target_GS))))

# 2. 打印Hub基因列表
print("Hub Genes:")
print(hub_genes)

# 3. 绘制散点图
p1 <- ggplot(ME9, aes(x = MM, y = !!sym(target_GS))) +
  geom_point(color = "grey80") +
  geom_point(data = hub_genes, color = "red") +
  geom_text_repel(data = hub_genes, aes(label = Gene), color = "red", size = 3, max.overlaps = 15) +
  labs(title = "Module Membership vs. Gene Significance", x = "MM", y = target_GS) +
  theme_minimal()

pdf("F:/workspace/CDH11/result_new/ME9_CDH11.pdf",4,4)
print(p1)
dev.off()
# 4. 绘制热图 (需要 datExpr 数据)
# hub_expr <- datExpr[, hub_genes$Gene]
# pheatmap(t(hub_expr), scale = "row", main = "Hub Genes Expression Heatmap")

# 5. 绘制Top Hub基因条形图
top_hub <- head(hub_genes, 15)
p2 <- ggplot(top_hub, aes(x = reorder(Gene, !!sym(target_GS)), y = !!sym(target_GS))) +
  geom_col(aes(fill = MM)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.85) +
  coord_flip() +
  labs(title = paste("Top Hub Genes for", target_GS), x = "Gene") +
  theme_minimal()

print(p2)




target_GS <- "GS_Stromal"
hub_genes <- ME5 %>%
  filter(abs(MM) > 0.3 & (!!sym(target_GS)) > 0.3) %>%
  arrange(desc(abs(!!sym(target_GS))))
print(data.frame(hub_genes[,1]))
target_GS <- "GS_CDH11"
hub_genes <- ME5 %>%
  filter(abs(MM) > 0.3 & (!!sym(target_GS)) > 0.3) %>%
  arrange(desc(abs(!!sym(target_GS))))
print(data.frame(hub_genes[,1]))







# 按基因显著性排序
greenyellow_top <- ME5[order(-abs(ME5$GS_Exhausted)), ]
grey_top <- ME9[order(-abs(ME9$GS_Stromal)), ]

# 取前30个Hub基因进行比较
top_30_gy <- head(greenyellow_top$Gene, 350)
top_30_grey <- head(grey_top$Gene, 350)

print("Greenyellow模块Top Hub基因:")
print(top_30_gy)
print("Grey模块Top Hub基因:")
print(top_30_grey)

# 对两个模块分别进行GO/KEGG富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
# 对两个模块分别进行KEGG富集分析
library(clusterProfiler)

# 将基因符号转换为Entrez ID
gy_entrez <- bitr(top_30_gy, fromType = "SYMBOL", 
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
grey_entrez <- bitr(top_30_grey, fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG富集分析
gy_kegg <- enrichKEGG(
  gene = gy_entrez$ENTREZID,
  organism = "hsa",  # 人类
  pvalueCutoff = 0.25,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"
)

grey_kegg <- enrichKEGG(
  gene = grey_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2, 
  pAdjustMethod = "BH"
)

# 查看结果
print("GreenYellow模块KEGG富集结果:")
print(head(gy_kegg, 10))

print("Grey模块KEGG富集结果:")
print(head(grey_kegg, 10))


# 对比可视化
compare_kegg <- compareCluster(
  geneCluster = list(
    GreenYellow = gy_entrez$ENTREZID,
    Grey = grey_entrez$ENTREZID
  ),
  fun = "enrichKEGG",
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# 绘制点图
kegg_dotplot <- dotplot(compare_kegg, 
                       showCategory = 8,
                       title = "KEGG Pathway Enrichment Comparison") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(kegg_dotplot)



# 分别绘制两个模块的KEGG富集图
library(enrichplot)

if(nrow(gy_kegg) > 0) {
  p1 <- dotplot(gy_kegg, showCategory = 10, title = "GreenYellow Module - KEGG Pathways") +
    theme(axis.text.y = element_text(size = 10))
  print(p1)
}

if(nrow(grey_kegg) > 0) {
  p2 <- dotplot(grey_kegg, showCategory = 10, title = "Grey Module - KEGG Pathways") +
    theme(axis.text.y = element_text(size = 10))
  print(p2)
}

# 使用patchwork组合图形
library(patchwork)
if(exists("p1") & exists("p2")) {
  combined_kegg <- p1 / p2 + 
    plot_annotation(title = "KEGG Pathway Enrichment in Two Modules",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  print(combined_kegg)
}