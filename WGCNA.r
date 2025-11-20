COADMatrix_fpkm <- assay(COADnaseqSE,"fpkm_unstrand")

ensembl_ids <- gsub("\\.[0-9]+$", "", row.names(COADMatrix_fpkm))


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

count_matrix_pc <- data.frame(COADMatrix_fpkm) %>%
  mutate(ensembl_id_no_version = gsub("\\.[0-9]+$", "", row.names(COADMatrix_fpkm))) %>%
  filter(ensembl_id_no_version %in% protein_coding_ensembl)
count_matrix_pc=count_matrix_pc[,-525]
print(paste("筛选后矩阵维度:", dim(count_matrix_pc)))



count_matrix_with_symbol <- count_matrix_pc

# 添加去版本化的Ensembl ID用于匹配
count_matrix_with_symbol$ensembl_id_no_version <- gsub("\\.[0-9]+$", "", row.names(count_matrix_pc))


gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
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
coad_fpkm=count_matrix_final
write.table(coad_fpkm,file="F:/workspace/CDH11/COAD_exp_fpkm.txt",sep="\t",quote=F)
library("estimate")
filterCommonGenes(input.f="F:/workspace/CDH11/COAD_exp_fpkm.txt", output.f="F:/workspace/CDH11/COAD_exp_estimate_score_fpkm.gct", id = "GeneSymbol")
estimateScore('F:/workspace/CDH11/COAD_exp_estimate_score_fpkm.gct','F:/workspace/CDH11/COAD_exp_purify_fpkm.gct' )
result <- read.table('F:/workspace/CDH11/COAD_exp_purify_fpkm.gct',header = T,skip = 2)
head(result)

stromal_score=as.vector(unlist(result[1,-c(1,2)]))
names(stromal_score)=colnames(result[1,-c(1,2)])

tumor_purity=as.vector(unlist(result[4,-c(1,2)]))
names(tumor_purity)=colnames(result[1,-c(1,2)])

sig_sam=names(tumor_purity)[which(tumor_purity>0.8)]
library(GSVA)

exhausted_gene=c("PDCD1","CTLA4","LAG3","HAVCR2","ENTPD1")
fibro_gene = c(
  "TGFB1",      # 转化生长因子β1，纤维化的核心驱动因子
  
  # 2. 胶原合成 - 细胞外基质主要成分
  "COL1A1",     # I型胶原α1链，最主要的纤维化胶原
  
  # 3. 肌成纤维细胞标志物 - 纤维化的效应细胞
  "ACTA2",      # α-平滑肌肌动蛋白，肌成纤维细胞的标志
  
  # 4. 基质 stiffness 调控
  "LOX",        # 赖氨酰氧化酶，介导胶原交联和基质硬化
  
  # 5. 纤维化特异性因子
  "CTGF",       # 结缔组织生长因子，TGFB1的下游效应因子
  
  # 6. 成纤维细胞活化标志物
  "FAP",        # 成纤维细胞活化蛋白，活化成纤维细胞的标志
  
  # 7. 生长因子受体
  "PDGFRB",     # 血小板衍生生长因子受体β，成纤维细胞增殖和迁移
  
  # 8. 转录调控因子
  "TWIST1",     # 上皮-间质转化的关键转录因子
  
  # 9. 基质金属蛋白酶抑制剂
  "TIMP1",      # 金属蛋白酶组织抑制剂1，抑制ECM降解
  
  # 10. 信号通路介质
  "SMAD3"       # TGFB信号通路的关键细胞内介质)"
  )
gene_sets <- list(Exhausted_Signature = exhausted_gene,Fibro_sig = fibro_gene)
expression_matrix = as.matrix(coad_fpkm)
sample_scores <- gsva(expression_matrix,
                         gene_sets,
                         method = "ssgsea",
                         kcdf = "Gaussian",
                         parallel.sz = 1)


tumor_samples <- grep("01A", colnames(expression_matrix), value = TRUE)
print(paste("找到", length(tumor_samples), "个肿瘤样本"))

exhausted_score=sample_scores[1,tumor_samples]
fibro_score=sample_scores[2,tumor_samples]
stromal_score=stromal_score[tumor_samples]
CDH11_exp=as.vector(unlist(expression_matrix["CDH11",tumor_samples]))
names(CDH11_exp)=tumor_samples

#cor.test(exhausted_score,stromal_score)

sample_names <- tumor_samples
trait_samples <- names(stromal_score) # 假设stromal_score是一个行名为样本的数据框
traits=data.frame(stromal_score=stromal_score,exhausted_score=exhausted_score,CDH11_exp=CDH11_exp,fibro_score=fibro_score)
# 找出共有的样本
common_samples <- intersect(sample_names, trait_samples)
caod_tpm_filtered <- expression_matrix[, common_samples]
traits_filtered <- traits[common_samples, ] # 确保traits是一个数据框，行名为样本

library(WGCNA)
filter_genes_for_wgcna <- function(expression_matrix) {
  # 方法1：基于表达水平和方差过滤
  filtered_genes <- function(expression_matrix) {
    # 计算基因在样本中的表达比例（例如：FPKM > 1的比例）
    gene_expr_proportion <- apply(expression_matrix, 1, function(x) sum(x > 1) / length(x))
    
    # 计算基因的表达方差
    gene_variance <- apply(expression_matrix, 1, var)
    
    # 设置过滤阈值
    keep_genes <- (gene_expr_proportion > 0.5) & (gene_variance > quantile(gene_variance, 0.25))
    
    cat("原始基因数量:", nrow(expression_matrix), "\n")
    cat("过滤后基因数量:", sum(keep_genes), "\n")
    cat("过滤比例:", round(1 - sum(keep_genes)/nrow(expression_matrix), 3) * 100, "%\n")
    
    return(expression_matrix[keep_genes, ])
  }
  
  # 方法2：使用WGCNA内置的好基因筛选
  filtered_genes_wgcna <- function(expression_matrix) {
    # 检测"好基因" - 在足够多的样本中有表达
    gsg <- goodSamplesGenes(expression_matrix, verbose = 3)
    
    if (!gsg$allOK) {
      # 如果存在有问题的基因或样本
      if (sum(!gsg$goodGenes) > 0) {
        print(paste("移除基因:", sum(!gsg$goodGenes)))
      }
      if (sum(!gsg$goodSamples) > 0) {
        print(paste("移除样本:", sum(!gsg$goodSamples)))
      }
      # 移除有问题的基因和样本
      expression_matrix <- expression_matrix[gsg$goodGenes, gsg$goodSamples]
    }
    
    return(expression_matrix)
  }
  
  # 执行过滤
  expr_filtered <- filtered_genes(expression_matrix)
  expr_filtered <- filtered_genes_wgcna(expr_filtered)
  
  return(expr_filtered)
}
# 最终确认性状数据框
exp_filter <- filter_genes_for_wgcna(caod_tpm_filtered)
datExpr <- t(exp_filter)
traits_filtered <- traits_filtered[colnames(exp_filter), ]

# 最终确认性状数据框
datTraits <- traits_filtered
# 1. 选择一组软阈值功率
powers <- c(1:20)

# 2. 调用网络拓扑分析函数
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# 3. 绘制结果图
par(mfrow = c(1, 2))
cex1 <- 0.9

# 拟合指数与软阈值的尺度独立性图
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") # 建议将R^2 > 0.9作为标准

# 平均连通性与软阈值图
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# 4. 选择第一个使R^2达到0.9以上的功率
softPower <- sft$powerEstimate
print(paste("Estimated soft threshold power:", softPower))

net <- blockwiseModules(datExpr,
                        power = softPower,
                        networkType = "signed", # 使用有符号网络，能更好区分正负相关
                        TOMType = "signed",
                        minModuleSize = 30,     # 模块最小基因数
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,  # 合并相似度高于0.75的模块
                        numericLabels = TRUE,   # 模块用数字表示
                        pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,       # 为节省空间，不保存TOM矩阵
                        verbose = 3)
table(net$colors)
# 3. 将数字标签转换为颜色
moduleColors <- labels2colors(net$colors)
unique_module_colors <- unique(moduleColors)[order(unique(moduleColors))]


# 1. 计算模块特征基因（Module Eigengene, ME）
MEs <- net$MEs

# 2. 计算模块-性状相关性
moduleTraitCor <- cor(MEs, datTraits[,-1], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))



# 设定p值阈值
# 设定相关系数绝对值阈值
cor_cutoff <- 0.2

# 找出与任意一个性状有较强相关性的模块
high_cor_modules <- which(apply(abs(moduleTraitCor), 1, function(x) any(x > cor_cutoff)))
yColors <- unique_module_colors[high_cor_modules]
# 后续筛选和绘图代码与策略一相同，只需替换变量名
moduleTraitCor_filtered <- moduleTraitCor[high_cor_modules, ]
moduleTraitPvalue_filtered <- moduleTraitPvalue[high_cor_modules, ]

# ... (其余代码与策略一相同)

# 更新文本矩阵（只包含显著模块）
textMatrix_filtered <- paste(signif(moduleTraitCor_filtered, 2), 
                             "\n(", 
                             signif(moduleTraitPvalue_filtered, 1), 
                             ")", sep = "")
dim(textMatrix_filtered) <- dim(moduleTraitCor_filtered)

pdf("F:/workspace/CDH11/result_new/WGCNA_heatmap_fpkm.pdf",4,7)

labeledHeatmap(Matrix = moduleTraitCor_filtered,
               xLabels = colnames(datTraits)[-1],
               yLabels = yColors,        # 使用颜色名称作为标签
               ySymbols = yColors,       # 符号也使用颜色名称
               colorLabels = TRUE,       # 关键参数：启用颜色标签
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_filtered,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships (Significant modules only)"))
dev.off()

labeledHeatmap(Matrix = moduleTraitCor_filtered,
               xLabels = colnames(datTraits)[-1],
               yLabels = names(high_cor_modules),        # 使用颜色名称作为标签
               ySymbols = names(high_cor_modules),       # 符号也使用颜色名称
               colorLabels = TRUE,       # 关键参数：启用颜色标签
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_filtered,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships (Significant modules only)"))

#pink ME3
#purple ME6
target_module_color <- "ME3" 
target_module_genes <- (moduleColors == "pink")

# 2. 计算基因显著性（Gene Significance, GS）：基因与性状的相关性
GS_stromal <- as.numeric(cor(datExpr, datTraits$stromal_score, use = "p"))
GS_exhausted <- as.numeric(cor(datExpr, datTraits$exhausted_score, use = "p"))
GS_CDH11 <- as.numeric(cor(datExpr, datTraits$CDH11_exp, use = "p"))

# 3. 计算模块成员关系（Module Membership, MM）：基因与模块特征基因的相关性
target_ME <- MEs[, target_module_color]
MM <- as.numeric(cor(datExpr, target_ME, use = "p"))

# 4. 绘制Gene Significance vs Module Membership散点图
# 这可以验证模块内的基因确实与外部性状高度相关
library(WGCNA)
par(mfrow = c(1, 3))
verboseScatterplot(MM, GS_stromal,
                   main = paste("MM vs. GS\n", target_module_color, "module"),
                   xlab = "Module Membership",
                   ylab = "Gene Significance (Stromal Score)",
                   abline = TRUE, abline.color = "red")
verboseScatterplot(MM, GS_exhausted,
                   main = paste("MM vs. GS\n", target_module_color, "module"),
                   xlab = "Module Membership",
                   ylab = "Gene Significance (Exhausted Score)",
                   abline = TRUE, abline.color = "red")
verboseScatterplot(MM, GS_CDH11,
                   main = paste("MM vs. GS\n", target_module_color, "module"),
                   xlab = "Module Membership",
                   ylab = "Gene Significance (CDH11 Exp)",
                   abline = TRUE, abline.color = "red")


gene_names <- colnames(datExpr)
target_module_gene_list <- gene_names[target_module_genes]
geneInfo <- data.frame(
    Gene = gene_names,
    ModuleColor = moduleColors,
    GS_Stromal = GS_stromal,
    GS_Exhausted = GS_exhausted,
    GS_CDH11 = GS_CDH11,
    MM = MM
)
library(dplyr)
ME3_target_geneInfo <- geneInfo[geneInfo$ModuleColor == "pink", ]
ME3_target_geneInfo <- ME3_target_geneInfo[order(-ME3_target_geneInfo$GS_Stromal), ]
ME6_target_geneInfo <- geneInfo[geneInfo$ModuleColor == "purple", ]
ME6_target_geneInfo <- ME6_target_geneInfo[order(-ME6_target_geneInfo$GS_Exhausted), ]

target_GS <- "GS_CDH11"
hub_genes <- ME3_target_geneInfo %>%
  filter(abs(MM) > 0.4 & (!!sym(target_GS)) > 0.4) %>%
  arrange(desc(abs(!!sym(target_GS))))
# 验证当前ME3_target_geneInfo中的GS值
# 重新计算这些基因的GS
ME3_hub=c("NDUFAB1","NDUFS3","MRPS7","ATP5PO","ATP5F1C")
ME6_hub=c("CDK1","CCNA2","EXO1","KIF11","CCNB1")
library(pheatmap)
pheatmap(ME3_target_geneInfo[which(ME3_target_geneInfo$Gene%in%ME3_hub),3:5])
pheatmap(ME6_target_geneInfo[which(ME6_target_geneInfo$Gene%in%ME6_hub),3:5])

current_genes <- ME3_target_geneInfo$Gene
GS_recalculated <- cor(datExpr[, current_genes], datTraits$GS_CDH11, use = "p")

# 比较新旧GS值
summary(ME3_target_geneInfo$GS_CDH11)
summary(GS_recalculated)

#ME3_target_geneInfo$Gene
# 1. 重新确认模块基因
ME3_genes <- ME3_target_geneInfo$Gene

# 2. 重新计算这些基因的GS
GS_ME3_genes <- cor(datExpr[, ME3_genes], datTraits$GS_CDH11, use = "p")

# 3. 检查新GS值的分布
summary(GS_ME3_genes)
hist(GS_ME3_genes, main = "Recalculated GS for ME3 genes")

# 4. 验证模块特征基因相关性
ME3_eigengene <- moduleEigengenes(datExpr[, ME3_genes])$eigengenes
cor(ME3_eigengene, traitData$CDH11, use = "p")





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


library(clusterProfiler)
library(org.Hs.eg.db)
me3_genes <- colnames(datExpr)[moduleColors == "pink"] # 请将 "steelblue" 替换为 ME11 对应的实际颜色me9 纤维化
me6_genes <- colnames(datExpr)[moduleColors == "purple"] # 请将 "midnightblue" 替换为 ME12 对应的实际颜色me5  耗竭
gene.df1 <- bitr(me3_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
kk <- enrichKEGG(gene = gene.df1$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.25)

gene.df2 <- bitr(me6_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
kk1 <- enrichKEGG(gene = gene.df2$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.25)

pdf("F:/workspace/CDH11/result_new/Module_enrich_pink_purple.pdf",6,6)
dotplot(kk,showCategory = 20)
dotplot(kk1,showCategory = 20)
dev.off()


go_enrichment_me3 <- enrichGO(gene          = gene.df1$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID", # 现在输入的是ENTREZID
                          ont           = "BP",      # "BP", "MF", "CC" 或 "ALL"
                          pAdjustMethod = "BH",       # 多重检验校正方法：Benjamini-Hochberg
                          pvalueCutoff  = 0.05,       # p值阈值
                          qvalueCutoff  = 0.2,        # q值阈值
                          readable      = TRUE)       # 将ENTREZID转换回Gene Symbol输出

# 查看结果摘要
#head(go_enrichment1, n = 20)[,2]

# 执行GO富集分析
go_enrichment_me6 <- enrichGO(gene          = gene.df2$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID", # 现在输入的是ENTREZID
                          ont           = "BP",      # "BP", "MF", "CC" 或 "ALL"
                          pAdjustMethod = "BH",       # 多重检验校正方法：Benjamini-Hochberg
                          pvalueCutoff  = 0.05,       # p值阈值
                          qvalueCutoff  = 0.2,        # q值阈值
                          readable      = TRUE)       # 将ENTREZID转换回Gene Symbol输出

# 查看结果摘要
head(go_enrichment2, n = 20)[,2]

ME3_hub=read.table("F:/workspace/CDH11/ME3_hub.txt",header=T)
ME6_hub=read.table("F:/workspace/CDH11/ME6_hub.txt",header=T)

ME3_hub_degree=seq(30,1,-1)
names(ME3_hub_degree)=ME3_hub[,1]

go_enrichment_me3_1=go_enrichment_me3[1:30,]
pa_item=go_enrichment_me3_1$Description
path_pva=go_enrichment_me3_1$qvalue
names(path_pva)=pa_item
gene_pa=matrix(rep(0,30*length(path_pva)),nrow=30,ncol=length(path_pva))
row.names(gene_pa)=ME3_hub[,1]
colnames(gene_pa)=pa_item
for(g in ME3_hub[,1]){
    pa_ls=pa_item[grep(g,go_enrichment_me3_1$geneID)]
    if(length(pa_ls)>0){
      gene_pa[g,pa_ls]=1
    }
}
gene_pa_sum=apply(gene_pa,1,sum)
gene_pa=gene_pa[-which(gene_pa_sum==0),]
# gene_pa_sum1=apply(gene_pa,2,sum)
# gene_pa=gene_pa[,-which(gene_pa_sum1==0)]

#基因富集弦图
# 加载必要的包
library(circlize)
library(dplyr)
library(RColorBrewer)

### 假设您已经有以下数据：
# enrichment_matrix: 22×22的富集矩阵（基因×通路）
# ME3_hub_degree: 基因的degree值向量
# path_pva: 通路富集p值向量

### 步骤1：数据准备和颜色映射 ###

# 提取基因和通路名称
enrichment_matrix=gene_pa
genes <- rownames(enrichment_matrix)
pathways <- colnames(enrichment_matrix)

# 为通路分配颜色
pathway_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(pathways))
names(pathway_colors) <- pathways

# 为基因分配颜色（使用通路颜色或统一颜色）
gene_colors <- rep("grey80", length(genes))
names(gene_colors) <- genes

# 合并颜色向量
all_colors <- c(gene_colors, pathway_colors)

# 创建degree的渐变色标
degree_min <- min(ME3_hub_degree)
degree_max <- max(ME3_hub_degree)
degree_color_fun <- colorRamp2(c(degree_min, degree_max), c("lightblue", "darkblue"))

# 创建p值的渐变色标（-log10转换使显著性更明显）
pval_min <- -log10(max(path_pva))  # 注意：p值越小越显著，所以用-max
pval_max <- -log10(min(path_pva))
pval_color_fun <- colorRamp2(c(pval_min, pval_max), c("lightgreen", "darkgreen"))

## 步骤2：绘制弦图 ###
pdf("F:/workspace/CDH11/result_new/ME3_Gene_Pathway_Chord_Diagram.pdf", width = 12, height = 10)
par(mar = c(1, 1, 3, 1))

# 初始化圆形坐标系
circos.par(gap.after = c(rep(1, length(genes)), rep(1, length(pathways))),
           start.degree = 90)

# 绘制基础弦图
chordDiagram(enrichment_matrix,
             grid.col = all_colors,
             transparency = 0.3,
             annotationTrack = "grid",
             preAllocateTracks = list(
               list(track.height = 0.03),  # 轨道1：基因名称
               list(track.height = 0.08),  # 轨道2：degree渐变色条
               list(track.height = 0.08)   # 轨道3：p值渐变色条
             ),
             directional = 0,
             link.sort = TRUE,
             link.decreasing = TRUE,
             link.lwd = 2,
             link.lty = 1,
             big.gap = 0)

### 步骤3：添加自定义轨道 ###

# 轨道1：添加基因名称（左侧）
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% genes) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                sector.index, 
                facing = "clockwise", 
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.6,
                col = "black")
  }
}, bg.border = NA)

# 轨道2：添加degree渐变色条（左侧基因）
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% genes) {
    degree_value <- ME3_hub_degree[sector.index]
    circos.rect(CELL_META$cell.xlim[1], CELL_META$ylim[1],
                CELL_META$cell.xlim[2], CELL_META$ylim[2],
                col = degree_color_fun(degree_value), 
                border = degree_color_fun(degree_value))
  }
}, bg.border = NA)

# 轨道3：添加p值渐变色条（右侧通路）
circos.track(track.index = 3, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% pathways) {
    pval_value <- -log10(path_pva[sector.index])  # 使用-log10(p值)
    circos.rect(CELL_META$cell.xlim[1], CELL_META$ylim[1],
                CELL_META$cell.xlim[2], CELL_META$ylim[2],
                col = pval_color_fun(pval_value), 
                border = pval_color_fun(pval_value))
  }
}, bg.border = NA)

# 清除圆周坐标系
circos.clear()

### 步骤4：创建图例 ###
par(mar = c(5, 2, 5, 12), xpd = TRUE)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     axes = FALSE, xlab = "", ylab = "", 
     main = "Gene-Pathway Enrichment Chord Diagram")

# 通路图例
legend("right", 
       legend = pathways,
       fill = pathway_colors,
       border = NA,
       title = "Pathways",
       bty = "n",
       cex = 0.6,
       inset = c(-0.25, 0))

# Degree图例
legend("bottomright",
       legend = c(paste("Low (", round(degree_min, 2), ")", sep = ""),
                  paste("High (", round(degree_max, 2), ")", sep = "")),
       fill = c("lightblue", "darkblue"),
       border = "black",
       title = "Gene Degree",
       bty = "n",
       cex = 0.7)

# P值显著性图例
pval_breaks <- c(min(path_pva), median(path_pva), max(path_pva))
legend("topright",
       legend = c(paste("More Sig (", sprintf("%.2e", min(path_pva)), ")", sep = ""),
                  paste("Less Sig (", sprintf("%.2e", max(path_pva)), ")", sep = "")),
       fill = c("darkgreen", "lightgreen"),
       border = "black",
       title = "Pathway P-value",
       bty = "n",
       cex = 0.7)

# 添加说明文字
text(0.1, 0.9, "Left: Genes\nColor bar shows degree", cex = 0.7, font = 2, adj = 0)
text(0.1, 0.1, "Right: Pathways\nColor bar shows -log10(p-value)", cex = 0.7, font = 2, adj = 0)

# 重置图形参数
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
dev.off()



#############
##ME6 xiantu
#############
ME6_hub_degree=seq(30,1,-1)
names(ME6_hub_degree)=ME6_hub[,1]

go_enrichment_me6_1=go_enrichment_me6[1:30,]
pa_item=go_enrichment_me6_1$Description
path_pva=go_enrichment_me6_1$qvalue
names(path_pva)=pa_item
gene_pa=matrix(rep(0,30*length(path_pva)),nrow=30,ncol=length(path_pva))
row.names(gene_pa)=ME6_hub[,1]
colnames(gene_pa)=pa_item
for(g in ME6_hub[,1]){
    pa_ls=pa_item[grep(g,go_enrichment_me6_1$geneID)]
    if(length(pa_ls)>0){
      gene_pa[g,pa_ls]=1
    }
}
gene_pa_sum=apply(gene_pa,1,sum)
gene_pa=gene_pa[-which(gene_pa_sum==0),]
# 提取基因和通路名称
enrichment_matrix=gene_pa
genes <- rownames(enrichment_matrix)
pathways <- colnames(enrichment_matrix)

# 为通路分配颜色
pathway_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(pathways))
names(pathway_colors) <- pathways

# 为基因分配颜色（使用通路颜色或统一颜色）
gene_colors <- rep("grey80", length(genes))
names(gene_colors) <- genes

# 合并颜色向量
all_colors <- c(gene_colors, pathway_colors)

# 创建degree的渐变色标
degree_min <- min(ME6_hub_degree)
degree_max <- max(ME6_hub_degree)
degree_color_fun <- colorRamp2(c(degree_min, degree_max), c("lightblue", "darkblue"))

# 创建p值的渐变色标（-log10转换使显著性更明显）
pval_min <- -log10(max(path_pva))  # 注意：p值越小越显著，所以用-max
pval_max <- -log10(min(path_pva))
pval_color_fun <- colorRamp2(c(pval_min, pval_max), c("lightgreen", "darkgreen"))

## 步骤2：绘制弦图 ###
pdf("F:/workspace/CDH11/result_new/ME6_Gene_Pathway_Chord_Diagram.pdf", width = 12, height = 10)
par(mar = c(1, 1, 3, 1))

# 初始化圆形坐标系
circos.par(gap.after = c(rep(1, length(genes)), rep(1, length(pathways))),
           start.degree = 90)

# 绘制基础弦图
chordDiagram(enrichment_matrix,
             grid.col = all_colors,
             transparency = 0.3,
             annotationTrack = "grid",
             preAllocateTracks = list(
               list(track.height = 0.03),  # 轨道1：基因名称
               list(track.height = 0.08),  # 轨道2：degree渐变色条
               list(track.height = 0.08)   # 轨道3：p值渐变色条
             ),
             directional = 0,
             link.sort = TRUE,
             link.decreasing = TRUE,
             link.lwd = 2,
             link.lty = 1,
             big.gap = 0)

### 步骤3：添加自定义轨道 ###

# 轨道1：添加基因名称（左侧）
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% genes) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                sector.index, 
                facing = "clockwise", 
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.6,
                col = "black")
  }
}, bg.border = NA)

# 轨道2：添加degree渐变色条（左侧基因）
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% genes) {
    degree_value <- ME6_hub_degree[sector.index]
    circos.rect(CELL_META$cell.xlim[1], CELL_META$ylim[1],
                CELL_META$cell.xlim[2], CELL_META$ylim[2],
                col = degree_color_fun(degree_value), 
                border = degree_color_fun(degree_value))
  }
}, bg.border = NA)

# 轨道3：添加p值渐变色条（右侧通路）
circos.track(track.index = 3, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  if(sector.index %in% pathways) {
    pval_value <- -log10(path_pva[sector.index])  # 使用-log10(p值)
    circos.rect(CELL_META$cell.xlim[1], CELL_META$ylim[1],
                CELL_META$cell.xlim[2], CELL_META$ylim[2],
                col = pval_color_fun(pval_value), 
                border = pval_color_fun(pval_value))
  }
}, bg.border = NA)

# 清除圆周坐标系
circos.clear()

### 步骤4：创建图例 ###
par(mar = c(5, 2, 5, 12), xpd = TRUE)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     axes = FALSE, xlab = "", ylab = "", 
     main = "Gene-Pathway Enrichment Chord Diagram")

# 通路图例
legend("right", 
       legend = pathways,
       fill = pathway_colors,
       border = NA,
       title = "Pathways",
       bty = "n",
       cex = 0.6,
       inset = c(-0.25, 0))

# Degree图例
legend("bottomright",
       legend = c(paste("Low (", round(degree_min, 2), ")", sep = ""),
                  paste("High (", round(degree_max, 2), ")", sep = "")),
       fill = c("lightblue", "darkblue"),
       border = "black",
       title = "Gene Degree",
       bty = "n",
       cex = 0.7)

# P值显著性图例
pval_breaks <- c(min(path_pva), median(path_pva), max(path_pva))
legend("topright",
       legend = c(paste("More Sig (", sprintf("%.2e", min(path_pva)), ")", sep = ""),
                  paste("Less Sig (", sprintf("%.2e", max(path_pva)), ")", sep = "")),
       fill = c("darkgreen", "lightgreen"),
       border = "black",
       title = "Pathway P-value",
       bty = "n",
       cex = 0.7)

# 添加说明文字
text(0.1, 0.9, "Left: Genes\nColor bar shows degree", cex = 0.7, font = 2, adj = 0)
text(0.1, 0.1, "Right: Pathways\nColor bar shows -log10(p-value)", cex = 0.7, font = 2, adj = 0)

# 重置图形参数
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
dev.off()


### 替代方案：使用自定义函数绘制渐变条 ###

# 绘制degree渐变条（垂直）
draw_colorbar_vertical <- function(x, y, width, height, col_fun, values, title) {
  n <- length(col_fun(values))
  y_seq <- seq(y, y + height, length.out = n + 1)
  
  for(i in 1:n) {
    rect(x, y_seq[i], x + width, y_seq[i+1], 
         col = col_fun(values)[i], border = NA)
  }
  rect(x, y, x + width, y + height, border = "black")
  
  # 添加标签
  text(x + width + 0.02, y, labels = round(min(values), 1), cex = 0.6, adj = c(0, 0.5))
  text(x + width + 0.02, y + height, labels = round(max(values), 1), cex = 0.6, adj = c(0, 0.5))
  text(x + width/2, y + height/2, labels = title, cex = 0.6, srt = 90, adj = c(0.5, 0.5))
}

# 绘制p-value渐变条（水平）
draw_colorbar_horizontal <- function(x, y, width, height, col_fun, values, title) {
  n <- length(col_fun(values))
  x_seq <- seq(x, x + width, length.out = n + 1)
  
  for(i in 1:n) {
    rect(x_seq[i], y, x_seq[i+1], y + height, 
         col = col_fun(values)[i], border = NA)
  }
  rect(x, y, x + width, y + height, border = "black")
  
  # 添加标签
  text(x, y - 0.03, labels = round(min(values), 1), cex = 0.6, adj = c(0.5, 0.5))
  text(x + width, y - 0.03, labels = round(max(values), 1), cex = 0.6, adj = c(0.5, 0.5))
  text(x + width/2, y - 0.06, labels = title, cex = 0.6, adj = c(0.5, 0.5))
}

# 使用自定义函数绘制图注
pdf("F:/workspace/CDH11/colorbar_legend.pdf",6,6)
par(mar = c(5, 2, 5, 12), xpd = TRUE)
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
     axes = FALSE, xlab = "", ylab = "", 
     main = "Gene-Pathway Enrichment Chord Diagram")

# 绘制degree渐变条
draw_colorbar_vertical(0.85, 0.7, 0.03, 0.2, 
                       degree_color_fun, 
                       seq(degree_min, degree_max, length.out = 100),
                       "Degree")

# 绘制p-value渐变条
draw_colorbar_horizontal(0.3, 0.15, 0.4, 0.03,
                         pval_color_fun,
                         seq(pval_min, pval_max, length.out = 100),
                         "-log10(qvalue)")

# 添加说明文字
text(0.1, 0.9, "Left: Genes\nColor bar shows degree", cex = 0.7, font = 2, adj = 0)
text(0.1, 0.1, "Right: Pathways\nColor bar shows -log10(qvalue)", cex = 0.7, font = 2, adj = 0)

# 重置图形参数
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)
dev.off()



save.image("F:/workspace/CDH11/result_new/WGCNA_CDH11_diff_survival.rda",compress = T,ascii = F)


load("F:/workspace/CDH11/result_new/WGCNA_CDH11_diff_survival.rda")

# 如果未安装，先安装包
# install.packages("VennDiagram")
# install.packages("RColorBrewer") # 用于提供漂亮的颜色
library(circlize)
library(dplyr) # 用于数据整理
library(VennDiagram)
library(RColorBrewer)
# 创建一个列表，这是函数要求的格式
gene_list <- list(
  ME3 = me3_genes,
  ME6 = me6_genes,
  UP = up_genes$gene
)

# 设置一个颜色组合
my_colors <- brewer.pal(3, "Set2") # 从Set2调色板中取3种颜色

# 绘制韦恩图
venn_plot <- venn.diagram(
  x = gene_list,
  filename = NULL, # 不直接保存为文件，先存为图形对象
  category.names = c("ME3 Genes", "ME6 Genes", "UP Genes"), # 设置类别名称
  fill = my_colors, # 填充颜色
  alpha = 0.6, # 透明度
  cat.cex = 1.2, # 类别名称字体大小
  cex = 1.5, # 区域内部数字字体大小
  cat.dist = 0.05, # 类别名称到圆圈的距离
  margin = 0.1,
  main = "Overlap of Gene Sets", # 标题
  main.cex = 1.5
)

# 在R的图形设备中显示图片
grid.draw(venn_plot)


data.frame(me3_genes)

data.frame(me6_genes)


