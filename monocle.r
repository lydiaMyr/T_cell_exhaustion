library("Matrix")
library("Seurat")
library("monocle3")
library("ggplot2")
load("F:/workspace/CDH11/seurat_stro_rename.rda")
# 预处理
cds_raw <- new_cell_data_set(
  expression_data = GetAssayData(seurat_stro_rename, assay = "RNA", slot = "counts"),
  cell_metadata = seurat_stro_rename@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(seurat_stro_rename), 
                            row.names = rownames(seurat_stro_rename))
)




# 重新预处理和降维，增加邻居数量
# cds <- preprocess_cds(cds, num_dim = 50)  # 增加PCA维度
# cds <- reduce_dimension(
#   cds,
#   preprocess_method = "PCA",
#   umap.n_neighbors = 30,  # 增加邻居数量，增强连接性
#   umap.min_dist = 0.1,    # 减小最小距离，让细胞更紧凑
#   umap.metric = "cosine"  # 尝试不同的距离度量
# )
# cds <- preprocess_cds(cds, 
#                      method = "PCA", 
#                      num_dim = 50, 
#                      norm_method = "log", 
#                      pseudo_count = 1)
# cds <- reduce_dimension(cds, 
#                        preprocess_method = "PCA", 
#                        reduction_method = "UMAP",  # 必须指定为 UMAP
#                        umap.min_dist = 0.1)

cds <- cluster_cells(cds_raw,
                    resolution = 1e-5,  # 低分辨率聚焦主轨迹
                    cluster_method = "leiden",
                    random_seed = 42)
# 学习轨迹图结构
cds <- learn_graph(cds, 
                  use_partition = FALSE, 
                  close_loop = FALSE,
                  learn_graph_control = list(ncenter = 100))

#根据umap图选择根节点
# 绘制基础轨迹图（未指定伪时间）
plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size = 1)
# 交互式选择根节点
cds <- order_cells(cds, reduction_method = "UMAP")

# 操作流程：
# 1. 运行上述代码后，R 会弹出 UMAP 图窗口
# 2. 点击轨迹起点附近的节点（如最早发育状态的细胞群）
# 3. 按 [Esc] 或关闭窗口完成选择
# 选择轨迹根节点（根据M1/M2标记表达）
# 轨迹图（按伪时间着色）
pdf("F:/workspace/TME_project_file/单细胞/pseudotime_plot.pdf",6.5,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1) +
  scale_color_viridis_c()
dev.off()


# 轨迹图（按伪时间着色）
plot_cells(cds,
           color_cells_by = "celltype",
           label_cell_groups = TRUE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1)

plot_genes_in_pseudotime(cds["CDH11"], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
pdf("F:/workspace/CDH11/CDH11_pseudotime_exp.pdf",5,3)
plot_genes_in_pseudotime(cds["CDH11"], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
# 绘制并应用颜色
pdf("F:/workspace/CDH11/result_new/cluster_pseudotime_exp.pdf",8,6)

plot_cells(cds, 
           color_cells_by = "celltype",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1) 
dev.off()

#耗竭和纤维化基因的表达
# T细胞耗竭相关基因
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TIM3", 
                      "TOX", "ENTPD1", "CD39", "BATF", "NR4A1", "NR4A2", 
                      "TCF7", "TCF1", "EOMES", "TBX21", "T-BET")

# 纤维化相关基因
fibrosis_genes <- c("ACTA2", "α-SMA", "COL1A1", "COL3A1", "FN1", "TGFB1", 
                   "TGFB2", "TGFB3", "MMP2", "MMP9", "TIMP1", "TIMP2", 
                   "POSTN", "FAP", "PDGFRA", "PDGFRB", "VIM", "S100A4")

# 过滤掉在数据中不存在的基因
available_genes <- rowData(cds)$gene_short_name
exhaustion_genes <- exhaustion_genes[exhaustion_genes %in% available_genes]
fibrosis_genes <- fibrosis_genes[fibrosis_genes %in% available_genes]

print(paste("找到的耗竭基因:", length(exhaustion_genes)))
print(paste("找到的纤维化基因:", length(fibrosis_genes)))

# 选择关键基因进行趋势分析
key_genes1 <- exhaustion_genes[1:6]
key_genes2 <- exhaustion_genes[7:12]
key_genes3 <- exhaustion_genes[13:17]

key_genes4 <- fibrosis_genes[1:6]
key_genes5 <- fibrosis_genes[7:12]
key_genes6 <- fibrosis_genes[13:17]


# 绘制基因表达随拟时序的变化
pdf("F:/workspace/CDH11/result_new/Exhausted_monocle.pdf",10,8)
plot_genes_in_pseudotime(cds[key_genes1, ],
  color_cells_by = "celltype",  # 按细胞类型着色
  min_expr = 0.1,
  nrow = 2,
  ncol = 3
) + theme(legend.position = "right")
plot_genes_in_pseudotime(cds[key_genes2, ],
  color_cells_by = "celltype",  # 按细胞类型着色
  min_expr = 0.1,
  nrow = 2,
  ncol = 3
) + theme(legend.position = "right")
plot_genes_in_pseudotime(cds[key_genes3, ],
  color_cells_by = "celltype",  # 按细胞类型着色
  min_expr = 0.1,
  nrow = 2,
  ncol = 3
) + theme(legend.position = "right")

dev.off()
pdf("F:/workspace/CDH11/result_new/Fibro_monocle.pdf",8,4)
plot_genes_in_pseudotime(cds[c("CDH11",key_genes4[1:3]), ],
  color_cells_by = "celltype",  # 按细胞类型着色
  min_expr = 0.1,
  nrow = 2,
  ncol = 2
) + theme(legend.position = "right")
dev.off()
# plot_genes_in_pseudotime(cds[key_genes5, ],
#   color_cells_by = "celltype",  # 按细胞类型着色
#   min_expr = 0.1,
#   nrow = 2,
#   ncol = 3
# ) + theme(legend.position = "right")
# plot_genes_in_pseudotime(cds[key_genes6, ],
#   color_cells_by = "celltype",  # 按细胞类型着色
#   min_expr = 0.1,
#   nrow = 2,
#   ncol = 3
# ) + theme(legend.position = "right")




m1_genes <- c("CD86", "CXCL9","IRF5")
m2_genes <- c("CD163", "CCL18","IRF4")
pdf("F:/workspace/TME_project_file/单细胞/M1_M2_pseudotime_exp.pdf",4,5)
plot_genes_in_pseudotime(cds[m1_genes, ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[m2_genes, ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)  
plot_genes_in_pseudotime(cds[c("IGJ","PTPRS","GZMB")], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
                     
pdf("F:/workspace/TME_project_file/单细胞/cellchat_gene_pseudotime.pdf",10,10)
plot_genes_in_pseudotime(cds[genes[1:5], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[6:10], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[11:15], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[16:20], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[21:28], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
#患者信息补充
cell_info = read.table("F:/workspace/TME_project_file/单细胞/single_cell_info",header=T,sep="\t")
sample_cell=table(cell_info[,5])
cell_sub = Idents(macro_cell_new)
library(dplyr)
cell_id_macro = cell_info%>%dplyr::filter(title%in%names(cell_sub))
sample_cell_count = table(cell_info[,5])
sample_response = unique(cell_info[,c(5,6)])
row.names(sample_response) = sample_response[,1]
new_df = data.frame(id=cell_id_macro[,"title"],cell_type=cell_sub[cell_id_macro[,"title"]],source=cell_id_macro[,5])
sub_count = table(new_df$cell_type,new_df$source)
sub_ratio = data.frame(apply(sub_count,1,function(x) x/sample_cell))

sub_ratio$response = sample_response[row.names(sub_ratio),2]
sub_ratio$id = row.names(sub_ratio)
sub_ratio$point = sapply(sub_ratio$id,function(x) as.vector(unlist(strsplit(x,"_")))[[1]])

library(ggpubr)
pdf("F:/workspace/TME_project_file/单细胞/sample_cell_boxplot.pdf",3.5,4)
for(i in seq(1,4)){
    df1 = data.frame(cell=sub_ratio[,i],group=sub_ratio[,5])
    df2 = data.frame(cell=sub_ratio[,i],group=sub_ratio[,7])
    p1 = ggplot(df1,aes(x=group,y=cell))+geom_violin(aes(fill=group,color=group))+geom_boxplot(width=0.1)+scale_color_brewer(palette="Pastel1")+scale_fill_brewer(palette="Pastel1",)+stat_compare_means(method="wilcox.test",label = "p.format",color="black")+theme_classic()+labs(x="",y="",title=colnames(sub_ratio)[i])+theme(axis.text = element_text(color="black",size=12),legend.position="none")
    p2 = ggplot(df2,aes(x=group,y=cell))+geom_violin(aes(fill=group,color=group))+geom_boxplot(width=0.1)+scale_color_brewer(palette="Pastel1")+scale_fill_brewer(palette="Pastel1",)+stat_compare_means(method="wilcox.test",label = "p.format",color="black")+theme_classic()+labs(x="",y="",title=colnames(sub_ratio)[i])+theme(axis.text = element_text(color="black",size=12),legend.position="none")
    print(p1)
    print(p2)
}
dev.off()



# 查看当前CDS对象中有多少个轨迹组件
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
print("组件数量：")
length(unique(partitions(cds)))