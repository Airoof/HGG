multi <- readRDS('./scMulti_7samples.rds')

##### rna umap #####
DefaultAssay(multi) = 'SCT'

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_colors <- c("Glial Cells" = "#CD5C5C",                 # 印度红
                 "Oligodendrocytes" = "#FF6347",            # 番茄红
                 "Pericytes" = "#FFD700",                   # 金色
                 "Endothelial Cells" = "#DAA520",           # 金棕色
                 "CNS-resident Myeloid Cells" = "#1E90FF",  # 道奇蓝
                 "T Cells" = "#4682B4",                     # 钢蓝
                 "Monocyte-derived Macrophages" = "#5F9EA0",# 军蓝
                 "B Cells" = "#87CEFA")                     # 浅天蓝

Idents(multi) <- multi$cellType
Idents(multi) <- factor(Idents(multi), levels = cell_order)

p <- DimPlot(multi, reduction = "umap", cols = cell_colors, pt.size = 0.01) + 
  theme_minimal()+
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.title = element_blank(), # 去掉坐标轴标题
    axis.text = element_blank(),  # 去掉坐标轴刻度
    axis.ticks = element_blank()  # 去掉坐标轴刻度线
  )+
  
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +
  guides(colour = guide_legend(override.aes = list(size = 4.5)))

ggsave('./allcell_UMAP.pdf',p,height = 6, width = 10)

##### atac umap #####
DefaultAssay(multi) = 'ATAC'

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_colors <- c("Glial Cells" = "#CD5C5C",                 # 印度红
                 "Oligodendrocytes" = "#FF6347",            # 番茄红
                 "Pericytes" = "#FFD700",                   # 金色
                 "Endothelial Cells" = "#DAA520",           # 金棕色
                 "CNS-resident Myeloid Cells" = "#1E90FF",  # 道奇蓝
                 "T Cells" = "#4682B4",                     # 钢蓝
                 "Monocyte-derived Macrophages" = "#5F9EA0",# 军蓝
                 "B Cells" = "#87CEFA")                     # 浅天蓝

Idents(multi) <- multi$cellType
Idents(multi) <- factor(Idents(multi), levels = cell_order)

p <- DimPlot(multi, reduction = "atac.umap", cols = cell_colors, pt.size = 0.01) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.title = element_blank(), # 去掉坐标轴标题
    axis.text = element_blank(),  # 去掉坐标轴刻度
    axis.ticks = element_blank(),  # 去掉坐标轴刻度线
    legend.position = "none" # 去掉图例
  )

ggsave('./allcell_atac_UMAP.pdf',p,height = 6, width = 7)

##### cluster umap #####
DefaultAssay(multi) = 'SCT'

p <- DimPlot(multi, group.by = "rna.seuratClusters", pt.size = 0.01,label = T,label.size = 5) + 
  theme_minimal()+
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.title = element_blank(), # 去掉坐标轴标题
    axis.text = element_blank(),  # 去掉坐标轴刻度
    axis.ticks = element_blank()  # 去掉坐标轴刻度线
  )+
  
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +
  guides(colour = guide_legend(override.aes = list(size = 4.5)))

ggsave('./allcell_rnacluster_UMAP.pdf',p,height = 8, width = 10)


p = DimPlot(multi, group.by = 'atac.seuratClusters',pt.size = 0.01,label=T, reduction = 'atac.umap',label.size = 5) + 
  theme_minimal()+
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.title = element_blank(), # 去掉坐标轴标题
    axis.text = element_blank(),  # 去掉坐标轴刻度
    axis.ticks = element_blank()  # 去掉坐标轴刻度线
  )+
  
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +
  guides(colour = guide_legend(override.aes = list(size = 4.5)))

ggsave('./allcell_ataccluster_UMAP.pdf',p,height = 8, width = 10)



##### infercnv #####
data<-multi@assays$RNA@counts

gene <- fread('/public/home/gzzhang/exp/commonDBfiles_AllStored/FurtherAnalysis/inferCNV/mart_export.txt',sep = '\t')
gene <- gene[gene$`Chromosome/scaffold name` %in% c(1:22,'X','Y'),]
gene <- gene[!duplicated(gene$`Gene name`),]
gene <- gene[gene$`Gene name`>0,]
gene <- gene[order(as.numeric(gene$`Chromosome/scaffold name`)),]
gene <- as.data.frame(gene)
rownames(gene) <- gene$`Gene name`
gene <- gene[,-1]

meta<-multi@meta.data
meta<-data.frame(meta[,'cellType'],row.names = rownames(meta))
meta$meta....cellType.. <- as.character(meta$meta....cellType..)

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data,gene_order_file=gene,
                                               annotations_file=meta,
                                               ref_group_names=c('B Cells','T Cells', 'Endothelial Cells', 'Oligodendrocytes'))
infercnv_obj<-infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="allCells_final", 
                            cluster_by_groups=TRUE, 
                            denoise=TRUE,
                            output_format = "pdf",
                            HMM=TRUE)

##### dot plot #####
### all type markers
DefaultAssay(multi) = 'SCT'
markers = c('LSAMP','CTNND2','LRP1B','PDGFRB','FBLN1','COL6A2','FLT1','CLDN5','VWF',
            'PLP1','PTGDS','MBP', 'P2RY12','CX3CR1','NAV3', 'HEXB','CD163','STAB1',
            'CD3D','THEMIS','CD247', 'MS4A1', 'FCRL1', 'CD79B')

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_order = rev(cell_order)
multi$tmp <- multi$cellType
multi$tmp <- factor(multi$tmp, levels = cell_order)

p <- DotPlot(multi, features = markers, group.by = 'tmp', dot.scale = 12) +  
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 16, face = "italic"), # 设置 x 轴字体大小为 20
    axis.text.y = element_text(size = 16),                        # 设置 y 轴字体大小为 20
    legend.text = element_text(size = 16),                        # 设置图例字体大小为 20
    legend.title = element_text(size = 16)                        # 设置图例标题字体大小为 20
  ) +
  scale_color_gradientn(
    colors = c("#4A75A3", "#DCEBE0", "#8B3626"),  # 蓝色, 白色, 红色
    values = scales::rescale(c(-1, 0, 2.5))       # -1对应蓝色, 0对应白色, 2.5对应红色
  )

p

ggsave('allCell_dotplot.pdf',p, height = 6, width = 16)

##### coverage plot #####
DefaultAssay(multi) = 'ATAC'
cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_colors <- c("Glial Cells" = "#CD5C5C",                 # 印度红
                 "Pericytes" = "#FFD700",                   # 金色
                 "Endothelial Cells" = "#DAA520",           # 金棕色
                 "Oligodendrocytes" = "#FF6347",            # 番茄红
                 "CNS-resident Myeloid Cells" = "#1E90FF",  # 道奇蓝
                 "Monocyte-derived Macrophages" = "#5F9EA0",# 军蓝
                 "T Cells" = "#4682B4",                     # 钢蓝
                 "B Cells" = "#87CEFA")                     # 浅天蓝

Idents(multi) <- multi$cellType
Idents(multi) <- factor(Idents(multi), levels = cell_order)
markers = c('LSAMP','CTNND2','LRP1B','PDGFRB','FBLN1','COL6A2','FLT1','CLDN5','VWF',
            'PLP1','PTGDS','MBP', 'P2RY12','CX3CR1','NAV3', 'HEXB','CD163','STAB1',
            'CD3D','THEMIS','CD247', 'MS4A1', 'FCRL1', 'CD79B')

for(i in 1:length(markers)){
  p = CoveragePlot(
    object = multi,
    region = markers[i],
    features = markers[i],
    expression.assay = "SCT",
    idents = unique(Idents(multi)),
    extend.upstream = 10000,
    extend.downstream = 10000,peaks = F
  )& scale_fill_manual(values = as.vector(cell_colors))
  
  ggsave(paste0('./coverage/',markers[i],'.pdf'),p,height = 8,width = 20)
}

###### cellType Percentage #####
gc()

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_colors <- c("Glial Cells" = "#CD5C5C",                 # 印度红
                 "Oligodendrocytes" = "#FF6347",            # 番茄红
                 "Pericytes" = "#FFD700",                   # 金色
                 "Endothelial Cells" = "#DAA520",           # 金棕色
                 "CNS-resident Myeloid Cells" = "#1E90FF",  # 道奇蓝
                 "T Cells" = "#4682B4",                     # 钢蓝
                 "Monocyte-derived Macrophages" = "#5F9EA0",# 军蓝
                 "B Cells" = "#87CEFA")                     # 浅天蓝

multi$cellType <- factor(multi$cellType, levels = cell_order)

df <- as.data.frame(prop.table(table(multi$cellType)))
colnames(df) <- c('CellType', 'Percentage')
df$Percentage <- df$Percentage * 100  # 转换为百分比
df$color <- cell_colors[as.character(df$CellType)]  
df <- df[order(-df$Percentage),]
labs <- paste0(df$CellType, " \n(", round(df$Percentage, 2), "%)")
df <- df[order(df$Percentage),]

pdf('./allCell_proportion.pdf',height = 8, width = 10)

pie(df$Percentage, labels = NA, init.angle = -90, col = df$color, border = "#3D3D3D")
legend("topright", legend = paste0(rev(df$CellType), " (", round(rev(df$Percentage), 2), "%)"), 
       fill = rev(df$color), cex = 1.5, inset = c(0,0))  # inset 可以控制图例在图形中的偏移

dev.off()

##### TF analysis - chromvar #####
motif_name = unlist(multi@assays$ATAC@motifs@motif.names)

DefaultAssay(multi) <- 'chromvar'

Idents(multi) = multi$cellType
cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")
Idents(multi) = factor(Idents(multi), levels = cell_order)
tf_diff <- FindAllMarkers(
  object = multi,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
tf_diff$TF_Name <- motif_name[tf_diff$gene]

tf_diff %>%
  group_by(cluster) %>%
  dplyr::filter(avg_diff > 2) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

multi@assays$chromvar@data[1:5,1:5]
multi@assays$chromvar@scale.data = multi@assays$chromvar@data

tmp = multi
expr_matrix <- tmp@assays$chromvar@data
rownames(expr_matrix) <- motif_name[rownames(expr_matrix)]

tmp@assays$chromvar@data <- expr_matrix
multi = tmp

### heatmap
tmp = subset(multi, downsample = 2000)
DoHeatmap(tmp, features = top5$TF_Name, slot = 'data')

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")
Idents(multi) = factor(Idents(multi), levels = cell_order)

averageExp_markers <- AverageExpression(multi, assays = 'chromvar',slot = 'data',
                                        features = top5$TF_Name)
averageExp_markers = averageExp_markers$chromvar
mat <- averageExp_markers

log_mat = mat
log_mat <- log_mat[c(1,2,5,7,8,10,11,13,14,15,16,17,20,21,24,22,23,25,26,30,32,31,33,34),]

p = pheatmap(log_mat, 
             scale = 'row',
             cluster_rows = FALSE, 
             cluster_cols = FALSE, 
             border_color = '#3D3D3D',
             fontsize = 14,        # 整体字体大小
             fontsize_row = 14,    # 行标签的字体大小
             fontsize_col = 14)    # 列标签的字体大小

ggsave('allCell_TF_heatmap.pdf',p,height = 16,width = 6)

### featureplot
p = FeaturePlot(multi, features = c('NFIB'),ncol = 1, cols = c('#FFFFF0','#8B0000')) + scale_colour_gradientn(colors = c('#FFFFF0', '#FFE4C4', '#FF6347', '#8B0000'), 
                                                                                                              values = c(0, 0.2, 0.4, 1)) 
p
ggsave('allCell_NFIB_dotplot.pdf',p, height = 10,width = 10)

p = FeaturePlot(multi, features = c('NFIC'),ncol = 1, cols = c('#FFFFF0','#8B0000')) + scale_colour_gradientn(colors = c('#FFFFF0', '#FFE4C4', '#FF6347', '#8B0000'),values = c(0, 0.2, 0.4, 1)) 
p
ggsave('allCell_NFIC_dotplot.pdf',p, height = 10,width = 10)

p = FeaturePlot(multi, features = c('RFX1'),ncol = 1, cols = c('#FFFFF0','#8B0000')) + scale_colour_gradientn(colors = c('#FFFFF0', '#FFE4C4', '#FF6347', '#8B0000'), values = c(0, 0.2, 0.4, 1)) 
p
ggsave('allCell_RFX1_dotplot.pdf',p, height = 10,width = 10)

p = FeaturePlot(multi, features = c('RFX2'),ncol = 1, cols = c('#FFFFF0','#8B0000')) + scale_colour_gradientn(colors = c('#FFFFF0', '#FFE4C4', '#FF6347', '#8B0000'), values = c(0, 0.2, 0.4, 1)) 
p
ggsave('allCell_RFX2_dotplot.pdf',p, height = 10,width = 10)

##### TF footprinting #####
### change fragments file path
DefaultAssay(multi.glial) = 'ATAC'

### calculate footprints
DefaultAssay(multi.glial) = 'ATAC'
Idents(multi.glial) <- multi.glial$cellType
cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")
Idents(multi.glial) <- factor(Idents(multi.glial), levels = cell_order)

motif_list <- multi.glial@assays$ATAC@motifs@motif.names
data_frame <- data.frame(ID = 'test_ID', Value = 'test_Value', stringsAsFactors = FALSE)
for (i in seq_along(motif_list)) {
  data_frame <- rbind(data_frame, c(ID = names(motif_list)[i], Value = motif_list[[i]]))
}
data_frame <- data_frame[-1,]

# gather the footprinting information for sets of motifs
multi.glial <- Footprint(
  object = multi.glial,
  motif.name = data_frame$Value,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
for(i in data_frame$Value){
  p <- PlotFootprint(multi.glial, features = i,idents = as.character(unique(Idents(multi.glial))))
  ggsave(paste0('allCells/',i,'.png'),p,width = 8,height = 4)
}
