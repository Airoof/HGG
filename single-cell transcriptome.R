###### umap #####
multi = readRDS('./scRNA_17samples.rds')

Idents(multi) = multi$cellType
table(Idents(multi))
multi = subset(multi, cellType == 'Neurons',invert=T)
cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes","CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_colors <- c("Glial Cells" = "#CD5C5C",                 # 印度红
                 "Oligodendrocytes" = "#FF6347",
                 "Pericytes" = "#FFD700",                   # 金色
                 "Endothelial Cells" = "#DAA520",           # 金棕色
                 "CNS-resident Myeloid Cells" = "#1E90FF",  # 道奇蓝
                 "T Cells" = "#4682B4",                     # 钢蓝
                 "Monocyte-derived Macrophages" = "#5F9EA0",# 军蓝
                 "B Cells" = "#87CEFA")                     # 浅天蓝

Idents(multi) <- multi$cellType
Idents(multi) <- factor(Idents(multi), levels = cell_order)

umap_bcc = as.data.frame(multi@reductions$umap@cell.embeddings)

p <- DimPlot(multi, reduction = "umap", cols = cell_colors)+scale_size_continuous(range = c(0.01))  + 
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

ggsave('./scRNA_UMAP.pdf',p,height = 6, width = 11)


##### dotplot #####
DefaultAssay(multi) = 'SCT'
multi <- PrepSCTFindMarkers(multi)

markers = c('LSAMP','CTNND2','LRP1B','ZIC1','PDGFRB','FBLN1','FLT1','CLDN5','VWF',
            'PLP1','PTGDS','MBP', 'P2RY12','CX3CR1','NAV3', 'HEXB','CD163','STAB1',
            'CD3D','THEMIS','CD247', 'MS4A1', 'FCRL1', 'CD79B')

cell_order <- c("Glial Cells", "Pericytes", "Endothelial Cells",
                "Oligodendrocytes", "CNS-resident Myeloid Cells","Monocyte-derived Macrophages",
                "T Cells", "B Cells")

cell_order = rev(cell_order)
multi$tmp <- multi$cellType
multi$tmp <- factor(multi$tmp, levels = cell_order)

p <- DotPlot(multi, features = markers, group.by = 'tmp', dot.scale = 25) +  
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

ggsave('./scRNA_dotplot.pdf',p, height = 9, width = 25)


##### vlnplot #####
### rna
tmp = readRDS('./scRNA_17samples.rds')

order = c('GBM-01',
          'GBM-03',
          'GBM-05',
          'GBM-06',
          'GBM-10',
          'GBM-12',
          'GBM-13',
          'GBM-15',
          'GBM-18',
          'NBT-02',
          'NBT-05',
          'G4A-01',
          'G4A-06',
          'G4A-09',
          'G4A-10',
          'NBT-01',
          'NBT-06'
)
table(tmp$tmp)
tmp$tmp = factor(tmp$tmp, levels = order)

base_colors <- brewer.pal(10, "Paired")
extra_colors <- darken(base_colors[1:7], amount = 0.2)  # 让新颜色更亮
colors <- c(base_colors, extra_colors)

p = VlnPlot(tmp, 
            features = c('nCount_RNA'), 
            group.by = "tmp", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('./scRNA_Vln_nCount.pdf', p,height = 6,width = 12 )

p = VlnPlot(tmp, 
            features = c('nFeature_RNA'), 
            group.by = "tmp", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('.scRNA_Vln_nFeature.pdf', p,height = 6,width = 12 )

p = VlnPlot(tmp, 
            features = c('percent.mt'), 
            group.by = "tmp", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('./scRNA_Vln_mt.pdf', p,height = 6,width = 12 )

### atac
tmp_combined <- readRDS(tmp_combined, './scMulti_7samples.rds')

tmp_combined$orig = factor(tmp_combined$orig, levels = c('GBM-15','GBM-18','G4A-09','G4A-10','NBT-02','NBT-05','NBT-06'))

colors <- brewer.pal(7, "Paired")

p = VlnPlot(tmp_combined, 
            features = c('nCount_ATAC'), 
            group.by = "orig", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('./scATAC_Vln_nCount.pdf', p,height = 4,width = 8 )

p = VlnPlot(tmp_combined, 
            features = c('nFeature_ATAC'), 
            group.by = "orig", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('./scATAC_Vln_nFeature.pdf', p,height = 4,width = 8 )

p = VlnPlot(tmp_combined, 
            features = c('TSS.enrichment'), 
            group.by = "orig", 
            pt.size = 0) + 
  scale_fill_manual(values = colors)
ggsave('./scATAC_Vln_tss.pdf', p,height = 4,width = 8 )


##### proportion #####
multi = readRDS('./scRNA_17samples.rds')

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

multi$cellType <- factor(multi$cellType, levels = cell_order)
multi$tissue =multi$idh

wtallPer <- as.data.frame(prop.table(table(subset(multi, tissue == 'IDH-wildtype')$cellType)))
wtallPer <- data.frame(wtallPer, Patient = 'IDH-wildtype')

mutallPer <- as.data.frame(prop.table(table(subset(multi, tissue == 'IDH-mutant')$cellType)))
mutallPer <- data.frame(mutallPer, Patient = 'IDH-mutant')

allPer <- rbind(wtallPer, mutallPer)
colnames(allPer) <- c('CellType', 'Percentage', 'Patient')
allPer$Patient = factor(allPer$Patient, levels = c('IDH-wildtype','IDH-mutant'))

p = ggplot(allPer,aes(x = Patient, y = Percentage, fill = CellType), base_family='Arail')+
  geom_bar(stat="identity", position = 'stack', width = 0.8)+
  theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1, ))+
  xlab('')+
  guides(fill = guide_legend(title = NULL))+
  theme(axis.title.y = element_text(face = "bold", colour = 'black', size = 17))+
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold"))+
  theme(legend.text=element_text( size = 14))+
  scale_fill_manual(values = cell_colors)
p
ggsave('./scRNA_proportion.pdf',p,width = 8, height = 10)


