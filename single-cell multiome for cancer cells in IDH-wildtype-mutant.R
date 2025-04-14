##### wt mut percentage #####
multi <- readRDS('./scMulti_7samples.rds')

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

wtallPer <- as.data.frame(prop.table(table(subset(multi, tissue == 'IDH-wildtype')$cellType)))
wtallPer <- data.frame(wtallPer, Patient = 'IDH-wildtype')
mutallPer <- as.data.frame(prop.table(table(subset(multi, tissue == 'IDH-mutant')$cellType)))
mutallPer <- data.frame(mutallPer, Patient = 'IDH-mutant')

allPer <- rbind(wtallPer, mutallPer)

colnames(allPer) <- c('CellType', 'Percentage', 'Patient')

p = ggplot(allPer,aes(x = Patient, y = Percentage, fill = CellType), base_family='Arail')+
  geom_bar(stat="identity", position = 'dodge', width = 0.8)+
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
ggsave('./allCell_IDH_proportion.pdf',p,width = 8, height = 10)

df = allPer
p = ggplot(df, aes(x = CellType, y = Percentage, fill = Patient)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "Cell Type Proportion by IDH Status", 
       x = "Cell Type", 
       y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +  # 自定义颜色
  scale_y_continuous(expand = c(0, 0.02)) +  # 调整 y 轴 0 点位置
  theme(panel.grid.major = element_blank(),   # 去掉主网格线
        panel.grid.minor = element_blank(),   # 去掉次网格线
        text = element_text(size = 14),       # 调整整体字体大小
        axis.title = element_text(size = 16), # 调整坐标轴标题字体大小
        axis.text = element_text(size = 12),  # 调整坐标轴文字字体大小
        legend.title = element_text(size = 14),  # 图例标题字体大小
        legend.text = element_text(size = 12),   # 图例内容字体大小
        plot.title = element_text(size = 18, hjust = 0.5))  # 图标题字体大小

ggsave('./IDH_cellProportion.pdf',p, height = 8, width = 10)

##### infercnv #####
all = readRDS('./scMulti_7samples.rds')

DefaultAssay(all) = 'RNA'
all$tmp = all$cellType_IDH
all <- subset(all, tmp == 'Glial Cells', invert = T)

data<-all@assays$RNA@counts

gene <- fread('/public/home/gzzhang/exp/commonDBfiles_AllStored/FurtherAnalysis/inferCNV/mart_export.txt',sep = '\t')
gene <- gene[gene$`Chromosome/scaffold name` %in% c(1:22,'X','Y'),]
gene <- gene[!duplicated(gene$`Gene name`),]
gene <- gene[gene$`Gene name`>0,]
gene <- gene[order(as.numeric(gene$`Chromosome/scaffold name`)),]
gene <- as.data.frame(gene)
rownames(gene) <- gene$`Gene name`
gene <- gene[,-1]

meta<-all@meta.data
meta<-data.frame(meta[,'tmp'],row.names = rownames(meta))
meta$meta....tmp.. <- as.character(meta$meta....tmp..)

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data,gene_order_file=gene,
                                               annotations_file=meta,
                                               ref_group_names=c('B Cells','T Cells', 'Endothelial Cells', 'Oligodendrocytes'))
infercnv_obj<-infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="cancerCells_IDH", 
                            cluster_by_groups=TRUE, 
                            denoise=TRUE,
                            output_format = "pdf",
                            HMM=TRUE)

##### chrom plot from infercnv #####
require(IdeoViz)
require(RColorBrewer) 
data(binned_multiSeries)
data(hg18_ideo)
data(binned_singleSeries)

hg38_ideo <- read.table('./inferCNV/hg38_ideogram.txt',sep = '\t')
colnames(hg38_ideo) <- colnames(hg18_ideo)

all <- readRDS('./scMulti_7samples.rds')
meta <- all@meta.data

meta$barcodes = rownames(meta)
meta <- meta[,c('cellType_IDH','barcodes')]
head(meta)

infercnv_tumor <- read.table('../inferCNV/infercnv_rnaClusters/infercnv.observations.txt', header = T)
infercnv_long <- as.data.frame(t(infercnv_tumor))
rownames(infercnv_long) <- gsub("\\.", "-", rownames(infercnv_long))
infercnv_long$Cell <- rownames(infercnv_long)

infercnv_long <- merge(infercnv_long, meta, by.x = "Cell", by.y = "row.names")
infercnv_melted <- infercnv_long %>%
  gather(Gene, Expression, -Cell, -cellType_IDH)
infercnv_melted$Expression <- as.numeric(as.character(infercnv_melted$Expression))
average_expression <- infercnv_melted %>%
  group_by(cellType_IDH, Gene) %>%
  dplyr::summarise(Average_Expression = mean(Expression, na.rm = TRUE))
average_expression_wide <- spread(average_expression, Gene, Average_Expression)
expr_cellType <- average_expression

gtf <- import('./inferCNV/hg38.refGene.gtf.gz')
gtf <- as.data.frame(gtf)

merged_data <- merge(expr_cellType, gtf, by.x = "Gene", by.y = "gene_name")
merged_data$Chromosome_Region <- paste0(merged_data$seqnames, ":", merged_data$start, "-", merged_data$end)
final_data <- merged_data[, c("cellType_IDH", "Chromosome_Region", "Average_Expression",'Gene','seqnames','start','end')]

normal_chromosomes <- c(paste0("chr", 1:22),'chrX','chrY','chrM')
final_data <- final_data[final_data$seqnames%in%normal_chromosomes,]
write.csv(final_data,'../inferCNV/final_data.csv')

### for IDH-mutant
mut_data <- final_data[final_data$cellType_IDH == 'IDH-mutant',]

create_intervals <- function(chr, max_position, interval_size = 1000000) {
  start_positions <- seq(1, max_position, by = interval_size)
  end_positions <- c((start_positions[-1] - 1), max_position)
  return(data.frame(chr = chr, start = start_positions, end = end_positions))
}
chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
all_intervals <- do.call(rbind, lapply(chromosomes, function(chr) {
  max_pos <- max(mut_data$end[mut_data$seqnames == chr])
  create_intervals(chr, max_pos)
}))

all_intervals$expression_mean <- 0

for (i in 1:nrow(all_intervals)) {
  interval <- all_intervals[i, ]
  
  genes_in_interval <- mut_data[mut_data$seqnames == interval$chr &
                                  mut_data$start >= interval$start &
                                  mut_data$end <= interval$end, ]
  
  if (nrow(genes_in_interval) > 0) {
    all_intervals$expression_mean[i] <- mean(genes_in_interval$Average_Expression)
  }
}
mut_intervals = all_intervals

### for IDH-wildtype
wt_data <- final_data[final_data$cellType_IDH == 'IDH-wildtype',]

create_intervals <- function(chr, max_position, interval_size = 1000000) {
  start_positions <- seq(1, max_position, by = interval_size)
  end_positions <- c((start_positions[-1] - 1), max_position)
  return(data.frame(chr = chr, start = start_positions, end = end_positions))
}
chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
all_intervals <- do.call(rbind, lapply(chromosomes, function(chr) {
  max_pos <- max(wt_data$end[wt_data$seqnames == chr])
  create_intervals(chr, max_pos)
}))

all_intervals$expression_mean <- 0

for (i in 1:nrow(all_intervals)) {
  interval <- all_intervals[i, ]
  
  # 找到位于该区间内的基因
  genes_in_interval <- wt_data[wt_data$seqnames == interval$chr &
                                 wt_data$start >= interval$start &
                                 wt_data$end <= interval$end, ]
  
  if (nrow(genes_in_interval) > 0) {
    all_intervals$expression_mean[i] <- mean(genes_in_interval$Average_Expression)
  }
}
wt_intervals = all_intervals

idh_intervals = cbind(mut_intervals, wt_intervals)
binned_multiSeries

gr_mut <- GRanges(seqnames = mut_intervals$chr,
                  ranges = IRanges(start = mut_intervals$start, end = mut_intervals$end),
                  IDH_mutant = mut_intervals$expression_mean)

gr_wt <- GRanges(seqnames = wt_intervals$chr,
                 ranges = IRanges(start = wt_intervals$start, end = wt_intervals$end),
                 IDH_wildtype = wt_intervals$expression_mean)

gr_idh <- GRanges(seqnames = mut_intervals$chr,
                  ranges = IRanges(start = mut_intervals$start, end = mut_intervals$end),
                  IDH_mutant = mut_intervals$expression_mean,IDH_wildtype = wt_intervals$expression_mean)

hg38_ideo

pdf('IDH_cnv_level.pdf',height = 6,width = 10)
plotOnIdeo(chrom=c('chr7','chr10'), # which chrom to plot?
           ideoTable=hg38_ideo, # ideogram name
           values_GR=gr_idh, # data goes here
           col = c("#1F77B4", "#FF7F0E"),
           value_cols=colnames(mcols(gr_idh)), # col to plot
           val_range=c(0.8,1.2), # set y-axis range
           ylab="array intensities",
           plot_title="Trendline example",
           lwd = 2)
dev.off()

##### scomatic #####
### run in local - prepare meta file
library(Seurat)
all = readRDS('./scMulti_7samples.rds')

meta = all@meta.data
meta = meta[,c('cellType_IDH','orig')]

meta$barcodes = rownames(meta)
colnames(meta) = c('Cell_type', 'orig', 'Index')

cell_counts <- table(meta$Cell_type)
valid_cell_types <- names(cell_counts[cell_counts > 500])
meta <- meta[meta$Cell_type%in%valid_cell_types,]
meta = meta[,c('Index', 'Cell_type', 'orig')]

meta$Cell_type <- gsub(" ", "_", meta$Cell_type)

write.table(meta,'meta_all.tsv', row.names = F, quote = F, col.names = T, sep = '\t')


### prepare splitting bams
python SComatic-main/scripts/SplitBam/SplitBamCellTypes.py \
--bam /public/home/gzzhang/exp/myself_gliomaMultiome_20220913_AllStored/dataInServer/bam/merged.bam \
--meta /public/home/gzzhang/exp/myself_gliomaMultiome_20220913_AllStored/dataInServer/analysis/scomatic/meta_all.tsv \
--id all \
--n_trim 5 \
--max_nM 5 \
--max_NH 1 \
--outdir ./out/bamSplit &
  
  ### collect base count information
  output_dir2=baseCount
mkdir -p $output_dir2

for bam in $(ls -d bamSplit/*bam);do

# Cell type
cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

# Temp folder
temp=$output_dir2/temp_${cell_type}
mkdir -p $temp

# Command line to submit to cluster
python ../SComatic-main/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
--ref /public/home/gzzhang/exp/commonDBfiles_AllStored/ref_genome/chr_hg38/hg38.fa \
--chrom all \
--out_folder $output_dir2 \
--min_bq 30 \
--tmp_dir $temp \
--nprocs 25

rm -rf $temp
done

### merge matrix

python ../SComatic-main/scripts/MergeCounts/MergeBaseCellCounts.py --tsv_folder baseCount \
--outfile ./baseCountMerged/all.BaseCellCounts.AllCellTypes.tsv


### detect somatic mutations
python SComatic-main/scripts/BaseCellCalling/BaseCellCalling.step1.py \
--infile ./out/baseCountMerged/all.BaseCellCounts.AllCellTypes.tsv \
--outfile ./out/variantCalling/all \
--ref /public/home/gzzhang/exp/commonDBfiles_AllStored/ref_genome/chr_hg38/hg38.fa

python SComatic-main/scripts/BaseCellCalling/BaseCellCalling.step2.py \
--infile ./out/variantCalling/all.calling.step1.tsv \
--outfile ./out/variantCalling/all \
--editing SComatic-main/RNAediting/AllEditingSites.hg38.txt \
--pon SComatic-main/PoNs/PoN.scRNAseq.hg38.tsv


### high quality snp
bedtools intersect -header \
-a all.calling.step2.tsv \
-b ../../SComatic-main/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed \
| awk '$1 ~ /^#/ || $6 == "PASS"' > all.calling.step2.pass.tsv

### high quality snp for each cell type
python ../../SComatic-main/scripts/GetCallableSites/GetAllCallableSites.py \
--infile all.calling.step1.tsv  \
--outfile ./all \
--max_cov 150 --min_cell_types 2

############## process snp data
summary = read.table('../scomatic/variantCalling/all.coverage_cell_count.per_chromosome.report.tsv', header = T)
table(summary$Cell_types)

summary_max <- summary[summary$Cov==150,]
summary_max <- summary_max[summary_max$CHROM!='chrM',]
summary_max <- summary_max[summary_max$Cell_types %in% c('IDH-mutant','IDH-wildtype'),c('Cell_types','CHROM','DP')]
summary_tmp = data.frame(Cell_types = c('IDH-mutant','IDH-wildtype'),CHROM = c('ALL','ALL'),DP = c(17368140,13533033))
summary_max <- rbind(summary_max, summary_tmp)

df <- data.frame(
  Cell_types = c(rep("IDH-mutant", 25), rep("IDH-wildtype", 25)),
  CHROM = rep(c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
                "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", 
                "chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "ALL"), 2),
  DP = c(1336971, 797554, 1113555, 798704, 469485, 648155, 602798, 343885, 638665, 
         400981, 425616, 1326393, 235151, 275648, 244590, 1711697, 780777, 959814, 
         872029, 1534590, 783997, 635868, 403755, 13260, 17368140, 
         1202393, 249874, 644045, 587102, 362954, 603097, 382220, 343825, 418645, 
         248266, 502412, 1089395, 282706, 170330, 179798, 1248756, 666143, 877626, 
         680739, 1382833, 505226, 508893, 377455, 3189, 13533033)
)

df <- df[!df$CHROM %in% c("chrX", "chrY"), ]

df$CHROM <- factor(df$CHROM, levels = c("ALL", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                                        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
                                        "chr18", "chr19", "chr20", "chr21", "chr22"))

p = ggplot(df, aes(x = CHROM, y = DP, fill = Cell_types)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +  # 调整 position_dodge 和 bar 的宽度
  theme_minimal() + 
  labs(title = "Comparison of Mutation Counts by Chromosome (Including ALL) in IDH-mutant and IDH-wildtype Cells",
       x = "Chromosome", 
       y = "Mutation Counts (DP)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +  # 自定义颜色
  scale_y_break(c(2.5e6, 10e6), scales = 0.5, space = 0.2) +  # 添加 y 轴断轴
  scale_y_continuous(labels = scales::comma) +  # 使用逗号分隔符让大数值更易读
  theme(axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())  # 去掉右侧y轴

ggsave('IDH_snp_summary.pdf',p,height = 8,width = 10)

##### snp - genes panel #####
header = read.table('./scomatic/variantCalling/header.txt',header = T)
snp = read.table('../scomatic/variantCalling/all.calling.step2.pass.tsv',header = F)
colnames(snp) = colnames(header)
snp_anno = read.table('../scomatic/variantCalling/annovar/all.variants.hg38_multianno.csv', header = T, sep = ',')
snp_final = cbind(snp, snp_anno[,1:7])

panel = read.table('../scomatic/gene_panel.txt')
snp_panel = snp_final[snp_final$Gene.refGene%in%panel$V1,]
snp_panel = snp_panel[,c(1:7,11:40)]
snp_panel = snp_panel[,c(1:11,21:30,36,37)]
snp_panel = snp_panel[,c(-14,-15,-16,-19,-20,-21)]
snp_panel <- snp_panel[snp_panel$Cell_types%in%c('IDH-mutant','IDH-wildtype'),]

split_data <- strsplit(as.character(snp_panel$IDH.mutant), "\\|")
snp_panel$IDH.mutant_allCellNumber <- sapply(split_data, function(x) x[2])
snp_panel$IDH.mutant_allCellStatus <- sapply(split_data, function(x) x[3])

split_data <- strsplit(as.character(snp_panel$IDH.wildtype), "\\|")
snp_panel$IDH.wildtype_allCellNumber <- sapply(split_data, function(x) x[2])
snp_panel$IDH.wildtype_allCellStatus <- sapply(split_data, function(x) x[3])

snp_panel$REF_seq <- snp_panel$REF
snp_panel$REF_seq = gsub('A','1',snp_panel$REF_seq)
snp_panel$REF_seq = gsub('C','2',snp_panel$REF_seq)
snp_panel$REF_seq = gsub('T','3',snp_panel$REF_seq)
snp_panel$REF_seq = gsub('G','4',snp_panel$REF_seq)
snp_panel$REF_seq <- as.numeric(snp_panel$REF_seq)

snp_panel$ALT_seq <- snp_panel$ALT
snp_panel$ALT_seq = gsub('A','1',snp_panel$ALT_seq)
snp_panel$ALT_seq = gsub('C','2',snp_panel$ALT_seq)
snp_panel$ALT_seq = gsub('T','3',snp_panel$ALT_seq)
snp_panel$ALT_seq = gsub('G','4',snp_panel$ALT_seq)
snp_panel$ALT_seq <- as.numeric(snp_panel$ALT_seq)


split_data <- strsplit(as.character(snp_panel$IDH.mutant_allCellStatus), ":")
snp_panel$IDH.mutant_REFCellNumber <- mapply(function(x, y) x[y], split_data, snp_panel$REF_seq)
split_data <- strsplit(as.character(snp_panel$IDH.mutant_allCellStatus), ":")
snp_panel$IDH.mutant_ALTCellNumber <- mapply(function(x, y) x[y], split_data, snp_panel$ALT_seq)

split_data <- strsplit(as.character(snp_panel$IDH.wildtype_allCellStatus), ":")
snp_panel$IDH.wildtype_REFCellNumber <- mapply(function(x, y) x[y], split_data, snp_panel$REF_seq)
split_data <- strsplit(as.character(snp_panel$IDH.wildtype_allCellStatus), ":")
snp_panel$IDH.wildtype_ALTCellNumber <- mapply(function(x, y) x[y], split_data, snp_panel$ALT_seq)

snp_target <- snp_panel[,c('Cell_types','Func.refGene','Gene.refGene',"IDH.mutant_allCellNumber",'IDH.wildtype_allCellNumber')]
snp_target <- cbind(snp_target, snp_panel[,24:27])
tmp = snp_target

snp_target = tmp

snp_target$IDH.mutant_REFCellNumber[is.na(snp_target$IDH.mutant_REFCellNumber)] <- 0
snp_target$IDH.mutant_ALTCellNumber[is.na(snp_target$IDH.mutant_ALTCellNumber)] <- 0

snp_target$IDH.mutant_allCellNumber <- as.numeric(snp_target$IDH.mutant_allCellNumber)
snp_target$IDH.mutant_allCellNumber[is.na(snp_target$IDH.mutant_allCellNumber)] <- 0.1
snp_target$IDH.mutant_ALTCellNumber <- as.numeric(snp_target$IDH.mutant_ALTCellNumber)
snp_target$IDH.mutant_ALTPercentage <- snp_target$IDH.mutant_ALTCellNumber/snp_target$IDH.mutant_allCellNumber

snp_target$IDH.wildtype_allCellNumber <- as.numeric(snp_target$IDH.wildtype_allCellNumber)
snp_target$IDH.wildtype_allCellNumber[is.na(snp_target$IDH.wildtype_allCellNumber)] <- 0.1
snp_target$IDH.wildtype_ALTCellNumber <- as.numeric(snp_target$IDH.wildtype_ALTCellNumber)
snp_target$IDH.wildtype_ALTPercentage <- snp_target$IDH.wildtype_ALTCellNumber/snp_target$IDH.wildtype_allCellNumber

snp_target_1 <- aggregate(IDH.mutant_ALTCellNumber ~ Gene.refGene, 
                          data = snp_target, 
                          FUN = sum, 
                          na.rm = TRUE)

snp_target_2 <- aggregate(IDH.mutant_allCellNumber ~ Gene.refGene, 
                          data = snp_target, 
                          FUN = sum, 
                          na.rm = TRUE)

snp_target_3 <- aggregate(IDH.wildtype_ALTCellNumber ~ Gene.refGene, 
                          data = snp_target, 
                          FUN = sum, 
                          na.rm = TRUE)

snp_target_4 <- aggregate(IDH.wildtype_allCellNumber ~ Gene.refGene, 
                          data = snp_target, 
                          FUN = sum, 
                          na.rm = TRUE)
snp_heatmap <- cbind(snp_target_1, snp_target_2, snp_target_3, snp_target_4)
rownames(snp_heatmap) <- snp_heatmap$Gene.refGene
snp_heatmap <- snp_heatmap[,c(-1,-3,-5,-7)]

snp_heatmap$IDH.mutant_ALTPercentage <- snp_heatmap$IDH.mutant_ALTCellNumber/snp_heatmap$IDH.mutant_allCellNumber
snp_heatmap$IDH.wildtype_ALTPercentage <- snp_heatmap$IDH.wildtype_ALTCellNumber/snp_heatmap$IDH.wildtype_allCellNumber

snp_heat_final <- snp_heatmap[,5:6]

pheatmap(snp_heat_final)

my_colors <- colorRampPalette(c("#3F1C4E", "#FFFF00"))(100)

my_colors_with_grey <- c("#D3D3D1", my_colors)

breaks <- c(0, seq(0.001, 0.35, length.out = 75), seq(0.36, 1, length.out = 26))
genes <- c('PTEN','TP53','NTRK2','FGFR2','NF1','ATRX','PIK3CA','EGFR','MTAP','PIK3R1','FGFR1','SMARCA4','IDH1','CDKN2A')
snp_heat_final <- snp_heat_final[genes,]

p = pheatmap(snp_heat_final, 
             color = my_colors_with_grey, 
             breaks = breaks, cluster_rows = F, cluster_cols = F,border_color = '#4D4D4D')

ggsave('IDH_genePanel_snp.pdf',p,height = 10, width = 2.5)


gc()

##### wt mut markers volcano #####
multi <- readRDS('../cancerCell_multi.rds')

Idents(multi) <- multi$tissue

DefaultAssay(multi) = 'SCT'
multi <- PrepSCTFindMarkers(multi)
marker <- FindAllMarkers(multi, only.pos = T, assay = 'SCT',logfc.threshold = 0, max.cells.per.ident = 2000)

### save table
DEG_up_wt <- subset(marker, avg_log2FC>0.5&p_val<0.001&cluster == 'IDH-wildtype')
DEG_up_mut <- subset(marker, avg_log2FC>0.5&p_val<0.001&cluster == 'IDH-mutant')

write.csv(DEG_up_wt, './table/wt_deg.csv')
write.csv(DEG_up_mut, './table/mut_deg.csv')
### volcano
marker <- marker %>%
  mutate(avg_log2FC = ifelse(cluster == "IDH-mutant", -avg_log2FC, avg_log2FC)) %>%
  mutate(difference = case_when(
    cluster == "IDH-wildtype" & avg_log2FC > 0.5 & p_val < 0.001 ~ 'IDH-wildtype',
    cluster == "IDH-mutant" & avg_log2FC < -0.5 & p_val < 0.001 ~ 'IDH-mutant',
    TRUE ~ 'non-significant'
  ))

marker <- marker %>%
  mutate(
    avg_log2FC = ifelse(avg_log2FC > 4, 4, ifelse(avg_log2FC < -4, -4, avg_log2FC)),
    p_val = ifelse(p_val < 1e-300, 1e-300, p_val)
  )

p = ggplot(marker, aes(x = avg_log2FC, y = -log10(p_val), color = difference)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("IDH-wildtype" = "red", "IDH-mutant" = "blue", "non-significant" = "gray")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed") +
  labs(title = "Volcano Plot of IDH-wildtype vs IDH-mutant",
       x = "log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal()

ggsave('IDH_volcano.pdf',p,height = 8,width = 10)

### add label
wt <- read.csv('./GO_wt.csv')

wt_term <- c('oxidative phosphorylation',
             'aerobic respiration',
             'ATP synthesis coupled electron transport',
             'energy derivation by oxidation of organic compounds',
             'respiratory electron transport chain',
             'ribonucleoside triphosphate metabolic process',
             'ATP metabolic process',
             'ATP biosynthetic process',
             'positive regulation of signal transduction by p53 class mediator',
             'regulation of cell-substrate adhesion',
             'epithelial cell migration',
             'epithelium migration',
             'regulation of epithelial cell migration',
             'tissue migration',
             'endothelial cell migration',
             'positive regulation of epithelial cell migration',
             'cell-substrate adhesion',
             'gliogenesis',
             'glial cell differentiation',
             'phosphatidylinositol 3-kinase/protein kinase B signal transduction',
             'regulation of stress-activated MAPK cascade',
             'neuronal stem cell population maintenance',
             'astrocyte development',
             'astrocyte differentiation',
             'regulation of mesenchymal cell proliferation',
             'mesenchymal cell differentiation'
)
go_gene_mapping <- list()

for (i in 1:nrow(wt)) {
  term <- wt$Description[i]
  genes <- unlist(strsplit(as.character(wt$geneID[i]), "/"))
  go_gene_mapping[[term]] <- genes
}

##### wt mut GO analysis #####
DEG_up_wt <- subset(marker, avg_log2FC>0.5&p_val<0.05&cluster == 'IDH-wildtype')

# BP
ego_up<-enrichGO(gene = rownames(DEG_up_wt),
                 OrgDb = 'org.Hs.eg.db',
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
res <- ego_up@result

res_filter <- res[res$pvalue<0.05,]
write.csv(res_filter, 'GO_wt.csv')

ego.simplify <- clusterProfiler::simplify(ego_up,cutoff = 0.7)
p = barplot(ego.simplify, showCategory = 50)
p = dotplot(ego.simplify, showCategory = 50)

DEG_up_mut <- subset(marker, avg_log2FC>0.5&p_val<0.05&cluster == 'IDH-mutant')

# BP
ego_up<-enrichGO(gene = rownames(DEG_up_mut),
                 OrgDb = 'org.Hs.eg.db',
                 keyType = 'SYMBOL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.05)
res <- ego_up@result

res_filter <- res[res$pvalue<0.05,]
write.csv(res_filter,'./GO_mut.csv')

ego.simplify <- clusterProfiler::simplify(ego_up,cutoff = 0.7)
p = barplot(ego.simplify, showCategory = 50)
p = dotplot(ego.simplify, showCategory = 50)

wt <- read.csv('./GO_wt.csv')

wt_term <- c('oxidative phosphorylation',
             'aerobic respiration',
             'ATP synthesis coupled electron transport',
             'energy derivation by oxidation of organic compounds',
             'respiratory electron transport chain',
             'ribonucleoside triphosphate metabolic process',
             'ATP metabolic process',
             'ATP biosynthetic process',
             'positive regulation of signal transduction by p53 class mediator',
             'regulation of cell-substrate adhesion',
             'epithelial cell migration',
             'epithelium migration',
             'regulation of epithelial cell migration',
             'tissue migration',
             'endothelial cell migration',
             'positive regulation of epithelial cell migration',
             'cell-substrate adhesion',
             'gliogenesis',
             'glial cell differentiation',
             'phosphatidylinositol 3-kinase/protein kinase B signal transduction',
             'regulation of stress-activated MAPK cascade',
             'neuronal stem cell population maintenance',
             'astrocyte development',
             'astrocyte differentiation',
             'regulation of mesenchymal cell proliferation',
             'mesenchymal cell differentiation'
)
wt <- wt[wt$Description%in%wt_term,]


mut <- read.csv('./GO_mut.csv')

mut_term <- c('modulation of chemical synaptic transmission',
              'regulation of trans-synaptic signaling',
              'synapse organization',
              'axonogenesis',
              'regulation of nervous system development',
              'dendrite development',
              'regulation of axonogenesis',
              'neuron projection extension',
              'dendrite morphogenesis',
              'regulation of neurogenesis',
              'neuron migration',
              'ERBB signaling pathway',
              'regulation of Wnt signaling pathway',
              'calcium-mediated signaling',
              'Wnt signaling pathway',
              'cell-cell signaling by wnt',
              'oligodendrocyte differentiation',
              'ERBB4 signaling pathway',
              'response to hydroperoxide',
              'oligodendrocyte development'
)
mut_term <- rev(mut_term)

rownames(mut) <- mut$Description
mut <- mut[mut_term,]

wt$updown <- 'IDH-wildtype'
mut$updown <- 'IDH-mutant'

graphTable <- rbind(wt,mut)


graphTable$logTransPvalue <- ifelse(graphTable$updown == 'IDH-wildtype',-log10(graphTable$pvalue),log10(graphTable$pvalue))

terms <- c(wt_term, mut_term)

graphTable$Description <- factor(terms, levels = rev(terms))

p = ggplot(graphTable, aes(x = logTransPvalue, y = Description, fill = updown)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#36648B', '#8B2323')) +
  xlab('t - value')+
  theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(colour = "black", size = 1, ))+
  theme(axis.title.x = element_text(face = "bold", vjust = 0.5, hjust = 0.55, size = 18))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  theme(axis.text.y = element_text(size = 10, color = "black"))+
  guides(fill=guide_legend(title=NULL, reverse = T))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.position=c(1, .5))+
  ylab('')

ggsave('./GO_barPlot.pdf',p,height = 15,width = 15)


##### survival plot #####

#### integrate micro and mRNA
clinical_data <- read.table("./public/TCGA_GBM/gbm_tcga/gbm_tcga/data_clinical_patient.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clinical_data <- clinical_data[, c("PATIENT_ID", "OS_MONTHS", "OS_STATUS")]
clinical_data$OS_STATUS <- ifelse(clinical_data$OS_STATUS == "1:DECEASED", 1, 0)
clinical_data <- clinical_data[complete.cases(clinical_data), ]
clinical_data <- clinical_data[!is.na(as.numeric(clinical_data$OS_MONTHS)),]
colnames(clinical_data) <- c("patient_id", "os_months", "os_status")



micro_data <- read.table("./public/TCGA_GBM/gbm_tcga/gbm_tcga/data_mrna_agilent_microarray.txt", 
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")

micro_data <- micro_data %>% group_by(Hugo_Symbol) %>% summarise_all(mean)
micro_data <- micro_data[!is.na(micro_data$Hugo_Symbol), ]
micro_data <- as.data.frame(micro_data)
row.names(micro_data) <- micro_data$Hugo_Symbol
micro_data <- micro_data[, -which(names(micro_data) == "Hugo_Symbol")]
micro_data <- micro_data[, -which(names(micro_data) == "Entrez_Gene_Id")]

cols_to_remove <- grep("\\.02$", colnames(micro_data))
colnames(micro_data) <- gsub('.','-',colnames(micro_data),fixed = T)
colnames(micro_data) <- sapply(colnames(micro_data), function(x) {
  return(substr(x, 1, nchar(x) - 3))
})
intersect(colnames(micro_data), clinical_data$patient_id)




mrna_data <- read.table("./public/TCGA_GBM/gbm_tcga/gbm_tcga/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
mrna_data <- mrna_data %>% group_by(Hugo_Symbol) %>% summarise_all(mean)
mrna_data <- mrna_data[!is.na(mrna_data$Hugo_Symbol), ]
mrna_data <- as.data.frame(mrna_data)
row.names(mrna_data) <- mrna_data$Hugo_Symbol
mrna_data <- mrna_data[, -which(names(mrna_data) == "Hugo_Symbol")]
mrna_data <- mrna_data[, -which(names(mrna_data) == "Entrez_Gene_Id")]
cols_to_remove <- grep("\\.02$", colnames(mrna_data))
colnames(mrna_data) <- gsub('.','-',colnames(mrna_data),fixed = T)
colnames(mrna_data) <- sapply(colnames(mrna_data), function(x) {
  return(substr(x, 1, nchar(x) - 3))
})
intersect(colnames(mrna_data), clinical_data$patient_id)

chip_data = micro_data
common_genes <- intersect(rownames(chip_data), rownames(mrna_data))
chip_data_common <- chip_data[common_genes, ]
mrna_data_common <- mrna_data[common_genes, ]
chip_data_normalized <- normalizeBetweenArrays(as.matrix(chip_data_common))
mrna_data_normalized <- normalizeBetweenArrays(as.matrix(mrna_data_common))
batch <- c(rep("chip", ncol(chip_data_normalized)), rep("mrna", ncol(mrna_data_normalized)))
combined_data <- cbind(chip_data_normalized, mrna_data_normalized)
combined_data_no_na <- combined_data[complete.cases(combined_data), ]
combat_data <- ComBat(dat = combined_data_no_na, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
integrated_data <- as.data.frame(combat_data)
head(integrated_data)
common_patients <- intersect(clinical_data$patient_id, rownames(t(integrated_data)))

merged_data <- merge(clinical_data[clinical_data$patient_id %in% common_patients, ], 
                     t(integrated_data)[common_patients, ], 
                     by.x = "patient_id", by.y = "row.names")
tmp = merged_data

gene_of_interest <- "VEGFA"
merged_data = tmp
merged_data <- merged_data[!is.na(merged_data[[gene_of_interest]]), ]
median_expression <- median(merged_data[[gene_of_interest]], na.rm = TRUE)
merged_data$expression_group <- ifelse(merged_data[[gene_of_interest]] > median_expression, "High", "Low")

surv_object <- Surv(time = as.numeric(merged_data$os_months), event = merged_data$os_status)
fit <- survfit(surv_object ~ expression_group, data = merged_data)
ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = paste("Survival Analysis of", gene_of_interest, "Expression"),
           xlab = "Time (Days)", ylab = "Survival Probability")
cox_fit <- coxph(surv_object ~ expression_group, data = merged_data)
summary(cox_fit)

### pathway score in wt
wt <- read.csv('./GO_wt.csv')

wt_term <- c('oxidative phosphorylation',
             'aerobic respiration',
             'ATP synthesis coupled electron transport',
             'energy derivation by oxidation of organic compounds',
             'respiratory electron transport chain',
             'ribonucleoside triphosphate metabolic process',
             'ATP metabolic process',
             'ATP biosynthetic process',
             'positive regulation of signal transduction by p53 class mediator',
             'regulation of cell-substrate adhesion',
             'epithelial cell migration',
             'epithelium migration',
             'regulation of epithelial cell migration',
             'tissue migration',
             'endothelial cell migration',
             'positive regulation of epithelial cell migration',
             'cell-substrate adhesion',
             'gliogenesis',
             'glial cell differentiation',
             'phosphatidylinositol 3-kinase/protein kinase B signal transduction',
             'regulation of stress-activated MAPK cascade',
             'neuronal stem cell population maintenance',
             'astrocyte development',
             'astrocyte differentiation',
             'regulation of mesenchymal cell proliferation',
             'mesenchymal cell differentiation'
)
wt = wt[wt$Description%in%wt_term,]

go_gene_mapping <- list()

for (i in 1:nrow(wt)) {
  term <- wt$Description[i]
  genes <- unlist(strsplit(as.character(wt$geneID[i]), "/"))
  go_gene_mapping[[term]] <- genes
}

names(go_gene_mapping[20]) = 'phosphatidylinositol 3-kinase_or_protein kinase B signal transduction'
for (num in c(1:19,21:26)){
  go_gene_mapping[num]
  oxphos_genes <- go_gene_mapping[[num]]
  
  valid_genes <- intersect(oxphos_genes, colnames(merged_data))
  oxphos_gene_data <- merged_data[, valid_genes, drop = FALSE]）
  merged_data$oxphos_score <- rowMeans(oxphos_gene_data, na.rm = TRUE)
  
  q1 <- quantile(merged_data$oxphos_score, 0.25, na.rm = TRUE)
  q3 <- quantile(merged_data$oxphos_score, 0.75, na.rm = TRUE)
  
  merged_data$oxphos_group <- cut(merged_data$oxphos_score,
                                  breaks = c(-Inf, q1, q3, Inf),
                                  labels = c("Low", "Intermediate", "High"))
  
  merged_data_extremes <- merged_data[merged_data$oxphos_group != "Intermediate", ]
  
  surv_object_extremes <- Surv(time = as.numeric(merged_data_extremes$os_months), 
                               event = merged_data_extremes$os_status)
  
  fit_extremes <- survfit(surv_object_extremes ~ oxphos_group, data = merged_data_extremes)

  tit = paste0('Survival Analysis of ',toTitleCase(names(go_gene_mapping[num])),' (Low vs High)')
  p = ggsurvplot(fit_extremes, data = merged_data_extremes, pval = TRUE, risk.table = TRUE,
                 title = tit,
                 xlab = "Time (days)", ylab = "Survival probability",
                 palette = c("oxphos_group=Low" = "blue", "oxphos_group=High" = "red"))
  p
  p = p$plot
  saveName = paste0('./survival/',toTitleCase(names(go_gene_mapping[num])),'.pdf')
  ggsave(saveName,p,height = 8,width = 10)
}
num = 20
go_gene_mapping[num]
oxphos_genes <- go_gene_mapping[[num]]

valid_genes <- intersect(oxphos_genes, colnames(merged_data))
oxphos_gene_data <- merged_data[, valid_genes, drop = FALSE]

merged_data$oxphos_score <- rowMeans(oxphos_gene_data, na.rm = TRUE)

q1 <- quantile(merged_data$oxphos_score, 0.25, na.rm = TRUE)
q3 <- quantile(merged_data$oxphos_score, 0.75, na.rm = TRUE)

merged_data$oxphos_group <- cut(merged_data$oxphos_score,
                                breaks = c(-Inf, q1, q3, Inf),
                                labels = c("Low", "Intermediate", "High"))

merged_data_extremes <- merged_data[merged_data$oxphos_group != "Intermediate", ]

surv_object_extremes <- Surv(time = as.numeric(merged_data_extremes$os_months), 
                             event = merged_data_extremes$os_status)

fit_extremes <- survfit(surv_object_extremes ~ oxphos_group, data = merged_data_extremes)

tit = paste0('Survival Analysis of ',toTitleCase(names(go_gene_mapping[num])),' (Low vs High)')
p = ggsurvplot(fit_extremes, data = merged_data_extremes, pval = TRUE, risk.table = TRUE,
               title = tit,
               xlab = "Time (days)", ylab = "Survival probability",
               palette = c("oxphos_group=Low" = "blue", "oxphos_group=High" = "red"))
p
p = p$plot
saveName = './survival/phosphatidylinositol 3-kinase_or_protein kinase B signal transduction.pdf'
ggsave(saveName,p,height = 8,width = 10)


### pathway score in mut
mut <- read.csv('./GO_mut.csv')
mut = mut[mut$qvalue<0.01,]

mut_term <- c('modulation of chemical synaptic transmission',
              'regulation of trans-synaptic signaling',
              'synapse organization',
              'axonogenesis',
              'regulation of nervous system development',
              'dendrite development',
              'regulation of axonogenesis',
              'neuron projection extension',
              'dendrite morphogenesis',
              'regulation of neurogenesis',
              'neuron migration',
              'ERBB signaling pathway',
              'regulation of Wnt signaling pathway',
              'calcium-mediated signaling',
              'Wnt signaling pathway',
              'cell-cell signaling by wnt',
              'oligodendrocyte differentiation',
              'ERBB4 signaling pathway',
              'response to hydroperoxide',
              'oligodendrocyte development'
)
mut = mut[mut$Description%in%mut_term,]
go_gene_mapping <- list()
for (i in 1:nrow(mut)) {
  term <- mut$Description[i]
  genes <- unlist(strsplit(as.character(mut$geneID[i]), "/"))
  go_gene_mapping[[term]] <- genes
}

for (num in 1:nrow(mut)){
  go_gene_mapping[num]
  oxphos_genes <- go_gene_mapping[[num]]
  
  valid_genes <- intersect(oxphos_genes, colnames(merged_data))
  oxphos_gene_data <- merged_data[, valid_genes, drop = FALSE]

  merged_data$oxphos_score <- rowMeans(oxphos_gene_data, na.rm = TRUE)
  
  q1 <- quantile(merged_data$oxphos_score, 0.25, na.rm = TRUE)
  q3 <- quantile(merged_data$oxphos_score, 0.75, na.rm = TRUE)
  
  merged_data$oxphos_group <- cut(merged_data$oxphos_score,
                                  breaks = c(-Inf, q1, q3, Inf),
                                  labels = c("Low", "Intermediate", "High"))
  
  merged_data_extremes <- merged_data[merged_data$oxphos_group != "Intermediate", ]
  
  surv_object_extremes <- Surv(time = as.numeric(merged_data_extremes$os_months), 
                               event = merged_data_extremes$os_status)
  
  fit_extremes <- survfit(surv_object_extremes ~ oxphos_group, data = merged_data_extremes)
  
  tit = paste0('Survival Analysis of ',toTitleCase(names(go_gene_mapping[num])),' (Low vs High)')
  p = ggsurvplot(fit_extremes, data = merged_data_extremes, pval = TRUE, risk.table = TRUE,
                 title = tit,
                 xlab = "Time (days)", ylab = "Survival probability",
                 palette = c("oxphos_group=Low" = "blue", "oxphos_group=High" = "red"))
  p
  p = p$plot
  saveName = paste0('./survival/mut/',toTitleCase(names(go_gene_mapping[num])),'.pdf')
  ggsave(saveName,p,height = 8,width = 10)
}
num = 20
go_gene_mapping[num]
oxphos_genes <- go_gene_mapping[[num]]

valid_genes <- intersect(oxphos_genes, colnames(merged_data))
oxphos_gene_data <- merged_data[, valid_genes, drop = FALSE]

merged_data$oxphos_score <- rowMeans(oxphos_gene_data, na.rm = TRUE)

q1 <- quantile(merged_data$oxphos_score, 0.25, na.rm = TRUE)
q3 <- quantile(merged_data$oxphos_score, 0.75, na.rm = TRUE)

merged_data$oxphos_group <- cut(merged_data$oxphos_score,
                                breaks = c(-Inf, q1, q3, Inf),
                                labels = c("Low", "Intermediate", "High"))

merged_data_extremes <- merged_data[merged_data$oxphos_group != "Intermediate", ]

surv_object_extremes <- Surv(time = as.numeric(merged_data_extremes$os_months), 
                             event = merged_data_extremes$os_status)

fit_extremes <- survfit(surv_object_extremes ~ oxphos_group, data = merged_data_extremes)

tit = paste0('Survival Analysis of ',toTitleCase(names(go_gene_mapping[num])),' (Low vs High)')
p = ggsurvplot(fit_extremes, data = merged_data_extremes, pval = TRUE, risk.table = TRUE,
               title = tit,
               xlab = "Time (days)", ylab = "Survival probability",
               palette = c("oxphos_group=Low" = "blue", "oxphos_group=High" = "red"))
p
p = p$plot
saveName = './survival/phosphatidylinositol 3-kinase_or_protein kinase B signal transduction.pdf'
ggsave(saveName,p,height = 8,width = 10)

##### tumor tissue wt vs mut methylation #####
mutations <- read.table("../public/TCGA_GBM/gbm_tcga/gbm_tcga/data_mutations.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")

idh1_mutations <- mutations[mutations$Hugo_Symbol == 'IDH1', c('Hugo_Symbol', 'Tumor_Sample_Barcode')]

idh1_mutations$Tumor_Sample_Barcode_new <- sapply(idh1_mutations$Tumor_Sample_Barcode, function(x) {
  return(substr(x, 1, 12))
})

clinical_data <- read.table("../public/TCGA_GBM/gbm_tcga/gbm_tcga/data_clinical_patient.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
merged_data <- merge(clinical_data, idh1_mutations, by.x = "PATIENT_ID", by.y = "Tumor_Sample_Barcode_new", all.x = TRUE)
merged_data$IDH1_status <- ifelse(is.na(merged_data$Hugo_Symbol), "Wildtype", "Mutated")

meth_data <- read.table("../public/TCGA_GBM/gbm_tcga/gbm_tcga/data_methylation_hm450.txt", 
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(meth_data) <- gsub("\\.", "-", colnames(meth_data))
sample_ids <- colnames(meth_data)[-c(1, 2)]
sample_ids_clean <- substr(sample_ids, 1, 12) 

meth_data_t <- as.data.frame(t(meth_data[, -c(1, 2)]))
colnames(meth_data_t) <- meth_data$Hugo_Symbol  
meth_data_t$Sample <- sample_ids_clean
merged_data <- merge(merged_data, meth_data_t, by.x = "PATIENT_ID", by.y = "Sample", all.x = TRUE)
head(merged_data)
table(merged_data$IDH1_status)
sum(!is.na(merged_data$IDH1))

merged_data$mean_methylation <- rowMeans(merged_data[, -c(1:41)], na.rm = TRUE)  # 假设前3列是非甲基化数据
t_test_result <- t.test(mean_methylation ~ IDH1_status, data = merged_data)
print(t_test_result)
table(merged_data$IDH1_status)
merged_data$IDH1_status = ifelse(merged_data$IDH1_status=='Wildtype','IDH-wildtype','IDH-mutant')
table(merged_data$IDH1_status)
merged_data$IDH1_status = factor(merged_data$IDH1_status, levels = c('IDH-wildtype','IDH-mutant'))

p = ggplot(merged_data, aes(x = IDH1_status, y = mean_methylation, fill = IDH1_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +  # 设置箱线宽度和颜色
  scale_fill_manual(values = c( "#FF7F0E","#1F77B4")) +  # 自定义颜色，粉色和蓝色
  labs(title = "Methylation Levels in IDH1 Mutant vs Wildtype GBM",
       x = "IDH1 Status",
       y = "Mean Methylation Level") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # 标题居中且加粗
    axis.title.x = element_text(size = 14, face = "bold"),  # X轴标签加粗
    axis.title.y = element_text(size = 14, face = "bold"),  # Y轴标签加粗
    axis.text = element_text(size = 12),  # 坐标轴刻度字体大小
    legend.position = "none",  # 去掉图例
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank()   # 移除次网格线
  ) +
  coord_cartesian(ylim = c(0.18, 0.3)) +  # 设置Y轴范围
  scale_y_continuous(breaks = seq(0.18, 0.3, by = 0.02))  # 设置Y轴间隔为0.1
p

ggsave('public_tissue_mutVSwt_methy.pdf',p,height = 6, width = 6)

##### single cell wt vs mut methylation #####
data <- read.table('../public/mutAndwt_scRNA+scMet/data/analysis_EPIC_betas.tsv', sep = '\t', header = TRUE, row.names = 1)
mut_samples <- data[, c('SM001','SM002','SM004','SM008','SM015','SM019')]
wt_samples <- data[, c('SM006', 'SM011', 'SM012', 'SM017', 'SM018')]）
mut_means <- rowMeans(mut_samples)
wt_means <- rowMeans(wt_samples)
mut_means <- mut_means[mut_means >= 0.25]
wt_means <- wt_means[wt_means >= 0.25]
plot_data <- data.frame(
  SampleGroup = rep(c("IDH-mutant", "IDH-wildtype"), c(length(mut_means), length(wt_means))),
  MethylationLevel = c(mut_means, wt_means)
)

table(plot_data$SampleGroup)
plot_data$SampleGroup = factor(plot_data$SampleGroup, levels = c('IDH-wildtype','IDH-mutant'))

p = ggplot(plot_data, aes(x = SampleGroup, y = MethylationLevel, fill = SampleGroup)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Comparison of Methylation Levels", y = "Methylation Level", x = "Sample Group") +
  scale_fill_manual(values = c("#FF7F0E","#1F77B4"))  # 将这里的颜色换成你选的颜色
p

ggsave('./public_sc_mutVSwt_methy.pdf',p, height = 8, width = 10)

##### wt vs mut TF numbers #####
data <- data.frame(
  Category = c("IDH-mutant TFs", "IDH-wildtype TFs"),
  Value = c(-21, 123)
)

p <- ggplot(data, aes(x = Category, y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_y_continuous(labels = abs) + # 去掉负号
  coord_flip() +
  theme_minimal(base_size = 15) + # 设置基础字体大小
  labs(title = "TFs",
       x = "",
       y = "Count") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)) + # 添加黑色框
  scale_fill_manual(values = c("IDH-mutant TFs" = "#1f77b4", "IDH-wildtype TFs" = "#ff7f0e")) +
  geom_text(aes(label = abs(Value)), 
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) # 添加竖线

ggsave('./enrichedTFinDiffIDHstatus.pdf',p,height = 4,width = 10)

##### TF footprinting #####

multi.glial = readRDS('../cancerCell_multi.rds')
DefaultAssay(multi.glial) = 'ATAC'

### calculate footprints
DefaultAssay(multi.glial) = 'ATAC'
table(multi.glial$tissue)
Idents(multi.glial) <- multi.glial$tissue
cell_order <- c('IDH-wildtype','IDH-mutant')
Idents(multi.glial) <- factor(Idents(multi.glial), levels = cell_order)

motif_list <- multi.glial@assays$ATAC@motifs@motif.names
data_frame <- data.frame(ID = 'test_ID', Value = 'test_Value', stringsAsFactors = FALSE)
for (i in seq_along(motif_list)) {
  data_frame <- rbind(data_frame, c(ID = names(motif_list)[i], Value = motif_list[[i]]))
}
data_frame <- data_frame[-1,]

data_frame = data_frame[data_frame$Value%in%c('CTCF','NFIC','STAT3','FOSL2','OLIG1','OLIG2','OLIG3','BHLHE23'),]

multi.glial <- Footprint(
  object = multi.glial,
  motif.name = data_frame$Value,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
CelltypeColors=c('IDH-wildtype'='#FF7F0E','IDH-mutant'='#1F77B4')

for(i in data_frame$Value){
  p <- PlotFootprint(multi.glial, features = i,idents = as.character(unique(Idents(multi.glial))),show.expected = F,label = F)& scale_colour_manual(values = CelltypeColors)
  ggsave(paste0('./footprint/new/',i,'.pdf'),p,width = 6,height = 4)
}
