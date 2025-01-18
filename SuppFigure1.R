#Import dependencies
library(tidyverse)
library(matrixStats)
library(ggstar)
library(ggrepel)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/SuppFigure1/"

#Supp. Figure 1A
##Read in CNV data (log2-transformed with a pseudocount of 1) (21Q3) (downloaded on 8/16/2024) and produce lung and breast cancer copy number dataset
cnv <- read.csv(paste0(input_dir, "CCLE_gene_cn_21Q3.csv"), row.names = 1, check.names = F)
colnames(cnv) <- gsub(" (.*)$", "", colnames(cnv))
cnv$DepMap_ID <- rownames(cnv)
cnv <- cnv[, c("DepMap_ID", "CKS1B")]
cnv$absolute_copies <- ((2^cnv$CKS1B) - 1)*2
meta <- read.csv(paste0(input_dir, "sample_info_21Q3.csv"))
lines <- meta[meta$primary_disease %in% c("Breast Cancer", "Lung Cancer"),]
cks1b <- cnv[cnv$DepMap_ID %in% lines$DepMap_ID,]
cks1b$CellLine <- meta$stripped_cell_line_name[match(cks1b$DepMap_ID, meta$DepMap_ID)]
###Smooth number of copies
#Ruleset: i) if number between previous_integer.8 and integer.2, it becomes integer.
#ii) if number between integer.2 and integer.8, it becomes integer.5.
cks1b$decimal <- cks1b$absolute_copies-floor(cks1b$absolute_copies)
cks1b$smoothed_copies <- ifelse(cks1b$decimal < 0.2, floor(cks1b$absolute_copies), 
                                ifelse(cks1b$decimal >= 0.2 & cks1b$decimal <= 0.8, floor(cks1b$absolute_copies) + 0.5, 
                                       ceiling(cks1b$absolute_copies)))
###Write CKS1B copy number
write.csv(cks1b, paste0(output_dir, "cks1b_copynumber_pancancer.csv"), row.names = F)

##Import pan-cancer dependency data 
##24Q2: Chronos-processed CRISPR data and metadata (downloaded on 8/16/2024)
crispr <- read.csv(paste0(input_dir, "CRISPRGeneEffect_24Q2.csv"), row.names = 1, check.names = F)
meta <- read.csv(paste0(input_dir, "Model_24Q2.csv"))
lines <- meta[meta$OncotreeLineage %in% c("Lung", "Breast"),]
nrow(lines)
##351 lung and breast cancer cell lines in the project
sum(rownames(crispr) %in% lines$ModelID)
##176 lung and breast cancer cell lines with dependency data
dataset.ach <- crispr[rownames(crispr) %in% lines$ModelID,]
rownames(dataset.ach) <- lines$StrippedCellLineName[match(rownames(dataset.ach), lines$ModelID)]
colnames(dataset.ach) <- gsub(" (.*)$", "", colnames(dataset.ach))
dataset.ach <- data.frame(t(dataset.ach), check.names = F)
nrow(dataset.ach)
##18,443 genes
write.csv(dataset.ach, paste0(output_dir, "dataset.ach_pancancer.csv"), row.names = T)

##Read in RNA-seq RPKM data (18Q2) and compute expression effect size
rpkm <- read.delim(paste0(input_dir, "CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct"), skip=2, check.names = F)
ccle_annotations <- read.delim(paste0(input_dir, "Cell_lines_annotations_20181226.txt"))
rpkm_lb <- rpkm[, c("Name", "Description", colnames(rpkm)[colnames(rpkm) %in% ccle_annotations$CCLE_ID[ccle_annotations$Site_Primary %in% c("lung", "breast")]])]
colnames(rpkm_lb) <- gsub("_LUNG|_BREAST", "", colnames(rpkm_lb))
rpkm_lb <- rpkm_lb[, colnames(rpkm_lb) %in% c(cks1b$CellLine, "Name", "Description")]
amplified_lb_rpkm <- colnames(rpkm_lb)[colnames(rpkm_lb) %in% cks1b$CellLine[cks1b$smoothed_copies > 2]]
length(amplified_lb_rpkm)
##151
nonamplified_lb_rpkm <- colnames(rpkm_lb)[colnames(rpkm_lb) %in% cks1b$CellLine[cks1b$smoothed_copies == 2]]
length(nonamplified_lb_rpkm)
##55
rpkm_lb$median.amp <- rowMedians(as.matrix(rpkm_lb[,amplified_lb_rpkm]), na.rm = T)
rpkm_lb$median.non <- rowMedians(as.matrix(rpkm_lb[,nonamplified_lb_rpkm]), na.rm = T)
rpkm_lb$exp.effect.size <- rpkm_lb$median.amp - rpkm_lb$median.non
rpkm_lb$ExpFC <- log2(rpkm_lb$median.amp/rpkm_lb$median.non)
final.just.exp.effect.size.lb <- rpkm_lb[,c("Description", "exp.effect.size", "ExpFC")]

##Call hits
amplified_lb <- as.character(cks1b$CellLine)[cks1b$smoothed_copies > 2]
amplified_lb <- amplified_lb[amplified_lb %in% colnames(dataset.ach)]
length(amplified_lb)
##104
nonamplified_lb <- as.character(cks1b$CellLine)[cks1b$smoothed_copies == 2]
nonamplified_lb <- nonamplified_lb[nonamplified_lb %in% colnames(dataset.ach)]
length(nonamplified_lb)
##43
dataset.ach$median.amp <- rowMedians(as.matrix(dataset.ach[, amplified_lb]), na.rm = T)
dataset.ach$median.non <- rowMedians(as.matrix(dataset.ach[, nonamplified_lb]), na.rm = T)
dataset.ach$effect.size <- dataset.ach$median.amp - dataset.ach$median.non
sum(is.na(dataset.ach$median.amp))
##0
sum(is.na(dataset.ach$effect.size))
##0
dataset.ach$exp.effect.size <- final.just.exp.effect.size.lb$exp.effect.size[match(rownames(dataset.ach), final.just.exp.effect.size.lb$Description)]
dataset.ach$ExpFC <- final.just.exp.effect.size.lb$ExpFC[match(rownames(dataset.ach), final.just.exp.effect.size.lb$Description)]
pc_ach.neg.more.neg <- subset(dataset.ach, dataset.ach$median.amp < 0 & dataset.ach$effect.size <0,)
##Only consider as hits genes upregulated in amplified cell lines
pc_ach.neg.more.neg$hit <- ifelse(pc_ach.neg.more.neg$exp.effect.size > 2, 1, 0)
##Annotate hits based on whether they came up in the MM screen analysis or not
out_high <- read.csv("results/Figure1/DependencyHits_chr1q.csv")
pc_ach.neg.more.neg$color_lab <- ifelse(pc_ach.neg.more.neg$ExpFC > 0, 1, 0)
pc_ach.neg.more.neg <- pc_ach.neg.more.neg[!is.na(pc_ach.neg.more.neg$ExpFC),]
pc_ach.neg.more.neg$confirmed <- ifelse(rownames(pc_ach.neg.more.neg) %in% out_high$Gene[out_high$DatasetCount > 1] & pc_ach.neg.more.neg$hit == "1", 1, 0)
pc_ach.neg.more.neg$shape_lab <- ifelse(pc_ach.neg.more.neg$confirmed == "1", 2, ifelse(pc_ach.neg.more.neg$confirmed == "0" & pc_ach.neg.more.neg$ExpFC > 0, 1, 0))
pc_ach.neg.more.neg$lab <- ifelse(rownames(pc_ach.neg.more.neg) %in% c("MCL1", "PIP5K1A"), as.character(rownames(pc_ach.neg.more.neg)), NA)
##Only keep genes located on chr1q
chr1q_genes <- read.csv(paste0(input_dir, "chr1q_genes.csv"), header = T)
chr1q_genes <- chr1q_genes[!duplicated(chr1q_genes$hg38.kgXref.geneSymbol), ]
chr1q_genes <- chr1q_genes[-grep("ENSG", chr1q_genes$hg38.kgXref.geneSymbol),]
pc_ach.neg.more.neg$Start <- chr1q_genes$hg38.knownGene.txStart[match(rownames(pc_ach.neg.more.neg), chr1q_genes$hg38.kgXref.geneSymbol)]
pc_ach.neg.more.neg <- pc_ach.neg.more.neg[!is.na(pc_ach.neg.more.neg$Start), ]
pc_ach.neg.more.neg$ExpFC <- ifelse(is.infinite(pc_ach.neg.more.neg$ExpFC) & pc_ach.neg.more.neg$ExpFC >0, max(pc_ach.neg.more.neg$ExpFC[is.finite(pc_ach.neg.more.neg$ExpFC)]),
                                    ifelse(is.infinite(pc_ach.neg.more.neg$ExpFC) & pc_ach.neg.more.neg$ExpFC <0, min(pc_ach.neg.more.neg$ExpFC[is.finite(pc_ach.neg.more.neg$ExpFC)]),
                                           as.numeric(as.character(pc_ach.neg.more.neg$ExpFC))))
write.csv(pc_ach.neg.more.neg, paste0(output_dir, "SuppFig1_SourceData.csv"), row.names = F)
##Plot
png(paste0(output_dir, "Dependency_PanCancerv2.png"), res=300, unit="in", width=6, height=5)
ggplot(pc_ach.neg.more.neg) + 
  geom_star(aes(ExpFC, -effect.size, fill = factor(shape_lab), size = factor(shape_lab), alpha = factor(shape_lab), starshape = factor(shape_lab)), color = "black") +
  scale_alpha_manual(values = c(0.7, 0.7, 1), guide = "none") +
  scale_starshape_manual(values = c(15, 15, 1), guide = "none") +
  scale_fill_manual(values = c("steelblue", "tomato2", "lemonchiffon"), guide = "none") +
  scale_size_manual(values = c(2, 2, 4), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.5) + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14), 
        legend.position = c(0.2, 0.8),
        plot.title = element_text(size = 16, hjust = 0.5)) +
  xlab(expression("Expression Log"[2]*" fold-change")) + 
  ylab("Differential Dependency") +
  xlim(min(pc_ach.neg.more.neg$ExpFC)*1.1, max(pc_ach.neg.more.neg$ExpFC)*1.1) +
  ylim(0, max(-pc_ach.neg.more.neg$effect.size)*1.1) +
  ggtitle("DepMap Lung and Breast Ca: 1q+ vs non-1q+") + 
  geom_text_repel(aes(ExpFC, -effect.size, label = lab), segment.alpha = 0.5, nudge_y = 0.01, nudge_x = 0.15, max.overlaps = Inf, fontface = "italic")
dev.off()

