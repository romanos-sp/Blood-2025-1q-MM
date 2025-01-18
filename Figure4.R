#Import dependencies
library(tidyverse)
library(ggpubr)
library(viridis)
library(cowplot)
library(matrixStats)
library(ggrepel)
library(qusage)
library(org.Hs.eg.db)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/Figure4/"

#Import single-cell RNA-seq data from patient samples
sc_df <- read.csv(paste0(input_dir, "amp1q_non1q_cell_obs.csv"), row.names = 1)

#Figure 4A
##Plot sample UMAP
a <- ggplot(sc_df) + 
  geom_point(aes(UMAP1, UMAP2, color = pid), alpha = 0.5, size = 1, stroke = 0, show.legend = F) +
  scale_color_manual(values = c("lightblue", "steelblue", "purple", "tomato2", "pink", "orange"), name = "Sample") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank())
b <- ggplot(sc_df) + 
  geom_point(aes(UMAP1, UMAP2, color = pid), alpha = 1) +
  scale_color_manual(values = c("lightblue", "steelblue", "purple", "tomato2", "pink", "orange"), name = "Sample") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_blank(), axis.ticks = element_blank())
b_leg <- get_legend(b)
png(paste0(output_dir, "Sample_UMAP.SC.png"), res=300, unit="in", width=6, height=5)
plot_grid(a, b_leg, ncol = 2, rel_widths = c(0.6, 0.1), align = "hv", axis = "tlbr")
dev.off()

#Figure 4B
#Plot 1q+-status UMAP
sc_df$amp1q_genotype <- factor(sc_df$amp1q_genotype, levels = c("Non-1q", "Amp1q"))
png(paste0(output_dir, "Amp1q_UMAP.SC.png"), res=300, unit="in", width=5, height=5)
ggplot(sc_df) + 
  geom_point(aes(UMAP1, UMAP2, color = amp1q_genotype), alpha = 0.5, size = 1, stroke = 0) +
  scale_color_manual(values = c("navajowhite1", "orangered4"), name = "Status", labels = c("Non-1q+", "1q+")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.position = c(0.2, 0.8))
dev.off()

#Figure 4C
#Read in CNV data and assign status
make_seg <- function(sampleid, joint, seg){
  joint$seglabel <- paste0(joint$seg, "_", joint$cnv_state)
  joint_clone <- aggregate(joint$p_cnv, list(joint$clone, joint$seglabel), median)
  colnames(joint_clone) <- c("Clone", "SegLabel", "Median")
  joint_clone$state <- gsub("^.*_", "", joint_clone$SegLabel)
  joint_clone$seg <- gsub("_.*$", "", joint_clone$SegLabel)
  joint_clone$state <- ifelse(joint_clone$Median <= 0.35, "neu", as.character(joint_clone$state))
  joint_clone$clone_seg <- paste0(joint_clone$Clone, "_", joint_clone$seg)
  final_seg <- seg[, c("CHROM", "seg", "seg_start", "seg_end", "cnv_state_post")]
  final_seg$clone <- 1
  final_seg2 <- final_seg
  final_seg2$clone <- 2
  final_seg <- rbind(final_seg, final_seg2)
  final_seg$clone_seg <- paste0(final_seg$clone, "_", final_seg$seg)
  final_seg$state <- joint_clone$state[match(final_seg$clone_seg, joint_clone$clone_seg)]
  final_seg$state[is.na(final_seg$state)] <- "neu"
  final_seg$ymin <- ifelse(final_seg$clone == "1", 0, 1)
  final_seg$ymax <- ifelse(final_seg$clone == "1", 1, 2)
  final_seg$state_n <- ifelse(final_seg$state == "neu", 0, ifelse(final_seg$state %in% c("bamp", "amp"), 1, ifelse(final_seg$state %in% c("del", "bdel"), -1, NA)))
  final_seg$Sample <- sampleid
  return(final_seg)
}

##P1_BM
clones <- read.delim(paste0(input_dir, "CNV/P1_BM_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P1_BM_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P1_BM_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Both clones are malignant and will be kept for visualization
P1_BM <- make_seg("P1_BM", joint = joint, seg = seg)

##P2_BM
clones <- read.delim(paste0(input_dir, "CNV/P2_BM_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P2_BM_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P2_BM_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Clone 1 is normal cells and will be discarded for visualization
joint <- joint[joint$clone %in% c("2", "3", "4"),]
##Clones 3 & 4 correspond to the 1q+ clones and will be merged for visualization
joint$clone <- ifelse(joint$clone %in% c("3", "4"), "1", "2")
P2_BM <- make_seg("P2_BM", joint = joint, seg = seg)

##P3_BM
clones <- read.delim(paste0(input_dir, "CNV/P3_BM_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P3_BM_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P3_BM_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Clone 1 is normal cells and will be discarded for visualization
joint <- joint[joint$clone %in% c("2", "3", "4"),]
##Clones 2 & 3 correspond to the non-1q+ clones and will be merged for visualization
joint$clone <- ifelse(joint$clone == "4", "1", "2")
P3_BM <- make_seg("P3_BM", joint = joint, seg = seg)

##P4_BM
clones <- read.delim(paste0(input_dir, "CNV/P4_BM_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P4_BM_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P4_BM_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Clone 1 is normal cells and will be discarded for visualization
joint <- joint[joint$clone %in% c("2", "3", "4"),]
##Clones 2 & 3 correspond to the non-1q+ clones and will be merged for visualization
joint$clone <- ifelse(joint$clone == "4", "1", "2")
P4_BM <- make_seg("P4_BM", joint = joint, seg = seg)

##P4_PB
clones <- read.delim(paste0(input_dir, "CNV/P4_PB_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P4_PB_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P4_PB_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Clone 1 is normal cells and will be discarded for visualization
joint <- joint[joint$clone %in% c("2", "3", "4"),]
##Clones 2 & 3 correspond to the non-1q+ clones and will be merged for visualization
joint$clone <- ifelse(joint$clone == "4", "1", "2")
P4_PB <- make_seg("P4_PB", joint = joint, seg = seg)

##P5_PB
clones <- read.delim(paste0(input_dir, "CNV/P5_PB_clone_post_2.tsv"))
joint <- read.delim(paste0(input_dir, "CNV/P5_PB_joint_post_2.tsv"))
seg <- read.delim(paste0(input_dir, "CNV/P5_PB_segs_consensus_2.tsv"))
joint$clone <- clones$clone_opt[match(joint$cell, clones$cell)]
##Clone 1 is normal cells and will be discarded for visualization
joint <- joint[joint$clone %in% c("2", "3"),]
##Clone 3 corresponds to the 1q+ clone and will be renamed clone "1" for visualization
joint$clone <- ifelse(joint$clone == "3", "1", "2")
P5_PB <- make_seg("P5_PB", joint = joint, seg = seg)

##Compile all segs
final_seg <- rbind(P1_BM, P2_BM, P3_BM, P4_BM, P4_PB, P5_PB)
final_seg$clone <- factor(final_seg$clone, levels = c("1", "2"))
final_seg$CHROM <- factor(final_seg$CHROM, levels = unique(final_seg$CHROM[order(final_seg$CHROM, decreasing = F)]))

png(paste0(output_dir, "PatientSC_CNV.png"), res=300, unit="in", width=11, height=5)
ggplot(final_seg, 
       aes(xmin = seg_start, xmax = seg_end, ymin = ymin, ymax = ymax)) +
  geom_rect(aes(fill = state_n), colour = NA, size=0) +
  geom_hline(yintercept = 1, size=0.3, color="darkgrey", alpha=0.5) +
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", mid = "#F7F7F7", midpoint = 0, na.value = "#F7F7F7") +
  labs(fill = "Copy number") +
  facet_grid(rows = vars(Sample),
             cols = vars(CHROM),
             scales = "free", 
             space = "free",
             switch = "y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill="#F7F7F7"),
        panel.spacing.y = unit(0.2, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(size=.3, fill = NA),
        strip.text.y.left = element_text(angle = 0, size = 7),
        strip.text.x.top = element_text(angle = 0, size = 7),
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0.02,0)) +
  scale_x_continuous(expand = c(0.05,0))
dev.off()

#Figure 4D
##Plot MCL1 violin plot
sc_df$pid <- factor(sc_df$pid, levels = c("P1_BM", "P2_BM", "P3_BM", "P4_BM", "P4_PB", "P5_PB"))
sc_df$amp1q_genotype2 <- ifelse(sc_df$amp1q_genotype %in% "Amp1q", "1q+",
                                ifelse(sc_df$amp1q_genotype %in% "Non-1q", "Non-1q+", as.character(sc_df$amp1q_genotype)))
sc_df$amp1q_genotype2 <- factor(sc_df$amp1q_genotype2, levels = c("Non-1q+", "1q+"))
stats <- compare_means(MCL1 ~ amp1q_genotype, group.by = "pid", data=sc_df, p.adjust.method = "BH")
stats$pid <- factor(stats$pid, levels = levels(sc_df$pid))
stats <- stats[order(stats$pid, decreasing = F),]
stats$xmin <- as.numeric(stats$pid) - 0.25
stats$xmax <- as.numeric(stats$pid) + 0.25
stats$padj <- ifelse(stats$p.adj < 2e-16, "q<2e-16", paste0("q=", signif(stats$p.adj, 2)))
png(paste0(output_dir, "MCL1_Violin.SC.png"), res=300, unit="in", width=6, height=5)
ggplot(sc_df) +
  geom_violin(aes(pid, MCL1, fill = amp1q_genotype), alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(pid, MCL1, fill = amp1q_genotype), alpha = 0.7, outlier.size = -1, show.legend = F, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("navajowhite1", "orangered4"), name = "Status", labels = c("Non-1q+", "1q+")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.x = element_text(color = "black", size = 14), 
        axis.text.y = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14)) +
  xlab("") +
  ylab("Scaled MCL1 expression (z-score)") + 
  geom_bracket(aes(xmin = xmin, xmax = xmax, label = padj), data = stats[stats$p.adj < 0.05,], y.position = 4.5, size = 0.4, label.size = 4)
dev.off()

#Figure 4E
##Plot PI3K UMAP
png(paste0(output_dir, "PI3K_UMAP.SC.png"), res=300, unit="in", width=5, height=5)
ggplot(sc_df) + 
  geom_point(aes(UMAP1, UMAP2, color = PI3K_AKT_mTOR_pathway), alpha = 0.5, size = 1, stroke = 0) +
  scale_color_viridis(name = "PI3K-AKT-mTOR") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = c(0.2, 0.8))
dev.off()

#Figure 4G
##Import results from single-cell RNA-seq DE between subclones with and without 1q+ 
sc <- read.csv(paste0(input_dir, "amp1q_single_cell_DE_n3_targets.csv"))

##Identify DepMap dependency hits that are differentially expressed between subclones
cnv <- read.csv(paste0(input_dir, "Final_chr1q_CN.csv"))
##DepMap CRISPR dataset
dataset.ach <- read.csv("results/Figure1/dataset.ach.csv", row.names = 1, check.names = F)
dataset.ach.neg.more.neg <- dataset.ach
ncol(dataset.ach)
##17
sum(!colnames(dataset.ach) %in% cnv$CellLine)
##5
colnames(dataset.ach)[!colnames(dataset.ach) %in% cnv$CellLine]
##[1] "L363"   "OCIMY7" "INA6"   "KMS34"  "KMM1"  
ach.amp <- colnames(dataset.ach)[colnames(dataset.ach) %in% cnv$CellLine[cnv$cnv > 2]]
length(ach.amp)
##11
ach.amp
##[1] "OPM2"     "LP1"      "KMS20"    "KMS27"    "KMS26"    "JJN3"     "KMS11"    "MM1S"    
##[9] "RPMI8226" "EJM"      "AMO1"  
ach.non <- colnames(dataset.ach)[colnames(dataset.ach) %in% cnv$CellLine[cnv$cnv == 2]]
length(ach.non)
##1
ach.non
##[1] "SKMM2"
dataset.ach.neg.more.neg$median.amp <- rowMedians(as.matrix(dataset.ach.neg.more.neg[, ach.amp]), na.rm = T)
dataset.ach.neg.more.neg$median.non <- rowMedians(as.matrix(dataset.ach.neg.more.neg[, ach.non]), na.rm = T)
dataset.ach.neg.more.neg$effect.size <- dataset.ach.neg.more.neg$median.amp - dataset.ach.neg.more.neg$median.non
sum(is.na(dataset.ach.neg.more.neg$median.amp))
##512
sum(is.na(dataset.ach.neg.more.neg$effect.size))
##656
dataset.ach.neg.more.neg <- dataset.ach.neg.more.neg[!is.na(dataset.ach.neg.more.neg$effect.size) & !is.na(dataset.ach.neg.more.neg$median.amp),]
dataset.ach.neg.more.neg <- dataset.ach.neg.more.neg[dataset.ach.neg.more.neg$median.amp < 0 & dataset.ach.neg.more.neg$effect.size < 0,]
intermediate <- dataset.ach.neg.more.neg
intermediate$Description <- rownames(intermediate)
gw.achilles.hits.sc <- as.character(intermediate$Description[intermediate$Description %in% sc$names])
write.table(gw.achilles.hits.sc, paste0(output_dir, "gw.achilles.hits.sc.txt"), sep = "\t")

##DepMap RNAi
dataset.rnai <- read.csv("results/Figure1/dataset.rnai.csv", row.names = 1, check.names = F)
dataset.rnai.neg.more.neg <- dataset.rnai
ncol(dataset.rnai)
##15
sum(!colnames(dataset.rnai) %in% cnv$CellLine)
##7
colnames(dataset.rnai)[!colnames(dataset.rnai) %in% cnv$CellLine]
##[1] "INA6"    "KMM1"    "KMS28BM" "KMS34"   "L363"    "OCIMY7"  "OPM1" 
rnai.amp <- colnames(dataset.rnai)[colnames(dataset.rnai) %in% cnv$CellLine[cnv$cnv > 2]]
length(rnai.amp)
##7
rnai.amp
##[1] "AMO1"     "JJN3"     "KMS12BM"  "KMS27"    "MM1S"     "OPM2"     "RPMI8226"
rnai.non <- colnames(dataset.rnai)[colnames(dataset.rnai) %in% cnv$CellLine[cnv$cnv == 2]]
length(rnai.non)
##1
rnai.non
##[1] "SKMM2"
dataset.rnai.neg.more.neg$median.amp <- rowMedians(as.matrix(dataset.rnai.neg.more.neg[,rnai.amp]), na.rm = T)
dataset.rnai.neg.more.neg$median.non <- rowMedians(as.matrix(dataset.rnai.neg.more.neg[,rnai.non]), na.rm = T)
dataset.rnai.neg.more.neg$effect.size <- dataset.rnai.neg.more.neg$median.amp - dataset.rnai.neg.more.neg$median.non
sum(is.na(dataset.rnai.neg.more.neg$median.amp))
##58
sum(is.na(dataset.rnai.neg.more.neg$effect.size))
##6173
dataset.rnai.neg.more.neg <- dataset.rnai.neg.more.neg[!is.na(dataset.rnai.neg.more.neg$effect.size) & !is.na(dataset.rnai.neg.more.neg$median.amp),]
dataset.rnai.neg.more.neg <- dataset.rnai.neg.more.neg[dataset.rnai.neg.more.neg$median.amp < 0 & dataset.rnai.neg.more.neg$effect.size < 0,]
intermediate <- dataset.rnai.neg.more.neg
intermediate$Description <- rownames(intermediate)
gw.rnai.hits.sc <- as.character(intermediate$Description[intermediate$Description %in% sc$names])
write.table(gw.rnai.hits.sc, paste0(output_dir, "gw.rnai.hits.sc.txt"), sep = "\t")

##Compile hits across the 2 genome-wide screens
vector.of.hits.final.sc <- c(unique(gw.rnai.hits.sc), unique(gw.achilles.hits.sc))
frequency.table.final.sc <- data.frame(table(vector.of.hits.final.sc))
colnames(frequency.table.final.sc) <- c("Hit", "Count")
write.csv(frequency.table.final.sc, paste0(output_dir, "GenomeWide_FrequencyTable.sc.csv"), row.names = F)
table(frequency.table.final.sc$Count)
##  1   2 
##255  61 
##Annotate genes that hit in both datasets with genomic coordinates
hits <- as.character(frequency.table.final.sc$Hit)[frequency.table.final.sc$Count>1]
res3 <- select(org.Hs.eg.db, keys=as.vector(hits), 
               columns = c("SYMBOL", "ENTREZID", "CHR"), keytype= "SYMBOL")
res3 <- res3[res3$CHR %in% c(1:22),]
res3$chr <- paste0("chr", res3$CHR)
res3 <- res3[!duplicated(res3$SYMBOL),]

##Visualize hits
frequency.table.final.sc$Ach <- ifelse(frequency.table.final.sc$Hit %in% gw.achilles.hits.sc, 1, 0)
frequency.table.final.sc$Rnai <- ifelse(frequency.table.final.sc$Hit %in% gw.rnai.hits.sc, 1, 0)
frequency.table.final = frequency.table.final.sc
out <- data.frame(Gene = NA, Dependency_Amp1q = NA, Dependency_Non1q = NA)
for(i in 1:nrow(frequency.table.final)){
  if(i%%100==0){
    print(paste0(i, "..."))
  }
  out[i, "Gene"] <- as.character(frequency.table.final$Hit)[i]
  if(frequency.table.final$Count[i]==1 & frequency.table.final$Ach[i]==1){
    out[i, "Dependency_Amp1q"] <- dataset.ach.neg.more.neg$median.amp[rownames(dataset.ach.neg.more.neg) %in% frequency.table.final$Hit[i]]
    out[i, "Dependency_Non1q"] <- dataset.ach.neg.more.neg$median.non[rownames(dataset.ach.neg.more.neg) %in% frequency.table.final$Hit[i]]
  } else if(frequency.table.final$Count[i]==1 & frequency.table.final$Rnai[i]==1){
    out[i, "Dependency_Amp1q"] <- dataset.rnai.neg.more.neg$median.amp[rownames(dataset.rnai.neg.more.neg) %in% frequency.table.final$Hit[i]]
    out[i, "Dependency_Non1q"] <- dataset.rnai.neg.more.neg$median.non[rownames(dataset.rnai.neg.more.neg) %in% frequency.table.final$Hit[i]]
  } else if(frequency.table.final$Count[i]==2 & frequency.table.final$Ach[i]==1 & frequency.table.final$Rnai[i]==1){
    diff1 <- dataset.ach.neg.more.neg$effect.size[rownames(dataset.ach.neg.more.neg) %in% frequency.table.final$Hit[i]]
    diff2 <- dataset.rnai.neg.more.neg$effect.size[rownames(dataset.rnai.neg.more.neg) %in% frequency.table.final$Hit[i]]
    if(abs(diff1) > abs(diff2)){
      out[i, "Dependency_Amp1q"] <- dataset.ach.neg.more.neg$median.amp[rownames(dataset.ach.neg.more.neg) %in% frequency.table.final$Hit[i]]
      out[i, "Dependency_Non1q"] <- dataset.ach.neg.more.neg$median.non[rownames(dataset.ach.neg.more.neg) %in% frequency.table.final$Hit[i]]
    } else {
      out[i, "Dependency_Amp1q"] <- dataset.rnai.neg.more.neg$median.amp[rownames(dataset.rnai.neg.more.neg) %in% frequency.table.final$Hit[i]]
      out[i, "Dependency_Non1q"] <- dataset.rnai.neg.more.neg$median.non[rownames(dataset.rnai.neg.more.neg) %in% frequency.table.final$Hit[i]]
    }
  } 
}
out$DatasetCount <- frequency.table.final.sc$Count[match(as.character(out$Gene), as.character(frequency.table.final.sc$Hit))]
out$DiffDep <- out$Dependency_Amp1q-out$Dependency_Non1q
write.csv(out, paste0(output_dir, "DependencyHit_Matrix.sc.csv"), row.names = F)

##Annotate hit table with differential dependency information
res3$LFC <- sc$median_lfc[match(res3$SYMBOL, sc$names)]
res3$DiffDep <- out$DiffDep[match(res3$SYMBOL, out$Gene)]
res3$DatasetCount <- out$DatasetCount[match(res3$SYMBOL, out$Gene)]
res3$color_lab <- factor(ifelse(res3$LFC > 0, 1, 0))

##Keep hits not located on chr1
plotdf <- res3[res3$chr != "chr1",]
nrow(plotdf)
##34
plotdf$lab <- ifelse(plotdf$DatasetCount == "2", as.character(plotdf$SYMBOL), NA)
write.csv(plotdf, paste0(output_dir, "Dependency_NonChr1Hits.csv"), row.names = F)

##Plot
png(paste0(output_dir, "Dependency_NonChr1.SC.png"), res=300, unit="in", width=6, height=5)
ggplot(plotdf) +
  geom_jitter(aes(LFC, -DiffDep, size = factor(DatasetCount), fill = color_lab), shape = 21, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.title = element_text(size = 14), 
        legend.position = c(0.2, 0.8)) +
  xlab(expression("Median within-tumor Log"[2]*" fold-change")) +
  ylab("Differential Dependency") +
  scale_size_manual(values = c(2, 3), name = "# of datasets") +
  scale_fill_manual(values = c("steelblue", "tomato2"), guide = "none") +
  geom_text_repel(aes(LFC, -DiffDep, label = lab), segment.alpha = 0.3, max.overlaps = Inf, fontface = "italic")
dev.off()

#Figure 4H
#Import HALLMARK pathways and perform hypergeometric tests for hits in 2 or more datasets
##HALLMARK pathways downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp on 4/21/2023
hallmark <- read.gmt(paste0(input_dir, "h.all.v2023.1.Hs.symbols.gmt"))
urn <- max(nrow(dataset.ach), nrow(dataset.rnai))
urn
##18,443
dataset_hits <- plotdf$SYMBOL
for(i in seq_along(hallmark)){
  print(i)
  pathway <- hallmark[[i]]
  if(i==1){
    out <- data.frame(Pathway = names(hallmark)[i], GenesInGeneSet = length(pathway), NumberOfGenesInOverlap = length(intersect(dataset_hits, pathway)), GenesInOverlap = paste(intersect(dataset_hits, pathway), collapse = ", "), GenesConsidered = urn, GenesInList = length(dataset_hits), Hypergeometric_p =  phyper(q = length(intersect(dataset_hits, pathway)) - 1, m = length(pathway), n = urn - length(pathway), k = length(dataset_hits), lower.tail = F))
  } else {
    out <- rbind(out, data.frame(Pathway = names(hallmark)[i], GenesInGeneSet = length(pathway), NumberOfGenesInOverlap = length(intersect(dataset_hits, pathway)), GenesInOverlap = paste(intersect(dataset_hits, pathway), collapse = ", "), GenesConsidered = urn, GenesInList = length(dataset_hits), Hypergeometric_p =  phyper(q = length(intersect(dataset_hits, pathway)) - 1, m = length(pathway), n = urn - length(pathway), k = length(dataset_hits), lower.tail = F)))
  }
}
out$FDR <- p.adjust(out$Hypergeometric_p, method = "BH")
write.csv(out, paste0(output_dir, "DependencyHits_nonchr1q_morethan2datasets_HALLMARK_Hypergeometric.csv"), row.names = F)

#Visualize results from pathway analysis (FDR < 0.1)
out <- read.csv(paste0(output_dir, "DependencyHits_nonchr1q_morethan2datasets_HALLMARK_Hypergeometric.csv"))
out_sig <- out[out$FDR < 0.1, ]
out_sig$Pathway <- factor(out_sig$Pathway, levels = out_sig$Pathway[order(out_sig$FDR, decreasing = T)])
out_sig$GeneSetName <- gsub("^HALLMARK_", "", out_sig$Pathway)
png(paste0(output_dir, "DependencyHits_nonchr1q_PathwayEnrichment.png"), res=300, unit="in", width=7, height=2)
ggplot(out_sig, aes(y=Pathway)) + 
  geom_bar(aes(x= GenesInGeneSet), position = "dodge", stat = "identity", fill = "skyblue1") +
  geom_bar(aes(x= NumberOfGenesInOverlap), position = "dodge", stat = "identity", fill = "darkorchid3") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(color = "black", size = 12), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 14)) + 
  xlab("Genes in Overlap/Genes in Gene Set") +
  ylab("") +
  geom_text(aes(x = GenesInGeneSet, label = GenesInGeneSet), hjust=1.2, color="black", size=5) +
  geom_text(aes(x = NumberOfGenesInOverlap, label = NumberOfGenesInOverlap), hjust=1.4, color="white", size=5) +
  geom_text(aes(x = GenesInGeneSet + 5, label = paste0("q=", signif(FDR, digits = 2))), hjust=0, color="black", size=3.5, fontface = "italic") +
  scale_y_discrete(breaks = out_sig$Pathway, labels = out_sig$GeneSetName, expand = c(0.04, 0.04)) +
  scale_x_log10(breaks = c(0, 10, 100, 200), limits = c(1, 500), expand = c(0.01, 0.01))
dev.off()

