#Import dependencies
library(tidyverse)
library(ggrepel)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/SuppFigure2/"

#Import IC50 data
ic50 <- read.csv("data/Drugs/IC50_Data_Final.csv", check.names = F)
rownames(ic50) <- ic50$Drug
ic50$Drug <- NULL
ic50 <- data.frame(t(ic50))
ic50$CellLine <- gsub("-", "", rownames(ic50))

#Read in RNA-seq RPKM data (18Q2) and annotate with MCL1 RNA levels
rpkm <- read.delim(paste0(input_dir, "CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct"), skip=2, check.names = F)
ccle_annotations <- read.delim(paste0(input_dir, "Cell_lines_annotations_20181226.txt"))
rpkm_mm <- rpkm[, c("Name", "Description", colnames(rpkm)[colnames(rpkm) %in% ccle_annotations$CCLE_ID[ccle_annotations$Hist_Subtype1=="plasma_cell_myeloma"]])]
###Remove KE97, HUNS1 which have been misannotated as MM, and KMS18/KMS21BM which presented genotyping issues in the DepMap dataset
rpkm_mm <- rpkm_mm[, -grep("KE97|HUNS1|KMS18|KMS21BM", colnames(rpkm_mm))]
colnames(rpkm_mm) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", colnames(rpkm_mm))
##Remove cell lines without chr1q copy number information
rpkm_mm_lines <- colnames(rpkm_mm)[!colnames(rpkm_mm) %in% c("Name", "Description")] 
cnv <- read.csv(paste0(input_dir, "Final_chr1q_CN.csv"))
rpkm_mm_lines[!rpkm_mm_lines %in% cnv$CellLine]
##[1]"KMM1"        "KMS28BM"     "KMS34"       "L363"     "U266B1"   
rpkm_mm <- rpkm_mm[, colnames(rpkm_mm) %in% c("Name", "Description", rpkm_mm_lines[rpkm_mm_lines %in% cnv$CellLine]
)]
all(colnames(rpkm_mm)[!colnames(rpkm_mm) %in% c("Name", "Description")] %in% cnv$CellLine)
##TRUE
gep <- rpkm_mm[!duplicated(rpkm_mm$Description),]
rownames(gep) <- gep$Description
gep$Name <- NULL
gep$Description <- NULL
gep <- data.frame(t(gep))
ic50$MCL1 <- gep$MCL1[match(rownames(ic50), rownames(gep))]

#Import BH3 profiling data and annotate with MCL1 dependency
bl <- read.csv(paste0(input_dir, "Baseline_BH3_Final.csv"))
bl_long <- gather(bl, "Sample", "CytoC", -Peptide)
bl_long$Cell_Line <- gsub("_.", "", bl_long$Sample)
bl_long$CL_Peptide <- paste0(bl_long$Cell_Line, ":", bl_long$Peptide)
bl_long_agg <- data.frame(aggregate(bl_long$CytoC, list(bl_long$CL_Peptide), mean))
colnames(bl_long_agg) <- c("CL_Peptide", "Mean_CytoC")
bl_long_agg2 <- data.frame(aggregate(bl_long$CytoC, list(bl_long$CL_Peptide), sd))
colnames(bl_long_agg2) <- c("CL_Peptide", "SD_CytoC")
bl_long_agg$SD_CytoC <- bl_long_agg2$SD_CytoC[match(bl_long_agg$CL_Peptide, bl_long_agg2$CL_Peptide)]
bl_long_agg$Cell_Line <- gsub(":.*$", "", bl_long_agg$CL_Peptide)
bl_long_agg$Peptide <- gsub("^.*:", "", bl_long_agg$CL_Peptide)
bl_long_agg$Concentration <- as.numeric(as.character(gsub("^.*_", "", bl_long_agg$Peptide)))
bl_long_agg$Peptide <- gsub("_", " (", bl_long_agg$Peptide)
bl_long_agg$Peptide <- paste0(bl_long_agg$Peptide, "uM)")
cnv <- read.csv(paste0(input_dir, "Final_chr1q_CN.csv"))
sum(bl_long_agg$Cell_Line %in% cnv$CellLine)/nrow(bl_long_agg)
##1
cnv_sub <- cnv[cnv$CellLine %in% bl_long_agg$Cell_Line,]
bl_long_agg$Cell_Line <- factor(bl_long_agg$Cell_Line, levels = cnv_sub$CellLine[order(cnv_sub$cnv, decreasing = F)])
bl_long_agg$Pep <- gsub(" \\(.*$", "", bl_long_agg$Peptide)
bl_long_agg$Pep <- factor(bl_long_agg$Pep, levels = c("A133", "ABT199", "FS1", "MS1", "HRKy", "PUMA", "BAD", "BIM"))
bl_long_agg$Peptide <- factor(bl_long_agg$Peptide, levels = c(unique(bl_long_agg$Peptide[order(bl_long_agg$Pep, bl_long_agg$Concentration, decreasing = F)])))
ms1 <- bl_long_agg[bl_long_agg$Peptide == "MS1 (5uM)",]
ic50$CytoC <- ms1$Mean_CytoC[match(ic50$CellLine, ms1$Cell_Line)]

#Import DepMap CRISPR screen data and annotate with DepMap MCL1 dependency 
dataset.ach <- read.csv("results/Figure1/dataset.ach.csv", row.names = 1, check.names = F)
mcl1_crispr <- data.frame(t(dataset.ach[rownames(dataset.ach) == "MCL1",]))
ic50$MCL1_CRISPR <- mcl1_crispr$MCL1[match(rownames(ic50), rownames(mcl1_crispr))]

#Import FISH data and annotate with chr1q copy number
fish <- read.csv(paste0(input_dir, "FISHData_Final.csv"), check.names = F)
ic50$FISH_1q <- fish$`1q22`[match(ic50$CellLine, fish$CellLine)]
write.csv(ic50, paste0(output_dir, "Annotated_IC50_Data.csv"), row.names = F)

#Plot association of IC50 with MCL1 RNA levels
###S63845
ic50 <- read.csv(paste0(output_dir, "Annotated_IC50_Data.csv"))
png(paste0(output_dir, "IC50_By_MCL1levels_S63845.png"), res = 300, units = "in", width = 5, height = 5)
ggplot(ic50) + 
  geom_smooth(aes(S63845, MCL1), method = "lm", fill = "orange", color = "steelblue") +
  geom_point(aes(S63845, MCL1), size = 3, fill = "tomato2", shape = 21) + 
  scale_x_log10() + 
  geom_text_repel(aes(S63845, MCL1, label = CellLine)) + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14)) + 
  ylab("MCL1 expression (RPKM)") + 
  xlab("IC50 w/ MCL1 inhibitor (S63845)") + 
  annotate("text", x=1, y = 300, label = paste0("r = ", signif(cor.test(log10(ic50$S63845), ic50$MCL1)$estimate, 3), ", p =", signif(cor.test(log10(ic50$S63845), ic50$MCL1)$p.value, 2)), fontface = "italic", size = 5)
dev.off()

###AZD5991
png(paste0(output_dir, "IC50_By_MCL1levels_AZD5991.png"), res = 300, units = "in", width = 5, height = 5)
ggplot(ic50) + 
  geom_smooth(aes(AZD5991, MCL1), method = "lm", fill = "orange", color = "steelblue") +
  geom_point(aes(AZD5991, MCL1), size = 3, fill = "tomato2", shape = 21) + 
  scale_x_log10() + 
  geom_text_repel(aes(AZD5991, MCL1, label = CellLine)) + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14)) + 
  ylab("MCL1 expression (RPKM)") + 
  xlab("IC50 w/ MCL1 inhibitor (AZD5991)") + 
  annotate("text", x=1, y = 300, label = paste0("r = ", signif(cor.test(log10(ic50$AZD5991), ic50$MCL1)$estimate, 2), ", p =", signif(cor.test(log10(ic50$AZD5991), ic50$MCL1)$p.value, 2)), fontface = "italic", size = 5)
dev.off()
