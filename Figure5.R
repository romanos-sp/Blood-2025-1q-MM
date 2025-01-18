#Import dependencies
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(viridis)

#Set input/output directory
input_dir <- "data/WGS/"
output_dir <- "results/Figure5/"

#Figure 5B
#Read in data from the "modeled_segments_tumor" GATK output
C1 <- read_tsv(paste0(input_dir, "1_KMS12BM-Nang_C1_HFWG5DRX2.1.aligned.duplicates_marked.modelFinal.seg"), comment = "@")
C2 <- read_tsv(paste0(input_dir, "2_KMS12BM-Nang_C2_HFWG5DRX2.2.aligned.duplicates_marked.modelFinal.seg"), comment = "@")

#Read in chromosome sizes
hg19.chrom.sizes <- read_tsv(paste0(input_dir, "hg19.chrom.sizes.txt"), col_names = c("CONTIG", "SIZE"))
hg19.chrom.sizes$Start=0
hg19.chrom.sizes$End=hg19.chrom.sizes$SIZE
hg19.chrom.levels <-  paste0("chr", c(1:22))
hg19.chr.levels <- as.character(1:22)
hg19.chrom.sizes <- hg19.chrom.sizes %>% 
  mutate(Chromosome=factor(CONTIG, levels=hg19.chrom.levels)) %>% 
  filter(!is.na(Chromosome)) |>
  mutate(chr=factor(str_remove(Chromosome, "chr"), levels=hg19.chr.levels))
chr.color <- rep(c("#1F78B4", "#FF7F00"), 11)
names(chr.color) <- hg19.chr.levels

#Read in cytoband info
hg19.cent <- read_tsv(paste0(input_dir, "cytoBand.hg19.txt"), 
                      show_col_types = FALSE, 
                      col_names = c("Chromosome", "Start", "End", "Cytoband", "Type"))
hg19.cent <- hg19.cent |>
  filter(Type=="acen") |>
  mutate(Chromosome=factor(Chromosome, levels=hg19.chrom.levels)) %>% 
  filter(!is.na(Chromosome)) |>
  mutate(chr=factor(str_remove(Chromosome, "chr"), levels=hg19.chr.levels))


# Load DATA ---------------------------------------------------------------
C1 <- C1 |> 
  mutate(Chromosome=factor(paste0("chr", CONTIG), levels=hg19.chrom.levels),
         Start=START, End=END,
         chr=factor(str_remove(Chromosome, "chr"), levels=hg19.chr.levels)) |>
  filter(!is.na(chr))

C2 <- C2 |> 
  mutate(Chromosome=factor(paste0("chr", CONTIG), levels=hg19.chrom.levels),
         Start=START, End=END,
         chr=factor(str_remove(Chromosome, "chr"), levels=hg19.chr.levels)) |>
  filter(!is.na(chr))


# MAKE BACKGROUND ---------------------------------------------------------
cn.background <- ggplot(hg19.chrom.sizes) +
  geom_rect(aes(xmin=Start, xmax=End, ymin=-Inf, ymax=Inf), color=NA, fill=NA) +
  geom_rect(data = hg19.cent, aes(xmin=Start, xmax=End, ymin=-Inf, ymax=Inf), color="lightgrey", fill="lightgrey") +
  facet_grid(~chr, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = c(0, 5e7, 1e8, 1.5e8, 2e8, 2.5e8)) +
  theme_bw(base_size = 6, base_line_size = 0.25, base_rect_size = 0.25) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(), 
        panel.spacing = unit(0, 'lines'), 
        strip.background = element_blank())

# C1 ----------------------------------------------------------------------
##Denoised copy ratio
cn.plot <- cn.background +
  geom_segment(data=C1, aes(y=2^LOG2_COPY_RATIO_POSTERIOR_50, yend=2^LOG2_COPY_RATIO_POSTERIOR_50, x=START, xend=END, color=chr),  linewidth=1, alpha=1, position = position_nudge(y=0.01)) +
  geom_hline(yintercept = c(0.6, 0.9, 1.2, 1.5, 1.8), linetype=2, color="#333333", alpha=0.5, linewidth=0.25) +
  scale_color_manual(values=chr.color) +
  scale_fill_manual(values=chr.color) +
  labs(x="Chromosome", y="Denoised\nCopy Ratio") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))
png(paste0(output_dir, "C1_denoised_ratio.png"), res = 600, units = "in", width = 5, height = 1)
cn.plot
dev.off()

##Allelic imbalance
cn.maf <- cn.background +
  geom_hline(yintercept = c(0.5), linetype=2, color="#333333", alpha=0.5, linewidth=0.25) +
  geom_rect(data=C1, aes(ymin=MINOR_ALLELE_FRACTION_POSTERIOR_10, ymax=MINOR_ALLELE_FRACTION_POSTERIOR_90, xmin = START, xmax=END, fill=chr), alpha=0.3, color="#333333", position = position_nudge(y=0.01), size=0.25) +
  geom_rect(data=C1, aes(ymin=1-MINOR_ALLELE_FRACTION_POSTERIOR_10, ymax=1-MINOR_ALLELE_FRACTION_POSTERIOR_90, xmin = START, xmax=END, fill=chr), alpha=0.3, color="#333333", position = position_nudge(y=0.01), size=0.25) +
  geom_segment(data=C1, aes(y=MINOR_ALLELE_FRACTION_POSTERIOR_50, yend=MINOR_ALLELE_FRACTION_POSTERIOR_50, x=START, xend=END,), color="#333333",  linewidth=0.25, alpha=1, position = position_nudge(y=0.01), size=0.25) +
  geom_segment(data=C1, aes(y=1-MINOR_ALLELE_FRACTION_POSTERIOR_50, yend=1-MINOR_ALLELE_FRACTION_POSTERIOR_50, x=START, xend=END), color="#333333", linewidth=0.25, alpha=1, position = position_nudge(y=0.01), size=0.25) +
  scale_color_manual(values=chr.color) +
  scale_fill_manual(values=chr.color) +
  scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
  labs(x="Chromosome", y="Minor Allele\nFrequency") +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"))
png(paste0(output_dir, "C1_minor_allele.png"), res = 600, units = "in", width = 5, height = 1)
cn.maf
dev.off()

##Put together
cn.ratio.maf <- cn.plot + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank()) + 
  cn.maf + 
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  plot_layout(ncol = 1)
png(paste0(output_dir, "C1_copy_ratio_minor_allele.png"), res = 600, units = "in", width = 5, height = 2)
cn.ratio.maf
dev.off()


# C2 ----------------------------------------------------------------------
##Denoised copy ratio
cn.plot <- cn.background +
  geom_segment(data=C2, aes(y=2^LOG2_COPY_RATIO_POSTERIOR_50, yend=2^LOG2_COPY_RATIO_POSTERIOR_50, x=START, xend=END, color=chr),  linewidth=1, alpha=1, position = position_nudge(y=0.01)) +
  geom_hline(yintercept = c(0.6, 0.9, 1.2, 1.5, 1.8), linetype=2, color="#333333", alpha=0.5, linewidth=0.25) +
  scale_color_manual(values=chr.color) +
  scale_fill_manual(values=chr.color) +
  labs(x="Chromosome", y="Denoised\nCopy Ratio") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))
png(paste0(output_dir, "C2_denoised_ratio.png"), res = 600, units = "in", width = 5, height = 1)
cn.plot
dev.off()

##Allelic imbalance
cn.maf <- cn.background +
  geom_hline(yintercept = c(0.5), linetype=2, color="#333333", alpha=0.5, linewidth=0.25) +
  geom_rect(data=C2, aes(ymin=MINOR_ALLELE_FRACTION_POSTERIOR_10, ymax=MINOR_ALLELE_FRACTION_POSTERIOR_90, xmin = START, xmax=END, fill=chr), alpha=0.3, color="#333333", position = position_nudge(y=0.01), size=0.25) +
  geom_rect(data=C2, aes(ymin=1-MINOR_ALLELE_FRACTION_POSTERIOR_10, ymax=1-MINOR_ALLELE_FRACTION_POSTERIOR_90, xmin = START, xmax=END, fill=chr), alpha=0.3, color="#333333", position = position_nudge(y=0.01), size=0.25) +
  geom_segment(data=C2, aes(y=MINOR_ALLELE_FRACTION_POSTERIOR_50, yend=MINOR_ALLELE_FRACTION_POSTERIOR_50, x=START, xend=END,), color="#333333",  linewidth=0.25, alpha=1, position = position_nudge(y=0.01), size=0.25) +
  geom_segment(data=C2, aes(y=1-MINOR_ALLELE_FRACTION_POSTERIOR_50, yend=1-MINOR_ALLELE_FRACTION_POSTERIOR_50, x=START, xend=END), color="#333333", linewidth=0.25, alpha=1, position = position_nudge(y=0.01), size=0.25) +
  scale_color_manual(values=chr.color) +
  scale_fill_manual(values=chr.color) +
  scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
  labs(x="Chromosome", y="Minor Allele\nFrequency") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))
png(paste0(output_dir, "C2_minor_allele.png"), res = 600, units = "in", width = 5, height = 1)
cn.maf
dev.off()

##Put together
cn.ratio.maf <- cn.plot + 
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank()) + 
  cn.maf + 
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  plot_layout(ncol = 1)
png(paste0(output_dir, "C2_copy_ratio_minor_allele.png"), res = 600, units = "in", width = 5, height = 2)
cn.ratio.maf
dev.off()

#Figure 5C
##Plot FISH results for KMS12BM clones
fish <- read.csv(paste0(input_dir, "../Clones_FISHData_Final.csv"), check.names = F)
fish$MCL1_Probe1 <- fish$`MCL1_WI2-2997P20`
fish$MCL1_Probe2 <- fish$`MCL1_RP11-54A4`
fish_long <- gather(fish[, c("CellLine", "MCL1_Probe1", "MCL1_Probe2", "1q22", "1q42.13")], "Assay", "CN", -CellLine)
fish_long$CellLine <- factor(fish_long$CellLine, levels = c("KMS12BMClone1", "KMS12BMClone2"))
fish_long$Assay <- factor(fish_long$Assay, levels = c("1q22", "MCL1_Probe1", "MCL1_Probe2", "1q42.13"))
png(paste0(output_dir, "Clones_FISHCopyNumber.png"), res=300, unit="in", width=4, height=3)
ggplot(fish_long) + 
  geom_tile(aes(x = CellLine, y = Assay, fill = CN)) + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black", angle = 65, hjust = 1),
        axis.title = element_blank()) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "lightgrey") + 
  geom_text(aes(x = CellLine, y = Assay, label = CN)) + 
  scale_y_discrete(labels = c("FISH: 1q22", "FISH: MCL1 Probe 1", "FISH: MCL1 Probe 2", "FISH: 1q42.13"), expand = c(0,0)) + 
  scale_x_discrete(labels = c("KMS12BM-Clone1", "KMS12BM-Clone2"), expand = c(0,0))
dev.off()

#Figure 5F
##Plot BH3 profiling results for the two clones
bl <- read.csv(paste0(input_dir, "../Clones_BH3_Final.csv"))
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
bl_long_agg$Cell_Line <- factor(bl_long_agg$Cell_Line, levels = c("KMS12BMClone1", "KMS12BMClone2"))
bl_long_agg$Pep <- gsub(" \\(.*$", "", bl_long_agg$Peptide)
##Remove BIM (0.01uM), as Clone1 has NA for this concentration
bl_long_agg <- bl_long_agg[!bl_long_agg$Peptide %in% "BIM (0.01uM)",]
bl_long_agg$Pep <- factor(bl_long_agg$Pep, levels = c("A133", "ABT199", "FS1", "MS1", "HRKy", "PUMA", "BAD", "BIM"))
bl_long_agg$Peptide <- factor(bl_long_agg$Peptide, levels = c(unique(bl_long_agg$Peptide[order(bl_long_agg$Pep, bl_long_agg$Concentration, decreasing = F)])))
png(paste0(output_dir, "Baseline_BH3.png"), res=300, unit="in", width=5, height=6)
ggplot(bl_long_agg) +
  geom_tile(aes(Cell_Line, Peptide, fill = Mean_CytoC)) +
  scale_fill_viridis(name = "Cytochrome C Loss (%)", na.value = "lightgrey") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(angle = 65, hjust = 0)) +
  xlab("") +
  ylab("") +
  scale_x_discrete(position = "top", expand = c(0,0), labels = c("KMS12BM-Clone 1", "KMS12BM-Clone 2")) +
  scale_y_discrete(expand = c(0,0))
dev.off()
