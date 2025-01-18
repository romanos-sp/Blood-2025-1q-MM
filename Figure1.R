#Import dependencies
library(tidyverse)
library(matrixStats)
library(qusage)
library(cowplot)
library(ggstar)
library(UpSetR)
library(circlize)
library(biomaRt)
library(ggrepel)
library(ggpubr)
library(viridis)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/Figure1/"

#Prepare data
##24Q2: Chronos-processed CRISPR data and metadata (downloaded on 8/16/2024)
crispr <- read.csv(paste0(input_dir, "CRISPRGeneEffect_24Q2.csv"), row.names = 1, check.names = F)
meta <- read.csv(paste0(input_dir, "Model_24Q2.csv"))
mm_lines <- meta[meta$OncotreeSubtype == "Plasma Cell Myeloma",]
nrow(mm_lines)
##34
##Remove non-myeloma cell lines (KE97, HUNS1, COLO775), and lines with genotyping issues (KMS18, KMS21BM)
mm_lines <- mm_lines[!mm_lines$StrippedCellLineName %in% c("KE97", "HUNS1", "COLO775", "KMS18", "KMS21BM"),]
nrow(mm_lines)
##29 MM cell lines in the project
sum(rownames(crispr) %in% mm_lines$ModelID)
##17 MM cell lines with dependency data
dataset.ach <- crispr[rownames(crispr) %in% mm_lines$ModelID,]
rownames(dataset.ach) <- mm_lines$StrippedCellLineName[match(rownames(dataset.ach), mm_lines$ModelID)]
colnames(dataset.ach) <- gsub(" (.*)$", "", colnames(dataset.ach))
dataset.ach <- data.frame(t(dataset.ach), check.names = F)
nrow(dataset.ach)
##18,443 genes
colnames(dataset.ach)
##[1] "OPM2"     "L363"     "LP1"      "SKMM2"    "KMS20"    "OCIMY7"   "INA6"     "KMS34"   
##[9] "KMS27"    "KMS26"    "JJN3"     "KMS11"    "MM1S"     "RPMI8226" "EJM"      "AMO1"    
##[17] "KMM1"   
write.csv(dataset.ach, paste0(output_dir, "dataset.ach.csv"), row.names = T)

##DEMETER2 V6: RNAi data and metadata (downloaded on 8/16/2024)
rnai <- read.csv(paste0(input_dir, "D2_combined_gene_dep_scores_V6.csv"), row.names = 1, check.names = F)
meta <- read.csv(paste0(input_dir, "sample_info_DEMETER2V6.csv"))
mm_lines <- meta[meta$disease %in% "multiple_myeloma",]
nrow(mm_lines)
##16
##Remove KMS18, due to genotyping issues in the dataset
mm_lines <- mm_lines[!mm_lines$CCLE_ID %in% "KMS18_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",]
nrow(mm_lines)
##15 MM cell lines in the project
sum(colnames(rnai) %in% mm_lines$CCLE_ID)
##15 MM cell lines with dependency data
dataset.rnai <- rnai[,colnames(rnai) %in% mm_lines$CCLE_ID]
colnames(dataset.rnai) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", colnames(dataset.rnai))
rownames(dataset.rnai) <- gsub(" (.*)$", "", rownames(dataset.rnai))
##For rows annotated with multiple genes (for example, FAM99B & FAM99A), keep the first one
rownames(dataset.rnai) <- gsub("&.*$", "", rownames(dataset.rnai))
nrow(dataset.rnai)
##17,309 genes
colnames(dataset.rnai)
##[1] "AMO1"     "INA6"     "JJN3"     "KMM1"     "KMS12BM"  "KMS27"    "KMS28BM"  "KMS34"   
##[9] "L363"     "MM1S"     "OCIMY7"   "OPM1"     "OPM2"     "RPMI8226" "SKMM2"   
write.csv(dataset.rnai, paste0(output_dir, "dataset.rnai.csv"), row.names = T)

#Determine chr1q copy number by merging data from re-analyzed CCLE WGS/WES and FISH
ccle <- read.csv("CCLE/CKS1B_CCLE.csv")
fish <- read.csv(paste0(input_dir, "FISHData_Final.csv"), check.names = F)
fish <- fish[, c("CellLine", "1q22")]
colnames(fish) <- c("CellLine", "FISH")
cnv <- merge(ccle, fish, by = "CellLine", all = T)
##Determine final chr1q copy number based on FISH or, if FISH testing not done, WES/WGS.
cnv$cnv <- ifelse(is.na(cnv$CN), cnv$FISH,
                  ifelse(is.na(cnv$FISH), cnv$CN, cnv$FISH))
write.csv(cnv, paste0(input_dir, "Final_chr1q_CN.csv"), row.names = F)

#Read in RNA-seq RPKM data (18Q2)
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
rpkm.non <- as.character(cnv$CellLine)[cnv$cnv=="2"]
rpkm.non <- rpkm.non[rpkm.non %in% colnames(rpkm_mm)]
length(rpkm.non)
##1
rpkm.non
##[1]"SKMM2"
rpkm.amp <- colnames(rpkm_mm)[!colnames(rpkm_mm) %in% c(rpkm.non, "Name", "Description")]
length(rpkm.amp)
##17
rpkm.amp
##[1] "AMO1"      "EJM"       "JJN3"      "KARPAS620" "KHM1B"     "KMS11"     "KMS12BM"  
##[8] "KMS20"     "KMS26"     "KMS27"     "LP1"       "MM1S"      "MOLP2"     "MOLP8"    
##[15] "NCIH929"   "OPM2"      "RPMI8226"  
##Compute expression-level effect sizes
rpkm_working.exp <- rpkm_mm
rpkm_working.exp$median.amp <- rowMedians(as.matrix(rpkm_working.exp[,rpkm.amp]), na.rm = T)
rpkm_working.exp$median.non <- rowMedians(as.matrix(rpkm_working.exp[,rpkm.non]), na.rm = T)
rpkm_working.exp$exp.effect.size <- rpkm_working.exp$median.amp - rpkm_working.exp$median.non
rpkm_working.exp$ExpFC <- log2(rpkm_working.exp$median.amp/rpkm_working.exp$median.non)
final.just.exp.effect.size <- rpkm_working.exp[, c("Description", "exp.effect.size", "ExpFC")]
final.just.exp.effect.size$Gene <- final.just.exp.effect.size$Description
final.just.exp.effect.size <- final.just.exp.effect.size[!duplicated(final.just.exp.effect.size$Description),]
summary(final.just.exp.effect.size$exp.effect.size)
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##-37346.23     -0.03      0.00     -2.31      0.00   2898.81   
quantile(final.just.exp.effect.size$exp.effect.size, probs = c(0.95))
##95% 
##1.955212  
sum(abs(final.just.exp.effect.size$exp.effect.size) > 2)
##6,868

#Plot expression thresholds
png(paste0(output_dir, "ExpressionThreshold.png"), res=300, unit="in", width=5, height=5)
ggplot(final.just.exp.effect.size) +
  geom_histogram(aes(exp.effect.size), bins = 30, fill = "lightblue", color = "black") +
  scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5), labels = c(-5, -2.5, 0, 2.5, 5), limits = c(-5, 5)) + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "steelblue") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "tomato2") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 12)) +
  xlab("Difference in Median Expression: Amp1q - Non-Amp1q") +
  ylab("# of genes")
dev.off()

#Determine hits for each dataset: a difference in median expression between HMCLs with and without 1q+ of > 2 or < -2, a negative median dependency score in 1q+ HMCLs and a negative difference in median dependency scores between HMCLs with and without 1q+
##DepMap RNAi dataset
dataset.rnai <- read.csv(paste0(output_dir, "dataset.rnai.csv"), row.names = 1, check.names = F)
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

##DepMap CRISPR dataset
dataset.ach <- read.csv(paste0(output_dir, "dataset.ach.csv"), row.names = 1, check.names = F)
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

##Unique 1q+ and non-1q+ HMCLs across the two datasets
length(unique(c(ach.amp, rnai.amp)))
##12
length(unique(c(ach.non, rnai.non)))
##1

##Filter for genes with an expression effect size of >2 or < -2
thres = 2
##DepMap CRISPR dataset
intermediate <- dataset.ach.neg.more.neg
intermediate$Gene <- rownames(intermediate)
integrated <- merge(intermediate, final.just.exp.effect.size, by = "Gene")
gw.achilles.hits.more.than.2 <- integrated$Gene[abs(integrated$exp.effect.size) > thres]
write.table(gw.achilles.hits.more.than.2, paste0(output_dir, "gw.achilles.hits.more.than.2.txt"), sep = "\t")

##DepMap RNAi dataset
intermediate <- dataset.rnai.neg.more.neg
intermediate$Gene <- rownames(intermediate)
integrated <- merge(intermediate, final.just.exp.effect.size, by = "Gene")
gw.rnai.hits.more.than.2 <- integrated$Gene[abs(integrated$exp.effect.size) > thres]
write.table(gw.rnai.hits.more.than.2, paste0(output_dir, "gw.rnai.hits.more.than.2.txt"), sep = "\t")

##Compile and write hits across the 2 genome-wide screens
vector.of.hits.final <- c(unique(gw.rnai.hits.more.than.2), unique(gw.achilles.hits.more.than.2))
frequency.table.final <- data.frame(table(vector.of.hits.final))
colnames(frequency.table.final) <- c("Hit", "Count")
write.csv(frequency.table.final, paste0(output_dir, "GenomeWide_FrequencyTable.csv"), row.names = F)
table(frequency.table.final$Count)
##1    2 
##2309  499 

#Plot bar plot of dataset frequency among hits
data <- data.frame(table(frequency.table.final$Count))
colnames(data) <- c("NDataset", "Count")
png(paste0(output_dir, "DependencyHits_Bar.png"), res=300, unit="in", width=5, height=5)
ggplot(data) +
  geom_bar(aes(factor(NDataset), Count, fill = factor(NDataset)), stat = "identity", show.legend = F, color = "black") +
  geom_text(aes(factor(NDataset), Count, label = Count), vjust = -1) +
  scale_fill_manual(values = c("lightblue", "steelblue", "dodgerblue4")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14)) +
  xlab("# of datasets") +
  ylab("# of hits") +
  ylim(0, 2400)
dev.off()

#Figure 1B
#Plot upsetR plot of dataset overlap among hits
vector.of.hits.final_list <- list("DepMap RNAi screen" = unique(gw.rnai.hits.more.than.2), "DepMap CRISPR screen" = unique(gw.achilles.hits.more.than.2))
png(paste0(output_dir, "DependencyHits_Upset.png"), res=300, unit="in", width=5, height=5)
upset(fromList(vector.of.hits.final_list), 
      order.by = "freq", 
      point.size = 2.5, 
      line.size = 1, 
      text.scale = c(2, 1.5, 1.5, 1, 1.5, 2), 
      sets.x.label = "# of hits/screen", 
      mainbar.y.label = "# of hits", 
      sets = c("DepMap RNAi screen", "DepMap CRISPR screen"), 
      keep.order = T, 
      sets.bar.color= "black", 
      main.bar.color = c("steelblue",  "orange", "tomato2"), 
      att.color = "black",
      matrix.color = "black", 
      shade.color = "pink", 
      group.by = "degree")
dev.off()

#Figure 1C
#Visualize positional enrichment across the genome for genes that hit in both datasets
hits <- as.character(frequency.table.final$Hit)[frequency.table.final$Count>1]
ensembl <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
which.max(c(nrow(dataset.ach), nrow(dataset.rnai)))
##1
res2 <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "chromosome_name", "start_position", "end_position", "band"),
              filters = "hgnc_symbol",
              values = rownames(dataset.ach), 
              mart = ensembl)
res2 <- res2[res2$chromosome_name %in% c(1:22),]
res2$chr <- paste0("chr", res2$chromosome_name)
res2$start <- res2$start_position
res2$end <- res2$end_position
res2$hit <- ifelse(res2$hgnc_symbol %in% hits, "Hit", "No Hit")
res2 <- res2 %>% group_by(chr) %>% mutate(bin = cut(start, 30))
res2$chr_bin <- paste0(res2$chr, ":", res2$bin)
write.csv(res2, paste0(output_dir, "Positional_Data.csv"), row.names = F)
res2 <- read.csv(paste0(output_dir, "Positional_Data.csv"))
freq <- data.frame(table(res2$chr_bin, res2$hit))
##Test segments with at least 10 hits for enrichment
hitfreq <- freq[freq$Var2=="Hit",]
bins <- as.character(hitfreq$Var1)[hitfreq$Freq>= 10]
for(i in seq_along(bins)){
  res2$tmp <- ifelse(res2$chr_bin == bins[i], bins[i], "Other")
  if(i==1){
    region_hyper <- data.frame("chrom_bin" = bins[i], "hyper_p" = phyper(q = sum(res2$chr_bin == bins[i] & res2$hit == "Hit") - 1, m = sum(res2$chr_bin == bins[i]), n = sum(res2$chr_bin != bins[i]), k = sum(res2$hit == "Hit"), lower.tail = F), "fisher_p" = fisher.test(table(res2$hit, res2$tmp))$p.value)
  } else {
    region_hyper <- rbind(region_hyper, data.frame("chrom_bin" = bins[i], "hyper_p" = phyper(q = sum(res2$chr_bin == bins[i] & res2$hit == "Hit") - 1, m = sum(res2$chr_bin == bins[i]), n = sum(res2$chr_bin != bins[i]), k = sum(res2$hit == "Hit"), lower.tail = F), "fisher_p" = fisher.test(table(res2$hit, res2$tmp))$p.value))
  }
}
region_hyper$hyper_FDR <- p.adjust(region_hyper$hyper_p, method = "BH")
region_hyper$fisher_FDR <- p.adjust(region_hyper$fisher_p, method = "BH")
write.csv(region_hyper, paste0(output_dir, "DependencyHits_RegionHypergeometric.csv"), row.names = F)
##Plot segments with an FDR < 0.1
region_hyper <- region_hyper[region_hyper$hyper_FDR < 0.1,]
region_hyper$start <- gsub("^chr.*\\(", "", region_hyper$chrom_bin)
region_hyper$start <- gsub(",.*$", "", region_hyper$start)
region_hyper$end <- gsub("^chr.*\\(.*,", "", region_hyper$chrom_bin)
region_hyper$end <- gsub("\\]$", "", region_hyper$end)
region_hyper$start <- as.numeric(as.character(region_hyper$start))
region_hyper$end <- as.numeric(as.character(region_hyper$end))
region_hyper$chr <- gsub("\\:.*$", "", region_hyper$chrom_bin)
region_hyper$value <- 1
bed <- region_hyper[, c("chr", "start", "end", "value")]
nohits_bed <- res2[res2$hit != "Hit", c("chr", "start", "end")]
hits_bed <- res2[res2$hit == "Hit", c("chr", "start", "end")]
png(paste0(output_dir, "Dependency_Circos.png"), res=300, unit="in", width=6, height=6)
par(cex = 1)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), species = "hg38")
circos.genomicDensity(hits_bed, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(nohits_bed, col = c("orange"), track.height = 0.1)
circos.genomicTrack(bed, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, col = "steelblue", border = NA)
                    }, ylim = c(0, 1), bg.border = F)
dev.off()

#Figure 1D
#Import HALLMARK pathways and perform hypergeometric tests for hits in both datasets
##HALLMARK pathways downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp on 4/21/2023
hallmark <- read.gmt(paste0(input_dir, "h.all.v2023.1.Hs.symbols.gmt"))
urn <- max(nrow(dataset.ach), nrow(dataset.rnai))
urn
##18,443
dataset_hits <- as.character(frequency.table.final$Hit)[frequency.table.final$Count>1]
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
write.csv(out, paste0(output_dir, "DependencyHits_morethan2datasets_HALLMARK_Hypergeometric.csv"), row.names = F)

##Visualize results from pathway analysis (FDR < 0.05) 
out <- read.csv(paste0(output_dir, "DependencyHits_morethan2datasets_HALLMARK_Hypergeometric.csv"))
out_sig <- out[out$FDR < 0.05, ]
nrow(out_sig)
##26
out_sig$Pathway <- factor(out_sig$Pathway, levels = out_sig$Pathway[order(out_sig$FDR, decreasing = T)])
out_sig$GeneSetName <- gsub("^HALLMARK_", "", out_sig$Pathway)
png(paste0(output_dir, "DependencyHits_PathwayEnrichment.png"), res=300, unit="in", width=7, height=5)
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
  geom_text(aes(x = GenesInGeneSet, label = GenesInGeneSet), hjust=1.2, color="black", size=3.5) +
  geom_text(aes(x = NumberOfGenesInOverlap, label = NumberOfGenesInOverlap), hjust=1.2, color="white", size=3.5) +
  geom_text(aes(x = GenesInGeneSet + 5, label = paste0("q=", signif(FDR, digits = 2))), hjust=0, color="black", size=3.5, fontface = "italic") +
  scale_y_discrete(breaks = out_sig$Pathway, labels = out_sig$GeneSetName, expand = c(0.04, 0.04)) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200), limits = c(0, 260), expand = c(0.01, 0.01))
dev.off()

#Figure 1E
##Combine the two datasets for visualization purposes
colnames(dataset.ach) <- paste0(colnames(dataset.ach),".1")
colnames(dataset.rnai) <- paste0(colnames(dataset.rnai),".2")
tmp <- merge(dataset.ach, dataset.rnai, by = "row.names", all = T)
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL
large.dataset <- tmp
write.csv(large.dataset, paste0(output_dir, "Dependency_LargeDataset.csv"), row.names = T)
##Annotate hits with dataset of origin
large.dataset <- read.csv(paste0(output_dir, "Dependency_LargeDataset.csv"), row.names = 1)
frequency.table.final <- read.csv(paste0(output_dir, "GenomeWide_FrequencyTable.csv"))
gw.achilles.hits.more.than.2 <- as.character(read.delim(paste0(output_dir, "gw.achilles.hits.more.than.2.txt"))$x)
gw.rnai.hits.more.than.2 <- as.character(read.delim(paste0(output_dir, "gw.rnai.hits.more.than.2.txt"))$x)
frequency.table.final$Ach <- ifelse(frequency.table.final$Hit %in% gw.achilles.hits.more.than.2, 1, 0)
frequency.table.final$Rnai <- ifelse(frequency.table.final$Hit %in% gw.rnai.hits.more.than.2, 1, 0)
##For each gene, record its median RPKM level in 1q+ vs non-1q+ HMCLs
##Because 5 hits (LSP1, KMT2B, IDS, IL3RA, ASMTL) are duplicated in the RPKM dataset, the maximum level of expression will be considered across duplicates
##For each gene:
##if it hit in 1 dataset, record its median dependency score in 1q+ vs non-1q+ HMCLs in that respective dataset
##if it hit in 2 datasets, record the median dependency scores from the dataset with the largest effect size
out <- data.frame(Gene = NA, Dependency_Amp1q = NA, Dependency_Non1q = NA, Expression_Amp1q = NA, Expression_Non1q = NA)
for(i in 1:nrow(frequency.table.final)){
  if(i%%100==0){
    print(paste0(i, "..."))
  }
  out[i, "Gene"] <- as.character(frequency.table.final$Hit)[i]
  out[i, "Expression_Amp1q"] <- max(rpkm_working.exp$median.amp[rpkm_working.exp$Description %in% frequency.table.final$Hit[i]], na.rm = T)
  out[i, "Expression_Non1q"] <- max(rpkm_working.exp$median.non[rpkm_working.exp$Description %in% frequency.table.final$Hit[i]], na.rm = T)
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
out$DatasetCount <- frequency.table.final$Count[match(as.character(out$Gene), as.character(frequency.table.final$Hit))]
write.csv(out, paste0(output_dir, "DependencyHit_Matrix.csv"), row.names = F)

##Import chr1q genes and cytobands from UCSC Table Browser with chr1:125,100,000-248,956,422 position and https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/
cytobands <- read.delim(paste0(input_dir, "cytoBandIdeo.txt"), header = F)
colnames(cytobands) <- c("Chrom", "Start", "End", "Cytoband", "Giemsa")
cytobands$Giemsa <- as.character(cytobands$Giemsa)
chr1q_cytobands <- cytobands[cytobands$Chrom=="chr1",]
chr1q_cytobands$Giemsa <- factor(chr1q_cytobands$Giemsa, levels = c("acen", "gvar", "gneg", "gpos25", "gpos50", "gpos75", "gpos100"))
chr1q_cytobands$Middle <- chr1q_cytobands$Start + ((chr1q_cytobands$End - chr1q_cytobands$Start)/2)
row_odd <- seq_len(nrow(chr1q_cytobands)) %% 2  
chr1q_cytobands$Cytoband_odd <- chr1q_cytobands$Cytoband
chr1q_cytobands$Cytoband_even <- chr1q_cytobands$Cytoband
chr1q_cytobands$Cytoband_even[!!row_odd] <- NA
chr1q_cytobands$Cytoband_odd[!row_odd] <- NA
chr1q_genes <- read.csv(paste0(input_dir, "chr1q_genes.csv"), header = T, check.names = F)
chr1q_genes <- chr1q_genes[!duplicated(chr1q_genes$hg38.kgXref.geneSymbol), ]
chr1q_genes <- chr1q_genes[-grep("ENSG", chr1q_genes$hg38.kgXref.geneSymbol),]
out$Start <- chr1q_genes$hg38.knownGene.txStart[match(out$Gene, chr1q_genes$hg38.kgXref.geneSymbol)]
out$chr <- chr1q_genes$hg38.knownGene.chrom[match(out$Gene, chr1q_genes$hg38.kgXref.geneSymbol)]
##Isolate hits located on chr1q and upregulated in expression (effect size > 2)
out_high <- out[((out$Expression_Amp1q-out$Expression_Non1q) > 2) & (!is.na(out$Start)), ]
##Annotate pan-essential genes
panessential <- read.csv(paste0(input_dir, "AchillesCommonEssentialControls.csv"))
panessential$gene <- gsub(" \\(.*$", "", panessential$Gene)
sum(out_high$Gene %in% panessential$gene)
##28
out_high$panessential <- ifelse(out_high$Gene %in% panessential$gene, "panessential", "no")
out_high$lab <- ifelse(out_high$DatasetCount > 1, as.character(out_high$Gene), NA)
write.csv(out_high, paste0(output_dir, "DependencyHits_chr1q.csv"), row.names = F)

##Plot
out_high <- read.csv(paste0(output_dir, "DependencyHits_chr1q.csv"))
out_high$DiffDep <- out_high$Dependency_Amp1q-out_high$Dependency_Non1q
a <- ggplot(out_high) +
  geom_point(aes(Start, -DiffDep, fill=factor(DatasetCount), size = factor(DatasetCount), alpha = factor(DatasetCount)), shape = 21) +
  geom_rug(aes(x=Start, y=NULL)) +
  scale_alpha_manual(values = c(0.5, 1, 1), guide = "none") +
  scale_size_manual(values = c(2, 4), name = "# datasets") + 
  theme(panel.background = element_blank(), panel.border = element_blank(), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 14), legend.position = "top", axis.line.y.left   = element_line(color = 'black'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("") +
  ylab("Differential Dependency") +
  scale_fill_manual(values = c("navajowhite2", "orangered3"), name = "# datasets") +
  geom_text_repel(aes(Start, -DiffDep, label = lab), segment.alpha = 0.3, max.overlaps = Inf, fontface = "italic") + 
  xlim(125100000, 248956422) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.3, 1.7), breaks = c(-0.3, 0, 0.5, 1, 1.5, 1.7), labels = c("", "0", "0.5", "1", "1.5", "")) 

b <- ggplot(chr1q_cytobands[chr1q_cytobands$Start >= 123400000,]) +
  geom_rect(aes(xmin = Start - 0.2, xmax = End + 0.2, ymin = 0.95, ymax = 1.05, fill=Giemsa), show.legend = F, color = "black") +
  scale_fill_manual(values = c("red", "#660099", "#FFFFFF", "#C0C0C0", "#808080", "#404040", "#000000")) +
  theme(panel.background = element_blank(), panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 14)) +
  scale_y_discrete(expand = c(-0.2, 0.2)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  geom_text(aes(x = Middle, y = 1.2, label = Cytoband_odd), size = 3) +
  geom_text(aes(x = Middle, y = 1.4, label = Cytoband_even), size = 3) +
  xlab("chr1q") +
  ylab("")

png(paste0(output_dir, "Dependency_chr1q_hits_v2.png"), res=300, unit="in", width=10, height=7)
plot_grid(a, NULL, b, align = "v", rel_heights = c(1, -0.04, 0.2), ncol = 1)
dev.off()

#Compute hypergeometric/Fisher's exact p-value for enrichment of 1q21-1q23
roi <- c(chr1q_cytobands$Start[chr1q_cytobands$Cytoband=="q21.1"], chr1q_cytobands$End[chr1q_cytobands$Cytoband=="q23.3"])
chr1q_genes$cytoband <- ifelse(chr1q_genes$hg38.knownGene.txStart >= roi[1] & chr1q_genes$hg38.knownGene.txStart <= roi[2], "1q21-1q23", "Other")
chr1q_genes$hit <- ifelse(chr1q_genes$hg38.kgXref.geneSymbol %in% out_high$Gene, "Hit", "No hit")
table(chr1q_genes$hit, chr1q_genes$cytoband)
##        1q21-1q23 Other
##Hit           84    58
##No hit       603  1075
phyper(q = sum(chr1q_genes$cytoband == "1q21-1q23" & chr1q_genes$hit == "Hit") - 1, m = sum(chr1q_genes$cytoband == "1q21-1q23"), n = sum(chr1q_genes$cytoband != "1q21-1q23"), k = sum(chr1q_genes$hit == "Hit"), lower.tail = F)
##5.928304e-08
fisher.test(table(chr1q_genes$hit, chr1q_genes$cytoband))
##p-value = 1.108e-07
##alternative hypothesis: true odds ratio is not equal to 1
##95 percent confidence interval:
##  1.796336 3.727796
##sample estimates:
##  odds ratio 
##2.580501 

#Figure 1F
#Compare chr1q copy number with MCL1 expression level
cnv <- read.csv(paste0(input_dir, "Final_chr1q_CN.csv"))
rpkm <- read.delim(paste0(input_dir, "CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct"), skip=2, check.names = F)
ccle_annotations <- read.delim(paste0(input_dir, "Cell_lines_annotations_20181226.txt"))
rpkm_mm <- rpkm[, c("Name", "Description", colnames(rpkm)[colnames(rpkm) %in% ccle_annotations$CCLE_ID[ccle_annotations$Hist_Subtype1=="plasma_cell_myeloma"]])]
###Remove KE97, HUNS1 which have been misannotated as MM, and KMS18/KMS21BM which presented genotyping issues in the DepMap datasetrpkm_mm <- rpkm_mm[, -grep("KE97|HUNS1|KMS18|KMS21BM", colnames(rpkm_mm))]
rpkm_mm <- rpkm_mm[, -grep("KE97|HUNS1|KMS18|KMS21BM", colnames(rpkm_mm))]
colnames(rpkm_mm) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", colnames(rpkm_mm))
##Remove cell lines without chr1q copy number information
rpkm_mm_lines <- colnames(rpkm_mm)[!colnames(rpkm_mm) %in% c("Name", "Description")] 
rpkm_mm_lines[!rpkm_mm_lines %in% cnv$CellLine]
##[1]"KMM1"        "KMS28BM"     "KMS34"       "L363"     "U266B1"   
rpkm_mm <- rpkm_mm[, colnames(rpkm_mm) %in% c("Name", "Description", rpkm_mm_lines[rpkm_mm_lines %in% cnv$CellLine]
)]
all(colnames(rpkm_mm)[!colnames(rpkm_mm) %in% c("Name", "Description")] %in% cnv$CellLine)
##TRUE
##Append MCL1 expression levels to copy number data and compare
data <- cnv
gep <- rpkm_mm[!duplicated(rpkm_mm$Description),]
rownames(gep) <- gep$Description
gep$Name <- NULL
gep$Description <- NULL
gep <- data.frame(t(gep))
data$MCL1 <- gep$MCL1[match(data$CellLine, rownames(gep))]
png(paste0(output_dir, "MCL1byChr1qCopyNumber.png"), res = 300, units = "in", width = 5, height = 5)
ggplot(data) + 
  geom_smooth(aes(cnv, MCL1), method = "lm", fill = "orange", color = "steelblue") +
  geom_point(aes(cnv, MCL1), size = 3, fill = "darkblue", shape = 21) + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14)) + 
  xlab("Chr1q copy number") + 
  ylab("DepMap MCL1 Expression") + 
  annotate("text", y = 1, x = 6.5, label = paste0("r = ", signif(cor.test(data$MCL1, data$cnv)$estimate, 2), ", p =", signif(cor.test(data$MCL1, data$cnv)$p.value, 2)), fontface = "italic", size = 5) + 
  geom_text_repel(aes(cnv, MCL1, label = CellLine), size = 3, segment.alpha = 0.2, segment.linetype = 2, segment.curvature = 1, max.overlaps = Inf, force = 2, fontface = "italic")
dev.off()

#Figure 1G
##Compare MCL1 dependency (DepMap CRISPR screens) to chr1q copy number
cnv <- read.csv(paste0(input_dir, "Final_chr1q_CN.csv"))
dataset.ach <- read.csv(paste0(output_dir, "dataset.ach.csv"), row.names = 1, check.names = F)
mcl1_crispr <- data.frame(t(dataset.ach[rownames(dataset.ach) == "MCL1",]))
cnv$MCL1 <- mcl1_crispr$MCL1[match(cnv$CellLine, rownames(mcl1_crispr))]
png(paste0(output_dir, "MCL1dependency_By_chr1qCN.png"), res = 300, units = "in", width = 5, height = 5)
ggplot(cnv) + 
  geom_smooth(aes(cnv, MCL1), method = "lm", fill = "orange", color = "steelblue") +
  geom_point(aes(cnv, MCL1), size = 3, fill = "tomato2", shape = 21) + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14)) + 
  xlab("Chr1q copy number") + 
  ylab("DepMap CRISPR MCL1 Dependency Score") + 
  annotate("text", y = -0.5, x = 6.5, label = paste0("r = ", signif(cor.test(cnv$MCL1, cnv$cnv)$estimate, 2), ", p =", signif(cor.test(cnv$MCL1, cnv$cnv)$p.value, 2)), fontface = "italic", size = 5) +
  geom_text_repel(aes(cnv, MCL1, label = CellLine), size = 3, segment.alpha = 0.2, segment.linetype = 2, segment.curvature = 1, max.overlaps = Inf, force = 2, fontface = "italic")
dev.off()

#Figure 1I
##Plot MCL1 copy number for all cell lines
fish <- read.csv(paste0(input_dir, "FISHData_Final.csv"), check.names = F)
fish$MCL1_Probe1 <- fish$`MCL1_WI2-2997P20`
fish$MCL1_Probe2 <- fish$`MCL1_RP11-54A4`
ccle <- read.csv("CCLE/CKS1B_CCLE.csv")
fish$WGS_WES <- ccle$CN[match(fish$CellLine, ccle$CellLine)]
fish_long <- gather(fish[, c("CellLine", "MCL1_Probe1", "MCL1_Probe2", "1q22", "WGS_WES")], "Assay", "CN", -CellLine)
fish_long$CellLine <- factor(fish_long$CellLine, levels = fish$CellLine[order(fish$`1q22`, decreasing = F)])
order <- fish$CellLine[order(fish$`1q22`, decreasing = F)]
png(paste0(output_dir, "MCL1_FISHCopyNumber.png"), res=300, unit="in", width=6, height=3)
ggplot(fish_long) + 
  geom_tile(aes(x = CellLine, y = Assay, fill = CN)) + 
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black", angle = 65, hjust = 1),
        axis.title = element_blank()) + 
  scale_fill_gradient(low = "white", high = "red", na.value = "lightgrey") + 
  geom_text(aes(x = CellLine, y = Assay, label = CN)) + 
  scale_y_discrete(labels = c("FISH: MCL1 Probe 1", "FISH: MCL1 Probe 2", "FISH: 1q22", "WGS/WES"), expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0))
dev.off()

#Figure 1J
#Import baseline data and plot heatmap
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
fish <- read.csv(paste0(input_dir, "FISHData_Final.csv"), check.names = F)
order <- fish$CellLine[order(fish$`1q22`, decreasing = F)]
bl_long_agg$Cell_Line <- factor(bl_long_agg$Cell_Line, levels = order)
bl_long_agg$Pep <- gsub(" \\(.*$", "", bl_long_agg$Peptide)
bl_long_agg$Pep <- factor(bl_long_agg$Pep, levels = c("A133", "ABT199", "FS1", "MS1", "HRKy", "PUMA", "BAD", "BIM"))
bl_long_agg$Peptide <- factor(bl_long_agg$Peptide, levels = c(unique(bl_long_agg$Peptide[order(bl_long_agg$Pep, bl_long_agg$Concentration, decreasing = F)])))
png(paste0(output_dir, "Baseline_BH3.png"), res=300, unit="in", width=8, height=6)
ggplot(bl_long_agg) +
  geom_tile(aes(Cell_Line, Peptide, fill = Mean_CytoC)) +
  scale_fill_viridis(name = "Cytochrome C Loss (%)") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(size = 14, color = "black"), axis.text.x = element_text(angle = 65, hjust = 0)) +
  xlab("") +
  ylab("") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

#Figure 1K
#Plot baseline boxplot for MS1 (5uM)
ms1 <- bl_long_agg[bl_long_agg$Peptide == "MS1 (5uM)", ]
png(paste0(output_dir, "Baseline_BH3_MS1_5uM.png"), res=300, unit="in", width=5, height=5)
ggplot(ms1) + 
  geom_bar(aes(Cell_Line, Mean_CytoC, fill = Cell_Line), stat = "identity", show.legend = F, color = "black") +
  geom_errorbar(aes(Cell_Line, ymin = Mean_CytoC-SD_CytoC, ymax = Mean_CytoC+SD_CytoC)) +
  geom_point(data = bl_long[bl_long$Peptide == "MS1_5",], aes(Cell_Line, CytoC, fill = Cell_Line), position = position_jitter(width = 0.2), shape =21, size = 4, alpha = 0.7, show.legend = F) +
  scale_fill_manual(values = c(c("skyblue3", "skyblue1", "brown4", "orangered4", "orangered3", "tomato3", "tomato1", "salmon", "salmon1"))) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black", angle = 65, hjust = 1),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 14)) +
  xlab("") +
  ylab("Cytochrome C Loss (%)") +
  ggtitle("MS1 (5uM)")
dev.off()

#Figure 1L
#Visualize hits from targeted shRNA screen in the 1q21-1q23 region and genome-wide CRISPR screen
##Identify hits in in-house genome-wide CRISPR screen
dataset.mr <- read.delim(paste0(input_dir, "MR.CRISPR.Screen.Ceres.Scores.txt"), row.names = 1)
colnames(dataset.mr) <- gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", colnames(dataset.mr))
dataset.mr.neg.more.neg <- dataset.mr
mr.amp <- "MM1S"
length(mr.amp)
##1
mr.non <- "KMS18"
length(mr.non)
##1
dataset.mr.neg.more.neg$median.amp <- rowMedians(as.matrix(dataset.mr.neg.more.neg[, mr.amp]), na.rm = T)
dataset.mr.neg.more.neg$median.non <- rowMedians(as.matrix(dataset.mr.neg.more.neg[, mr.non]), na.rm = T)
dataset.mr.neg.more.neg$effect.size <- dataset.mr.neg.more.neg$median.amp - dataset.mr.neg.more.neg$median.non
sum(is.na(dataset.mr.neg.more.neg$median.amp))
##136
sum(is.na(dataset.mr.neg.more.neg$effect.size))
##136
dataset.mr.neg.more.neg <- dataset.mr.neg.more.neg[!is.na(dataset.mr.neg.more.neg$effect.size) & !is.na(dataset.mr.neg.more.neg$median.amp),]
dataset.mr.neg.more.neg <- dataset.mr.neg.more.neg[dataset.mr.neg.more.neg$median.amp < 0 & dataset.mr.neg.more.neg$effect.size < 0,]
##Filter for upregulation in expression only (effect size > 2)
thres = 2
intermediate <- dataset.mr.neg.more.neg
intermediate$Gene <- rownames(intermediate)
integrated <- merge(intermediate, final.just.exp.effect.size, by = "Gene")
gw.mr.hits.more.than.2 <- integrated$Gene[integrated$exp.effect.size > thres]
write.table(gw.mr.hits.more.than.2, paste0(output_dir, "gw.mr.hits.more.than.2.txt"), sep = "\t")
##Identify hits in targeted shRNA screen (1q21-1q23)
dataset.t <- read.delim(paste0(input_dir, "Targeted_shRNA_dataset.txt"))
t.amp <- c("NCIH929", "MM1S", "OPM2")
length(t.amp)
##3
t.non <- c("KMS18")
length(t.non)
##1
dataset.t.neg.more.neg <- dataset.t
dataset.t.neg.more.neg$median.amp <- rowMedians(as.matrix(dataset.t.neg.more.neg[,t.amp]), na.rm = T)
dataset.t.neg.more.neg$median.non <- rowMedians(as.matrix(dataset.t.neg.more.neg[,t.non]), na.rm = T)
dataset.t.neg.more.neg$effect.size <- dataset.t.neg.more.neg$median.amp - dataset.t.neg.more.neg$median.non
sum(is.na(dataset.t.neg.more.neg$median.amp))
##1
sum(is.na(dataset.t.neg.more.neg$effect.size))
##1
dataset.t.neg.more.neg <- dataset.t.neg.more.neg[!is.na(dataset.t.neg.more.neg$effect.size) & !is.na(dataset.t.neg.more.neg$median.amp),]
dataset.t.neg.more.neg <- dataset.t.neg.more.neg[dataset.t.neg.more.neg$median.amp < 0 & dataset.t.neg.more.neg$effect.size < 0,]
dataset.t.neg.more.neg <- dataset.t.neg.more.neg[!is.na(dataset.t.neg.more.neg$hgnc_symbol),]
dataset.t.neg.more.neg$exp.effect.size <- final.just.exp.effect.size$exp.effect.size[match(dataset.t.neg.more.neg$hgnc_symbol, final.just.exp.effect.size$Gene)]
dataset.t.neg.more.neg$ExpFC <- final.just.exp.effect.size$ExpFC[match(dataset.t.neg.more.neg$hgnc_symbol, final.just.exp.effect.size$Description)]
dataset.t.neg.more.neg$ExpFC <- ifelse(is.infinite(dataset.t.neg.more.neg$ExpFC), max(dataset.t.neg.more.neg$ExpFC[is.finite(dataset.t.neg.more.neg$ExpFC)], na.rm = T), dataset.t.neg.more.neg$ExpFC)
dataset.t.neg.more.neg$hit <- ifelse(dataset.t.neg.more.neg$exp.effect.size > 2, 1, 0)
targeted.shrna.hits.more.than.2 <- as.character(dataset.t.neg.more.neg$hgnc_symbol)[dataset.t.neg.more.neg$hit %in% "1"]
write.table(targeted.shrna.hits.more.than.2, paste0(output_dir, "targeted.shrna.hits.more.than.2.txt"), sep = "\t")
dataset.t.neg.more.neg$confirmed <- ifelse(dataset.t.neg.more.neg$hgnc_symbol %in% gw.mr.hits.more.than.2, 1, 0)
dataset.t.neg.more.neg$shape_lab <- ifelse(dataset.t.neg.more.neg$confirmed %in% "1", 2, ifelse(dataset.t.neg.more.neg$confirmed %in% "0" & dataset.t.neg.more.neg$ExpFC > 0, 1, 0))
dataset.t.neg.more.neg$lab <- ifelse(dataset.t.neg.more.neg$confirmed %in% "1", as.character(dataset.t.neg.more.neg$hgnc_symbol), NA)
sum(dataset.t.neg.more.neg$confirmed %in% 1)
##13
##Plot scatterplot of hits
png(paste0(output_dir, "Dependency_TargetedshRNAScreen.png"), res=300, unit="in", width=6, height=5)
ggplot(dataset.t.neg.more.neg) + 
  geom_star(aes(ExpFC, -effect.size, fill = factor(shape_lab), size = factor(shape_lab), alpha = factor(shape_lab), starshape = factor(shape_lab)), color = "black") +
  scale_alpha_manual(values = c(0.7, 0.7, 1), guide = "none") +
  scale_starshape_manual(values = c(15, 15, 1), guide = "none") +
  scale_fill_manual(values = c("steelblue", "tomato2", "lemonchiffon"), guide = "none") +
  scale_size_manual(values = c(2, 2, 4), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.5) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size = 14), legend.position = c(0.2, 0.8)) +
  xlab(expression("Expression Log"[2]*" fold-change")) + 
  ylab("Differential Dependency") +
  xlim(-5, 5) +
  geom_text_repel(aes(ExpFC, -effect.size, label = lab), segment.alpha = 0.5, max.overlaps = Inf, fontface = "italic")
dev.off()