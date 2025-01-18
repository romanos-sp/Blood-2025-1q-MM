#Import dependencies
library(tidyverse)
library(data.table)
library(matrixStats)
library(stringr)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/Figure2/"

################Drug repurposing screen###########
#Figure 2A
##Import data and compound annotation from the Drug Repurposing Screen (annotation downloaded from https://clue.io/repurposing#download-data, version: 05/16/2018)
data <- read.csv(paste0(input_dir, "drug_dataset.csv"), na.strings = "", check.names = F)
annotation <- read.delim(paste0(input_dir, "repurposing_samples_20180516.txt"), skip = 9, header = T)
##Annotate data with compound information
data2 <- merge(data, annotation, by.x = "Drug", by.y = "broad_id", all.x = T)
data2$controls <- ifelse(is.na(data2$Drug), "NTC", 
                         ifelse(data2$Drug=="DMSO", "DMSO", 
                                ifelse(data2$Drug=="BTZ", "BTZ", 
                                       "Compounds")))
##Compare viability between MCL1-dep-high and low 1q+ HMCL
##Filter out controls and unannotated compounds
data2.noctrl <- data2[data2$controls=="Compounds" & !is.na(data2$pert_iname),]
##Select HMCLs
rownames(data2.noctrl) <- data2.noctrl$Drug
rownames(data2.noctrl) <- gsub("\\.", "-", rownames(data2.noctrl))
data2.noctrl <- data2.noctrl[,c("KMS18", "KMS12BM", "MM1S", "H929")]
colnames(data2.noctrl)[4] <- "NCIH929"
amp1q <- c("KMS12BM", "MM1S", "NCIH929")
length(amp1q)
##3
non1q <- "KMS18"
length(non1q)
##1
##Identify drugs for which: (i) the median viability for MCL1-dep-high 1q+ HMCLs is < 75%
##and (ii) the median viability for MCL1-dep-high 1q+ HMCLs is smaller than 75% of the score of the MCL1-dep-low 1q+ HMCL
data2.noctrl$Amp1q <- rowMedians(as.matrix(data2.noctrl[, colnames(data2.noctrl) %in% amp1q]), na.rm = T)
data2.noctrl$Non <- rowMedians(as.matrix(data2.noctrl[, colnames(data2.noctrl) %in% non1q]), na.rm = T)
data2.noctrl$Diff <- data2.noctrl$Amp1q-data2.noctrl$Non
data2.noctrl$hit <- ifelse(data2.noctrl$Amp1q < 0.75*data2.noctrl$Non & data2.noctrl$Diff < 0 & data2.noctrl$Amp1q < 75, 1, 0)
data2.noctrl$pert_iname <- data2$pert_iname[match(rownames(data2.noctrl), data2$Drug)]
data2.noctrl$vendor <- data2$vendor[match(rownames(data2.noctrl), data2$Drug)]
data2.noctrl$vendor_name <- data2$vendor_name[match(rownames(data2.noctrl), data2$Drug)]
data2.noctrl$catalog_no <- data2$catalog_no[match(rownames(data2.noctrl), data2$Drug)]
table(data2.noctrl$hit)
##0    1 
##4515  470 
##Annotate drugs with latest annotation (3/24/2020) from: https://clue.io/repurposing#download-data
drug_annotation <- read.delim(paste0(input_dir, "repurposing_drugs_20200324.txt"), skip = 9, header = T, na.strings = "")
data2.noctrl$moa <- drug_annotation$moa[match(data2.noctrl$pert_iname, drug_annotation$pert_iname)]
data2.noctrl$clinical_phase <- drug_annotation$clinical_phase[match(data2.noctrl$pert_iname, drug_annotation$pert_iname)]
data2.noctrl$target <- drug_annotation$target[match(data2.noctrl$pert_iname, drug_annotation$pert_iname)]
write.csv(data2.noctrl, paste0(output_dir, "DrugRepurposingScreen_hits.csv"), row.names = T)
##Visualize distributions for hits/no hits
data_long <- gather(data2.noctrl[, c("KMS18", "KMS12BM", "NCIH929", "MM1S", "pert_iname", "hit")], "Cell_Line", "Viability", -pert_iname, -hit)
data_long$Cell_Line <- as.character(data_long$Cell_Line)
data_long$Cell_Line <- factor(data_long$Cell_Line, levels = c("KMS18", "NCIH929", "MM1S", "KMS12BM"))
data_long$hit_lab <- ifelse(data_long$hit=="1", "Hits", "Other Compounds")
data_long$fill_lab <- ifelse(data_long$hit_lab == "Hits" & data_long$Cell_Line == "KMS18", "Hits_non1q", ifelse(data_long$hit_lab == "Hits" & data_long$Cell_Line != "KMS18", "Hits_1q", "No hits"))
png(paste0(output_dir, "DrugRepurposingScreen_Dist.png"), res=300, unit="in", width=6, height=5)
ggplot(data_long) + 
  geom_histogram(aes(Viability, fill=fill_lab), show.legend = F, bins = 30) +
  scale_fill_manual(values = c("tomato2", "steelblue",  "orange")) +
  scale_y_log10() +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 9), axis.title = element_text(size=14), strip.background = element_blank(), strip.text = element_text(size = 11, face = "italic")) +
  ylab("Count") +
  xlab("(%) Viability") +
  facet_grid(Cell_Line~hit_lab)
dev.off()

#Figure 2B
#Compare enrichment of different pathways with > 5 hits between hits and no hits
data2.noctrl$moa2 <- ifelse(data2.noctrl$moa %in% c("PI3K inhibitor", "mTOR inhibitor", "mTOR inhibitor|PI3K inhibitor"), "PI3K/mTOR inhibitor", as.character(data2.noctrl$moa))
freq <- data.frame(table(data2.noctrl$moa2, data2.noctrl$hit))
colnames(freq) <- c("moa", "hit", "count")
freq <- spread(freq, "hit", "count")
freq <- freq[freq$`1`>5,]
freq$p.value <- NA
freq$OR <- NA
for(i in 1:nrow(freq)){
  freq$p.value[i] <- fisher.test(data.frame(c(freq$`1`[i], freq$`0`[i]), c(sum(data2.noctrl$hit == "1", na.rm = T) - freq$`1`[i], sum(data2.noctrl$hit == "0", na.rm = T) - freq$`0`[i])))$p.value
  freq$OR[i] <- fisher.test(data.frame(c(freq$`1`[i], freq$`0`[i]), c(sum(data2.noctrl$hit == "1", na.rm = T) - freq$`1`[i], sum(data2.noctrl$hit == "0", na.rm = T) - freq$`0`[i])))$estimate
}
freq$FDR <- p.adjust(freq$p.value, method = "BH")
freq_plot <- freq[order(freq$`1`, decreasing = T),]
freq_plot$MoA <- c("PI3K/mTORi", "HDACi", "TOPi", "BRDi", "Calcium CB", "RRi", "EGFRi", "NFkBi")
freq_plot$MoA <- factor(freq_plot$MoA, levels = freq_plot$MoA[order(freq_plot$`1`, decreasing = F)])
freq_plot$lab <- ifelse(freq_plot$FDR < 0.05, paste0("q=", as.character(signif(freq_plot$FDR, digits = 1))), "NS")
freq_plot$lab <- paste0("n=", freq_plot$`1`, ", ", freq_plot$lab)
png(paste0(output_dir, "DrugRepurposingScreen_MoA_Bar.png"), res=300, unit="in", width=6, height=6)
ggplot(freq_plot) +
  geom_bar(aes(`1`, MoA, fill=OR), stat = "identity", color = "black") +
  geom_text(aes(`1`, MoA, label = lab), hjust = -0.1) +
  scale_fill_gradient(low = "navajowhite1", high = "orangered3") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size=14), legend.position = c(0.9, 0.2)) +
  xlab("# of compounds") +
  ylab("") +
  scale_x_continuous(breaks = c(0, 10, 20, 30), labels = c("0", "10", "20", "30"), limits = c(0, 40))
dev.off()

#Figure 2C
#Plot heatmap of PI3K inhibitors
hits <- data2.noctrl[data2.noctrl$hit==1,]
pi3k <- as.character(hits$pert_iname[hits$moa %in% "PI3K inhibitor"])
hits$FC <- hits$Amp1q/hits$KMS18
plot_df <- hits[hits$pert_iname %in% c(pi3k), c("KMS18", "KMS12BM", "MM1S", "NCIH929", "pert_iname", "Diff", "FC")]
plot_df <- gather(plot_df, "Cell_Line", "Viability", -pert_iname, -Diff, -FC)
plot_df$Cell_Line <- factor(plot_df$Cell_Line, levels = c("KMS18",  "KMS12BM", "NCIH929", "MM1S"))
plot_df$pert_iname <- factor(plot_df$pert_iname, levels = unique(plot_df$pert_iname[order(plot_df$FC, decreasing = F)]))
png(paste0(output_dir, "DrugRepurposingScreen_PI3Ki.png"), res=300, unit="in", width=6, height=5)
ggplot(plot_df) + 
  geom_tile(aes(y=pert_iname, x=Cell_Line, fill=Viability)) + 
  scale_fill_gradient(low = "navyblue", high="skyblue", name = "(%) Viability") + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size=14)) +
  ylab("Compound") +
  xlab("") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

################CMap-guided drug screen###########
#Figure 2D
##Import CMap data
cmap <- read.csv(paste0(input_dir, "CMap_screen.csv"))
##Select HMCLs
rownames(cmap) <- cmap$Compounds
cmap$Compounds <- NULL
amp1q <- c("KMS11", "NCIH929", "OPM2")
length(amp1q)
##3
non1q <- c("KMS18")
length(non1q)
##1
##Identify hits (see above for approach)
cmap$Amp1q <- rowMedians(as.matrix(cmap[, colnames(cmap) %in% amp1q]), na.rm = T)
cmap$Non1q <- rowMedians(as.matrix(cmap[, colnames(cmap) %in% non1q]), na.rm = T)
cmap$Diff <- cmap$Amp1q-cmap$Non1q
cmap$hit <- ifelse(cmap$Amp1q < 0.75*cmap$Non1q & cmap$Diff < 0 & cmap$Amp1q < 75, 1, 0)
table(cmap$hit)
##0   1 
##155  47
cmap$FC <- cmap$Amp1q/cmap$Non1q
rownames(cmap) <- gsub("\\.", "-", rownames(cmap))
write.csv(cmap, paste0(output_dir, "CMapScreen_hits.csv"), row.names = T)
##Annotate compounds
cmap_hits <- read.csv(paste0(output_dir, "CMapScreen_hits.csv"), row.names = 1)
cmap_annotation <- read.csv(paste0(input_dir, "CMap_Annotation.csv"), row.names = 1)
cmap_hits$Vendor_ID <- rownames(cmap_hits)
cmap_hits$Compound <- gsub("^.* ", "", rownames(cmap_hits))
cmap_hits <- merge(cmap_hits, cmap_annotation, by = "Compound", all.x = T)
rownames(cmap_hits) <- cmap_hits$Vendor_ID
write.csv(cmap_hits, paste0(output_dir, "CMapScreen_Annotatedhits.csv"), row.names = T)
##Visualize distributions for hits/no hits
cmap_long <- gather(cmap_hits[, c("KMS18", "KMS11", "NCIH929", "OPM2", "Compound", "hit")], "Cell_Line", "Viability", -Compound, -hit)
cmap_long$Cell_Line <- as.character(cmap_long$Cell_Line)
cmap_long$Cell_Line <- factor(cmap_long$Cell_Line, levels = c("KMS18", "NCIH929", "OPM2", "KMS11"))
cmap_long$hit_lab <- ifelse(cmap_long$hit=="1", "Hits", "Other Compounds")
cmap_long$fill_lab <- ifelse(cmap_long$hit_lab == "Hits" & cmap_long$Cell_Line %in% non1q, "Hits_non1q", ifelse(cmap_long$hit_lab == "Hits" & cmap_long$Cell_Line %in% amp1q, "Hits_1q", "No hits"))
png(paste0(output_dir, "CMapScreen_Dist.png"), res=300, unit="in", width=6, height=5)
ggplot(cmap_long[cmap_long$Compound != "Bortezomib", ]) + 
  geom_histogram(aes(Viability, fill=fill_lab), show.legend = F, bins = 20) +
  scale_fill_manual(values = c("tomato2", "steelblue",  "orange")) +
  scale_y_log10() +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 9), axis.title = element_text(size=14), strip.background = element_blank(), strip.text = element_text(size = 11, face = "italic")) +
  ylab("Count") +
  xlab("(%) Viability") +
  facet_grid(Cell_Line~hit_lab)
dev.off()

#Figure 2E
#Compare enrichment of different pathways with > 2 hits between hits and no hits
cmap_hits <- read.csv(paste0(output_dir, "CMapScreen_Annotatedhits.csv"), row.names = 1)
##unify PI3K/mTOR inhibitors: "PI3K inhibitor", "mTOR inhibitor", "mTOR inhibitor|PI3K inhibitor"
cmap_hits$moa2 <- ifelse(cmap_hits$moa %in% c("PI3K ", "mTOR", "DNA-PJ, PI3K, mTOR", "PI3K, mTOR"), "PI3K/mTOR inhibitor", as.character(cmap_hits$moa))
freq <- data.frame(table(cmap_hits$moa2, cmap_hits$hit))
colnames(freq) <- c("moa", "hit", "count")
freq <- spread(freq, "hit", "count")
freq <- freq[freq$`1`>1,]
freq$p.value <- NA
freq$OR <- NA
for(i in 1:nrow(freq)){
  freq$p.value[i] <- fisher.test(data.frame(c(freq$`1`[i], freq$`0`[i]), c(sum(cmap_hits$hit == "1", na.rm = T) - freq$`1`[i], sum(cmap_hits$hit == "0", na.rm = T) - freq$`0`[i])))$p.value
  freq$OR[i] <- fisher.test(data.frame(c(freq$`1`[i], freq$`0`[i]), c(sum(cmap_hits$hit == "1", na.rm = T) - freq$`1`[i], sum(cmap_hits$hit == "0", na.rm = T) - freq$`0`[i])))$estimate
}
freq$FDR <- p.adjust(freq$p.value, method = "BH")
freq_plot <- freq[order(freq$`1`, decreasing = T),]
freq_plot$MoA <- c("PI3K/mTORi", "CDKi", "TOPi", "antimetabolite", "AURKi", "HDACi")
freq_plot$MoA <- factor(freq_plot$MoA, levels = freq_plot$MoA[order(freq_plot$`1`, decreasing = F)])
freq_plot$lab <- ifelse(freq_plot$FDR < 0.05, paste0("q=", as.character(signif(freq_plot$FDR, digits = 1))), "NS")
freq_plot$lab <- paste0("n=", freq_plot$`1`, ", ", freq_plot$lab)
png(paste0(output_dir, "CMapScreen_MoA_Bar.png"), res=300, unit="in", width=6, height=6)
ggplot(freq_plot) +
  geom_bar(aes(`1`, MoA, fill=OR), stat = "identity", color = "black") +
  geom_text(aes(`1`, MoA, label = lab), hjust = -0.1) +
  scale_fill_gradient(low = "navajowhite1", high = "orangered3") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size=14), legend.position = c(0.9, 0.2)) +
  xlab("# of compounds") +
  ylab("") +
  scale_x_continuous(breaks = c(0, 5, 10), labels = c("0", "5", "10"), limits = c(0, 13))
dev.off()

#Figure 2F
##Plot heatmap of PI3K inhibitors that hit
pi3k <- rownames(cmap_hits)[grep("PI3K", cmap_hits$moa)]
plot_df <- cmap_hits[rownames(cmap_hits) %in% c(pi3k), c("KMS18", "KMS11", "OPM2", "NCIH929", "FC", "Compound", "Vendor_ID", "hit")]
plot_df <- plot_df[plot_df$hit %in% "1",]
plot_df <- gather(plot_df, "Cell_Line", "Viability", -Compound, -FC, -Vendor_ID, -hit)
plot_df$Cell_Line <- factor(plot_df$Cell_Line, levels = c("KMS18", "KMS11", "OPM2", "NCIH929"))
plot_df$Compound <- ifelse(plot_df$Compound %in% c("528113", "440206"), paste0("EMD-", plot_df$Compound), as.character(plot_df$Compound))
plot_df$Compound <- factor(plot_df$Compound, levels = unique(plot_df$Compound[order(plot_df$FC, decreasing = F)]))
png(paste0(output_dir, "CMapScreen_PI3Ki.png"), res=300, unit="in", width=6, height=5)
ggplot(plot_df) + 
  geom_tile(aes(y=Compound, x=Cell_Line, fill=Viability)) + 
  scale_fill_gradient(low = "navyblue", high="skyblue", name = "(%) Viability") + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="black"), axis.text = element_text(color = "black", size = 12), axis.title = element_text(size=14)) +
  ylab("Compound") +
  xlab("") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()
