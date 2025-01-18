#Import dependencies
library(tidyverse)
library(ggpubr)
library(cowplot)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/Figure6/"

#Set seed
set.seed(1234)

#Import downsampled read data and determine cell cycle phases
data <- read.csv(paste0(input_dir, "prolif_mito_after_ds_table_18Krpc.csv"))

##Apply new quality filter
data$newfilter <- ifelse(data$nCount_RNA < 10000 | data$nFeature_RNA < 2500 | data$nFeature_RNA > 4500 | data$percent.mito > 10, "remove", "keep")
data <- data[data$newfilter %in% "keep",]

##Determine new Phases
data$newphase <- ifelse(data$S.Score > 0, "S", 
                        ifelse(data$G2M.Score > -0.2, "G2M", "G1"))
data$newphase <- ifelse(data$newphase != "G1" & data$S.Score <= 0.5*data$G2M.Score, "G2M", as.character(data$newphase))
data$newphase <- factor(data$newphase, levels = c("G1", "S", "G2M"))
write.csv(data, paste0(output_dir, "CellLine_SingleCellData.csv"), row.names = F)
          
#Figure 6B
##UMAP embedding
data <- read.csv(paste0(output_dir, "CellLine_SingleCellData.csv"), check.names = F)
data$Condition <- factor(data$Condition, levels = c("DMSO", "MCL1i", "PI3Ki", "COMBO"))
data$cell_line <- factor(data$cell_line, levels = c("KMS12BM", "KMS11"))
data$newphase <- factor(data$newphase, levels = c("G1", "S", "G2M"))
a <- ggplot(data[data$cell_line == "KMS12BM",]) + 
  geom_point(aes(UMAP_1, UMAP_2, color = newphase), stroke = 1, size = 1) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black")) + 
  scale_color_manual(values = c("lightblue", "tomato2", "orange"), name = "Phase") + 
  xlab("") + 
  ylab("") + 
  facet_grid(cell_line ~ Condition, scales = "free")
a_leg <- get_legend(a)
a <- ggplot(data[data$cell_line == "KMS12BM",]) + 
  geom_point(aes(UMAP_1, UMAP_2, color = newphase), stroke = 0, size = 1, show.legend = F) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black")) + 
  scale_color_manual(values = c("lightblue", "tomato2", "orange"), name = "Phase") + 
  xlab("") + 
  ylab("") + 
  facet_grid(cell_line ~ Condition, scales = "free")
b <- ggplot(data[data$cell_line == "KMS11",]) + 
  geom_point(aes(UMAP_1, UMAP_2, color = newphase), stroke = 0, size = 1, show.legend = F) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, color = "black")) + 
  scale_color_manual(values = c("lightblue", "tomato2", "orange"), name = "Phase") + 
  xlab("UMAP1") + 
  ylab("UMAP2") + 
  ylim(-5, 4.5) +
  facet_grid(cell_line ~ Condition, scales = "free")
plots <- plot_grid(a, b, ncol = 1, align = "vh")
png(paste0(output_dir, "CellLine_UMAP.png"), res = 300, units = "in", width = 8, height = 5)
plot_grid(plots, a_leg, ncol = 2, rel_widths = c(0.8, 0.1))
dev.off()

#Figure 6C
#Compare cycling phase proportions between conditions by bootstrapping
data <- read.csv(paste0(output_dir, "CellLine_SingleCellData.csv"), check.names = F)
data$Condition <- factor(data$Condition, levels = c("DMSO", "MCL1i", "PI3Ki", "COMBO"))
data$cell_line <- factor(data$cell_line, levels = c("KMS12BM", "KMS11"))
data$newphase <- factor(data$newphase, levels = c("G1", "S", "G2M"))
##Bootstrap p-values for KMS12BM
kms12bm <- data[data$cell_line %in% "KMS12BM",]
conditions <- c("COMBO", "MCL1i", "PI3Ki")
phases <- c("G1", "G2M", "S")
iter <- 10000
nullg1 <- c()
nulls <- c()
nullg2m <- c()
#Sample from both conditions iter number of times to create null distribution 
for(c in 1:length(conditions)){
  condition <- conditions[c]
  print(paste0("Condition: ", condition))
  data2 <- kms12bm[kms12bm$Condition %in% c(condition, "DMSO"),]
  for(i in 1:iter){
    if(i%%1000==0){
      print(paste0("iteration ", i))
    }
    c1 <- data2[sample(x = nrow(data2), size = sum(data2$Condition==condition)),]
    c2 <- data2[sample(x = nrow(data2), size = sum(data2$Condition=="DMSO")),]
    p1 <- table(c1$newphase)/nrow(c1)
    p2 <- table(c2$newphase)/nrow(c2)
    diff <- p1 - p2
    nullg1[i] <- diff[1]
    nullg2m[i] <- diff[3]
    nulls[i] <- diff[2]
  }
  #Calculate empirical p-value
  c1 <- data2[data2$Condition==condition,]
  c2 <- data[data2$Condition=="DMSO",]
  p1 <- table(c1$newphase)/nrow(c1)
  p2 <- table(c2$newphase)/nrow(c2)
  diff <- p1 - p2
  if(c==1){
    out <- data.frame("Condition" = rep(condition, 3), "Phase" = c("G1", "G2M", "S"), "Diff" = as.vector(unname(diff)), "p" = c(sum(abs(nullg1) >= abs(diff[1]))/iter, sum(abs(nullg2m) >= abs(diff[3]))/iter, sum(abs(nulls) >= abs(diff[2]))/iter))
  } else {
    out <- rbind(out, data.frame("Condition" = rep(condition, 3), "Phase" = c("G1", "G2M", "S"), "Diff" = as.vector(unname(diff)), "p" = c(sum(abs(nullg1) >= abs(diff[1]))/iter, sum(abs(nullg2m) >= abs(diff[3]))/iter, sum(abs(nulls) >= abs(diff[2]))/iter)))
  }
}
kms12bm_p <- out
write.csv(kms12bm_p, paste0(output_dir, "KMS12BM_CellCycle_EmpiricalP.csv"), row.names = F)

##Bootstrap p-values for KMS11
kms11 <- data[data$cell_line=="KMS11",]
conditions <- c("COMBO", "MCL1i", "PI3Ki")
phases <- c("G1", "G2M", "S")
iter <- 10000
nullg1 <- c()
nulls <- c()
nullg2m <- c()
#Sample from both conditions iter number of times to create null distribution 
for(c in 1:length(conditions)){
  condition <- conditions[c]
  print(paste0("Condition: ", condition))
  data2 <- kms11[kms11$Condition %in% c(condition, "DMSO"),]
  for(i in 1:iter){
    if(i%%1000==0){
      print(paste0("iteration ", i))
    }
    c1 <- data2[sample(x = nrow(data2), size = sum(data2$Condition==condition)),]
    c2 <- data2[sample(x = nrow(data2), size = sum(data2$Condition=="DMSO")),]
    p1 <- table(c1$newphase)/nrow(c1)
    p2 <- table(c2$newphase)/nrow(c2)
    diff <- p1 - p2
    nullg1[i] <- diff[1]
    nullg2m[i] <- diff[3]
    nulls[i] <- diff[2]
  }
  #Calculate empirical p-value
  c1 <- data2[data2$Condition==condition,]
  c2 <- data2[data2$Condition=="DMSO",]
  p1 <- table(c1$newphase)/nrow(c1)
  p2 <- table(c2$newphase)/nrow(c2)
  diff <- p1 - p2
  if(c==1){
    out <- data.frame("Condition" = rep(condition, 3), "Phase" = c("G1", "G2M", "S"), "Diff" = as.vector(unname(diff)), "p" = c(sum(abs(nullg1) >= abs(diff[1]))/iter, sum(abs(nullg2m) >= abs(diff[3]))/iter, sum(abs(nulls) >= abs(diff[2]))/iter))
  } else {
    out <- rbind(out, data.frame("Condition" = rep(condition, 3), "Phase" = c("G1", "G2M", "S"), "Diff" = as.vector(unname(diff)), "p" = c(sum(abs(nullg1) >= abs(diff[1]))/iter, sum(abs(nullg2m) >= abs(diff[3]))/iter, sum(abs(nulls) >= abs(diff[2]))/iter)))
  }
}
kms11_p <- out
write.csv(kms11_p, paste0(output_dir, "KMS11_CellCycle_EmpiricalP.csv"), row.names = F)

##Compute 95% confidence intervals for Phase frequency per cell line
##KMS12BM
subsample_factor <- 0.7
iter <- 10000
conditions <- c("DMSO", "MCL1i", "PI3Ki", "COMBO")
for(c in 1:length(conditions)){
  condition <- conditions[c]
  print(paste0("Condition: ", condition))
  data2 <- kms12bm[kms12bm$Condition==condition,]
  dists <- c()
  distg1 <- c()
  distg2m <- c()
  for(i in 1:iter){
    if(i%%1000==0){
      print(paste0("iteration ", i))
    }
    sample_df <- data2[sample(x = nrow(data2), size = round(subsample_factor*nrow(data2))), ]
    prop_df <- table(sample_df$newphase)/nrow(sample_df)
    distg1[i] <- unname(prop_df)[1]
    distg2m[i] <- unname(prop_df)[3]
    dists[i] <- unname(prop_df)[2]
  }
  if(c==1){
    out <- data.frame("Condition" = rep(condition,3), "Phase" = c("G1", "G2M", "S"), "Mean" = c(mean(distg1, na.rm = T), mean(distg2m, na.rm = T), mean(dists, na.rm = T)), "Median" = c(median(distg1, na.rm = T), median(distg2m, na.rm = T), median(dists, na.rm = T)), "SD" = c(sd(distg1, na.rm = T), sd(distg2m, na.rm = T), sd(dists, na.rm = T)), "Low" = c(quantile(distg1, 0.025), quantile(distg2m, 0.025), quantile(dists, 0.025)), "High" = c(quantile(distg1, 0.975), quantile(distg2m, 0.975), quantile(dists, 0.975)))
  } else {
    out <- rbind(out, data.frame("Condition" = rep(condition,3), "Phase" = c("G1", "G2M", "S"), "Mean" = c(mean(distg1, na.rm = T), mean(distg2m, na.rm = T), mean(dists, na.rm = T)), "Median" = c(median(distg1, na.rm = T), median(distg2m, na.rm = T), median(dists, na.rm = T)), "SD" = c(sd(distg1, na.rm = T), sd(distg2m, na.rm = T), sd(dists, na.rm = T)), "Low" = c(quantile(distg1, 0.025), quantile(distg2m, 0.025), quantile(dists, 0.025)), "High" = c(quantile(distg1, 0.975), quantile(distg2m, 0.975), quantile(dists, 0.975))))
  }
}
kms12bm_conf <- out
write.csv(kms12bm_conf, paste0(output_dir, "KMS12BM_CellCycle_ConfidenceInterval.csv"), row.names = F)

##KMS11
subsample_factor <- 0.7
iter <- 10000
conditions <- c("DMSO", "MCL1i", "PI3Ki", "COMBO")
for(c in 1:length(conditions)){
  condition <- conditions[c]
  print(paste0("Condition: ", condition))
  data2 <- kms11[kms11$Condition==condition,]
  dists <- c()
  distg1 <- c()
  distg2m <- c()
  for(i in 1:iter){
    if(i%%1000==0){
      print(paste0("iteration ", i))
    }
    sample_df <- data2[sample(x = nrow(data2), size = round(subsample_factor*nrow(data2))), ]
    prop_df <- table(sample_df$newphase)/nrow(sample_df)
    distg1[i] <- unname(prop_df)[1]
    distg2m[i] <- unname(prop_df)[3]
    dists[i] <- unname(prop_df)[2]
  }
  if(c==1){
    out <- data.frame("Condition" = rep(condition,3), "Phase" = c("G1", "G2M", "S"), "Mean" = c(mean(distg1, na.rm = T), mean(distg2m, na.rm = T), mean(dists, na.rm = T)), "Median" = c(median(distg1, na.rm = T), median(distg2m, na.rm = T), median(dists, na.rm = T)), "SD" = c(sd(distg1, na.rm = T), sd(distg2m, na.rm = T), sd(dists, na.rm = T)), "Low" = c(quantile(distg1, 0.025), quantile(distg2m, 0.025), quantile(dists, 0.025)), "High" = c(quantile(distg1, 0.975), quantile(distg2m, 0.975), quantile(dists, 0.975)))
  } else {
    out <- rbind(out, data.frame("Condition" = rep(condition,3), "Phase" = c("G1", "G2M", "S"), "Mean" = c(mean(distg1, na.rm = T), mean(distg2m, na.rm = T), mean(dists, na.rm = T)), "Median" = c(median(distg1, na.rm = T), median(distg2m, na.rm = T), median(dists, na.rm = T)), "SD" = c(sd(distg1, na.rm = T), sd(distg2m, na.rm = T), sd(dists, na.rm = T)), "Low" = c(quantile(distg1, 0.025), quantile(distg2m, 0.025), quantile(dists, 0.025)), "High" = c(quantile(distg1, 0.975), quantile(distg2m, 0.975), quantile(dists, 0.975))))
  }
}
kms11_conf <- out
write.csv(kms11_conf, paste0(output_dir, "KMS11_CellCycle_ConfidenceInterval.csv"), row.names = F)

###Combine results and plot
kms11_conf <- read.csv(paste0(output_dir, "KMS11_CellCycle_ConfidenceInterval.csv"))
kms12bm_conf <- read.csv(paste0(output_dir, "KMS12BM_CellCycle_ConfidenceInterval.csv"))
kms11_p <- read.csv(paste0(output_dir, "KMS11_CellCycle_EmpiricalP.csv"))
kms11_p$p <- ifelse(kms11_p$p == 0, 1e-04, as.numeric(as.character(kms11_p$p)))
kms11_p$p.adj <- p.adjust(kms11_p$p, method = "BH")
kms12bm_p <- read.csv(paste0(output_dir, "KMS12BM_CellCycle_EmpiricalP.csv"))
kms12bm_p$p <- ifelse(kms12bm_p$p == 0, 1e-04, as.numeric(as.character(kms12bm_p$p)))
kms12bm_p$p.adj <- p.adjust(kms12bm_p$p, method = "BH")
kms11_conf$Cell_Line <- "KMS11"
kms12bm_conf$Cell_Line <- "KMS12BM"
kms11_p$Cell_Line <- "KMS11"
kms12bm_p$Cell_Line <- "KMS12BM"
data <- rbind(kms11_conf, kms12bm_conf)
data$Phase <- factor(data$Phase, levels = c("G1", "S", "G2M"))
stats <- rbind(kms11_p, kms12bm_p)
stats$p.sign <- ifelse(stats$p.adj < 0.05 & stats$p.adj >= 0.01, "*", 
                       ifelse(stats$p.adj < 0.01 & stats$p.adj >=0.001, "**", 
                              ifelse(stats$p.adj < 0.001, "***", "")))
stats$y.pos <- 0
for(i in 1:nrow(stats)){
  vec <- data[(data$Cell_Line==as.character(stats$Cell_Line[i])) & (data$Phase==as.character(stats$Phase[i])), "Mean"]
  stats$y.pos[i] <- max(vec, na.rm = T)
}
cell_lines <- c("KMS11", "KMS12BM")
for(c in 1:length(cell_lines)){
  for(p in 1:length(phases)){
    stats$y.pos[stats$Cell_Line==cell_lines[c] & stats$Phase==phases[p] & stats$Condition=="PI3Ki"] <- stats$y.pos[stats$Cell_Line==cell_lines[c] & stats$Phase==phases[p] & stats$Condition=="MCL1i"] + 0.05
    stats$y.pos[stats$Cell_Line==cell_lines[c] & stats$Phase==phases[p] & stats$Condition=="COMBO"] <- stats$y.pos[stats$Cell_Line==cell_lines[c] & stats$Phase==phases[p] & stats$Condition=="MCL1i"] + 0.1
  }
}
stats$control <- "DMSO"
data$Phase <- factor(data$Phase, levels = c("G1", "S", "G2M"))
data$Condition <- factor(data$Condition, levels = c("DMSO", "MCL1i", "PI3Ki", "COMBO"))
data$Cell_Line <- factor(data$Cell_Line, levels = c("KMS12BM", "KMS11"))
stats$Condition <- factor(stats$Condition, levels = c("DMSO", "MCL1i", "PI3Ki", "COMBO"))
stats$Cell_Line <- factor(stats$Cell_Line, levels = c("KMS12BM", "KMS11"))
stats$Phase <- factor(stats$Phase, levels = c("G1", "S", "G2M"))

png(paste0(output_dir, "Downsampled_CellCycle_ErrorBars_Pvalues.png"), res=300, unit="in", width=7, height=5)
ggplot(data, aes(Condition, Mean*100, fill=Phase)) + 
  geom_bar(stat = "identity", show.legend = F) +
  geom_errorbar(aes(ymin=Mean*100 - 2.5*SD*100, ymax=Mean*100 + 2.5*SD*100), width=0.05, position="identity") +
  scale_fill_manual(values = c("lightblue", "tomato2", "orange")) + 
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color="Black"), axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14, color="black"), strip.background = element_blank(), strip.text = element_text(size = 12), axis.text.x = element_text(angle=65, hjust = 1)) +
  xlab("") +
  ylab("% of cells") +
  geom_bracket(aes(xmin = control, xmax = Condition, label = p.sign, y.position = y.pos*100+6), data = stats[stats$p < 0.05,], size = 0.3, label.size = 4, vjust = 0.6) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0,100)) +
  facet_grid(Cell_Line ~ Phase)
dev.off()

#Figure 6D
##Plot DE results of MCL1i vs DMSO
data <- read.csv(paste0(input_dir, "CellLine_DE_Results2.csv"))
data$mcl1_diff <- data$MCL1ivsDMSO_KMS11_Log2FC-data$MCL1ivsDMSO_KMS12BM_Log2FC
data$mcl1_lab <- ifelse(data$MCL1ivsDMSO_KMS11_Log2FC < -0.2 & data$MCL1ivsDMSO_KMS12BM_Log2FC < -0.2 & abs(data$mcl1_diff) < 0.2, 1, ifelse(data$MCL1ivsDMSO_KMS11_Log2FC > 0.2 & data$MCL1ivsDMSO_KMS12BM_Log2FC > 0.2 & abs(data$mcl1_diff) < 0.2, 1, 0))
data$mcl1_label <- ifelse(data$mcl1_lab=="1", as.character(data$Gene), NA)

png(paste0(output_dir, "MCL1i_DE.png"), res=300, unit="in", width=6, height=5)
ggplot(data) +
  annotate(geom = "rect", xmin=min(data$MCL1ivsDMSO_KMS11_Log2FC, na.rm = T) - 0.01, xmax=0, ymin=min(data$MCL1ivsDMSO_KMS12BM_Log2FC, na.rm = T) - 0.01, ymax=0, alpha=0.1, fill="blue") +
  annotate(geom = "rect", xmin=0, xmax=max(data$MCL1ivsDMSO_KMS11_Log2FC, na.rm = T) + 0.01, ymin=0, ymax=max(data$MCL1ivsDMSO_KMS12BM_Log2FC, na.rm = T) + 0.01, alpha=0.1, fill="red") +
  geom_point(aes(MCL1ivsDMSO_KMS11_Log2FC, MCL1ivsDMSO_KMS12BM_Log2FC, size = factor(mcl1_lab), alpha = factor(mcl1_lab), color = factor(mcl1_lab)), show.legend = F) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  scale_color_manual(values = c("gray45", "steelblue")) + 
  geom_text_repel(aes(MCL1ivsDMSO_KMS11_Log2FC, MCL1ivsDMSO_KMS12BM_Log2FC, label = mcl1_label), segment.alpha = 0.2, size = 3, fontface = 3, max.overlaps = Inf, force = 2) +
  xlab(expression("KMS11: MCL1i vs DMSO Log"[2]*"FC")) +
  ylab(expression("KMS12BM: MCL1i vs DMSO Log"[2]*"FC")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="Black"), 
        axis.text = element_text(size=12, color="black"), 
        axis.title = element_text(size=14, color="black"))
dev.off()

#Figure 6E
##Plot DE results of PI3Ki vs DMSO
data <- read.csv(paste0(input_dir, "CellLine_DE_Results2.csv"))
data$pi3k_diff <- data$PI3KivsDMSO_KMS11_Log2FC-data$PI3KivsDMSO_KMS12BM_Log2FC
data$pi3k_lab <- ifelse(data$PI3KivsDMSO_KMS11_Log2FC < -0.15 & data$PI3KivsDMSO_KMS12BM_Log2FC < -0.15 & abs(data$pi3k_diff) < 0.2, 1, 
                        ifelse(data$PI3KivsDMSO_KMS11_Log2FC > 0.15 & data$PI3KivsDMSO_KMS12BM_Log2FC > 0.15 & abs(data$pi3k_diff) < 0.2, 1, 0))
data$pi3k_label <- ifelse(data$pi3k_lab=="1", as.character(data$Gene), NA)
data2 <- data[data$PI3KivsDMSO_KMS11_Log2FC < 0.5 & data$PI3KivsDMSO_KMS12BM_Log2FC < 0.5, ]
png(paste0(output_dir, "PI3Ki_DE.png"), res=300, unit="in", width=6, height=5)
ggplot(data2) +
  annotate(geom = "rect", xmin=min(data2$PI3KivsDMSO_KMS11_Log2FC, na.rm = T) - 0.01, xmax=0, ymin=min(data2$PI3KivsDMSO_KMS12BM_Log2FC, na.rm = T) - 0.01, ymax=0, alpha=0.1, fill="steelblue") +
  annotate(geom = "rect", xmin=0, xmax=max(data2$PI3KivsDMSO_KMS11_Log2FC, na.rm = T) + 0.01, ymin=0, ymax=max(data2$PI3KivsDMSO_KMS12BM_Log2FC, na.rm = T) + 0.01, alpha=0.1, fill="orange") +
  geom_point(aes(PI3KivsDMSO_KMS11_Log2FC, PI3KivsDMSO_KMS12BM_Log2FC, size = factor(pi3k_lab), alpha = factor(pi3k_lab), color = factor(pi3k_lab)), show.legend = F) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  scale_color_manual(values = c("gray45", "steelblue")) + 
  geom_text_repel(aes(PI3KivsDMSO_KMS11_Log2FC, PI3KivsDMSO_KMS12BM_Log2FC, label = pi3k_label), segment.alpha = 0.2, size = 3, fontface = 3, max.overlaps = Inf, force = 2) +
  xlab(expression("KMS11: PI3Ki vs DMSO Log"[2]*"FC")) +
  ylab(expression("KMS12BM: PI3Ki vs DMSO Log"[2]*"FC")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="Black"), 
        axis.text = element_text(size=12, color="black"), 
        axis.title = element_text(size=14, color="black"))
dev.off()

#Figure 6F
##Plot dynamic BH3 profiling figure
dyn <- read.csv(paste0(input_dir, "Dynamic_BH3.csv"))
dyn_long <- gather(dyn, "Sample", "CytoC", -Peptide)
dyn_long$Cell_Line <- gsub("_.*$", "", dyn_long$Sample)
dyn_long$Concentration <- as.numeric(as.character(gsub("^.*_", "", dyn_long$Peptide)))
dyn_long$Peptide <- gsub("_", " (", dyn_long$Peptide)
dyn_long$Peptide <- paste0(dyn_long$Peptide, "uM)")
dyn_long$Pep <- gsub(" \\(.*$", "", dyn_long$Peptide)
dyn_long$Pep <- factor(dyn_long$Pep, levels = c("A133", "ABT199", "FS1", "MS1", "HRKy", "PUMA", "BAD", "BIM"))
dyn_long$Peptide <- factor(dyn_long$Peptide, levels = c(unique(dyn_long$Peptide[order(dyn_long$Pep, dyn_long$Concentration, decreasing = F)])))
dyn_long$Dose <- gsub("KMS11_", "", dyn_long$Sample)
dyn_long$Dose <- gsub("_.*$", "", dyn_long$Dose)
dyn_long$Dose <- paste0(dyn_long$Dose, "uM")
dyn_long$Dose <- paste0("AZD8186 ", dyn_long$Dose)
dyn_long$Dose <- factor(dyn_long$Dose, levels = c("AZD8186 2.5uM", "AZD8186 5uM"))
png(paste0(output_dir, "Dynamic_BH3.png"), res=300, unit="in", width=7, height=5)
ggplot(dyn_long[dyn_long$Dose == "AZD8186 5uM",]) +
  geom_point(aes(Peptide, CytoC, color = Dose), show.legend = F) +
  geom_bar(aes(Peptide, CytoC, fill = Dose), stat = "summary", fun.y = "mean", show.legend = F) +
  scale_fill_manual(values = c("orange", "tomato2")) +
  scale_color_manual(values = c("orange", "tomato2")) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.text.x = element_text (angle = 65, hjust = 1), 
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 14), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14, color = "black")) +
  xlab("") +
  ylab("Delta Cytochrome C Loss (%)") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -10, linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(.~Dose, ncol = 1)
dev.off()  

#Figure 6G
##Plot DE results of Combo vs DMSO
data <- read.csv(paste0(input_dir, "CellLine_DE_Results2.csv"))
data$combo_diff <- data$COMBOvsDMSO_KMS11_Log2FC-data$COMBOvsDMSO_KMS12BM_Log2FC
data$combo_lab <- ifelse(data$COMBOvsDMSO_KMS11_Log2FC < -0.2 & data$COMBOvsDMSO_KMS12BM_Log2FC < -0.2 & abs(data$combo_diff) < 0.25, 1, ifelse(data$COMBOvsDMSO_KMS11_Log2FC > 0.2 & data$COMBOvsDMSO_KMS12BM_Log2FC > 0.2 & abs(data$combo_diff) < 0.25, 1, 0))
data$combo_label <- ifelse(data$combo_lab=="1", as.character(data$Gene), NA)
data2 <- data
data2 <- data[data$COMBOvsDMSO_KMS11_Log2FC < 0.65 & data$COMBOvsDMSO_KMS12BM_Log2FC < 0.65, ]
##Remove label for AS transcript
data2$combo_label <- as.character(data2$combo_label)
data2$combo_label[data2$combo_label=="VPS9D1.AS1"] <- NA
png(paste0(output_dir, "Combo_DE.png"), res=300, unit="in", width=6, height=5)
ggplot(data2) +
  annotate(geom = "rect", xmin=min(data2$COMBOvsDMSO_KMS11_Log2FC, na.rm = T) - 0.01, xmax=0, ymin=min(data2$COMBOvsDMSO_KMS12BM_Log2FC, na.rm = T) - 0.01, ymax=0, alpha=0.1, fill="lightblue") +
  annotate(geom = "rect", xmin=0, xmax=max(data2$COMBOvsDMSO_KMS11_Log2FC, na.rm = T) + 0.01, ymin=0, ymax=max(data2$COMBOvsDMSO_KMS12BM_Log2FC, na.rm = T) + 0.01, alpha=0.1, fill="yellow") +
  geom_point(aes(COMBOvsDMSO_KMS11_Log2FC, COMBOvsDMSO_KMS12BM_Log2FC, size = factor(combo_lab), alpha = factor(combo_lab), color = factor(combo_lab)), show.legend = F) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.05, 1)) +
  scale_color_manual(values = c("gray45", "steelblue")) + 
  geom_text_repel(aes(COMBOvsDMSO_KMS11_Log2FC, COMBOvsDMSO_KMS12BM_Log2FC, label = combo_label), segment.alpha = 0.2, size = 3, fontface = 3, max.overlaps = Inf, force = 2) +
  xlab(expression("KMS11: Combo vs DMSO Log"[2]*"FC")) +
  ylab(expression("KMS12BM: Combo vs DMSO Log"[2]*"FC")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill=NA, color="Black"), 
        axis.text = element_text(size=12, color="black"), 
        axis.title = element_text(size=14, color="black"))
dev.off()


