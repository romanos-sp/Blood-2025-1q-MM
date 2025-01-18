#Import dependencies
library(tidyverse)
library(ggrepel)

#Set input/output directories
input_dir <- "data/"
output_dir <- "results/Figure3/"

#Summarize sensitivity to inhibitors
data <- read.csv("results/SuppFigure2/Annotated_IC50_Data.csv")
data$Status <- factor(ifelse(data$CellLine %in% c("SKMM2"), "MCL1-dep-low non-1q+",
                                  ifelse(data$CellLine %in% c("OCIMY5"), "MCL1-dep-high non-1q+",
                                         ifelse(data$CellLine %in% c("KMS18", "EJM"), "MCL1-dep-low 1q+", "MCL1-dep-high 1q+"))), 
                           levels = c("MCL1-dep-low non-1q+", "MCL1-dep-high non-1q+", "MCL1-dep-low 1q+", "MCL1-dep-high 1q+"))
median(data$AZD5991[data$Status == "MCL1-dep-high 1q+"])
##[1] 0.04709
median(data$S63845[data$Status == "MCL1-dep-high 1q+"])
##[1] 0.0149
median(data$AZD8186[data$Status == "MCL1-dep-high 1q+"])
##[1] 4.79
median(data$AZD5991[data$Status == "MCL1-dep-low non-1q+"])
##[1] 28.8
median(data$S63845[data$Status == "MCL1-dep-low non-1q+"])
##[1] 27.03
median(data$AZD8186[data$Status == "MCL1-dep-low non-1q+"])
##[1] 56.87
median(data$AZD5991[data$Status == "MCL1-dep-low 1q+"])
##[1] 0.9605
median(data$S63845[data$Status == "MCL1-dep-low 1q+"])
##[1] 0.4168
median(data$AZD8186[data$Status == "MCL1-dep-low 1q+"])
##[1] 21.175
median(data$AZD5991[data$Status == "MCL1-dep-high non-1q+"])
##[1] 0.1345
median(data$S63845[data$Status == "MCL1-dep-high non-1q+"])
##[1] 0.1197
median(data$AZD8186[data$Status == "MCL1-dep-high non-1q+"])
##[1] 16.74

#Figure 3B
##Plot association of MCL1 dependency with IC50s against S63845 and AZD5991
data <- read.csv("results/SuppFigure2/Annotated_IC50_Data.csv")
data_long <- gather(data[, c("CytoC", "S63845", "AZD5991", "CellLine")], "Drug", "IC50", -CytoC, -CellLine)
data_long$Drug <- paste0("MCL1i: ", data_long$Drug)
data_long$Drug <- factor(data_long$Drug, levels = c("MCL1i: S63845", "MCL1i: AZD5991"))
stats <- data.frame("Drug" = c("MCL1i: S63845", "MCL1i: AZD5991"), 
                    "x" = c(median(data$S63845)*2.5, median(data$AZD5991)*2.5),
                    "r" = c(signif(cor.test(log10(data$S63845), data$CytoC)$estimate, 2), signif(cor.test(log10(data$AZD5991), data$CytoC)$estimate, 2)),
                    "p" = c(signif(cor.test(log10(data$S63845), data$CytoC)$p.value, 2), signif(cor.test(log10(data$AZD5991), data$CytoC)$p.value, 2)))
stats$Drug <- factor(stats$Drug, levels = c("MCL1i: S63845", "MCL1i: AZD5991"))
data_long$Status <- factor(ifelse(data_long$CellLine %in% c("SKMM2"), "MCL1-dep-low non-1q+",
                           ifelse(data_long$CellLine %in% c("OCIMY5"), "MCL1-dep-high non-1q+",
                                  ifelse(data_long$CellLine %in% c("KMS18", "EJM"), "MCL1-dep-low 1q+", "MCL1-dep-high 1q+"))), 
                           levels = c("MCL1-dep-low non-1q+", "MCL1-dep-high non-1q+", "MCL1-dep-low 1q+", "MCL1-dep-high 1q+"))
png(paste0(output_dir, "IC50_By_CytoC.png"), res = 300, units = "in", width = 8, height = 6)
ggplot(data_long) + 
  geom_smooth(aes(IC50, CytoC), method = "lm", fill = "orange", color = "steelblue") +
  geom_point(aes(IC50, CytoC, fill = Status), size = 3, shape = 21) + 
  scale_x_log10() + 
  scale_fill_manual(values = c("steelblue", "lightblue", "orange", "tomato2")) + 
  geom_text_repel(aes(IC50, CytoC, label = CellLine)) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 14),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        legend.position = "top") + 
  ylab("Cytochrome C Loss (%)") + 
  xlab("IC50") + 
  geom_text(data = stats, aes(x=x, y = 80, label = paste0("r=", r, ", p=", p)), fontface = "italic", size = 5) + 
  facet_wrap(.~Drug, scales = "free")
dev.off()