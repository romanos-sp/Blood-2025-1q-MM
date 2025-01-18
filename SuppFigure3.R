#Import dependencies
library(tidyverse)
library(ggrepel)

#Set input/output directory
input_dir <- "data/"
output_dir <- "results/SuppFigure3/"

#Supp. Figure 3A
##Plot quantified KMS11 treatment WB results
dens <- read.csv(paste0(input_dir, "Densitometry/KMS11_Densitometry.csv"))
dens$Condition <- gsub("KMS11 ", "", dens$SAMPLE)
dens$Condition <- factor(dens$Condition, levels = c("DMSO", "MCL1i", "PI3Ki", "COMBO"))
png(paste0(output_dir, "KMS11_Treatment_Densitometry.png"), res = 300, units = "in", width = 6, height = 5)
ggplot(dens) + 
  geom_bar(aes(Condition, dens$Protein_RelativeToVehicle, fill = Protein), stat = "identity", show.legend = F) + 
  theme(axis.text.x = element_text(angle = 65, hjust = 1, color = 'black', size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "italic")) +
  xlab("") + 
  ylab("Normalized protein level relative to DMSO") + 
  facet_grid(.~Protein)
dev.off()