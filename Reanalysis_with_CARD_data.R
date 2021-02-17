# Written by Julio Avelar-Barragan. Last updates: 7/3/20
library(ggplot2)
library(reshape2)
library(cowplot)
library(tidyverse)
library(vegan)

setwd("/media/julio/Storage/zr2514_report/Sarah_steele_collaboration_JA/")

rgi_annotations <- read.delim("/media/julio/Storage/zr2514_report/Sarah_steele_collaboration_JA/rgi_annotations.txt")
contig_table <- read.delim("/media/julio/Storage/zr2514_report/Sarah_steele_collaboration_JA/contig_table.txt", check.names = F)
metadata <- read.delim("Final_Mapping_File for microbiomeanalyst(copy).txt", row.names = 1)
annotated_counts <- merge(rgi_annotations, contig_table, by.x = "row.names", "#NAME")
annotated_counts <- aggregate(annotated_counts[,27:ncol(annotated_counts)], by = list(Best_Hit_ARO=annotated_counts$Best_Hit_ARO), FUN = sum)
#write.table(t(annotated_counts), "annotated_counts_rpkm.csv", quote = F, col.names = FALSE, sep = ",")
rpkm_table <- read.csv("annotated_counts_rpkm.csv", check.names = F)

#melt for stacked bar plot
Abx_melt <- melt(rpkm_table) %>% merge(metadata, ., by.x = "Alt_ID2", by.y = "Best_Hit_ARO")
Abx_melt <- mutate(Abx_melt,variable = str_replace(variable, "Bifidobacterium adolescentis rpoB mutants conferring resistance to rifampicin", "rpoB"))
Abx_melt$variable <- as.character(Abx_melt$variable)

#Take top 10 abx rest genes.
top_abx <- group_by(Abx_melt, variable) %>% summarise(., top_abx_tmp = sum(value)) %>% arrange(., desc(top_abx_tmp)) %>% slice(., 1:10)
high_abundance <- split(top_abx$variable, 1:NROW(top_abx))

#Change non top hits to other.
Abx_melt$variable[Abx_melt$variable %in% high_abundance != "TRUE"] <- "Other"
Abx_melt <- Abx_melt[order(Abx_melt$variable),] #Re order
Abx_melt <- rbind(Abx_melt[!(Abx_melt$variable == "Other"),],Abx_melt[(Abx_melt$variable == "Other"),]) #Move other to bottom
Abx_melt$variable <- factor(Abx_melt$variable, levels = unique(Abx_melt$variable)) #Fix the order

sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947", "gray")

#Bar plots
plot1 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "two"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "none", strip.background = element_rect(fill="lightblue")) +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1500)

plot2 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "one"), 
                  aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = "Abundance - RPKM") +
  theme(legend.position = "none", strip.background = element_rect(fill="lightblue"))+
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1500)

plot3 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "three"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank(), strip.background = element_rect(fill="lightblue")) + 
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1500)

plot4 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "four"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank(), strip.background = element_rect(fill="lightblue")) + 
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1500)


plot1.2 <- plot_grid(plot2, plot1, plot3, plot4,  nrow = 1, rel_widths = c(1.4,2.2,1.4,1.01))
plot1.2

ggsave(filename = "figure_1.2_CARD_data.svg", plot = plot1.2, width = 18, height = 9) #Write to file

### Statistical analysis ###
perma_df <- rpkm_table %>% merge(metadata, ., by.x = "Alt_ID2", by.y = "Best_Hit_ARO")

perma1 <- perma_df[perma_df$Timepoint == 1,]
adonis(formula = perma1[,12:ncol(perma1)] ~ perma1$Treatment, data = perma1, method = "bray", permutations = 9999, parallel = 0)
#No Significant differences between donors and pre-FMT, but its close. Maybe if there were more samples.

perma2 <- perma_df[perma_df$Timepoint %in% c(1,5) & !perma_df$Treatment == "Donors",]
adonis(formula = perma2[,12:ncol(perma2)] ~ as.factor(perma2$Individual) + as.factor(perma2$Timepoint), data = perma2, method = "bray", permutations = 9999, parallel = 0)
#Time point four was the most significant, at around 0.2, Individual explained around 60-80% in each case and was always significant.

###NMDS
NMDS <- metaMDS(perma_df[,12:ncol(perma_df)])
NMDS2 <- as.data.frame(NMDS$points) %>% cbind(., perma_df) %>% arrange(., Individual, Timepoint)

NMDS_plot <- ggplot(data = NMDS2) +
  aes(x = MDS1, y = MDS2) + 
  theme_classic() +
  geom_path(aes(linetype = as.factor(Individual), group = as.factor(Individual)), show.legend = F) +
  geom_point(aes(pch = Specific_donor, fill = as.factor(Timepoint)), size = 7, alpha = 0.7) +
  scale_shape_manual(values = c(21,22,23,24), name = "Donor", labels = c("3","1.1","2","1.2")) +
  geom_text(label = as.factor(NMDS2$Individual), size = 3) +
  guides(color = F, fill = guide_legend(override.aes = list(shape = 21))) +
  scale_fill_manual(name = "Time point", values=c("#DAF7A6", "#FFC300", "#FF5733", "#C70039", "#900C3F"), 
                    labels = c("Pre-FMT", "1 Week", "1 Month", "3 Months", "6 Months"))
NMDS_plot
ggsave("ARG_NMDS.svg", plot = NMDS_plot, device = "svg", units = "in", dpi = 1000, height = 4.5, width = 6)
