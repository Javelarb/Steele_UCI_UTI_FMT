# Written by Julio Avelar-Barragan. Last updates: 7/3/20
library(ggplot2)
library(reshape2)
library(cowplot)
library(tidyverse)

setwd("/media/julio/Storage/zr2514_report/Sarah_steele_collaboration_JA/") #Set directory
metadata2 <- read.delim("ABX_metadata2.txt") #Minor tweaks to allow for plotting donor on same line

#### Need to get data in correct shape first. ####
abx_genes <-read.delim("abx_res_merged.txt", header = F) #Load combined file made in bash
abx_genes_unique <- as.data.frame(unique(abx_genes$V1)) #Only want the unique genes to merge them all into a single table.
abx_genes_unique <- as.data.frame(sort(abx_genes_unique$`unique(abx_genes$V1)`)) #Sort entries for proper cbind later

colnames(abx_genes_unique) <- "Gene" #Rename column

files <- list.files(path = "/media/julio/Storage/zr2514_report/Antibiotic_Resistance_Results", pattern = "*.txt", full.names = T) #Load in all Abx summary files from every sample.
files <- files[files != "Antibiotic_Resistance_Results/abx_res_merged.txt"] #Exclude abx_res_merged.txt if its in there.

file_tbl <- sapply(files, read.delim, simplify=F) #Load in all Abx summary files from every sample.

abx_table <- abx_genes_unique
#Merge them all together
for (i in 1:24) {
  abx_merge <- merge(abx_genes_unique,file_tbl[[i]], by = "Gene", all = T)
  abx_merge <- abx_merge[,(colnames(abx_merge) %in% c("Gene", "Relative_abundance"))]
  colnames(abx_merge) <- c("Gene", names(file_tbl[i]))
  abx_table <- cbind(abx_table, abx_merge[2])
}

hits_table <- abx_genes_unique
#Get hits instead of rel. ab.
for (i in 1:24) {
  hits <- merge(abx_genes_unique,file_tbl[[i]], by = "Gene", all = T)
  hits <- hits[,(colnames(hits) %in% c("Gene", "Hits"))]
  colnames(hits) <- c("Gene", names(file_tbl[i]))
  hits_table <- cbind(hits_table, hits[2])
}

#Get gene length for every gene.
Length_table <- abx_genes_unique
for (i in 1:24) {
  hits <- merge(abx_genes_unique,file_tbl[[i]], by = "Gene", all = T)
  hits <- hits[,(colnames(hits) %in% c("Gene", "Average_Gene_Length"))]
  colnames(hits) <- c("Gene", names(file_tbl[i]))
  Length_table <- cbind(Length_table, hits[2])
}

#Replace NA with 0
abx_table[is.na(abx_table)] <- 0
hits_table[is.na(hits_table)] <- 0
Length_table[is.na(Length_table)] <- 0

#Write tables, uncomment if wish to rewrite, but you will have to rename the samples again using a bash command.

#write.table(t(abx_table), "/media/julio/Storage/zr2514_report/Abx_summary_all_samples.csv", quote = F, sep = ",", col.names = F)
#write.table(t(Length_table), "/media/julio/Storage/zr2514_report/Gene_lengths.csv", quote = F, sep = ",", col.names = F)
#write.table(t(hits_table), "/media/julio/Storage/zr2514_report/Abx_summary_hits_all_samples.csv", quote = F, sep = ",", col.names = F)

#Reading in tables fixes an issue of having entries as characters instead of numbers.
Abx_summary_all_samples <- read.csv("Abx_summary_all_samples.csv", check.names = F, header = T)
hits_table <- read.csv("Abx_summary_hits_all_samples.csv", check.names = F, header = T)
Length_table <- read.csv("Gene_lengths.csv", check.names = F, header = T)
#NOTE: MISLABEL BY ZYMO. 612487.FMT.R6.V3 SHOULD BE R9 INSTEAD OF R6.

#Remove col 9 is the first donor sample, which apparantly didnt have abx resistance?
Abx_summary_all_samples <- Abx_summary_all_samples[-9,]
hits_table <- hits_table[-9,]
Length_table <- Length_table[-9,]

#Next normalize by RPKM. Need a scaling factor to normalize raw hits across samples
Read_counts <- read.delim("Read_counts.txt") #File from zymo
Read_counts$scaler <- Read_counts$rawseq/1000000 #Divide to get 'reads per million'

#Merge hits table with reads table to divide by per million scaler.
rpkm_table <- merge(hits_table, Read_counts, by.x = "Gene", by.y = "UniqueLabel") 
rpkm_table <- rpkm_table[,2:(NCOL(rpkm_table)-3)]/rpkm_table$scaler
rownames(rpkm_table) <- hits_table$Gene

rpkm_table <- rpkm_table/Length_table[,-1] #Normalize by dividing with gene length.
rpkm_table[is.na(rpkm_table)] <- 0 #replace na with 0.

#### Bar plot: credit to Andrew Oliver ####
#Andrew's metadata to allow for plot parity. 
metadata3 <- read.delim("Final_Mapping_File for microbiomeanalyst(copy).txt", row.names = 1)

#Converty RPKM to relative abundances
rel_ab_rpkm <- as.data.frame(rowSums(rpkm_table)) %>% merge(., metadata3, by= "row.names")
rel_ab_rpkm2 <- merge(rpkm_table, rel_ab_rpkm[,1:2], by.x = "row.names", by.y = "Row.names")
rel_ab_rpkm <- rel_ab_rpkm2[,2:(ncol(rel_ab_rpkm2)-1)]/rel_ab_rpkm2$`rowSums(rpkm_table)`
rel_ab_rpkm$Row.names <- rel_ab_rpkm2$Row.names

#Melt to allow for stacked barplot.
Abx_melt <- melt(rel_ab_rpkm) %>% merge(metadata3, ., by.x = "row.names", by.y = "Row.names")
Abx_melt$variable <- as.character(Abx_melt$variable)

#Take the top 13 antibiotic gene hits
top_abx <- group_by(Abx_melt, variable) %>% summarise(., top_abx_tmp = sum(value)) %>% arrange(., desc(top_abx_tmp)) %>% slice(., 1:10)
high_abundance <- split(top_abx$variable, 1:NROW(top_abx))

#Change things that arent in top into Other
Abx_melt$variable[Abx_melt$variable %in% high_abundance != "TRUE"] <- "Other"
Abx_melt <- Abx_melt[order(Abx_melt$variable),] #Reorder
Abx_melt <- rbind(Abx_melt[!(Abx_melt$variable == "Other"),],Abx_melt[(Abx_melt$variable == "Other"),]) #Move other to bottom
Abx_melt$variable <- factor(Abx_melt$variable, levels = unique(Abx_melt$variable)) #Keep the order fixed

#Color palette
sarah_color <- c("#003f5c", "#665191", "#d45087", "#ff7c43","#ffa600", "#7F0A57", "#CD9ABB", "#39A9AB", "#71CFC5", "#007947", "gray")

#Barplots
donors1 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "two"), 
                  aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))+
  scale_y_discrete(expand = c(0,0))

donors2 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "one"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))+
  scale_y_discrete(expand = c(0,0))

donors3 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "three"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  scale_y_discrete(expand = c(0,0))

donors4 <- ggplot(data = subset(Abx_melt, Abx_melt$Specific_donor == "four"), aes(x = as.factor(Timepoint), weight = value, fill = variable)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  scale_y_discrete(expand = c(0,0))
donors1

plot1.1 <- plot_grid(donors2, donors1, donors3, donors4, labels = NULL, nrow =1, rel_widths = c(1,2,1.1,2.35)) #plot as single
plot1.1
ggsave(filename = "figure_1.1_relab_Zymo_data.svg", plot = plot1.1, width = 21, height = 9) #Write to file

#### Abundance instead of relative abundance stacked barplot ####
#melt for stacked bar plot
Abx_melt2 <- melt(as.matrix(rpkm_table)) %>% merge(metadata3, ., by.x = "row.names", by.y = "Var1")
Abx_melt2$Var2 <- as.character(Abx_melt2$Var2)

#Take top 13 abx rest genes.
top_abx2 <- group_by(Abx_melt2, Var2) %>% summarise(., top_abx_tmp = sum(value)) %>% arrange(., desc(top_abx_tmp)) %>% slice(., 1:10)
high_abundance2 <- split(top_abx2$Var2, 1:NROW(top_abx2))

#Change non top hits to other.
Abx_melt2$Var2[Abx_melt2$Var2 %in% high_abundance2 != "TRUE"] <- "Other"
Abx_melt2 <- Abx_melt2[order(Abx_melt2$Var2),] #Re order
Abx_melt2 <- rbind(Abx_melt2[!(Abx_melt2$Var2 == "Other"),],Abx_melt2[(Abx_melt2$Var2 == "Other"),]) #Move other to bottom
Abx_melt2$Var2 <- factor(Abx_melt2$Var2, levels = unique(Abx_melt2$Var2)) #Fix the order

#Bar plots
donors5 <- ggplot(data = subset(Abx_melt2, Abx_melt2$Specific_donor == "two"), aes(x = as.factor(Timepoint), weight = value, fill = Var2)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1.2)

donors6 <- ggplot(data = subset(Abx_melt2, Abx_melt2$Specific_donor == "one"), 
                  aes(x = as.factor(Timepoint), weight = value, fill = Var2)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = 'Abundance - RPKM') +
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1.2)

donors7 <- ggplot(data = subset(Abx_melt2, Abx_melt2$Specific_donor == "three"), aes(x = as.factor(Timepoint), weight = value, fill = Var2)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1.2)

donors8 <- ggplot(data = subset(Abx_melt2, Abx_melt2$Specific_donor == "four"), aes(x = as.factor(Timepoint), weight = value, fill = Var2)) +
  geom_bar(color = "black") +
  theme_bw(base_size = 16) + facet_grid(. ~  Individual, space = "free", scales = "free") + 
  scale_fill_manual(values = sarah_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank()) +
  labs(x = NULL,
       y = NULL, fill = "Gene") +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm")) +
  ylim(0,1.2)

plot1.2 <- plot_grid(donors6, donors5, donors7, donors8, labels = NULL, nrow =1, rel_widths = c(1,1.9,1.1,2.45)) #plot as one plot.
plot1.2

ggsave(filename = "figure_1.2_ab_Zymo_data.svg", plot = plot1.2, width = 21, height = 9) #Write to file

#### Relative abundances of bacteria that carry these ABX resistance genes ####
#Load in OTU table from Andrews folder.
midas <- as.data.frame(t(read.csv("Bacteria_Composition_Summary.txt", check.names = FALSE, sep = "\t", row.names = 1, header = F)))
midas[,-1] <- apply(midas[,-1], c(1,2) ,as.numeric)
midas[,-1] <- midas[,-1]*100 #multiply by 100 since abx relative abundance is in percentage.

#See which taxa are responsible for the high abundance genes.
species <- NULL
for (i in 1:24) {
  assign(paste0("species", i), file_tbl[[i]][file_tbl[[i]][,1] %in% top_abx2$Var2[1:6],] %>% .[-c(1:2,12),1:2] %>% mutate(Species = str_replace(Species, "  Spp.", "")) %>% mutate(Species = str_replace(Species, "fragilis ", "fragilis")) %>% mutate(Species = str_replace(Species, "ruminicola ", "ruminicola"))) 
  species <- rbind(species, get(paste0("species", i)))
}
species <- unique(species) %>% .[!(.[,2] %in% c("Terrabacteria group ", "uncultured bacterium", "Terrabacteria group")),]

#Grab only the microbes that are present in species variable.
for (i in 1:length(paste0(unique(species$Species)))) {
assign(paste0("midas_", i), as.data.frame(midas[,grepl(paste0(unique(species$Species)[[i]]), colnames(midas))]))
}

#Grab the sum of the relative abundance for that taxonomy and place into one table.
midas_list <- list(midas_1, midas_2, midas_3, midas_4, midas_5, midas_6, midas_7) %>% 
  lapply(.,rowSums) %>% 
  set_names(paste0(unique(species$Species))) %>% 
  bind_rows(., .id = "column_ID")

#Grab Sample IDs 
midas_list$column_ID <-midas$Taxa

#Grab antibiotics now
top_taxa <- Abx_summary_all_samples[,colnames(Abx_summary_all_samples) %in% species$Gene]
rownames(top_taxa) <- Abx_summary_all_samples$Gene

#merge metadata with taxonomy since they have different labeling scheme.
top_taxa_plus_abx <- merge(metadata2,midas_list, by.x = "Alt_ID", by.y = "column_ID")

#merge antibiotics with taxonomy.
top_taxa_plus_abx2 <- merge(top_taxa_plus_abx,top_taxa, by.x = "X.NAME", by.y = "row.names")
top_taxa_plus_abx2  <- melt(top_taxa_plus_abx2[,-(2:8)], id.vars = "X.NAME") %>% merge(.,metadata2, by = "X.NAME") #Melt for barplot
top_taxa_plus_abx2 <- top_taxa_plus_abx2[order(top_taxa_plus_abx2$variable),]
top_taxa_plus_abx2 <- top_taxa_plus_abx2[!(top_taxa_plus_abx2$Individual == 3),] #Exclude donor 3 otherwise he will show up in plot

top_taxa_plus_abx2$variable <- as.character(top_taxa_plus_abx2$variable)
top_taxa_plus_abx2$variable[top_taxa_plus_abx2$variable == "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)"] <- "ErmF"
top_taxa_plus_abx2$variable[top_taxa_plus_abx2$variable == "CfxA family class A broad-spectrum beta-lactamase"] <- "CfxA"
top_taxa_plus_abx2$variable[top_taxa_plus_abx2$variable == "macrolide efflux MFS transporter Mef(En2)"] <- "MefE"
top_taxa_plus_abx2$variable[top_taxa_plus_abx2$variable == "tetracycline resistance ribosomal protection protein Tet(Q)"] <- "TetQ"

species$Gene <- as.character(species$Gene)
species$Gene[species$Gene == "23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(F)"] <- "ErmF"
species$Gene[species$Gene == "CfxA family class A broad-spectrum beta-lactamase"] <- "CfxA"
species$Gene[species$Gene == "macrolide efflux MFS transporter Mef(En2)"] <- "MefE"
species$Gene[species$Gene == "tetracycline resistance ribosomal protection protein Tet(Q)"] <- "TetQ"

#For defining line types
legend_key <- c("Microbe" = "solid", "Antibiotic resistance gene" = "dashed")

#For loop to generate 9 figures for each taxa vs antibiotic
for (i in 1:length(species$Gene)) {
  assign(paste0("df_", i), top_taxa_plus_abx2[top_taxa_plus_abx2$variable %in% c(paste0(species[i,1]), paste0(species[i,2])),])
  assign(paste0("fig_", i), ggplot() +
    geom_line(data = get(paste0("df_", i))[get(paste0("df_", i))[,2] == paste0(species[i,1]),], aes(x = as.factor(Timepoint), y = as.numeric(value), group = variable, linetype = "Microbe"), color = "black") +
    geom_line(data = get(paste0("df_", i))[get(paste0("df_", i))[,2] == paste0(species[i,2]),], aes(x = as.factor(Timepoint), y = as.numeric(value),  group = variable, linetype = "Antibiotic resistance gene"), color = "black") +
    geom_point(data = get(paste0("df_", i)), aes(x = as.factor(Timepoint), y = as.numeric(value), color = Treatment), size = 2) +
    facet_grid(Individual ~ .) +
    theme_bw() +
    theme(plot.title = element_text(size=9), legend.position = "none", axis.title = element_blank()) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Rel. ab. of Antibiotic resistance gene")) +
    scale_linetype_manual(values = legend_key) +
    labs(title = paste(paste0(species[i,2]), "vs.", paste0(species[i,1])), x = "Time point", y = "Rel. ab. of microbe", linetype = "Line type"))
}

#Redo any plot without legend hiding parameter, then grab its legend to plot in cowplot.
Plot_legend <- get_legend(ggplot() +
  geom_line(data = get(paste0("df_", i))[get(paste0("df_", i))[,2] == paste0(species[i,1]),], aes(x = as.factor(Timepoint), y = as.numeric(value), group = variable, linetype = "Microbe"), color = "black") +
  geom_line(data = get(paste0("df_", i))[get(paste0("df_", i))[,2] == paste0(species[i,2]),], aes(x = as.factor(Timepoint), y = as.numeric(value),  group = variable, linetype = "Antibiotic resistance gene"), color = "black") +
  geom_point(data = get(paste0("df_", i)), aes(x = as.factor(Timepoint), y = as.numeric(value), color = Treatment), size = 2) +
  facet_grid(Individual ~ .) +
  theme_bw() +
  theme(plot.title = element_text(size=9)) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Rel. ab. of Antibiotic resistance gene")) +
  scale_linetype_manual(values = legend_key) +
  labs(title = paste(paste0(species[i,2]), "vs.", paste0(species[i,1])), x = "Time point", y = "Relative abundance of microbe", linetype = "Line type"))

#plot all figures together.
plot2 <- plot_grid(fig_1, fig_2, fig_3, fig_4, fig_5, fig_6, fig_7, fig_8, fig_9, nrow = 3, ncol = 4, rel_widths = c(1,1,1,1))
plot2

ggsave(filename = "figure_2_arg_tracking_Zymo_data.svg", plot = plot2, width = 9, height = 7.5) #Write to file
ggsave(filename = "figure_2_legend.svg", plot = Plot_legend, width = 3, height = 2)
