##Master's Thesis
##João Lopes
set.seed(1)
setwd("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Transcriptomics")

##Load a R package
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(seqinr)
library(data.table)
library(tidyverse)
library(tximport)
library(readxl)
library(tidyr)
library(multiGSEA)
windowsFonts(Calibri = windowsFont("Arial"))

##Create folders
dir.create("Blast_Drerio")
dir.create("Blast_Drerio/Files")
dir.create("Blast_Drerio/Files/Tables")
dir.create("Blast_Drerio/Plots")
dir.create("Blast_Drerio/Plots/Dotplot")

dir.create("Blast_SwissProt")
dir.create("Blast_SwissProt/Files")
dir.create("Blast_SwissProt/Files/Tables")
dir.create("Blast_SwissProt/Plots")
dir.create("Blast_SwissProt/Plots/Heatmap")
dir.create("Blast_SwissProt/Plots/PieCharts")

dir.create("Trinity")
dir.create("Trinity/Files")
dir.create("Trinity/Files/Tables")
dir.create("Trinity/Plots")
dir.create("Trinity/Plots/Volcano plots")
dir.create("Trinity/Plots/Heatmap")
dir.create("Trinity/Plots/PieCharts")

##Directories
##Folders
pathFiles_Trinity <- "Trinity/Files/" 
pathTables_Trinity <- "Trinity/Files/Tables/" 
pathVolcano_Trinity <- "Trinity/Plots/Volcano plots/"
pathHeatmap_Trinity <- "Trinity/Plots/Heatmap/"
pathPieCharts_Trinity <- "Trinity/Plots/PieCharts/"

pathFiles_BSwiss <- "Blast_SwissProt/Files/" 
pathTables_BSwiss <- "Blast_SwissProt/Files/Tables/" 
pathPieCharts_BSwiss <- "Blast_SwissProt/Plots/PieCharts/"
pathHeatmap_BSwiss <- "Blast_SwissProt/Plots/Heatmap/"

pathFiles_BDrerio <- "Blast_Drerio/Files/" 
pathTables_BDreio <- "Blast_Drerio/Files/Tables/"
pathDotplot <- "Blast_Drerio/Plots/Dotplot/"

pathBothOmicsPathways <- "C:/Users/jdpl2/OneDrive/Ambiente de Trabalho/Mestrado/2º Ano/Both Omics/"

#######################################Functions#######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

relativeExpression<-function(treatment1, treatment2, Design, fit){
  contrast<-paste(treatment1, "-", treatment2, sep = "")
  treatment1vtreatment2<-makeContrasts(contrasts = contrast, levels = Design)
  print(treatment1vtreatment2)
  lrt<-glmLRT(fit, contrast = treatment1vtreatment2)
  print(topTags(lrt))
  return(lrt)
}

volcanoPlot<-function(pathFiles_Trinity, colour, legend, title, contrast, i){
  load(paste(pathFiles_Trinity, "lrt.RData", sep = ""))
  print(colnames(DEG_trinity))
  signif<--log10(DEG_trinity[,paste("Pvalue-", contrast, sep = "")])
  plot(DEG_trinity[,paste("logFC-", contrast, sep = "")], signif, pch = ".", xlab = expression(log[2] * FC), ylab = expression(log[10] * FDR), main = NULL, 
       panel.first={
         points(0, 0, pch = 16, cex = 1e6, col = "grey95")
         grid(col = "white", lty = 1)
       }
  )
  points(DEG_trinity[which(DEG_trinity[contrast] == 1), paste("logFC-", contrast, sep = "")], -log10(DEG_trinity[which(DEG_trinity[contrast] == 1), paste("Pvalue-", contrast, sep = "")]), pch = 20, cex = 1.5, col = colour[1,])
  points(DEG_trinity[which(DEG_trinity[contrast] == -1), paste("logFC-", contrast, sep = "")], -log10(DEG_trinity[which(DEG_trinity[contrast] == -1), paste("Pvalue-", contrast, sep = "")]), pch = 20, cex = 1.5, col = colour[2,])
}

heatmapGraph<-function(colCutoff, rowCutoff, data, path){
  heatData<-as.matrix(data[,c(2,6,10,14,18)])
  rownames(heatData)<-data$Description
  head(heatData)
  clustFunction <- function(x) hclust(x, method="complete")
  distFunction <- function(x) dist(x,method="euclidean")
  #clusterint (sidebars)
  cCol<-colorRampPalette(brewer.pal(9, "Set1"))
  colFit<-clustFunction(distFunction(t(heatData)))
  cClusters<-cutree(colFit,h=colCutoff)
  cHeight<-length(unique(as.vector(cClusters)))
  colHeight = cCol(cHeight)
  cRow<-colorRampPalette(brewer.pal(9,"Set1"))
  rowFit<-clustFunction(distFunction(heatData))
  rClusters<-cutree(rowFit,h=rowCutoff)
  rHeight<-length(unique(as.vector(rClusters)));
  rowHeight = cRow(rHeight)
  xDimension<-25
  yDimension<-25
  windows(xDimension, yDimension)
  par(family="Arial")
  #Heatmap
  heatColour<-colorRampPalette(c("#fff04a", "#f28e2b", "#c23e40"))(n = 100)
  heatmap.2(heatData,
            hclust=clustFunction, 
            distfun=distFunction, 
            ColSideColors=colHeight[cClusters],
            RowSideColors=rowHeight[rClusters],
            density.info="none",
            col=heatColour, 
            trace="none", 
            scale="row", 
            cexCol = 0.8,
            lhei = c(1,5),
            lwid = c(1.5,5), 
            offsetRow = 0, 
            offsetCol = 0,
            srtCol = 45, 
            key.title = NA, 
            margins = c(8,13),
            labRow = "")
  recordPlot()
  dev.print(tiff, paste(path, "Heatmap", ".tif", sep =""), height = 15, width = 15,  units = 'cm', res=600)
}

pieChart <- function(data, i,type,contrasts,colour){
  if (type == "transcripts"){
    row_data <- data[i, ]
    long_data <- data.frame(
      Type = c("Overexpressed transcripts", "Underexpressed transcripts"),
      Count = c(row_data$`Overexpressed in treatment 2`, row_data$`Underexpressed in treatment 2`)
    )
    p <- ggplot(long_data, aes(x = "", y = Count, fill = Type)) +
      geom_bar(width = 1, stat = "identity", color = "black", linewidth = 1,linetype = "solid") +
      geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "white", size = 12) + # Add numbers inside the pie chart
      coord_polar(theta = "y") +
      theme_void() +
      scale_fill_manual(name = NULL, values = c("Overexpressed transcripts" = colour[,i][1], "Underexpressed transcripts" = colour[,i][2])) + 
      labs(title = NULL) +
      theme(legend.position = "bottom", legend.spacing.y = unit(0.1, "cm"),
            legend.text = element_text(size = 15))
    ggsave(paste(pathPieCharts_Trinity, "PieChart-", row_data[1], ".svg", sep =""), plot = p, height = 15, width = 15,  units = 'cm', dpi=600)
    ggsave(paste(pathPieCharts_Trinity, "PieChart-", row_data[1], ".tif", sep =""), plot = p, height = 15, width = 15,  units = 'cm', dpi=600)
  }
  else if (type == "ORFs"){
    row_data <- data[i, ]
    long_data <- data.frame(
      Type = c("Overexpressed ORFs", "Underexpressed ORFs"),
      Count = c(row_data$`Overexpressed in treatment 2`, row_data$`Underexpressed in treatment 2`)
    )
    p <- ggplot(long_data, aes(x = "", y = Count, fill = Type)) +
      geom_bar(width = 1, stat = "identity", color = "black", linewidth = 1, linetype = "dashed") +
      geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "white", size = 12) +
      coord_polar(theta = "y") +
      theme_void() +
      scale_fill_manual(name = NULL, values = c("Overexpressed ORFs" = colour[,i][1], "Underexpressed ORFs" = colour[,i][2])) +
      labs(title = NULL) +
      theme(legend.position = "bottom", legend.spacing.y = unit(0.1, "cm"),
            legend.text = element_text(size = 15))
    ggsave(paste(pathPieCharts_BSwiss, "PieChart-", row_data[1], ".svg", sep =""), plot = p, height = 15, width = 15,  units = 'cm', dpi=600)
    ggsave(paste(pathPieCharts_BSwiss, "PieChart-", row_data[1], ".tif", sep =""), plot = p, height = 15, width = 15,  units = 'cm', dpi=600)
  }
}

dotplot <- function (data,path) {
  o <- ggplot(data, aes(x = contrast, y = pathway)) +
    geom_point(aes(size = 1-padj, color = NES, fill = NES), shape = 21, stroke = 0.5) +  
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_size_continuous(range = c(3, 15)) +
    labs(x = "Contrast", y = "Pathway", 
         size = "1-p.adjust", 
         color = "NES",
         fill = "NES") +  
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),  
          legend.position = "bottom",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          panel.grid.major = element_line(linewidth = 1.2, color = "grey90"),  
          panel.grid.minor = element_line(linewidth = 0.8, color = "grey90"),  
          panel.background = element_rect(fill = "grey", color = NA))  
  ggsave(filename = paste(path, "Dotplot.svg", sep = ""), plot = o, device = "svg", width = 15, height = 15)
  ggsave(filename = paste(path, "Dotplot.tiff", sep = ""), plot = o, device = "tiff", width = 15, height = 15)
}

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Functions#######################################



##################################Objects/Descriptions#################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

titleCTL_10vMHW2_10<-expression("CTL day 10 vs. MHW2 day 10")
titleCTL_25vMHW1_25<-expression("CTL day 25 vs. MHW1 day 25")
titleCTL_25vMHW2_25<-expression("CTL day 25 vs. MHW2 day 25")
titleMHW1_25vMHW2_25<-expression("MHW1 day 25 vs. MHW2 day 25")
titleMHW2_10vMHW2_25<-expression("MHW2 day 10 vs. MHW2 day 25")

legend<-data.frame(CTL_10vMHW2_10 = c("Overexpressed in the control group in day 10", "Overexpressed in the MHW2 group in day 10"),
                   CTL_25vMHW1_25 = c("Overexpressed in the control group in day 25", "Overexpressed in the MHW1 group in day 25"),
                   CTL_25vMHW2_25 = c("Overexpressed in the control group in day 25","Overexpressed in MHW2 group in day 25"),
                   MHW1_25vMHW2_25 = c("Overexpressed in the MHW1 group in day 25","Overexpressed in MHW2 group in day 25"),
                   MHW2_10vMHW2_25 = c("Overexpressed in the MHW2 group in day 10","Overexpressed in MHW2 group in day 25"))
contrasts <- c("CTL_10vMHW2_10","CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25","MHW2_10vMHW2_25")

colour<-data.frame(CTL_10vMHW2_10 = c("#FFA500", "#0000FF"), 
                   CTL_25vMHW1_25 = c("#FF0000", "#006400"),
                   CTL_25vMHW2_25 = c("#B8860B", "#800080"),
                   MHW1_25vMHW2_25 = c("#654321", "#1E90FF"),
                   MHW2_10vMHW2_25 = c("#4B0082", "#696969")) 

sampleName <- c(
  "T13MHW2Pm1 MD25(M)",
  "T7MHW2Pm1 MD25(M)",
  "T3MHW2Pm1 MD10(M)",
  "T6MHW1Pm1 MD25(M)",
  "T1CTLPm1 MD10(M)",
  "T10CTLPm1 MD10(M)",
  "T1CTLPm1 MD25(M)",
  "T4CTLPm1 MD25(M)",
  "T12MHW1Pm1 MD25(M)",
  "T2MHW1Pm1 MD25(M)",
  "T10CTLPm1 MD25(M)",
  "T15CTLPm1 MD10(M)",
  "T3MHW2Pm1 MD25(M)",
  "T7MHW2Pm1 MD10(M)",
  "T5MHW2Pm1 MD10(M)"
)
sampleID <- c("RNA_11_N3066",
               "RNA_13_N3067",
               "RNA_14_N3058",
               "RNA_17_N3064",
               "RNA_18_N3054",
               "RNA_19_N3055",
               "RNA_20_N3060",
               "RNA_21_N3061",
               "RNA_3_N3062",
               "RNA_4_N3063",
               "RNA_5_N3059",
               "RNA_6_N3053",
               "RNA_7_N3065",
               "RNA_8_N3056",
               "RNA_9_N3057")
Day <- c(
  "25",
  "25",
  "10",
  "25",
  "10",
  "10",
  "25",
  "25",
  "25",
  "25",
  "25",
  "10",
  "25",
  "10",
  "10"
)
Treatment <- c(
  "MHW2",
  "MHW2",
  "MHW2",
  "MHW1",
  "CTL",
  "CTL",
  "CTL",
  "CTL",
  "MHW1",
  "MHW1",
  "CTL",
  "CTL",
  "MHW2",
  "MHW2",
  "MHW2"
)
Levels <- c(
  "MHW2_25",
  "MHW2_25",
  "MHW2_10",
  "MHW1_25",
  "CTL_10",
  "CTL_10",
  "CTL_25",
  "CTL_25",
  "MHW1_25",
  "MHW1_25",
  "CTL_25",
  "CTL_10",
  "MHW2_25",
  "MHW2_10",
  "MHW2_10"
)

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
##################################Objects/Descriptions#################################



#######################################Files/Data######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

#Kallisto folder with the abundance files
Kallistofolder<-"C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Outputs\\Server\\Kallisto_50\\abundance\\"
#########################################

#Trinity file
Trinityfile<-"C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Outputs\\Server\\Trinity\\Trinity.Trinity.fasta"
#############

#Blast with swissprot
Blast <- read.csv("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Outputs\\Server\\Blastp\\Blastp_for_analysis\\Pmicrops_Blastp_SwissProt.csv",sep =";",header = F)
#####################

#Blast with only drerio matches
Blast_drerio <- read.csv("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Outputs\\Server\\Blast_Drerio\\Pmicrops_Blastp_Drerio.csv",sep =";",header = F)
#############

#Transdecoder file
TransDecoder <- fread("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Outputs\\Server\\TransDecoder\\Trinity.Trinity.fasta.transdecoder.pep",header=FALSE)
##################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Files/Data######################################



####################Merge the counts with the respective cDNA sequence#################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Import the abundance results and estimate the counts
samples<-dir(Kallistofolder)
file<-c(paste(sep = "", Kallistofolder, samples, "\\", "abundance.h5"))
names(file)=sampleID
print(names(file))
Tsv<-tximport(file, type = "kallisto", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
print(head(Tsv$counts))
print(dim(Tsv$counts))
print(colnames(Tsv$counts))
######################################################
  
##Read the fasta file to get the ids and sequences
Fa<-read.fasta(Trinityfile, as.string = TRUE)
print(head(Fa))
TxID<-getName(Fa)
TxSEQ<-unlist(getSequence(Fa, as.string = TRUE))
TxSEQ<-as.data.frame(cbind(TxID, TxSEQ))
print(head(TxSEQ))
colnames(TxSEQ)[1]<-"ID"
print(nrow(TxSEQ))
Counts<-Tsv$counts
Counts<-cbind(as.data.frame(row.names(Counts)), Counts)
colnames(Counts)[1]<-"ID"
##################################################

##Create a matrix with  blast ids, sequence and counts
colnames(Blast) <- c("TrinityID","Accession", "%ID", "Aligment_Length", "Mismatches", "Gap_Openings", "Start_of_Aligment_TrinityID", "End_of_Aligment_TrinityID", "Start_of_Aligment_Accession", "End_of_Aligment_Accession", "Evalue", "Bit_Score","Accession","-2","Description")
Blast <- subset(Blast, grepl("^TRINITY", Blast[,1]))
Blast$ID<-unlist(sapply(Blast$TrinityID, function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
blast_trinity_merge<-merge(Blast,Counts, by = "ID")
blast_trinity_merge <- blast_trinity_merge[, -c(2,4:17)]
blast_trinity_merge$Accession <- str_extract(blast_trinity_merge$Accession, "(?<=\\|)[^|]+(?=\\|)")
Full_blast<-merge(blast_trinity_merge, TxSEQ, by="ID")
write.table(Full_blast, paste(pathBothOmicsPathways, "genes_ID.csv", sep = ""),sep = ";",col.names = TRUE, row.names = FALSE)
######################################################

##Create a matrix with  trinity ids, sequence and counts
Full_trinity <- merge(TxSEQ,Counts, by = "ID")
colnames(Full_trinity)<-c("ID", "sequence", sampleID)
########################################################

##Create a materix with blast_drerio ids, sequence and counts
colnames(Blast_drerio) <- c("TrinityID","Accession", "%ID", "Aligment_Length", "Mismatches", "Gap_Openings", "Start_of_Aligment_TrinityID", "End_of_Aligment_TrinityID", "Start_of_Aligment_Accession", "End_of_Aligment_Accession", "Evalue", "Bit_Score","Accession","-2","Description")
Blast_drerio <- subset(Blast_drerio, grepl("^TRINITY", Blast_drerio[,1]))
Blast_drerio$ID<-unlist(sapply(Blast_drerio$TrinityID, function(x) unlist(strsplit(x, ".", fixed = TRUE))[1]))
blast_drerio_trinity_merge<-merge(Blast_drerio,Counts, by = "ID")
blast_drerio_trinity_merge <- blast_drerio_trinity_merge[, -c(1,2,4:17)]
sum(duplicated(blast_drerio_trinity_merge$Accession))

##Duplicate with higher mean
w <- c()
for (i in c(1:length(unique(blast_drerio_trinity_merge$Accession)))){
  dups <- which(blast_drerio_trinity_merge$Accession %in% unique(blast_drerio_trinity_merge$Accession)[i])
  v <- c()
  for (j in dups){
    v <- c(v,mean(as.numeric(blast_drerio_trinity_merge[-c(1,2)][j,])))
  }
  w = c(w,dups[which(v == max(v))[1]])
}
blast_drerio_trinity_merge <- blast_drerio_trinity_merge[w,]
############################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
####################Merge the counts with the respective cDNA sequence#################



################################Statistics (with trinity)##############################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create a DGEList object from a table of counts 
Data_trinity <-
  DGEList(counts = Full_trinity[, 3:17],
          group = Levels,
          genes = Full_trinity[, 1])
################################################

##Filter out lowly expressed genes
keep_trinity <- filterByExpr(Data_trinity)
Data_trinity <- Data_trinity[keep_trinity, , keep.lib.sizes=FALSE]
##################################

##Calculate the normalization factors
Data_trinity <- calcNormFactors(Data_trinity)
#####################################

##Create the generalized linear models for comparing the expression levels from the different treatments
Design_trinity <- model.matrix( ~ 0 + Levels, data = Data_trinity$samples)
colnames(Design_trinity) <- levels(Data_trinity$samples$group)
Data_trinity <- estimateDisp(Data_trinity, Design_trinity)
fit_trinity <- glmFit(Data_trinity, Design_trinity)
########################################################################################################

##Expression levels between the treatments
CTL_10vMHW2_10_trinity<-relativeExpression(colnames(Design_trinity)[1],colnames(Design_trinity)[4],Design_trinity,fit_trinity)
MHW2_10vMHW2_25_trinity<-relativeExpression(colnames(Design_trinity)[4],colnames(Design_trinity)[5],Design_trinity,fit_trinity)
CTL_25vMHW1_25_trinity<-relativeExpression(colnames(Design_trinity)[2],colnames(Design_trinity)[3],Design_trinity,fit_trinity)
CTL_25vMHW2_25_trinity<-relativeExpression(colnames(Design_trinity)[2],colnames(Design_trinity)[5],Design_trinity,fit_trinity)
MHW1_25vMHW2_25_trinity<-relativeExpression(colnames(Design_trinity)[3],colnames(Design_trinity)[5],Design_trinity,fit_trinity)
##########################################

##Merge the logFC,logCPM,LR and pvalue with the trinity ids
Results_trinity<-cbind(Data_trinity$genes, CTL_10vMHW2_10_trinity$table, MHW2_10vMHW2_25_trinity$table, CTL_25vMHW1_25_trinity$table, CTL_25vMHW2_25_trinity$table, MHW1_25vMHW2_25_trinity$table)
colnames(Results_trinity)<-c(
  "ID",
  "logFC-CTL_10vMHW2_10",
  "logCPM-CTL_10vMHW2_10",
  "LR-CTL_10vMHW2_10",
  "Pvalue-CTL_10vMHW2_10",
  "logFC-MHW2_10vMHW2_25",
  "logCPM-MHW2_10vMHW2_25",
  "LR-MHW2_10vMHW2_25",
  "Pvalue-MHW2_10vMHW2_25",
  "logFC-CTL_25vMHW1_25",
  "logCPM-CTL_25vMHW1_25",
  "LR-CTL_25vMHW1_25",
  "Pvalue-CTL_25vMHW1_25",
  "logFC-CTL_25vMHW2_25",
  "logCPM-CTL_25vMHW2_25",
  "LR-CTL_25vMHW2_25",
  "Pvalue-CTL_25vMHW2_25",
  "logFC-MHW1_25vMHW2_25",
  "logCPM-MHW1_25vMHW2_25",
  "LR-MHW1_25vMHW2_25",
  "Pvalue-MHW1_25vMHW2_25"
)
head(Results_trinity)
write.table(Results_trinity,paste(pathFiles_Trinity,"Results.csv",sep =""), sep=";",col.names=NA)
#######################################

##Calculating the nº of DEGs in general 
expressionTable_trinity<-decideTests(Results_trinity[,grepl("Pvalue", colnames(Results_trinity))],coefficients = Results_trinity[,grepl("logFC", colnames(Results_trinity))], adjust.method = "fdr", lfc = 1.5)
expressionTable_trinity<-as.data.frame(expressionTable_trinity)
head(expressionTable_trinity)
colnames(expressionTable_trinity)<-c("CTL_10vMHW2_10", "MHW2_10vMHW2_25", "CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25")
DEG_trinity<-cbind(Results_trinity,expressionTable_trinity)
degTable_trinity<-DEG_trinity[which(abs(DEG_trinity$CTL_10vMHW2_10) == 1 | abs(DEG_trinity$MHW2_10vMHW2_25) == 1 | abs(DEG_trinity$CTL_25vMHW1_25) == 1 | abs(DEG_trinity$CTL_25vMHW2_25) == 1 | abs(DEG_trinity$MHW1_25vMHW2_25) == 1),]
head(degTable_trinity)
nrow(degTable_trinity)
write.table(degTable_trinity,paste(pathFiles_Trinity,"DEG.csv",sep =""),sep=";",col.names=NA)
save(CTL_10vMHW2_10_trinity, MHW2_10vMHW2_25_trinity, CTL_25vMHW1_25_trinity, CTL_25vMHW2_25_trinity, MHW1_25vMHW2_25_trinity, file = "Trinity\\Files\\lrt.RData")
#######################################

##Nº of DEGs per treatment
DEG_CTL_10vMHW2_10_trinity<-sum(degTable_trinity$CTL_10vMHW2_10 != 0);DEG_CTL_10vMHW2_10_trinity
DEG_MHW2_10vMHW2_25_trinity<-sum(degTable_trinity$MHW2_10vMHW2_25 != 0);DEG_MHW2_10vMHW2_25_trinity
DEG_CTL_25vMHW2_25_trinity<-sum(degTable_trinity$CTL_25vMHW2_25 != 0);DEG_CTL_25vMHW2_25_trinity
DEG_CTL_25vMHW1_25_trinity<-sum(degTable_trinity$CTL_25vMHW1_25 != 0);DEG_CTL_25vMHW1_25_trinity
DEG_MHW1_25vMHW2_25_trinity<-sum(degTable_trinity$MHW1_25vMHW2_25 != 0);DEG_MHW1_25vMHW2_25_trinity
##########################

##Nº of overexpressed genes in treatment 2
oedeg_CTL_10vMHW2_10_trinity<-sum(degTable_trinity$CTL_10vMHW2_10 == -1);oedeg_CTL_10vMHW2_10_trinity
oedeg_MHW2_10vMHW2_25_trinity<-sum(degTable_trinity$MHW2_10vMHW2_25 == -1);oedeg_MHW2_10vMHW2_25_trinity
oedeg_CTL_25vMHW2_25_trinity<-sum(degTable_trinity$CTL_25vMHW2_25 == -1);oedeg_CTL_25vMHW2_25_trinity
oedeg_CTL_25vMHW1_25_trinity<-sum(degTable_trinity$CTL_25vMHW1_25 == -1);oedeg_CTL_25vMHW1_25_trinity
oedeg_MHW1_25vMHW2_25_trinity<-sum(degTable_trinity$MHW1_25vMHW2_25 == -1);oedeg_MHW1_25vMHW2_25_trinity
#############################################

##Nº of underexpressed genes in treatment 
uedeg_CTL_10vMHW2_10_trinity<-sum(degTable_trinity$CTL_10vMHW2_10 == 1);uedeg_CTL_10vMHW2_10_trinity
uedeg_MHW2_10vMHW2_25_trinity<-sum(degTable_trinity$MHW2_10vMHW2_25 == 1);uedeg_MHW2_10vMHW2_25_trinity
uedeg_CTL_25vMHW2_25_trinity<-sum(degTable_trinity$CTL_25vMHW2_25 == 1);uedeg_CTL_25vMHW2_25_trinity
uedeg_CTL_25vMHW1_25_trinity<-sum(degTable_trinity$CTL_25vMHW1_25 == 1);uedeg_CTL_25vMHW1_25_trinity
uedeg_MHW1_25vMHW2_25_trinity<-sum(degTable_trinity$MHW1_25vMHW2_25 == 1);uedeg_MHW1_25vMHW2_25_trinity
############################################

##Table with total differential expressed transcripts, underexpressed transcripts and overexpressed transcripts per treatment
deg_per_treatment_trinity <- data.frame(Treatments = contrasts, 
                                No_of_Genes = c(DEG_CTL_10vMHW2_10_trinity,DEG_CTL_25vMHW1_25_trinity,DEG_CTL_25vMHW2_25_trinity,DEG_MHW1_25vMHW2_25_trinity,DEG_MHW2_10vMHW2_25_trinity),
                                No_of_underexpressed_Genes = c(uedeg_CTL_10vMHW2_10_trinity,uedeg_CTL_25vMHW1_25_trinity,uedeg_CTL_25vMHW2_25_trinity,uedeg_MHW1_25vMHW2_25_trinity,uedeg_MHW2_10vMHW2_25_trinity),
                                No_of_overexpressed_Genes = c(oedeg_CTL_10vMHW2_10_trinity,oedeg_CTL_25vMHW1_25_trinity,oedeg_CTL_25vMHW2_25_trinity,oedeg_MHW1_25vMHW2_25_trinity,oedeg_MHW2_10vMHW2_25_trinity))
colnames(deg_per_treatment_trinity) <- c("Treatments","transcripts with differential expression","Overexpressed in treatment 2","Underexpressed in treatment 2")
write.table(deg_per_treatment_trinity,paste(pathTables_Trinity,"DEG_per_treatment_trinity.csv",sep =""),sep=";",col.names=NA)
########################################################################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
###############################Statistics (with trinity)###############################



################################Statistics (with blast)################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create a DGEList object from a table of counts 
Data_blast <-
  DGEList(counts = Full_blast[, 3:17],
          group = Levels,
          genes = Full_blast[, 2])
################################################

##Filter out lowly expressed genes
keep_blast <- filterByExpr(Data_blast)
Data_blast <- Data_blast[keep_blast, , keep.lib.sizes=FALSE]
##################################

##Calculate the normalization factors
Data_blast <- calcNormFactors(Data_blast)
#####################################

##Create the generalized linear models for comparing the expression levels from the different treatments
Design_blast <- model.matrix( ~ 0 + Levels, data = Data_blast$samples)
colnames(Design_blast) <- levels(Data_blast$samples$group)
Data_blast <- estimateDisp(Data_blast, Design_blast)
fit_blast <- glmFit(Data_blast, Design_blast)
########################################################################################################

##Expression levels between the treatments
CTL_10vMHW2_10_blast<-relativeExpression(colnames(Design_blast)[1],colnames(Design_blast)[4],Design_blast,fit_blast)
MHW2_10vMHW2_25_blast<-relativeExpression(colnames(Design_blast)[4],colnames(Design_blast)[5],Design_blast,fit_blast)
CTL_25vMHW1_25_blast<-relativeExpression(colnames(Design_blast)[2],colnames(Design_blast)[3],Design_blast,fit_blast)
CTL_25vMHW2_25_blast<-relativeExpression(colnames(Design_blast)[2],colnames(Design_blast)[5],Design_blast,fit_blast)
MHW1_25vMHW2_25_blast<-relativeExpression(colnames(Design_blast)[3],colnames(Design_blast)[5],Design_blast,fit_blast)
##########################################

##Merge the logFC,logCPM,LR and pvalue with the trinity ids
Results_blast<-cbind(Data_blast$genes, CTL_10vMHW2_10_blast$table, MHW2_10vMHW2_25_blast$table, CTL_25vMHW1_25_blast$table, CTL_25vMHW2_25_blast$table, MHW1_25vMHW2_25_blast$table)
colnames(Results_blast)<-c(
  "ID",
  "logFC-CTL_10vMHW2_10",
  "logCPM-CTL_10vMHW2_10",
  "LR-CTL_10vMHW2_10",
  "Pvalue-CTL_10vMHW2_10",
  "logFC-MHW2_10vMHW2_25",
  "logCPM-MHW2_10vMHW2_25",
  "LR-MHW2_10vMHW2_25",
  "Pvalue-MHW2_10vMHW2_25",
  "logFC-CTL_25vMHW1_25",
  "logCPM-CTL_25vMHW1_25",
  "LR-CTL_25vMHW1_25",
  "Pvalue-CTL_25vMHW1_25",
  "logFC-CTL_25vMHW2_25",
  "logCPM-CTL_25vMHW2_25",
  "LR-CTL_25vMHW2_25",
  "Pvalue-CTL_25vMHW2_25",
  "logFC-MHW1_25vMHW2_25",
  "logCPM-MHW1_25vMHW2_25",
  "LR-MHW1_25vMHW2_25",
  "Pvalue-MHW1_25vMHW2_25"
)
head(Results_blast)
write.table(Results_blast,paste(pathFiles_BSwiss,"Results.csv",sep =""), sep=";",col.names=NA)
#######################################

##Calculating the nº of DEGs in general 
expressionTable_blast<-decideTests(Results_blast[,grepl("Pvalue", colnames(Results_blast))],coefficients = Results_blast[,grepl("logFC", colnames(Results_blast))], adjust.method = "fdr", lfc = 1.5)
expressionTable_blast<-as.data.frame(expressionTable_blast)
head(expressionTable_blast)
colnames(expressionTable_blast)<-c("CTL_10vMHW2_10", "MHW2_10vMHW2_25", "CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25")
DEG_blast<-cbind(Results_blast,expressionTable_blast)
degTable_blast<-DEG_blast[which(abs(DEG_blast$CTL_10vMHW2_10) == 1 | abs(DEG_blast$MHW2_10vMHW2_25) == 1 | abs(DEG_blast$CTL_25vMHW1_25) == 1 | abs(DEG_blast$CTL_25vMHW2_25) == 1 | abs(DEG_blast$MHW1_25vMHW2_25) == 1),]
head(degTable_blast)
nrow(degTable_blast)
write.table(degTable_blast,paste(pathBothOmicsPathways,"DEG.csv",sep =""),sep=";",col.names=NA)
save(CTL_10vMHW2_10_blast, MHW2_10vMHW2_25_blast, CTL_25vMHW1_25_blast, CTL_25vMHW2_25_blast, MHW1_25vMHW2_25_blast, file = "Blast_SwissProt\\Files\\lrt.RData")
#######################################

##Nº of DEGs per treatment
DEG_CTL_10vMHW2_10_blast<-sum(degTable_blast$CTL_10vMHW2_10 != 0);DEG_CTL_10vMHW2_10_blast
DEG_MHW2_10vMHW2_25_blast<-sum(degTable_blast$MHW2_10vMHW2_25 != 0);DEG_MHW2_10vMHW2_25_blast
DEG_CTL_25vMHW2_25_blast<-sum(degTable_blast$CTL_25vMHW2_25 != 0);DEG_CTL_25vMHW2_25_blast
DEG_CTL_25vMHW1_25_blast<-sum(degTable_blast$CTL_25vMHW1_25 != 0);DEG_CTL_25vMHW1_25_blast
DEG_MHW1_25vMHW2_25_blast<-sum(degTable_blast$MHW1_25vMHW2_25 != 0);DEG_MHW1_25vMHW2_25_blast
##########################

##Nº of overexpressed ORFs in treatment 2
oedeg_CTL_10vMHW2_10_blast<-sum(degTable_blast$CTL_10vMHW2_10 == -1);oedeg_CTL_10vMHW2_10_blast
oedeg_MHW2_10vMHW2_25_blast<-sum(degTable_blast$MHW2_10vMHW2_25 == -1);oedeg_MHW2_10vMHW2_25_blast
oedeg_CTL_25vMHW2_25_blast<-sum(degTable_blast$CTL_25vMHW2_25 == -1);oedeg_CTL_25vMHW2_25_blast
oedeg_CTL_25vMHW1_25_blast<-sum(degTable_blast$CTL_25vMHW1_25 == -1);oedeg_CTL_25vMHW1_25_blast
oedeg_MHW1_25vMHW2_25_blast<-sum(degTable_blast$MHW1_25vMHW2_25 == -1);oedeg_MHW1_25vMHW2_25_blast
#############################################

##Nº of underexpressed ORFs in treatment 2 
uedeg_CTL_10vMHW2_10_blast<-sum(degTable_blast$CTL_10vMHW2_10 == 1);uedeg_CTL_10vMHW2_10_blast
uedeg_MHW2_10vMHW2_25_blast<-sum(degTable_blast$MHW2_10vMHW2_25 == 1);uedeg_MHW2_10vMHW2_25_blast
uedeg_CTL_25vMHW2_25_blast<-sum(degTable_blast$CTL_25vMHW2_25 == 1);uedeg_CTL_25vMHW2_25_blast
uedeg_CTL_25vMHW1_25_blast<-sum(degTable_blast$CTL_25vMHW1_25 == 1);uedeg_CTL_25vMHW1_25_blast
uedeg_MHW1_25vMHW2_25_blast<-sum(degTable_blast$MHW1_25vMHW2_25 == 1);uedeg_MHW1_25vMHW2_25_blast
############################################

##Table with total diferantial expressed ORFs, underexpressed ORFs and overexpressed ORFs per treatment
deg_per_treatment_blast <- data.frame(Treatments = contrasts, 
                                      No_of_Genes = c(DEG_CTL_10vMHW2_10_blast,DEG_CTL_25vMHW1_25_blast,DEG_CTL_25vMHW2_25_blast,DEG_MHW1_25vMHW2_25_blast,DEG_MHW2_10vMHW2_25_blast),
                                      No_of_overexpressed_Genes = c(oedeg_CTL_10vMHW2_10_blast,oedeg_CTL_25vMHW1_25_blast,oedeg_CTL_25vMHW2_25_blast,oedeg_MHW1_25vMHW2_25_blast,oedeg_MHW2_10vMHW2_25_blast),
                                      No_of_underexpressed_Genes = c(uedeg_CTL_10vMHW2_10_blast,uedeg_CTL_25vMHW1_25_blast,uedeg_CTL_25vMHW2_25_blast,uedeg_MHW1_25vMHW2_25_blast,uedeg_MHW2_10vMHW2_25_blast))
colnames(deg_per_treatment_blast) <- c("Treatments","transcripts with differential expression","Overexpressed in treatment 2","Underexpressed in treatment 2")
write.table(deg_per_treatment_blast,paste(pathTables_BSwiss,"DEG_per_treatment_blast.csv",sep =""),sep=";",col.names=NA)
########################################################################################

##Table with nº of transcripts, nº of ORFs, nº of differential expressed transcripts and nº of ORFs with differential expression
transdecoder_n_transcripts <- grepl("^>", TransDecoder$V1)
transdecoder_filtered <- subset(TransDecoder, transdecoder_n_transcripts)
p_ORFs <- round((nrow(transdecoder_filtered)*100)/nrow(TxSEQ),3)
p_annotated_ORFs <- round((nrow(Blast)*100)/nrow(TxSEQ),3)
p_DETranscripts <- round((nrow(degTable_trinity)*100)/nrow(TxSEQ),3)
p_DEORFs <- round((nrow(degTable_blast)*100)/nrow(TxSEQ),3)
percentages <- c(100,p_ORFs,p_annotated_ORFs,p_DETranscripts,p_DEORFs)
percentages <- paste0(percentages, "%")
transcripts_of_interest <- data.frame(Stage = c(1,2,3,4,5),
                                      Objective = c("Transcriptome assembly","Transcripts with coding regions","Functional annotation","Differentially-Expressed Transcripts","ORFs with Differential Expression"),
                                      Specie = c(nrow(TxSEQ),nrow(transdecoder_filtered),nrow(Blast),nrow(degTable_trinity),nrow(degTable_blast)),
                                      Percentagens = percentages)
colnames(transcripts_of_interest) <- c("Analytical stage","Stage objective","Transcripts","Percentages")
write.table(transcripts_of_interest,paste(pathTables_BSwiss,"transcripts_of_interest.csv",sep =""),sep=";",col.names=NA)
################################################################################################################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
################################Statistics (with blast)################################



#############################Statistics (with blast_drerio)############################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create a DGEList object from a table of counts 
Data_blast_drerio <-
  DGEList(counts = blast_drerio_trinity_merge[, 2:16],
          group = Levels,
          genes = blast_drerio_trinity_merge[, 1])
################################################

##Filter out lowly expressed genes
keep_blast_drerio <- filterByExpr(Data_blast_drerio)
Data_blast_drerio <- Data_blast_drerio[keep_blast_drerio, , keep.lib.sizes=FALSE]
##################################

##Calculate the normalization factors
Data_blast_drerio <- calcNormFactors(Data_blast_drerio)
#####################################

##Create the generalized linear models for comparing the expression levels from the different treatments
Design_blast_drerio <- model.matrix( ~ 0 + Levels, data = Data_blast_drerio$samples)
colnames(Design_blast_drerio) <- levels(Data_blast_drerio$samples$group)
Data_blast_drerio <- estimateDisp(Data_blast_drerio, Design_blast_drerio)
fit_blast_drerio <- glmFit(Data_blast_drerio, Design_blast_drerio)
########################################################################################################

##Expression levels between the treatments
CTL_10vMHW2_10_blast_drerio<-relativeExpression(colnames(Design_blast_drerio)[1],colnames(Design_blast_drerio)[4],Design_blast_drerio,fit_blast_drerio)
MHW2_10vMHW2_25_blast_drerio<-relativeExpression(colnames(Design_blast_drerio)[4],colnames(Design_blast_drerio)[5],Design_blast_drerio,fit_blast_drerio)
CTL_25vMHW1_25_blast_drerio<-relativeExpression(colnames(Design_blast_drerio)[2],colnames(Design_blast_drerio)[3],Design_blast_drerio,fit_blast_drerio)
CTL_25vMHW2_25_blast_drerio<-relativeExpression(colnames(Design_blast_drerio)[2],colnames(Design_blast_drerio)[5],Design_blast_drerio,fit_blast_drerio)
MHW1_25vMHW2_25_blast_drerio<-relativeExpression(colnames(Design_blast_drerio)[3],colnames(Design_blast_drerio)[5],Design_blast_drerio,fit_blast_drerio)
##########################################

##Merge the logFC,logCPM,LR and pvalue with the trinity ids
Results_blast_drerio<-cbind(Data_blast_drerio$genes, CTL_10vMHW2_10_blast_drerio$table, MHW2_10vMHW2_25_blast_drerio$table, CTL_25vMHW1_25_blast_drerio$table, CTL_25vMHW2_25_blast_drerio$table, MHW1_25vMHW2_25_blast_drerio$table)
colnames(Results_blast_drerio)<-c(
  "ID",
  "logFC-CTL_10vMHW2_10",
  "logCPM-CTL_10vMHW2_10",
  "LR-CTL_10vMHW2_10",
  "Pvalue-CTL_10vMHW2_10",
  "logFC-MHW2_10vMHW2_25",
  "logCPM-MHW2_10vMHW2_25",
  "LR-MHW2_10vMHW2_25",
  "Pvalue-MHW2_10vMHW2_25",
  "logFC-CTL_25vMHW1_25",
  "logCPM-CTL_25vMHW1_25",
  "LR-CTL_25vMHW1_25",
  "Pvalue-CTL_25vMHW1_25",
  "logFC-CTL_25vMHW2_25",
  "logCPM-CTL_25vMHW2_25",
  "LR-CTL_25vMHW2_25",
  "Pvalue-CTL_25vMHW2_25",
  "logFC-MHW1_25vMHW2_25",
  "logCPM-MHW1_25vMHW2_25",
  "LR-MHW1_25vMHW2_25",
  "Pvalue-MHW1_25vMHW2_25"
)
head(Results_blast_drerio)
write.table(Results_blast_drerio,paste(pathFiles_BDrerio,"Results.csv",sep =""), sep=";",col.names=NA)
#######################################

##Calculating the nº of DEGs in general 
expressionTable_blast_drerio<-decideTests(Results_blast_drerio[,grepl("Pvalue", colnames(Results_blast_drerio))],coefficients = Results_blast_drerio[,grepl("logFC", colnames(Results_blast_drerio))], adjust.method = "fdr", lfc = 1.5)
expressionTable_blast_drerio<-as.data.frame(expressionTable_blast_drerio)
head(expressionTable_blast_drerio)
colnames(expressionTable_blast_drerio)<-c("CTL_10vMHW2_10", "MHW2_10vMHW2_25", "CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25")
DEG_blast_drerio<-cbind(Results_blast_drerio,expressionTable_blast_drerio)
degTable_blast_drerio<-DEG_blast_drerio[which(abs(DEG_blast_drerio$CTL_10vMHW2_10) == 1 | abs(DEG_blast_drerio$MHW2_10vMHW2_25) == 1 | abs(DEG_blast_drerio$CTL_25vMHW1_25) == 1 | abs(DEG_blast_drerio$CTL_25vMHW2_25) == 1 | abs(DEG_blast_drerio$MHW1_25vMHW2_25) == 1),]
head(degTable_blast_drerio)
nrow(degTable_blast_drerio)
write.table(degTable_blast_drerio,paste(pathFiles_BDrerio,"DEG.csv",sep =""),sep=";",col.names=NA)
save(CTL_10vMHW2_10_blast_drerio, MHW2_10vMHW2_25_blast_drerio, CTL_25vMHW1_25_blast_drerio, CTL_25vMHW2_25_blast_drerio, MHW1_25vMHW2_25_blast_drerio, file = "Blast_Drerio\\Files\\lrt.RData")
#######################################

##Nº of DEGs per treatment
DEG_drerio_CTL_10vMHW2_10<-sum(degTable_blast_drerio$CTL_10vMHW2_10 != 0);DEG_drerio_CTL_10vMHW2_10
DEG_drerio_MHW2_10vMHW2_25<-sum(degTable_blast_drerio$MHW2_10vMHW2_25 != 0);DEG_drerio_MHW2_10vMHW2_25
DEG_drerio_CTL_25vMHW2_25<-sum(degTable_blast_drerio$CTL_25vMHW2_25 != 0);DEG_drerio_CTL_25vMHW2_25
DEG_drerio_CTL_25vMHW1_25<-sum(degTable_blast_drerio$CTL_25vMHW1_25 != 0);DEG_drerio_CTL_25vMHW1_25
DEG_drerio_MHW1_25vMHW2_25<-sum(degTable_blast_drerio$MHW1_25vMHW2_25 != 0);DEG_drerio_MHW1_25vMHW2_25
##########################

##Nº of underexpressed genes per treatment 
uedeg_drerio_CTL_10vMHW2_10<-sum(degTable_blast_drerio$CTL_10vMHW2_10 == -1);uedeg_drerio_CTL_10vMHW2_10
uedeg_drerio_MHW2_10vMHW2_25<-sum(degTable_blast_drerio$MHW2_10vMHW2_25 == -1);uedeg_drerio_MHW2_10vMHW2_25
uedeg_drerio_CTL_25vMHW2_25<-sum(degTable_blast_drerio$CTL_25vMHW2_25 == -1);uedeg_drerio_CTL_25vMHW2_25
uedeg_drerio_CTL_25vMHW1_25<-sum(degTable_blast_drerio$CTL_25vMHW1_25 == -1);uedeg_drerio_CTL_25vMHW1_25
uedeg_drerio_MHW1_25vMHW2_25<-sum(degTable_blast_drerio$MHW1_25vMHW2_25 == -1);uedeg_drerio_MHW1_25vMHW2_25
#############################################

##Nº of overexpressed genes per treatment 
oedeg_drerio_CTL_10vMHW2_10<-sum(degTable_blast_drerio$CTL_10vMHW2_10 == 1);oedeg_drerio_CTL_10vMHW2_10
oedeg_drerio_MHW2_10vMHW2_25<-sum(degTable_blast_drerio$MHW2_10vMHW2_25 == 1);oedeg_drerio_MHW2_10vMHW2_25
oedeg_drerio_CTL_25vMHW2_25<-sum(degTable_blast_drerio$CTL_25vMHW2_25 == 1);oedeg_drerio_CTL_25vMHW2_25
oedeg_drerio_CTL_25vMHW1_25<-sum(degTable_blast_drerio$CTL_25vMHW1_25 == 1);oedeg_drerio_CTL_25vMHW1_25
oedeg_drerio_MHW1_25vMHW2_25<-sum(degTable_blast_drerio$MHW1_25vMHW2_25 == 1);oedeg_drerio_MHW1_25vMHW2_25
############################################

##Table with total DEG, underexpressed genes and overexpressed genes per treatment
deg_drerio_per_treatment <- data.frame(Treatments = contrasts, 
                                No_of_Genes = c(DEG_drerio_CTL_10vMHW2_10,DEG_drerio_CTL_25vMHW1_25,DEG_drerio_CTL_25vMHW2_25,DEG_drerio_MHW1_25vMHW2_25,DEG_drerio_MHW2_10vMHW2_25),
                                No_of_underexpressed_Genes = c(uedeg_drerio_CTL_10vMHW2_10,uedeg_drerio_CTL_25vMHW1_25,uedeg_drerio_CTL_25vMHW2_25,uedeg_drerio_MHW1_25vMHW2_25,uedeg_drerio_MHW2_10vMHW2_25),
                                No_of_overexpressed_Genes = c(oedeg_drerio_CTL_10vMHW2_10,oedeg_drerio_CTL_25vMHW1_25,oedeg_drerio_CTL_25vMHW2_25,oedeg_drerio_MHW1_25vMHW2_25,oedeg_drerio_MHW2_10vMHW2_25))
colnames(deg_drerio_per_treatment) <- c("Treatments","Nº of proteins","Nº of underexpressed proteins","Nº of overexpressed proteins")
write.table(deg_drerio_per_treatment,paste(pathTables_BSwiss,"DEG_per_treatment.csv",sep =""),sep=";",col.names=NA)
########################################################################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#############################Statistics (with blast_drerio)############################



####################################Pathway analysis###################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

DEG_pathway <- DEG_blast_drerio[-c(22:26)]
DEG_pathway$ID <- str_extract(DEG_pathway$ID, "\\|([^|]+)\\|") %>%
  str_replace_all("\\|", "")
write.table(DEG_pathway, paste(pathBothOmicsPathways, "DEG_pathway.csv", sep = ""),sep = ";",col.names = TRUE, row.names = FALSE)

##Create data frames with Accession, logFC and Pvalue for the pathway analysis
allg_CTL_10vMHW2_10<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC-CTL_10vMHW2_10`,DEG_pathway$`Pvalue-CTL_10vMHW2_10`)
names(allg_CTL_10vMHW2_10) <- c("Accession","logFC","Pvalue")

allg_CTL_25vMHW1_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC-CTL_25vMHW1_25`,DEG_pathway$`Pvalue-CTL_25vMHW1_25`)
names(allg_CTL_25vMHW1_25) <- c("Accession","logFC","Pvalue")

allg_CTL_25vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC-CTL_25vMHW2_25`,DEG_pathway$`Pvalue-CTL_25vMHW2_25`)
names(allg_CTL_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allg_MHW1_25vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC-MHW1_25vMHW2_25`,DEG_pathway$`Pvalue-MHW1_25vMHW2_25`)
names(allg_MHW1_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allg_MHW2_10vMHW2_25<- data.frame(DEG_pathway$ID,DEG_pathway$`logFC-MHW2_10vMHW2_25`,DEG_pathway$`Pvalue-MHW2_10vMHW2_25`)
names(allg_MHW2_10vMHW2_25) <- c("Accession","logFC","Pvalue")
############################################################################

##Create a data structure
layers <- c("transcriptome")
odataCTL_10vMHW2_10 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW1_25 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW1_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW2_10vMHW2_25 <- initOmicsDataStructure(layer=layers)
#########################

##Add transcriptome layer
odataCTL_10vMHW2_10$transcriptome <- rankFeatures(allg_CTL_10vMHW2_10$logFC,allg_CTL_10vMHW2_10$Pvalue)
names(odataCTL_10vMHW2_10$transcriptome) <- allg_CTL_10vMHW2_10$Accession
odataCTL_10vMHW2_10$transcriptome <- sort(odataCTL_10vMHW2_10$transcriptome)
head(odataCTL_10vMHW2_10$transcriptome)

odataCTL_25vMHW1_25$transcriptome <- rankFeatures(allg_CTL_25vMHW1_25$logFC,allg_CTL_25vMHW1_25$Pvalue)
names(odataCTL_25vMHW1_25$transcriptome) <- allg_CTL_25vMHW1_25$Accession
odataCTL_25vMHW1_25$transcriptome <- sort(odataCTL_25vMHW1_25$transcriptome)
head(odataCTL_25vMHW1_25$transcriptome)

odataCTL_25vMHW2_25$transcriptome <- rankFeatures(allg_CTL_25vMHW2_25$logFC,allg_CTL_25vMHW2_25$Pvalue)
names(odataCTL_25vMHW2_25$transcriptome) <- allg_CTL_25vMHW2_25$Accession
odataCTL_25vMHW2_25$transcriptome <- sort(odataCTL_25vMHW2_25$transcriptome)
head(odataCTL_25vMHW2_25$transcriptome)

odataMHW1_25vMHW2_25$transcriptome <- rankFeatures(allg_MHW1_25vMHW2_25$logFC,allg_MHW1_25vMHW2_25$Pvalue)
names(odataMHW1_25vMHW2_25$transcriptome) <- allg_MHW1_25vMHW2_25$Accession
odataMHW1_25vMHW2_25$transcriptome <- sort(odataMHW1_25vMHW2_25$transcriptome)
head(odataMHW1_25vMHW2_25$transcriptome)

odataMHW2_10vMHW2_25$transcriptome <- rankFeatures(allg_MHW2_10vMHW2_25$logFC,allg_MHW2_10vMHW2_25$Pvalue)
names(odataMHW2_10vMHW2_25$transcriptome) <- allg_MHW2_10vMHW2_25$Accession
odataMHW2_10vMHW2_25$transcriptome <- sort(odataMHW2_10vMHW2_25$transcriptome)
head(odataMHW2_10vMHW2_25$transcriptome)
####################

##Select the databases we want to query and download pathway definitions
databases <- c("kegg")
pathways <- getMultiOmicsFeatures(dbs = databases, layer = layers,
                                  returnTranscriptome = "UNIPROT",
                                  organism = "drerio",
                                  useLocal =  FALSE)
pathways_short <- lapply(names(pathways), function(name) {
  head(pathways[[name]], 2)
})
names(pathways_short) <- names(pathways)
pathways_short
pathways$transcriptome[8]
########################################################################

##Run the pathway enrichment
enrichment_scoresCTL_10vMHW2_10 <- multiGSEA(pathways,odataCTL_10vMHW2_10)
Tenrichment_scoresCTL_10vMHW2_10 <- as.data.frame(enrichment_scoresCTL_10vMHW2_10$transcriptome)
Tenrichment_scoresCTL_10vMHW2_10$leadingEdge <- sapply(Tenrichment_scoresCTL_10vMHW2_10$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_10vMHW2_10 <- Tenrichment_scoresCTL_10vMHW2_10[order(Tenrichment_scoresCTL_10vMHW2_10$padj),]
Tenrichment_scoresCTL_10vMHW2_10$contrast <- rep("CTL_10vMHW2_10")
top_10_esCTL_10vMHW2_10 <- head(Tenrichment_scoresCTL_10vMHW2_10,10)
write.table(top_10_esCTL_10vMHW2_10,paste(pathTables_BDreio,"top10_enrichment_scoresCTL_10vMHW2_10.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW1_25 <- multiGSEA(pathways,odataCTL_25vMHW1_25)
Tenrichment_scoresCTL_25vMHW1_25 <- as.data.frame(enrichment_scoresCTL_25vMHW1_25$transcriptome)
Tenrichment_scoresCTL_25vMHW1_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW1_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW1_25 <- Tenrichment_scoresCTL_25vMHW1_25[order(Tenrichment_scoresCTL_25vMHW1_25$padj),]
Tenrichment_scoresCTL_25vMHW1_25$contrast <- rep("CTL_25vMHW1_25")
top_10_esCTL_25vMHW1_25 <- head(Tenrichment_scoresCTL_25vMHW1_25,10)
write.table(top_10_esCTL_25vMHW1_25,paste(pathTables_BDreio,"top10_enrichment_scoresCTL_25vMHW1_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW2_25 <- multiGSEA(pathways,odataCTL_25vMHW2_25)
Tenrichment_scoresCTL_25vMHW2_25 <- as.data.frame(enrichment_scoresCTL_25vMHW2_25$transcriptome)
Tenrichment_scoresCTL_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW2_25 <- Tenrichment_scoresCTL_25vMHW2_25[order(Tenrichment_scoresCTL_25vMHW2_25$padj),]
Tenrichment_scoresCTL_25vMHW2_25$contrast <- rep("CTL_25vMHW2_25")
top_10_esCTL_25vMHW2_25 <- head(Tenrichment_scoresCTL_25vMHW2_25,10)
write.table(top_10_esCTL_25vMHW2_25,paste(pathTables_BDreio,"top10_enrichment_scoresCTL_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW1_25vMHW2_25 <- multiGSEA(pathways,odataMHW1_25vMHW2_25)
Tenrichment_scoresMHW1_25vMHW2_25 <- as.data.frame(enrichment_scoresMHW1_25vMHW2_25$transcriptome)
Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW1_25vMHW2_25 <- Tenrichment_scoresMHW1_25vMHW2_25[order(Tenrichment_scoresMHW1_25vMHW2_25$padj),]
Tenrichment_scoresMHW1_25vMHW2_25$contrast <- rep("MHW1_25vMHW2_25")
top_10_esMHW1_25vMHW2_25 <- head(Tenrichment_scoresMHW1_25vMHW2_25,10)
write.table(top_10_esMHW1_25vMHW2_25,paste(pathTables_BDreio,"top10_enrichment_scoresMHW1_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW2_10vMHW2_25 <- multiGSEA(pathways,odataMHW2_10vMHW2_25)
Tenrichment_scoresMHW2_10vMHW2_25 <- as.data.frame(enrichment_scoresMHW2_10vMHW2_25$transcriptome)
Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW2_10vMHW2_25 <- Tenrichment_scoresMHW2_10vMHW2_25[order(Tenrichment_scoresMHW2_10vMHW2_25$padj),]
Tenrichment_scoresMHW2_10vMHW2_25$contrast <- rep("MHW2_10vMHW2_25")
top_10_esMHW2_10vMHW2_25 <- head(Tenrichment_scoresMHW2_10vMHW2_25,10)
write.table(top_10_esMHW2_10vMHW2_25,paste(pathTables_BDreio,"top10_enrichment_scoresMHW2_10vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)
############################

##Create dataset with all enrichment scores
combined_df <- rbind(Tenrichment_scoresCTL_10vMHW2_10,Tenrichment_scoresCTL_25vMHW1_25,Tenrichment_scoresCTL_25vMHW2_25,Tenrichment_scoresMHW1_25vMHW2_25,Tenrichment_scoresMHW2_10vMHW2_25)
combined_df_top10 <- rbind(top_10_esCTL_10vMHW2_10,top_10_esCTL_25vMHW1_25,top_10_esCTL_25vMHW2_25,top_10_esMHW1_25vMHW2_25,top_10_esMHW2_10vMHW2_25)
combined_df_top10$pathway <- sub("^\\(KEGG\\) ", "", combined_df_top10$pathway)
###########################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
####################################Pathway analysis###################################



##########################################Plots########################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Volcano Plot
par(mfrow = c(1,1), mar = c(5,5,2,2), family = "Arial")
for(i in 1:length(contrasts)){
  volcanoPlot(pathFiles_Trinity, colour[contrasts[i]], legend[contrasts[i]], get(paste("title", contrasts[i], sep = "")), contrasts[i], i)
  dev.print(tiff, paste(pathVolcano_Trinity, "VolcanoPlot-", contrasts[i], ".tif", sep =""), height = 15, width = 15,  units = 'cm', res=600)
}
##############

##Heatmap
heatmapGraph(100, 20,degTable_trinity,pathHeatmap_Trinity)
heatmapGraph(60,24,degTable_blast_drerio,pathHeatmap_BSwiss)
#########

##Pie charts
for(i in 1:length(contrasts)){
  pieChart(deg_per_treatment_blast,i,"ORFs",contrasts,colour)
  pieChart(deg_per_treatment_trinity,i,"transcripts",contrasts,colour)
}
############

#Dotplot
dotplot(combined_df_top10,pathDotplot)
########

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
##########################################Plots########################################