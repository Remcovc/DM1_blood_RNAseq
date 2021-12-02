
# Original script made by R van Cruchten to visualize differential gene 
# expression after CBT in the OPTIMISTIC intervention group
# Slight modifications by D van As for publication purposes
# Last changes applied on 20/11/21


###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)
library("heatmap3")
library("cluster")
library("RColorBrewer")
library("gridExtra") #V 2.3

###########################
## Loading relevant data ##
###########################
# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "samples.RDATA")
# load CBT mixed effects model
load("CBT_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)


##########################################################################################
## Generate dataframe with relevant p-values, FDR-values and effect sizes for CBT model ##
##########################################################################################

CBT_df <- data.frame(lmer_fit_values["CBT"])
df <- data.frame(hgnc_symbol = CBT_df$CBT.hgnc_symbol,
                 ENSG = CBT_df$CBT.ENSG,
                 p.value = CBT_df$CBT.p.val,
                 FDR = CBT_df$CBT.FDR,
                 effect = CBT_df$CBT.Estimate)

################
# volcano plot #
################
A <- ggplot(df, aes(x = effect, y = -log10(p.value)))+
  geom_point(size= 1, col = ifelse(df$FDR < 0.05, "black","grey")) +
  geom_label_repel(label = ifelse(df$ENSG %in%  df$ENSG[order(df$p.value)][1:4], df$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 5, segment.size = 0.25, fontface ="italic") +
  xlab("CBT effect size") +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 6), 
                     breaks = seq(0, 7, 1),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-1.25, 1.25), 
                     breaks = c(-1.25, -0.75, -0.25, 0.25, 0.75, 1.25)) +
  labs(tag ="A") +
  theme(legend.position = "none", 
        aspect.ratio = 1.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major.x = element_line(size = 0.25, color = "grey"),
        panel.grid.major.y = element_line(size = 0.25, color = "grey"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color ="black", size = 14),
        axis.title.y = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold"))


##############
# Gene plots #
##############

gene_plots <- list()
#select genes with lowest pvalues
for (ENSG_ID in df$ENSG[order(df$p.value)][1:4]){
    
  df2 <- data.frame(samples[,"Visit"], v[["E"]][ENSG_ID,], samples[,"PatientID"])
  names(df2) <- c("Visit","counts","PatientID")
  df2 <- df2[order(df2$Visit),]
  
  gene_plots[[ENSG_ID]] <- ggplot(df2,  aes_string(x="Visit", y="counts",group="PatientID")) +    
    ggtitle(paste(hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id==ENSG_ID]))+
    geom_point(col = ifelse(df2$Visit == "V2", "blue", "red")) +
    geom_path(col= "black") + 
    ylab("logCPM") +
    xlab("") +
    scale_x_discrete(expand = c(0.1,0), labels = c("Baseline","10 Mo. CBT"))+ 
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      aspect.ratio = 1,
      panel.background = element_rect(fill="white"),
      panel.grid.major.x = element_line(size = 0.25, color = "grey"),
      panel.grid.major.y = element_line(size = 0.25, color = "grey"),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(color ="black", size = 14),
      axis.title.y = element_text(color = "black", size = 14),
      plot.title = element_text(color = "black", size = 14, face="italic"),
      plot.margin = margin(c(0.05,0.8,0.05,0.05), unit="cm"),
      plot.tag = element_text(color ="black", size= 20, face="bold")
      )
}
# Arrange plots and save 
plots <- plot_grid(plotlist = gene_plots, nrow = 1, scale=0.9)


#####################################
## Cast subplots into one and save ##
#####################################
# Heatmap will be added separetaly

gs <- list(A, plots)
lay <- rbind(c(1,1,3,3),
             c(1,1,3,3),
             c(1,1,3,3),
             c(2,2,2,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="Fig2_CBT_genes.jpeg", height = 10, width = 10, dpi = 600 )




######################################################
# Heat map of CBT induced changes in gene expression #
######################################################
counts <- v$E
V2_counts <- counts[,grepl("V2", colnames(counts))]
colnames(V2_counts) <- gsub(x = colnames(V2_counts), pattern = "_V2", replacement ="")
V2_counts <- V2_counts[,order(colnames(V2_counts))]

V4_counts <- counts[,grepl("V4", colnames(counts))]
colnames(V4_counts) <- gsub(x = colnames(V4_counts), pattern = "_V4", replacement ="")
V4_counts <- V4_counts[,order(colnames(V4_counts))]

delta_counts <- V4_counts - V2_counts

# Subset significant hits
load("CBT_coef.RDATA")
CBT <- data.frame(lmer_fit_values[1])
sig_genes <- CBT$CBT.ENSG[CBT$CBT.FDR < 0.05]
delta_counts_sig <- delta_counts[rownames(delta_counts) %in% sig_genes,]

hmcol <- colorRampPalette(brewer.pal(6, "YlOrRd"))(6)

reds <- colorRampPalette(rev(brewer.pal(9, "Reds")))(9)
blues <- colorRampPalette(brewer.pal(9, "Blues"))(9)
hmcol <- c(reds, blues)

# heatmap saved through console with 800*1200 resolution
heatmap3(delta_counts_sig, scale="row", 
         col=hmcol, hclustfun=hclust,
         showRowDendro = FALSE, labRow = FALSE, labCol = FALSE,
         main="",
         margin=c(2,2))










