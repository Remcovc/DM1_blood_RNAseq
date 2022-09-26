
# Original script made by R van Cruchten to visualize differential gene 
# expression after CBT in the OPTIMISTIC intervention group
# Modifications by D van As for publication purposes
# Last changes applied on 19/07/22

###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library("ggplot2")      #V 3.3.5
library("ggrepel")      #V 0.9.1
library("cowplot")      #V 1.1.1
library("scales")       #V 1.1.1
library("heatmap3")     #V 1.1.9
library("cluster")      #V 2.1.2
library("RColorBrewer") #V 1.1.2
library("gridExtra")    #V 2.3
library("psych")        #V 2.1.9
library("grid")         #V 4.0.4 
library("corrplot")     #V 0.92

###########################
## Loading relevant data ##
###########################
# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/for_paper_counts/metadata/samples.RDATA")
# load CBT mixed effects model
load("Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/for_paper_counts/CBT_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/for_paper_counts/v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("Z:/Analysis_projects/RvC_recognition_rnaseq_diff_expr/for_paper_counts/metadata/ENSG_geneSymbol.txt", sep =",", header=TRUE)


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
  geom_label_repel(label = ifelse(df$ENSG %in%  df$ENSG[order(df$p.value)][1:4], df$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 6, segment.size = 0.25, fontface ="italic") +
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
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title.x = element_text(color ="black", size = 18, face="bold"),
        axis.title.y = element_text(color = "black", size = 18, face="bold"),
        plot.title = element_text(color = "black", size = 18),
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
  df2$time <- df2$Visit
  df2$time[df2$time == "V2"] <- 0
  df2$time[df2$time == "V4"] <- 1
  df2$time <- as.numeric(df2$time)
  pcor <- corr.test(df2$time, df2$counts, method="pearson")
  
  
  gene_plots[[ENSG_ID]] <- ggplot(df2,  aes_string(x="Visit", y="counts",group="PatientID")) +    
    ggtitle(paste(hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id==ENSG_ID]))+
    geom_point(col = ifelse(df2$Visit == "V2", "blue", "red")) +
    geom_path(col= "black") + 
    ylab("logCPM") +
    xlab("") +
    scale_x_discrete(expand = c(0.1,0), labels = c("Baseline","10 Mo. CBT"))+ 
    scale_y_continuous(labels = label_number(accuracy = 0.1), 
                       limits = c(min(df2$counts), max(df2$counts)+0.6)) +
    annotation_custom(grobTree(textGrob(
      paste0("Rho = ", round(pcor$r, 2)), 
      x=0.05, y=0.92, just = "left",
      gp=gpar(fontsize=16))))+
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      aspect.ratio = 1,
      panel.background = element_rect(fill="white"),
      panel.grid.major.x = element_line(size = 0.25, color = "grey"),
      panel.grid.major.y = element_line(size = 0.25, color = "grey"),
      axis.text.x = element_text(color = "black", size = 16, face="bold"),
      axis.text.y = element_text(color = "black", size = 16),
      axis.title.x = element_text(color ="black", size = 18),
      axis.title.y = element_text(color = "black", size = 16, face="bold"),
      plot.title = element_text(color = "black", size = 18, face="italic"),
      plot.margin = margin(c(0.05,0.8,0.05,0.05), unit="cm"),
      plot.tag = element_text(color ="black", size= 20, face="bold")
      )
}
# Arrange plots and save 
plots <- plot_grid(plotlist = gene_plots, nrow = 1, scale=0.9)

#####################################
## Cast subplots into one and save ##
#####################################
# Heatmap will be added separately

gs <- list(A, plots)
lay <- rbind(c(1,1,3,3,3),
             c(1,1,3,3,3),
             c(1,1,3,3,3),
             c(2,2,2,2,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, 
       file ="Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/Publication_GitHub/Revision updates/Figures/Fig2_CBT_genes.png", 
       height = 10, width = 11, 
       dpi = 1200,
       device ="png")


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
load("Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/for_paper_counts/CBT_coef.RDATA")
CBT <- data.frame(lmer_fit_values[1])
sig_genes <- CBT$CBT.ENSG[CBT$CBT.FDR < 0.05]
delta_counts_sig <- delta_counts[rownames(delta_counts) %in% sig_genes,]

hmcol <- colorRampPalette(brewer.pal(6, "YlOrRd"))(6)

reds <- colorRampPalette(rev(brewer.pal(9, "Reds")))(9)
blues <- colorRampPalette(brewer.pal(9, "Blues"))(9)
hmcol <- c(reds, blues)

# create heatmap and save in high res (800*1200)
map <- heatmap3(delta_counts_sig, scale="row", 
                col=hmcol, hclustfun=hclust,
                showRowDendro = FALSE, labRow = FALSE, labCol = FALSE,
                main="",
                margin=c(2,2))

png("Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/Publication_GitHub/Revision updates/Figures/Fig2_Heatmap.png",
    width=400,
    height=600,
    res=300,
    units='mm')

heatmap3(delta_counts_sig, scale="row", 
         col=hmcol, hclustfun=hclust,
         showRowDendro = FALSE, labRow = FALSE, labCol = FALSE,
         main="",
         margin=c(2,2))
dev.off()




###############################
## Patient response addition ##
###############################

## Sort patient data by IDs
samples <- samples[order(samples$PatientID),]

## Split baseline and 10 month assessments
dfV2 <- samples[samples$Visit == "V2",]
dfV4 <- samples[samples$Visit == "V4",]

table(colnames(dfV2) == colnames(dfV4))
table(dfV2$PatientID == dfV4$PatientID)

## Calculate deltas and store in matrix
dDM1ActivC <- as.numeric(dfV4$DM1ActivC) - as.numeric(dfV2$DM1ActivC)
dSMWT <- as.numeric(dfV4$SMWT) - as.numeric(dfV2$SMWT)
CRS <- as.numeric(dfV2$scaled_mean_outcome)
ddf <- rbind(dDM1ActivC, dSMWT, CRS)
colnames(ddf) <- dfV2$PatientID

## Reorder patients based on gene expression dendogram
c_ord <- colnames(delta_counts_sig)[map$colInd]
ddf <- ddf[,c_ord]

## Scale dDM1ActivC and dSMWT
ddf[1,] <- scale(ddf[1,], center=F)
ddf[2,] <- scale(ddf[2,], center=F)

## Visualize results & save

png("Z:/Analysis_projects/DvA_recognition_rnaseq_diff_expr/Publication_GitHub/Revision updates/Figures/Fig2_Heatmap_cor.png",
    width=600,
    height=200,
    res=300,
    units='mm')

corrplot(ddf,
         is.corr = F,
         col = hmcol,
         method="color",
         tl.pos = 'l',
         cl.pos = "n")
dev.off()

