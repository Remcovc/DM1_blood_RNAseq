# Original script made by R van Cruchten to visualize the association of
# DM1-Activ-c with gene expression
# Slight modifications by D van As for publication purposes and updated results
# Last changes applied on 24/11/21


###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(scales)


###########################
## Loading relevant data ##
###########################
# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "samples.RDATA")
# load relevant mixed effect model data
load("outcome_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)


##########################################################################################
## Generate dataframe with relevant p-values, FDR-values and effect sizes for CBT model ##
##########################################################################################

DM1_df <- data.frame(outcome_coef["DM1ActivC"])

df <- data.frame(hgnc_symbol = DM1_df$DM1ActivC.DM1ActivC.hgnc_symbol,
                 ENSG = DM1_df$DM1ActivC.DM1ActivC.ENSG,
                 p.value = DM1_df$DM1ActivC.DM1ActivC.p.val,
                 FDR = DM1_df$DM1ActivC.DM1ActivC.FDR,
                 effect = DM1_df$DM1ActivC.DM1ActivC.Estimate)

##################
## volcano plot ##
##################

A <- ggplot(df, aes(x = effect, y = -log10(p.value)))+
  geom_point(size= 1, col = ifelse(df$FDR < 0.05, "black","grey")) +
  geom_label_repel(label = ifelse(df$ENSG %in%  df$ENSG[order(df$p.value)][1:4], df$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 5, segment.size = 0.25, fontface="italic") +
  xlab("DM1-Activ-c effect size") +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 4), 
                     breaks = seq(0, 4, 1),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-0.05, 0.05), 
                     breaks = seq(-0.05, 0.05, 0.025)) +
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
        plot.title = element_text(color = "black", size = 16),
        plot.tag = element_text(color ="black", size= 20, face="bold"),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"))


################
## Gene plots ##
################

gene_plots <- list()
#select genes with lowest pvalues
for (ENSG_ID in df$ENSG[order(df$p.value)][1:4]){
  
  df2 <- data.frame(samples[,"Visit"],v[["E"]][ENSG_ID,], samples[,"PatientID"], samples[,"DM1ActivC"])
  names(df2) <- c("Visit","counts","PatientID","DM1ActivC")
  df2 <- df2[order(df2$Visit),]
  
  gene_plots[[ENSG_ID]] <- ggplot(df2, aes_string(x="DM1ActivC", y="counts")) + ggtitle(hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id==ENSG_ID])+
    xlab("DM1-Activ-c score") +
    ylab("logCPM") +
    geom_smooth(method ="lm", formula =  y ~ x, se=F, col ="black") +
    geom_point(colour= ifelse(df2$Visit == "V2", "blue","red")) +
    scale_x_continuous(limits = c(20, 100),
                       breaks=seq(20 , 100, 20),
                       sec.axis = dup_axis(breaks = derive(), labels = NULL, name = NULL))+
    scale_y_continuous(sec.axis = dup_axis(breaks = derive(), labels = NULL, name = NULL))+
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      aspect.ratio = 0.4,
      panel.background = element_rect(fill="white"),
      panel.grid.major.x = element_line(size = 0.25, color = "grey"),
      panel.grid.major.y = element_line(size = 0.25, color = "grey"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 14),
      axis.title.x = element_text(color = "black", size = 14),
      plot.title = element_text(color = "black", size = 16, face="italic"),
      plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"))
  
}

plots <- plot_grid(plotlist = gene_plots, ncol = 1,
                   labels= c("B", "", "",""), label_size = 20, label_y = 0.92)

## Cast subplots into one and save
gs <- list(A, plots)
lay <- rbind(c(1,1,2),
             c(1,1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="FigS2_DM1ActivC.jpeg", height = 10, width = 10, dpi = 600 )







