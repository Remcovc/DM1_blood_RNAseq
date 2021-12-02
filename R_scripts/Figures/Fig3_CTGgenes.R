# Original script made by R van Cruchten to visualize the association of
# the CTG repeat length with gene expression
# Slight modifications by D van As for publication purposes and updated results
# Last changes applied on 24/11/21


###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)     #V 3.3.3
library(ggrepel)     #V 0.9.1
library(cowplot)     #V 1.1.1
library(scales)      #V 1.1.1
library("gridExtra") #V 2.3


###########################
## Loading relevant data ##
###########################
# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "samples.RDATA")
# load CTG mixed effects model
load("CTG_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)


##########################################################################################
## Generate dataframe with relevant p-values, FDR-values and effect sizes for CBT model ##
##########################################################################################

CTG_df <- data.frame(lmer_fit_values["CTG"])
df <- data.frame(hgnc_symbol = CTG_df$CTG.hgnc_symbol,
                 ENSG = CTG_df$CTG.ENSG,
                 p.value = CTG_df$CTG.p.val,
                 FDR = CTG_df$CTG.FDR,
                 effect = CTG_df$CTG.Estimate)


## Volcano plot
A <- ggplot(df, aes(x = effect*100, y = -log10(p.value)))+
  geom_point(size= 1, col = ifelse(df$FDR < 0.05, "black","grey")) +
  geom_label_repel(label = ifelse(df$ENSG %in%  df$ENSG[order(df$p.value)][1:4], df$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 5, segment.size = 0.25, fontface ="italic") +
  xlab("CTG repeat effect size (per 100 CTGs)") +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 6), 
                     breaks = seq(0, 7, 1),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-0.5, 0.5), 
                     breaks = seq(-0.5, 0.5, 0.25)) +
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
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
        plot.tag = element_text(color ="black", size= 20, face="bold")
  )


## Example plots

gene_plots <- list()
#select genes with lowest pvalues
for (ENSG_ID in df$ENSG[order(df$p.value)][1:4]){
  
  outcome <- "V2Mode"
  # for V2Mode, create a scatter plot with expression level vs CTG-repeat length
  df2 <- data.frame(as.numeric(samples[,outcome][!is.na(samples[,outcome])]), v$E[ENSG_ID,], 
                    samples[,"Visit"][!is.na(samples[,outcome])], samples[,"PatientID"][!is.na(samples[,outcome])])
  names(df2) <- c(outcome,"counts","Visit","PatientID")
  df2 <- df2[order(df2$Visit),]
  
  gene_plots[[ENSG_ID]] <- ggplot(df2, aes_string(x="V2Mode", y="counts")) +
    ggtitle(paste(hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id %in% ENSG_ID])) +
    xlab("CTG repeat length") +
    ylab("logCPM") +
    geom_point(aes_string(x="V2Mode", y="counts", group="PatientID"),
               size = 1, col=ifelse(df2$Visit == "V2", "blue","red")) +
    geom_smooth(method ="lm", formula =  y ~ x, se=F, col ="black") +
    scale_x_continuous(breaks=seq(0,1000,250), 
                       limits = c(0,1000)) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      aspect.ratio = 0.4,
      panel.background = element_rect(fill="white"),
      panel.grid.major.x = element_line(size = 0.25, color = "grey"),
      panel.grid.major.y = element_line(size = 0.25, color = "grey"),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(color = "black", size = 14),
      axis.title.x = element_text(color = "black", size=14),
      plot.title = element_text(color = "black", size = 16, face="italic"),
      plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
    )
}

plots <- plot_grid(plotlist = gene_plots, ncol = 1,
                   labels= c("B", "", "",""), label_size = 20, label_y = 0.92)

## Cast subplots into one and save
gs <- list(A, plots)
lay <- rbind(c(1,1,2),
             c(1,1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="Fig3_CTG_volc_genes.jpeg", height = 10, width = 10, dpi = 600 )




