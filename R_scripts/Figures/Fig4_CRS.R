# Original script made by R van Cruchten to visualize the association of
# the Compound Response scores with gene expression
# Slight modifications by D van As for publication purposes and updated results
# Last changes applied on 22/09/22


###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)     #V 3.3.3
library(ggrepel)     #V 0.9.1
library(cowplot)     #V 1.1.1
library(scales)      #V 1.1.1
library("gridExtra") #V 2.3
library("psych")     #V 4.0.5
library("grid")      #V 4.0.4 

###########################
## Loading relevant data ##
###########################
# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "samples.RDATA")
# load Compound Response Scores mixed effects model
load("scaled_response_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)


##########################################################################################
## Generate dataframe with relevant p-values, FDR-values and effect sizes for CBT model ##
##########################################################################################

CRS_df <- data.frame(lmer_fit_values["scaled_response"])
df <- data.frame(hgnc_symbol = CRS_df$scaled_response.hgnc_symbol,
                 ENSG = CRS_df$scaled_response.ENSG,
                 p.value = CRS_df$scaled_response.p.val,
                 FDR = CRS_df$scaled_response.FDR,
                 effect = CRS_df$scaled_response.Estimate)


## Volcano plot
A <- ggplot(df, aes(x = effect, y = -log10(p.value)))+
  geom_point(size= 1, col = ifelse(df$FDR < 0.05, "black","grey")) +
  geom_label_repel(label = ifelse(df$ENSG %in%  df$ENSG[order(df$p.value)][1:4], df$hgnc_symbol, ""), max.overlaps = 100,  force = 25, color = "black", size = 6, segment.size = 0.25, fontface="italic") +
  xlab("Compound Response effect size") +
  ylab("-10log(p-value)") +
  scale_y_continuous(limits = c(0, 8), 
                     breaks = seq(0, 8, 1),
                     sec.axis = dup_axis(breaks = derive(), labels = derive(), name = NULL))+
  scale_x_continuous(limits = c(-3, 3), 
                     breaks = seq(-3, 3, 1)) +
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
        plot.tag = element_text(color ="black", size= 20, face="bold"),
        plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm")
  )


## calculate change in expression all genes between baseline and 10M CBT
dE <- list()
for (patient in unique(samples$PatientID)){
  dE[[patient]]  <- v$E[, colnames(v$E) %in% paste0(patient, "_V4")] - v$E[, colnames(v$E) %in% paste0(patient, "_V2")]
} 
dE <- do.call(cbind, dE)
colnames(dE) <- unique(samples$PatientID)

gene_plots <- list()
#select genes with lowest pvalues

for (ENSG_ID in df$ENSG[order(df$p.value)][1:4]){
  df2 <- data.frame(dCounts = dE[ENSG_ID,],
                    response = as.numeric(unique(samples[,"scaled_mean_outcome"])))
  pcor <- corr.test(df2$response, df2$dCounts, method="pearson")
  
  gene_plots[[ENSG_ID]] <- ggplot(df2,  aes_string(x="response", y="dCounts")) +  ggtitle(paste(hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id==ENSG_ID]))+
    geom_point() +
    xlab("Compound Response Score") +
    ylab("Delta logCPM") +
    geom_hline(yintercept=0, col = "black") +
    geom_vline(xintercept=0, col = "black") +
    geom_smooth(method ="lm", formula =  y ~ x, se=F, col ="black") +
    scale_x_continuous(lim = c(-0.5, 1.5),
                       sec.axis = dup_axis(breaks = derive(), labels = NULL, name = NULL))+
    scale_y_continuous(labels = label_number(accuracy = 0.1), 
                       limits = c(-1.5, 1.75)) +
    annotation_custom(grobTree(textGrob(
      paste0("Rho = ", round(pcor$r, 2)), 
      x=0.55, y=0.90, just = "left",
      gp=gpar(fontsize=16))))+
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      aspect.ratio = 0.4,
      panel.background = element_rect(fill="white"),
      panel.grid.major.x = element_line(size = 0.25, color = "grey"),
      panel.grid.major.y = element_line(size = 0.25, color = "grey"),
      axis.text.x = element_text(color = "black", size = 16),
      axis.text.y = element_text(color = "black", size = 16),
      axis.title.y = element_text(color = "black", size = 18, face="bold"),
      axis.title.x = element_text(color = "black", size = 18, face="bold"),
      plot.title = element_text(color = "black", size = 18, face="italic"),
      plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm")
      )
}

plots <- plot_grid(plotlist = gene_plots, ncol = 1,
                   labels= c("B", "", "",""), label_size = 20, label_y = 0.92)

gs <- list(A, plots)
lay <- rbind(c(1,1,2),
             c(1,1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, 
       file ="Fig4_CRS_volc_genes.png", 
       height = 13, width = 13, 
       dpi = 1200,
       device ="png")




