# Original script made by R van Cruchten to visualize the association of
# the CTG repeat length with splice exclusion
# Slight modifications by D van As for publication purposes
# Last changes applied on 06/07/22

###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)     #V 3.3.3
library(ggrepel)     #V 0.9.1
library(cowplot)     #V 1.1.1
library(scales)      #V 1.1.1
library(stringr)     #V 1.4.0
library(gridExtra)   #V 2.3
library(readr)       #V 1.4.0
library(psych)       #V 2.1.9
library(grid)        #V 4.0.4


##################
## Loading data ##
##################

# a table with all samples and their metadata (generated in tableS3_metadata.rmd)
load(file = "samples.RDATA")
# Fitted models for rMATS splice exclusion data as generated in rMATS_SE.rmd 
load(file = "lmfit_rMATS_CTGrepeat.RData")
# extracted values from fits in previous line generated in rMATS_SE.rmd 
load(file = "lmfit_rMATS_CTGrepeat_values.RData")
# raw rMATS output
load(file = "rMATSSE.MATS.JC_filtered.RDATA")

#first load the raw rMATS output into the list
SE.MATS.JC <- list()
SE.MATS.JC[["raw_output"]] <- data.frame(read_delim("SE.MATS.JC.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
df_r <- data.frame(SE.MATS.JC$raw_output)
colnames(df_r)[1] <- "ID"

######################
## Data preparation ##
######################

#cleanup exon names by splitting the label in the unique exon ID and an identifier that indicates the gene and length of exon
lm_fit_values[["V2Mode_SE"]][,"exon_ID"] <- do.call(rbind, str_split(lm_fit_values[["V2Mode_SE"]][,"exon"], "_", 2))[,1]
lm_fit_values[["V2Mode_SE"]][,"exon_short"] <- do.call(rbind, str_split(lm_fit_values[["V2Mode_SE"]][,"exon"], "_", 2))[,2]
lm_fit_values[["V2Mode_SE"]][,"exon_short"] <- gsub("_nt", "nt", lm_fit_values[["V2Mode_SE"]][,"exon_short"])

# create dataframe with effects and significances
df <- data.frame(p.value = lm_fit_values[["V2Mode_SE"]][,"p.val"], 
                 FDR = lm_fit_values[["V2Mode_SE"]][,"FDR"], 
                 effect = lm_fit_values[["V2Mode_SE"]][,"Estimate"]) 


##################
## Volcano plot ##
##################

A <- ggplot(df, aes(x = effect*100, y = -log10(p.value)))+
  geom_label_repel(label = ifelse(lm_fit_values[["V2Mode_SE"]][,"exon"] %in%  
                                    lm_fit_values[["V2Mode_SE"]][,"exon"][order(lm_fit_values[["V2Mode_SE"]][,"p.val"])][1:4], 
                                    lm_fit_values[["V2Mode_SE"]][,"exon_short"], ""), 
                                    max.overlaps = 100, force = 25, color = "black", size = 5, fontface ="italic",
                                    segment.size = 0.25) +
  geom_point(size= 1, col = ifelse(df$FDR < 0.05, "black","grey"))+
  xlab("CTG repeat effect size (per 100 CTGs)")+
  ylab("-10log(p-value)")+
  labs(tag ="A") +
  theme(legend.position = "none", aspect.ratio = 4/3)  +
  scale_x_continuous(limits = c(-0.1, 0.1), 
                     breaks = seq(-0.1, 0.1, 0.025)) +
  scale_y_continuous(limits = c(0, 7), 
                     breaks = seq(0, 7, 1)) +
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
        plot.tag = element_text(color ="black", size= 20, face="bold"))


###################
## Example plots ##
###################


gene_plots <- list()
# exons with lowest pvalues
for (exon_ID in lm_fit_values[["V2Mode_SE"]][,"exon"][order(lm_fit_values[["V2Mode_SE"]][,"p.val"])][1:4]){
  
  # Obtain genomic coorindates of exon
  ID <- substr(exon_ID, 1, 5)
  coord <- paste(df_r$chr[df_r$ID == ID], ":", df_r$exonStart_0base[df_r$ID == ID], "-", df_r$exonEnd[df_r$ID == ID], sep="")
  
  outcome <- "V2Mode"
  # for V2Mode, create a scatter plot with expression level vs CTG-repeat length. V2Mode= modal CTG repeat length at V2
  df2 <- data.frame(as.numeric(samples[,outcome][!is.na(samples[,outcome])]), SE.MATS.JC_filtered[["IncLevel1"]][exon_ID,], samples[,"Visit"][!is.na(samples[,outcome])], samples[,"PatientID"][!is.na(samples[,outcome])])
  names(df2) <- c(outcome,"counts","Visit","PatientID")
  df2 <- df2[order(df2$Visit),]
  pcor <- corr.test(df2$V2Mode, df2$counts, method="pearson")
  pcor$p <- round(pcor$p, 4)
  if (pcor$p < 0.0001){
    pcor$p <- "< 0.0001"
  }
  
  gene_plots[[exon_ID]] <- ggplot(df2, aes_string(x=outcome, y="counts",group="PatientID")) +
    ggtitle(paste(lm_fit_values[["V2Mode_SE"]][lm_fit_values[["V2Mode_SE"]][,"exon"] == exon_ID,"exon_short"])) +
    labs(subtitle = coord) +
    xlab("CTG repeat length") +
    ylab("PSI") +
    geom_point(size = 1, col=ifelse(df2$Visit == "V2", "blue","red")) +
    scale_x_continuous(breaks=seq(0,1000,250), 
                       limits = c(0,1000)) +
    scale_y_continuous(breaks=seq(0, 1, 0.2), 
                       limits = c(0,1),
                       labels = label_number(accuracy = 0.1)) +
    annotation_custom(grobTree(textGrob(
        paste0("Rho = ", round(pcor$r,3)), 
        x=0.05, y=0.22, just = "left")))+
    annotation_custom(grobTree(textGrob(
        paste0("p = ", pcor$p), 
        x=0.05, y=0.1, just = "left")))+
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
      plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"))
  
}

plots <- plot_grid(plotlist = gene_plots, ncol = 1,
                   labels= c("B", "", "",""), label_size = 20, label_y = 0.96)

## Cast subplots into one and save
gs <- list(A, plots)
lay <- rbind(c(1,1,2),
             c(1,1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="FigS4_CTG_PSI.jpeg", height = 10, width = 10, dpi = 600 )








