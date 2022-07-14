# Original script made by R van Cruchten to visualize the association of
# the identified CTG repeat effect with the DM1 affect as assessed in tissues
# of other cohorts

# Slight modifications by D van As for publication purposes
# Last changes applied on 6/7/2022


###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)    #V 3.3.5
library(readxl)     #V 1.3.1
library(tidyverse)  #V 1.3.1
library(cowplot)    #V 1.1.1
library(stats)      #V 4.0.4
library(gplots)     #V 3.1.1
library(grid)       #V 4.0.4
library(psych)      #V 2.1.9


######################
## Loading datasets ##
######################

#Load a table with hgcn symbols and corresponding Ensembl geneIDs
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)

# load CTG fit results
load("CTG_coef.RDATA")
CTG_fit <- data.frame(lmer_fit_values[[2]])

# model values from wilcoxon tests in Expression_analyses_other_papers.rmd
Sznajder_blood <- read.csv(file = "wilcox_DM1_GSE138691_values.csv", sep = ",")
# model values from wilcoxon tests in Expression_analyses_other_papers.rmd
Charlet_heart <- read.csv(file = "wilcox_DM1_003_values.csv", sep = ",")
# model values from wilcoxon tests in Expression_analyses_other_papers.rmd
Otero_brain <- read.csv(file = "wilcox_DM1_041_values.csv", sep = ",")
# Table S5 in Wang et al HMG 2019 (DMseq)
DMseq_MBNLinf <- read_excel("Kallisto_table_s5_ddy432.xlsx")
# model values from wilcoxon tests in Expression_analyses_other_papers.rmd 
DMseq_tibialis <- read.csv(file = "wilcox_DM1_007_values.csv", sep = ",")
#Table EV10 from signorelli et al 2020  
signorelli_DMDblood <- full_join(read_excel("Spitali_DMD.xlsx", sheet = "1_body_measurements"), read_excel("Spitali_DMD.xlsx", sheet = "2_physical_tests"), by = "genes") 


######################
## Data preparation ##
######################

# add ensembl gene IDs to DMseq MBNLinf and muscle strength expression table
DMseq_MBNLinf <- inner_join(hgnc_symbol, DMseq_MBNLinf, by = c("hgnc_symbol" = "Gene"))
colnames(DMseq_MBNLinf)[colnames(DMseq_MBNLinf) == "ensembl_gene_id"] <- "ENSG"

# add ensembl gene IDs to DMseq wilcoxon test table
DMseq_tibialis <- inner_join(hgnc_symbol, DMseq_tibialis, by = c("hgnc_symbol" = "hgnc_symbol"))
colnames(DMseq_tibialis)[colnames(DMseq_tibialis) == "ensembl_gene_id"] <- "ENSG"

# fix column names of merged Signorelli table
colnames(signorelli_DMDblood) <- gsub("y", "DMDtests", colnames(signorelli_DMDblood))
colnames(signorelli_DMDblood) <- gsub("x", "DMDbody", colnames(signorelli_DMDblood))

# add ensembl gene IDs
signorelli_DMDblood <- inner_join(hgnc_symbol, signorelli_DMDblood, by = c("ensembl_gene_id" = "genes"))
colnames(signorelli_DMDblood)[colnames(signorelli_DMDblood) == "ensembl_gene_id"] <- "ENSG"

#### change column name of data from this study to facilitate labeling later 
#### colnames(model_values_S2) <- gsub("CTGRepeatLength","CTG_repeat_length", colnames(model_values_S2))

#make list of all separate relevant data sets
studies <- list("Sznajder_blood" = Sznajder_blood, 
                "Charlet_heart" = Charlet_heart,
                "Otero_brain" = Otero_brain,
                "DMseq_tibialis" = DMseq_tibialis, 
                "DMseq_inferred_MBNL_activity" = data.frame(DMseq_MBNLinf[, grep("Correlation to MBNLinf|hgnc", colnames(DMseq_MBNLinf))]),
                "DMseq_muscle_strength" = data.frame(DMseq_MBNLinf[, grep("Correlation to Strength|hgnc", colnames(DMseq_MBNLinf))]),
                "DMD_physical_tests"       = data.frame(signorelli_DMDblood[, grep("DMDtests|hgnc", colnames(signorelli_DMDblood))]),
                "DMD_body_measurements"      = data.frame(signorelli_DMDblood[, grep("DMDbody|hgnc", colnames(signorelli_DMDblood))])
)


###################
## Scatter plots ##
###################

# Remove rows with duplicated gene IDs
for (study in names(studies)){
  studies[[study]] <- studies[[study]][!duplicated(studies[[study]][,"hgnc_symbol"]),]
}

# find genes that are measured in all studies
measured_genes <- intersect(intersect(intersect(intersect(intersect(
  CTG_fit$hgnc_symbol, 
  studies[["Sznajder_blood"]][,"hgnc_symbol"]),
  studies[["Charlet_heart"]][,"hgnc_symbol"]),
  studies[["Otero_brain"]][,"hgnc_symbol"]),
  studies[["DMseq_tibialis"]][,"hgnc_symbol"]),
  studies[["DMD_physical_tests"]][,"hgnc_symbol"])

# extract info for only genes measured in all studies
for (study in names(studies)){
  studies[[study]] <- studies[[study]][studies[[study]][,"hgnc_symbol"] %in% measured_genes,]
}

# add data generated in this study to the list of studies to compare
studies[["CTGRepeat"]] <- CTG_fit[CTG_fit$hgnc_symbol %in% measured_genes,]

# create plot for each of the external studies vs this study
study_plots <- list()
for (study in names(studies)[!names(studies) == "CTGRepeat"]){
  
  #create dataframe with all DM1 vs expression data for both studies. Multiply CTG effect by 100 to obtain effect per 100 CTGs
  df <- data.frame(cbind(studies[["CTGRepeat"]][,"Estimate"]*100, studies[[study]][, grep("meandiff|Estimate|Correlation|logFC", colnames(studies[[study]]))]))
  
  colnames(df) <- c("ReCognitION", study)
  pcor <- corr.test(df[,"ReCognitION"],df[,study], method="pearson")
  pcor$p <- round(pcor$p, 4)
  if (pcor$p < 0.0001){
    pcor$p <- "< 0.0001"
  }
  
  study_plots[[study]] <- ggplot(df,  aes_string(x="ReCognitION", y = study)) + 
    geom_point(cex = 0.5)+
    ggtitle(paste(gsub("_", " ", study))) +
    theme(legend.position = "none", aspect.ratio = 1/1) +
    xlab("Repeat length effect ReCognitION") + 
    ylab(ifelse(study == "DMseq_inferred_MBNL_activity", "Correlation to MBNLinf", 
                ifelse(study == "DMseq_muscle_strength", "Correlation to Strength", 
                       ifelse(grepl("DMD", study), "logFC", paste("Expr. DM1 - ctrl.(logCPM)"))))) +
    geom_hline(yintercept=0, col = "red")+
    geom_vline(xintercept=0, col = "red") +
    annotation_custom(grobTree(textGrob(
      paste0("Rho = ", round(pcor$r,3)), 
      x=0.05, y=0.95, just = "left")))+
    annotation_custom(grobTree(textGrob(
      paste0("p = ", pcor$p), 
      x=0.05, y=0.85, just = "left")))+
    scale_x_continuous(labels = scales::comma(seq(-0.4, 0.4, 0.1), accuracy = 0.1), 
                       breaks = seq(-0.4, 0.4, 0.1), 
                       lim = c(-0.4, 0.4),
                       sec.axis = dup_axis(breaks = derive(), labels = NULL, name = NULL))+
    scale_y_continuous(lim = c(-max(abs(df[, study])), max(abs(df[, study]))),
                       sec.axis = dup_axis(breaks = derive(), labels = NULL, name = NULL))+
    theme_classic() +
    theme(panel.grid.major.x = element_line(size = 0.25, color = "grey"))+
    theme(panel.grid.major.y = element_line(size = 0.25, color = "grey")
    )
}

study_plots <- plot_grid(plotlist = study_plots, ncol = 2, labels="AUTO", label_size = 20)

ggsave(study_plots, file = "FigS3_Otherstudies.jpeg", height = 15, width = 10, dpi = 600 )
























