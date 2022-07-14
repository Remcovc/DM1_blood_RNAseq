
## Script generated on 6/7/22 by D van As, last update on 14/7/22
## Purpose: Explore the correlation of gene expression changes with changes
## in certain outcome measures

###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library("psych")        #V 2.1.9
library("ggplot2")      #V 3.3.5
library("psych")        #V 2.1.9


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
# results of the CTG fits
load("CTG_coef.RDATA")
CTG_fit_new <- data.frame(lmer_fit_values[[2]])
CTG_hits <- CTG_fit_new$ENSG[CTG_fit_new$FDR < 0.05]
# results of the scaled response fits
load("scaled_response_coef.RDATA")
new_fit <- data.frame(lmer_fit_values[[2]])
model_hits <- new_fit$ENSG[new_fit$FDR < 0.05]

#################################
## Calculate delta-DM1-Activ-c ##
#################################

## Sort patient data by IDs
samples <- samples[order(samples$PatientID),]

## Split baseline and 10 month assessments
dfV2 <- samples[samples$Visit == "V2",]
dfV4 <- samples[samples$Visit == "V4",]

table(colnames(dfV2) == colnames(dfV4))
table(dfV2$PatientID == dfV4$PatientID)

## Calculate deltas and store in matrix
dDM1ActivC <- as.numeric(dfV4$DM1ActivC) - as.numeric(dfV2$DM1ActivC)
names(dDM1ActivC) <- gsub("_V4", "", rownames(dfV4))
CRS <- as.numeric(dfV2$scaled_mean_outcome)
names(CRS) <- gsub("_V2", "", rownames(dfV2))
CTG <- as.numeric(dfV2$V2Mode)
names(CTG) <- gsub("_V2", "", rownames(dfV2))


############################################################################
## Calculate changes in gene expression and obtain DMPK expression values ##
############################################################################

counts <- v$E
V2_counts <- counts[,grepl("V2", colnames(counts))]
colnames(V2_counts) <- gsub(x = colnames(V2_counts), pattern = "_V2", replacement ="")
V2_counts <- V2_counts[,order(colnames(V2_counts))]

V4_counts <- counts[,grepl("V4", colnames(counts))]
colnames(V4_counts) <- gsub(x = colnames(V4_counts), pattern = "_V4", replacement ="")
V4_counts <- V4_counts[,order(colnames(V4_counts))]

table(colnames(V2_counts) == colnames(V4_counts))
table(rownames(V2_counts) == rownames(V4_counts))
delta_counts <- V4_counts - V2_counts


#########################################################
## Calculate pearson rho and p-value for specific hits ## 
#########################################################

table(colnames(delta_counts) == names(dDM1ActivC))
table(colnames(delta_counts) == names(CRS))
table(colnames(V2_counts) == names(CTG))

#DNAJB12
ENSG <- CTG_fit_new$ENSG[CTG_fit_new$hgnc_symbol=="DNAJB12"]
pcor <- corr.test(delta_counts[rownames(delta_counts) == ENSG], CRS, method="pearson")
pcor <- corr.test(V2_counts[rownames(V2_counts) == ENSG], CTG, method="pearson")

#HDAC5
ENSG <- CTG_fit_new$ENSG[CTG_fit_new$hgnc_symbol=="HDAC5"]
pcor <- corr.test(delta_counts[rownames(delta_counts) == ENSG], CRS, method="pearson")
pcor <- corr.test(V2_counts[rownames(V2_counts) == ENSG], CTG, method="pearson")

#TRIM8
ENSG <- CTG_fit_new$ENSG[CTG_fit_new$hgnc_symbol=="TRIM8"]
pcor <- corr.test(delta_counts[rownames(delta_counts) == ENSG], CRS, method="pearson")
pcor <- corr.test(V2_counts[rownames(V2_counts) == ENSG], CTG, method="pearson")

#ZNF22
ENSG <- CTG_fit_new$ENSG[CTG_fit_new$hgnc_symbol=="ZNF22"]
pcor <- corr.test(delta_counts[rownames(delta_counts) == ENSG], CRS, method="pearson")
pcor <- corr.test(V2_counts[rownames(V2_counts) == ENSG], CTG, method="pearson")


##############
## 97 genes ##
##############
## Here we calculate the pearson correlation for changes in DM1-Activ-C and 
## Compound Response with changes in the 97 candidate biomarkers
table(colnames(delta_counts) == names(dDM1ActivC))
table(colnames(delta_counts) == names(CRS))

DM1_rho <- c()
DM1_p <- c()
CRS_rho <- c()
CRS_p <- c()

gene_hits <- intersect(model_hits, CTG_hits)

for (ENSG in gene_hits){
  pcor_DM1 <- corr.test(delta_counts[rownames(delta_counts) == ENSG], dDM1ActivC, method="pearson")
  DM1_rho <- append(DM1_rho, pcor_DM1$r)
  DM1_p <- append(DM1_p, pcor_DM1$p)

  pcor_CRS <- corr.test(delta_counts[rownames(delta_counts) == ENSG], CRS, method="pearson")
  CRS_rho <- append(CRS_rho, pcor_CRS$r)
  CRS_p <- append(CRS_p, pcor_CRS$p)
}

cor_df <- as.data.frame(cbind(DM1_rho, DM1_p, CRS_rho, CRS_p))
colnames(cor_df) <- c("delta-DM1-Activ-c-rho", "delta-DM1-Activ-c-p",
                      "Compound-Response-rho", "Compound-Response-p")
rownames(cor_df) <- gene_hits
nrow(cor_df[cor_df$CRS_p < 0.05,])
nrow(cor_df[cor_df$DM1_p < 0.05,])

hist(cor_df$CRS_rho, n=30)
hist(cor_df$DM1_rho, n=30)

## Store results of the 97 genes
write.table(cor_df, file = paste("Pcor_97.xlsx", sep=""))



