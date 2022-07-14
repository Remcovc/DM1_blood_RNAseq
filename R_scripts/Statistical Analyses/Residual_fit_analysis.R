
## Script generated on 11/19/21 to study the impact of the correlation between
## CTG and Compound Response estimates on figure 5

###############
## Libraries ##
###############

knitr::opts_chunk$set(echo = TRUE)
library(limma)     #V 3.46.0
library(tidyverse) #V 1.3.1
library(edgeR)     #V 3.32.1
library(lmerTest)  #V 3.1.3


######################
## Loading datasets ##
######################

#Load a table with hgcn symbols and corresponding Ensembl geneIDs, and remove duplicate entries
hgnc_symbol <- read.table("ENSG_geneSymbol2.txt", sep =",", header=TRUE)
load(file = "outcome_types.RData")
load(file = "samples.RDATA")


######################################################
## Fitting Compound Response on CTG model residuals ##
######################################################

## load CTG models
load("CTG_fits.RDATA")

## Obtain residuals from CTG models
res_list <- list()
for (gene in names(lmer_fit)){
  res_list[[gene]] <- residuals(lmer_fit[[gene]])
}

## Fit Compound Response variable

# design matrix
design <- model.matrix(~samples[, "Visit"] * as.numeric(samples[, "scaled_mean_outcome"]))
design <- design[,-3]
colnames(design) <- c("(Intercept)", "CBT", "scaled_response")

# simple model
lmer_fit <- list()
for (gene in names(res_list)){
  lmer_fit[[gene]] <- lm(res_list[[gene]] ~ design[,"scaled_response"])
}

# extract regression coefficients and p-values
lmer_fit_coefficients <- lapply(names(lmer_fit), function(x) summary(lmer_fit[[x]])[["coefficients"]])
names(lmer_fit_coefficients) <- names(lmer_fit)

lmer_fit_values <- list()
lmer_fit_values[["scaled_response"]] <- do.call(rbind, lapply(names(lmer_fit_coefficients), function(x){lmer_fit_coefficients[[x]][2,]}))

# store relevant results per predictor in list, add gene names and FDR correction
for (fit in names(lmer_fit_values)){
  lmer_fit_values[[fit]] <- cbind(names(lmer_fit_coefficients),
                                  hgnc_symbol$hgnc_symbol[hgnc_symbol$ensembl_gene_id %in% names(lmer_fit_coefficients)],
                                  lmer_fit_values[[fit]], 
                                  p.adjust(lmer_fit_values[[fit]][, "Pr(>|t|)"], method = "fdr"))
  colnames(lmer_fit_values[[fit]]) <-c("ENSG", "hgnc_symbol", "Estimate", "Std.Error", "t_value", "p.val", "FDR")
}

# convert certain rows to numeric
lmer_fit_values[["scaled_response"]] <- data.frame(lmer_fit_values[["scaled_response"]])
lmer_fit_values[["scaled_response"]][,3:7] <- as.numeric(unlist(lmer_fit_values[["scaled_response"]][,3:7]))

# temporarily store and save results
save(lmer_fit_values, file = "CTG_res_CR_fit_values.RDATA")
save(lmer_fit, file = "CTG_res_CR_fit.RDATA")



