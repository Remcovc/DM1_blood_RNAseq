## Script created on 16/07/2022
## Title: Outcome types and directions of improvement
## Author: Remco van Cruchten

###############
## Libraries ##
###############

library(tidyverse)

######################
## Loading datasets ##
######################

# this file can be found in the Github and contains '-1' or '1' for each outcome measure to correct for measures that decrease upon clinical improvement
load(file = "corrections_outcomes.RData")

########################################
## Create list with outcomes per type ##
########################################

# Generate a list with strings of similar outcome measures, which can be queried to order the measures by type.
outcome_types <- list()
outcome_types[["Compound_scores"]]      <- c("DM1ActivC","MDHI","INQOLQolScore", "ASBQ","ICQ", "IMQ", "CISactivity")
outcome_types[["Physical"]]             <- c("SMWT", "MeanENMO" ,"PreBORG")
outcome_types[["Fatigue"]]              <- c("FDSS","CISFatigue", "JFCS")
outcome_types[["Cognitive"]]            <- c("TMT","StroopInterference")
outcome_types[["Pain_depr,_soc,_att."]] <- c("McGillPain","BDIFs", "SSLDScore", "SSLNScore", "SSLIScore","AEScScore","SES28")

###############################################
## Calculate scaled improvements per patient ##
###############################################

# Create list of change in outcome measures by looping over all outcome measures for all patients
dOutcomes <- list()
for (outcome in c(unlist(outcome_types))){
  for (patient in c(unique(samples$PatientID))){
    dOutcomes[[outcome]][[patient]] <- as.numeric(samples[paste0(patient,"_V4"), outcome]) - as.numeric(samples[paste0(patient,"_V2"), outcome])
  }  
  # Combine into one dataframe per outcome measure
  dOutcomes[[outcome]] <- do.call(cbind, dOutcomes[[outcome]])
}

# Create dataframe of all outcome measures and patients
dOutcomes <- data.frame(do.call(rbind, dOutcomes))
rownames(dOutcomes) <- outcome_measures

# Correct outcomes for change in 'positive' direction by looping over all changes in outcome measures and multiply by 1 or -1
dOutcomes_corrected <- list()
for (outcome in rownames(dOutcomes)){
  dOutcomes_corrected[[outcome]] <-  as.numeric(corrections[outcome]*dOutcomes[outcome,])
}

# Turn into dataframe and order
dOutcomes_corrected <- do.call(rbind, dOutcomes_corrected)
colnames(dOutcomes_corrected) <- colnames(dOutcomes)
dOutcomes_corrected <-  dOutcomes_corrected[,order(colnames(dOutcomes_corrected))]
