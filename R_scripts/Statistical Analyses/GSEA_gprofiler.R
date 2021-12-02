
## Script generated on 3/9/21, last change applied on 14/10/21
## Purpose: Gene set enrichment analysis for the different identified gene sets
## in the ReCognitION RNA-seq analyses: CBT-associated, CTG-associated, 
## ClinicalResponse-associated, Overlap_97_associated


###############
## Libraries ##
###############

library("gprofiler2")
library("readxl")

############################
## Load relevant datasets ##
############################

# CBT effect
CBT <- read_excel("CBT_effect.xlsx")
CBT <- data.frame(CBT)

# CTG effect
CTG <- read_excel("CTG_effects.xlsx")
CTG <- data.frame(CTG)
CTG <- CTG[,-c(3:5)]

# Compound Response effect
CR <- read_excel("Compound_Response.xlsx")
CR <- data.frame(CR)
CR <- CR[,-c(3:5)]


#################################
## Obtain ordered gene vectors ##
#################################

# Vector of background genes = all genes expressed in blood samples
background_vec <- CBT$hgnc_symbol

# Top 500 sig CBT genes ordered by absolute regression coefficient, decreasing
CBT_ordered <- CBT[order(CBT$CBT_p.val),]
CBT_subset <- CBT_ordered[1:500,]
CBT_subset <- CBT_subset[order(abs(CBT_subset$CBT_effect), decreasing = T),]
CBT_vec <- CBT_subset$hgnc_symbol

# Top 500 sig CTG genes ordered by absolute regression coefficients, decreasing
CTG_ordered <- CTG[order(CTG$CTG_p.val),]
CTG_subset <- CTG_ordered[1:500,]
CTG_subset <- CTG_subset[order(abs(CTG_subset$CTG_Effect), decreasing = T),]
CTG_vec <- CTG_subset$hgnc_symbol

# Top 500 sig Compound Response genes ordered by abs gene reg coef, decr
CR_ordered <- CR[order(CR$Compound_Response_p.val),]
CR_subset <- CR_ordered[1:500,]
CR_subset <- CR_subset[order(abs(CR_subset$Compound_Response_Effect), decreasing = T),]
CR_vec <- CR_subset$hgnc_symbol

# Top overlapping genes between CTG and Compound Response (FDR <= 10%), no weighted ORA
CTG_CR_vec <- CTG$hgnc_symbol[CTG$CTG_FDR < 0.1 & CR$Compound_Response_FDR < 0.1] 



##################################
## Over representation analyses ##
##################################

## Parameters
# ordered gene vector as input [except for overlap CTG_CR analysis]
# only estimate over representation
# default multiple testing correction
# domain scope: custom_bg = expressed genes in blood samples
# sources: limited to Wikipathways "WP"

# CBT ORA
CBT_ORA <- gost(query = CBT_vec,
                organism ="hsapiens", ordered_query = TRUE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "custom", custom_bg = background_vec,
                numeric_ns ="", sources = c("WP"), as_short_link = FALSE)
CBT_ORA$result$query <- "CBT_ORA"

# CTG ORA
CTG_ORA <- gost(query = CTG_vec,
                organism ="hsapiens", ordered_query = TRUE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "custom", custom_bg = background_vec,
                numeric_ns ="", sources = c("WP"), as_short_link = FALSE)
CTG_ORA$result$query <- "CTG_ORA"

# Compound Response ORA
CR_ORA <- gost(query = CR_vec,
                organism ="hsapiens", ordered_query = TRUE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "custom", custom_bg = background_vec,
                numeric_ns ="", sources = c("WP"), as_short_link = FALSE)
CR_ORA$result$query <- "CR_ORA"

# Overlapping sig results of CTG and CR ORA - not-weighted
CTG_CR_ORA <- gost(query = CTG_CR_vec,
                organism ="hsapiens", ordered_query = FALSE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "custom", custom_bg = background_vec,
                numeric_ns ="", sources = c("WP"), as_short_link = FALSE)
CTG_CR_ORA$result$query <- "CTG_CR_ORA"

################################################
## Combine, select and store relevant results ##
################################################

ORA_combined <- rbind(CBT_ORA$result, CTG_ORA$result, CTG_CR_ORA$result)
ORA_combined <- ORA_combined[,c(1, 3, 4, 6, 9, 11, 16)]
ORA_combined$p_value <- round(ORA_combined$p_value, 3)
colnames(ORA_combined)[colnames(ORA_combined) == "p_value"] <- "Adjusted p-value"

write.csv(ORA_combined, 
          "ORA_results.csv",
          row.names = F)







