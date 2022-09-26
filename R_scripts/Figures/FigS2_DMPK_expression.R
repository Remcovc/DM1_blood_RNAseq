
# Visualization of changes in DMPK expression before and after CBT
# Script generated on 22/09/22


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
load(file = "samples.RDATA")
# load CBT mixed effects model
load("CBT_coef.RDATA")
#this Voom object features counts with a cutoff of 50 based on the Visit (CBT) design. From Mixed_model_gene_expression_DVA.R
load(file = "v_visit.RDATA")
# df with all ENSG en hgnc symbols
hgnc_symbol <- read.table("ENSG_geneSymbol.txt", sep =",", header=TRUE)


#######################################
## Check significance of DMPK result ##
#######################################

lmer_fit_values$CBT[lmer_fit_values$CBT$hgnc_symbol=="DMPK",]

##########################################################################################
## Generate dataframe with relevant p-values, FDR-values and effect sizes for CBT model ##
##########################################################################################

CBT_df <- data.frame(lmer_fit_values["CBT"])
df <- data.frame(hgnc_symbol = CBT_df$CBT.hgnc_symbol,
                 ENSG = CBT_df$CBT.ENSG,
                 p.value = CBT_df$CBT.p.val,
                 FDR = CBT_df$CBT.FDR,
                 effect = CBT_df$CBT.Estimate)

############################################
## Calculate delta-gene expression values ##
############################################

counts <- v$E
V2_counts <- counts[,grepl("V2", colnames(counts))]
colnames(V2_counts) <- gsub(x = colnames(V2_counts), pattern = "_V2", replacement ="")
V2_counts <- V2_counts[,order(colnames(V2_counts))]

V4_counts <- counts[,grepl("V4", colnames(counts))]
colnames(V4_counts) <- gsub(x = colnames(V4_counts), pattern = "_V4", replacement ="")
V4_counts <- V4_counts[,order(colnames(V4_counts))]

table(colnames(V2_counts) == colnames(V4_counts))
delta_counts <- V4_counts - V2_counts

###############################################
## Check association of DMPK with CTG repeat ##
###############################################

# Obtain relevant CTG values
CTG <- as.numeric(samples$V2Mode[samples$Visit=="V2"])
names(CTG) <- samples$PatientID[samples$Visit=="V2"]
CTG <- CTG[colnames(delta_counts)]
table(colnames(delta_counts) == names(CTG))

# Calculate pearson correlations and store in dataframe
ENSG <- hgnc_symbol$ensembl_gene_id[hgnc_symbol$hgnc_symbol=="DMPK"]
pcor <- corr.test(delta_counts[rownames(delta_counts) == ENSG,], CTG, method="pearson")
df2 <- data.frame(cbind(CTG = CTG, DMPK = delta_counts[rownames(delta_counts) == ENSG,]))

plot1 <- ggplot(df2, aes_string(x="CTG", y="DMPK")) +
  ggtitle("DMPK") +
  xlab("CTG repeat length") +
  ylab("delta logCPM") +
  geom_point(size = 1) +
  geom_smooth(method ="lm", formula =  y ~ x, se=F, col ="black") +
  scale_x_continuous(breaks=seq(0,1000,250), 
                     limits = c(0,1000)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1), 
                     limits = c(min(df2$DMPK), max(df2$DMPK)+0.2)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor$r, 2)), 
    x=0.05, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  labs(tag ="B") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    aspect.ratio = 1,
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 18, face="bold"),
    axis.title.x = element_text(color = "black", size=18, face="bold"),
    plot.title = element_text(color = "black", size = 18, face="italic"),
    plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm",),
    plot.tag = element_text(color ="black", size= 20, face="bold")
  )


####################
# DMPK change plot #
####################

ENSG_ID <- df$ENSG[df$hgnc_symbol == "DMPK"]

df2 <- data.frame(samples[,"Visit"], v[["E"]][ENSG_ID,], samples[,"PatientID"])
names(df2) <- c("Visit","counts","PatientID")
df2 <- df2[order(df2$Visit),]

df2$time <- df2$Visit
df2$time[df2$time == "V2"] <- 0
df2$time[df2$time == "V4"] <- 1
df2$time <- as.numeric(df2$time)
pcor <- corr.test(df2$time, df2$counts, method="pearson")
  
plot2 <- ggplot(df2,  aes_string(x="Visit", y="counts",group="PatientID")) +    
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
    x=0.05, y=0.95, just = "left",
    gp=gpar(fontsize=16))))+
  labs(tag ="A") +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    aspect.ratio = 1,
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(color = "black", size = 16, face="bold"),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color ="black", size = 18, face="bold"),
    axis.title.y = element_text(color = "black", size = 18, face="bold"),
    plot.title = element_text(color = "black", size = 18, face="italic"),
    plot.margin = margin(c(0.05,0.8,0.05,0.05), unit="cm"),
    plot.tag = element_text(color ="black", size= 20, face="bold")
    )

############################
## Combine plots and save ##
############################

gs <- list(plot2, plot1)
lay <- rbind(c(1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, 
       file ="FigS2_DMPK.png", 
       height = 5, width = 10, 
       dpi = 1200,
       device ="png")


