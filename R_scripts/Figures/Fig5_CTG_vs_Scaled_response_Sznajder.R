# Original script made by R van Cruchten to visualize the association of
# CTG repeat effect with Compound Response effect sizes
# Slight modifications by D van As for publication purposes and updated results
# Last changes applied on 18/11/21: added Sznajder comparison


###############
## Libraries ##
###############

library(VennDiagram) #V 1.6.20
library(corrplot)    #V 0.89
library(ggplot2)     #V 3.3.3
library("gridExtra") #V 2.3

########################
## Load relevant data ##
########################

load("CTG_coef.RDATA")
CTG_fit_new <- data.frame(lmer_fit_values[[2]])
CTG_hits <- CTG_fit_new$ENSG[CTG_fit_new$FDR < 0.05]

load("scaled_response_coef.RDATA")
new_fit <- data.frame(lmer_fit_values[[2]])
model_hits <- new_fit$ENSG[new_fit$FDR < 0.05]

# model values from wilcoxon tests in GSE138691_gene_expression.rmd
Sznajder_blood <- read.csv(file = "wilcox_DM1_GSE138691_values.csv", 
                           sep = ",")

hit_list <- list("Compound response associated" = model_hits, "CTG repeat associated" = CTG_hits)


########################################
## Venn diagram to illustrate overlap ##
########################################

VennDiagram::venn.diagram(hit_list, 
                          cex = 0.75,
                          category.names = names(hit_list), 
                          cat.dist = 0.03,
                          cat.pos = c(0,0),
                          cat.cex = 0.75,
                          rotation.degree = 180,
                          fill = c("red","blue"), 
                          lwd = 1,
                          resolution =  600,
                          imagetype = "tiff",
                          filename ="Fig5A_venn.jpeg")

###################
## Scatter_plots ##
###################

#select all model data for the CTG-repeat and response for all genes related to either CTG repeat or response
Repeat <- CTG_fit_new
Response <- new_fit

# create dataframe with Repeat and Response, where repeat is multiplied by 100 for better interpretability (=effect per 100 CTGs). 
# Also, indicate which genes are significant in both sets

df <- data.frame(
  Repeat = Repeat[,"Estimate"]*100,
  Response = Response[,"Estimate"],
  both = ifelse(Response[,"FDR"] < 0.05 & Repeat[,"FDR"] < 0.05, "yes","no"))

# reorder so that genes significant in both are plotted on front of the others
df <- df[order(df$both),]

B <- ggplot(df,  aes_string(x="Repeat", y="Response")) + 
  geom_point(col = ifelse(df$both == "yes", "purple4", "grey70"), cex = ifelse(df$both == "yes", 1, 0.5))+
  ggtitle("")+
  xlab("CTG repeat effect size") + 
  ylab("Compound response effect size") +
  geom_hline(yintercept=0, col = "red")+
  geom_vline(xintercept=0, col = "red")+
  scale_y_continuous(breaks = seq(-2.5, 2.5, 0.5), 
                     limits = c(-2.5, 2.5)) +
  scale_x_continuous(breaks = seq(-0.25, 0.25, 0.05), 
                     limits = c(-0.25,0.25)) +
  labs(tag ="B") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    aspect.ratio = 1,
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 16),
    plot.tag = element_text(color ="black", size= 20, face="bold"),
    plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
    text = element_text(face="bold"),
  )



######################
## Data preparation ##
######################

#add ensembl IDs as separate column
Sznajder_blood[,"ENSG"] <- rownames(Sznajder_blood)

#remove genes not measurered in Sznajders experiment
CTG_fit_new <- CTG_fit_new[CTG_fit_new$ENSG %in% Sznajder_blood$ENSG,]
new_fit <- new_fit[new_fit$ENSG %in% Sznajder_blood$ENSG,]
Sznajder_blood <- Sznajder_blood[Sznajder_blood$ENSG %in% new_fit$ENSG,]

##################
## Scatter_plot ##
##################

# create dataframe with DM1eff and response 
# Also, indicate which genes belong to the selection of CTG & Clinical response association

df <- data.frame(
  DM1eff = Sznajder_blood[,"meandiff"],
  Response = new_fit[,"Estimate"],
  both = ifelse(CTG_fit_new[,"FDR"] < 0.05 & new_fit[,"FDR"] < 0.05, "yes","no"))

# reorder so that genes significant in both are plotted on front of the others
df <- df[order(df$both),]

C <- ggplot(df,  aes_string(x="DM1eff", y="Response")) + 
  geom_point(col = ifelse(df$both == "yes", "purple4", "grey70"), cex = ifelse(df$both == "yes", 1, 0.5))+
  ggtitle("")+
  xlab("DM1 effect size") + 
  ylab("Compound response effect size") +
  geom_hline(yintercept=0, col = "red")+
  geom_vline(xintercept=0, col = "red")+
  scale_y_continuous(breaks = seq(-2.5, 2.5, 0.5), 
                     limits = c(-2.5, 2.5)) +
  scale_x_continuous(breaks = seq(-1.5, 2.5, 0.5), 
                     limits = c(-1.5,1.5)) +
  labs(tag ="C") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    aspect.ratio = 1,
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 16),
    plot.margin = margin(c(0.05,0.05,0.05,0.05), unit="cm"),
    plot.tag = element_text(color ="black", size= 20, face="bold"),
    text = element_text(face="bold"),
  )

## Combine plots B and C, A is added separately
gs <- list(B, C)
lay <- rbind(c(1,2),
             c(1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="Fig5_B_C.jpeg", height = 5, width = 15, dpi = 600 )





