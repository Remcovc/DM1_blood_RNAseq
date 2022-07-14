
## Script generated on 19/11/21 to study the impact of the correlation between
## CTG and Compound Response estimates on figure 5
## Last update on 6/7/22

###############
## Libraries ##
###############

library(ggplot2)   #V 3.3.3
library(readxl)    #V 1.3.1
library(tidyverse) #V 1.3.1
library(cowplot)   #V 1.1.1
library(stats)     #V 4.0.4
library(gplots)    #V 3.1.1
library(grid)      #V 4.0.4
library("gridExtra") #V 2.3
library(psych)     #V 2.1.9

########################
## Load relevant data ##
########################

load("CTG_coef.RDATA")
CTG_fit_new <- data.frame(lmer_fit_values[[2]])
CTG_hits <- CTG_fit_new$ENSG[CTG_fit_new$FDR < 0.05]

load("scaled_response_coef.RDATA")
new_fit <- data.frame(lmer_fit_values[[2]])
model_hits <- new_fit$ENSG[new_fit$FDR < 0.05]

hit_list <- list("Compound response associated" = model_hits, "CTG repeat associated" = CTG_hits)

#Load CR fit on CTG res
load("CTG_res_CR_fit_values.RDATA")
CR_res <- data.frame(lmer_fit_values)


#################################################
## Plot old regression coefficient against new ##
#################################################
df2 <- data.frame("CR_effect" = new_fit$Estimate, 
"CR_res_effect" = CR_res$scaled_response.Estimate)
pcor <- corr.test(df2$CR_effect, df2$CR_res_effect, method="pearson")
pcor$p <- round(pcor$p, 4)
if (pcor$p < 0.0001){
  pcor$p <- "< 0.0001"
}

A <- ggplot(df2, aes(x=df2[,1], y=df2[,2])) +
  xlab("Compound response effect on gene expression") + 
  ylab("Compound response effect on residuals") +
  labs(tag ="A") +
  xlim(-2.5,2.5) +
  ylim(-2.5,2.5) +
  geom_point() +
  geom_abline(slope=1, intercept=0, col="red") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor$r, 2)), 
    x=0.05, y=0.95, just = "left")))+
  annotation_custom(grobTree(textGrob(
    paste0("p = ", pcor$p), 
    x=0.05, y=0.9, just = "left")))+
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
    plot.margin = margin(c(0.05,0.2,0.05,0.05), unit="cm"),
    text = element_text(face="bold"),
    plot.tag = element_text(color ="black", size= 20, face="bold")
  )


###########################
## Recreate scatter plot ##
###########################

#select all model data for the CTG-repeat and response for all genes related to either CTG repeat or response
Repeat <- CTG_fit_new
Response <- new_fit

# create dataframe with Repeat and Response, where repeat is multiplied by 100 for better interpretability (=effect per 100 CTGs). 
# Also, indicate which genes are significant in both sets

df <- data.frame(
  Repeat = Repeat[,"Estimate"]*100,
  Response = CR_res$scaled_response.Estimate, #in contrast to figure 5, estimates of residual fit are used
  both = ifelse(Response[,"FDR"] < 0.05 & Repeat[,"FDR"] < 0.05, "yes","no"))
pcor <- corr.test(df$Repeat[df$both == "yes"], df$Response[df$both=="yes"], method="pearson")
pcor$p <- round(pcor$p, 4)
if (pcor$p < 0.0001){
  pcor$p <- "< 0.0001"
}

# reorder so that genes significant in both are plotted on front of the others
df <- df[order(df$both),]

B <- ggplot(df,  aes_string(x="Repeat", y="Response")) + 
  geom_point(col = ifelse(df$both == "yes", "purple4", "grey70"), cex = ifelse(df$both == "yes", 1, 0.5))+
  xlab("CTG repeat effect size") + 
  ylab("Compound response effect size on CTG fit residuals") +
  geom_hline(yintercept=0, col = "red")+
  geom_vline(xintercept=0, col = "red")+
  scale_y_continuous(breaks = seq(-2.5, 2.5, 0.5), 
                     limits = c(-2.5, 2.5)) +
  scale_x_continuous(breaks = seq(-0.25, 0.25, 0.1), 
                     limits = c(-0.25,0.25)) +
  annotation_custom(grobTree(textGrob(
    paste0("Rho = ", round(pcor$r, 2)), 
    x=0.05, y=0.95, just = "left",
    gp=gpar(fontsize=14))))+
  annotation_custom(grobTree(textGrob(
    paste0("p = ", pcor$p), 
    x=0.05, y=0.9, just = "left",
    gp=gpar(fontsize=14))))+
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
    plot.margin = margin(c(0.05,0.05,0.05,0.05), unit="cm"),
    text = element_text(face="bold"),
  )


## Combine plots A and B
gs <- list(A, B)
lay <- rbind(c(1,2),
             c(1,2))
plot <- arrangeGrob(grobs = gs, layout_matrix=lay)
ggsave(plot, file ="FigS6_CR_residuals.jpeg", height = 6, width = 12, dpi = 600 )






