
# Plot generated to visualize the patient characteristics age, sex and CTG
# of the 27 patients which were selected for RNA-sequencing
# Script generated on 6/8/21, last modifications applied on 26/09/22

###############
## Libraries ##
###############

library(ggplot2)     #V 3.3.3
library("gridExtra") #V 2.3

#######################
## Load patient data ##
#######################

load(file = "samples.RDATA")

df <- samples[samples$Visit=="V2",]
df$V2Mode <- as.numeric(df$V2Mode)

# Age distribution per sex

A <- ggplot(df, aes_string(x= "Sex", y="AgeBaseline",fill="Sex")) +
  geom_boxplot(lwd=1) +
  geom_jitter(shape=16, position = position_jitter(0)) +
  xlab("") + 
  ylab("Age at baseline") +
  labs(tag ="A") +
  theme(
    aspect.ratio = 1.5,
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    plot.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 16, face="bold"),
    axis.text.x = element_text(color = "black", size = 14, face="bold"),
    axis.text.y = element_text(color = "black", size = 14),
    plot.tag = element_text(color ="black", size= 20, face="bold")
  )

# CTG distribution per sex

B <- ggplot(df, aes_string(x="Sex", y="V2Mode",fill="Sex")) +
  geom_boxplot(lwd=1) +
  geom_jitter(shape=16, position = position_jitter(0)) +
  xlab("") + 
  ylab("CTG Repeat length") +
  labs(tag ="B") +
  theme(
    aspect.ratio = 1.5,
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    plot.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 16, face="bold"),
    axis.text.x = element_text(color = "black", size = 14, face="bold"),
    axis.text.y = element_text(color = "black", size = 14),
    plot.tag = element_text(color ="black", size= 20, face="bold")
  )

# Delta DM1-Activ-c

df2 <- samples[samples$Visit == "V4",]
df <- df[order(df$PatientID),]
df2 <- df2[order(df2$PatientID),]
df$ddm1a <- df2$DM1ActivC - df$DM1ActivC

C <- ggplot(df, aes_string(x= "Sex", y="ddm1a",fill="Sex")) +
  geom_boxplot(lwd=1) +
  geom_jitter(shape=16, position = position_jitter(0)) +
  xlab("") + 
  ylab("Delta DM1-Activ-c") +
  labs(tag ="C") +
  theme(
    aspect.ratio = 1.5,
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill="white"),
    panel.grid.major.x = element_line(size = 0.25, color = "grey"),
    panel.grid.major.y = element_line(size = 0.25, color = "grey"),
    plot.margin = margin(c(0.1,0.1,0.1,0.1), unit="cm"),
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 16, face="bold"),
    axis.text.x = element_text(color = "black", size = 14, face="bold"),
    axis.text.y = element_text(color = "black", size = 14),
    plot.tag = element_text(color ="black", size= 20, face="bold")
  )

# Cast subplots into one
plot <- arrangeGrob(A, B, C, ncol = 3, nrow=1)

ggsave("FigS1_pat_char.png", 
       dpi = 1200, height = 7, 
       width = 15, 
       plot,
       device="png")



