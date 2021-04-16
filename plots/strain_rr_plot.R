install.packages("ggplot2")
install.packages("matrixStats")
install.packages("robustbase")
library(ggplot2)
library(tidyr)
library(matrixStats)
library(robustbase)

# Calculates mean/median and SD/IQR for the radial strains. Then plots these with the ranges as a ribbon. 

## Input:
## - radial_strainrate: Radial strain data in 50 phases for the scatter plot with density contours in Fig.1.
# radial_strainrate <- read.csv("Radial_strain_39k_50phases.csv", header=T)


me_SR<-colMeans(as.matrix(radial_strainrate), na.rm=T)
mean_SR<-as.data.frame(cbind(c(1:50),me_SR))
colnames(mean_SR)<-c("Phase","mean")


sd_SR<-colSds(as.matrix(radial_strainrate), na.rm=T)

median_SR<-colMedians(as.matrix(radial_strainrate), na.rm=T)


iq<-as.data.frame(lapply(radial_strainrate,quantile, probs=c(0.25,0.75), na.rm=T))

iqr_SR<-t(iq)

mean_SR$sd<-sd_SR;mean_SR$median<-median_SR;mean_SR$Q1<-iqr_SR[,1];mean_SR$Q3<-iqr_SR[,2]

theme_update(
  panel.border = element_blank(),
  axis.line.y = element_line(colour="black", size=.7),
  axis.line.x = element_blank(),
  axis.text = element_text(colour="black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.text.y  = element_text(size=14),
  axis.title.x  = element_text(size=20, vjust=0.9, face="plain"),
  axis.title.y  = element_text(size=20, face = "plain", vjust=0.9, angle = 90),
  axis.line = element_line(size = 1.2, linetype = "solid"),
  axis.ticks.y = element_line(colour="black"),
  axis.ticks.x = element_blank()
) 

# Scatterplot with density contours and marginal density plots
ggplot(mean_SR, aes(x=Phase, y=median)) + 
  geom_ribbon(aes(ymin=Q1, ymax=Q3), fill="steelblue2", alpha=0.3) +
  ylab("Radial strain")+
  geom_line(colour="black") 


#print size 10.53 x 4.19

