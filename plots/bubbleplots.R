
install.packages("dplyr")
install.packages("ggplot2")
install.packages("plotly")
install.packages("forcats")
install.packages("tidyverse")
install.packages("hrbrthemes")

library(dplyr)
library(ggplot2)
library(plotly)
library(forcats)
library(tidyverse)
library(hrbrthemes)

# load beta coefficient and corrected pvalues from multiple linear regression analysis
# from https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/phenotype_analysis/data
multiple_analysis <- read.table("multiple_analysis.txt", header = TRUE)
# first column contains the groups of imaging phenotypes, second column contains the non-imaging phenotypes,
# third column contains the -log10(P), forth column contains the names of the Imaging phenotypes, 
# fifth sixth column contains the beta coefficients and the absolute beta coefficients between the imaging 
# and the non-imaging phenotypes.

# Plot for Beta coefficient

p<-multiple_analysis %>%
  mutate(Imaging_group=factor(Imaging_group,levels=c("Strains", "Strain rates","AAo - DAo","LA",
                                   "LV", "RA","RV","Not significant")))%>%
  mutate(Nonimaging=factor(Nonimaging,levels=unique(Nonimaging)))%>%
  mutate(Imaging=factor(Imaging,levels=unique(Imaging)))%>%
  # prepare text for tooltip
  mutate(text = paste("Imaging phenotypes: ", Imaging, "\nNon-Imaging phenotypes: ", Nonimaging, 
                      "\nBeta coefficinet: ", round(Beta,4),
                      "\n-log(p): ", round(LogP,2), sep=""))%>%
  
  ggplot( aes(x=Nonimaging, y=Beta,color=Imaging_group, size=LogP, text=text)) +
  geom_point(alpha=0.8) +
  scale_color_manual(values = c("Strains"="red3","Strain rates"="orangered",
                                "AAo - DAo"="maroon4", "RV"="springgreen4","RA"="yellowgreen", "LV"="dodgerblue4",
                                "LA"="turquoise2", "Not significant" = "lightgrey"))+
  geom_hline(yintercept =0, color="black", size=0.2)+
  scale_size(range = c(1, 10)) +
  theme_ipsum() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour="black"),
        axis.text.x  = element_text(angle=45,size=14),
        axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.line.x = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid"))+
  xlab("")+ylab("")

ggplotly(p, tooltip="text")%>%layout(yaxis=list(title=TeX("\\text{Beta coefficient}")))%>%
  config(mathjax = "cdn")

# Plot for -log10(P)

mt_Bh<-median(-log10(t_BH)) # get the median value form t_BH calculated from phenotype_statistical_analysis.R
ytick<-c(-6,mt_Bh,50,100,150,200,250,300)

p<-multiple_analysis %>%
  mutate(Imaging_group=factor(Imaging_group,levels=c("Strains", "Strain rates","AAo - DAo","LA",
                                   "LV", "RA","RV","Not significant")))%>%
  mutate(Nonimaging=factor(Nonimaging,levels=unique(Nonimaging)))%>%
  mutate(Imaging=factor(Imaging,levels=unique(Imaging)))%>%
  # prepare text for tooltip
  mutate(text = paste("Imaging phenotypes: ", Imaging, "\nNon-Imaging phenotypes: ", Nonimaging, 
                      "\nBeta coefficinet: ", round(Beta,4),
                      "\n-log(p): ", round(LogP,2), sep=""))%>%
  
  ggplot( aes(x=Nonimaging, y=LogP,color=Imaging_group, size=AbsBeta, text=text)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = mt_Bh, linetype="dashed", color="black", size=0.3)+
  scale_y_continuous(breaks = ytick,labels=c("0","BH","50","100","150","200","250","300"))+
  scale_color_manual(values = c("Strains"="red3","Strain rates"="orangered",
                                "AAo - DAo"="maroon4", "RV"="springgreen4","RA"="yellowgreen", "LV"="dodgerblue4",
                                "LA"="turquoise2", "Not significant" = "lightgrey"))+
  scale_size(range = c(1, 10)) +
  theme_ipsum() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour="black"),
        axis.text.x  = element_text(angle=45,size=14),
        axis.text.y  = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_line(size = 1),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.line = element_line(size = 0.5, linetype = "solid"))+
  xlab("")+ylab("")

ggplotly(p, tooltip="text")%>%layout(yaxis=list(title=TeX("-log_{10}(P)")))%>%
  config(mathjax = "cdn")

# END
