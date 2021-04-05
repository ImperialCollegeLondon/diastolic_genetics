install.packages("plyr")
install.packages("data.table")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("ggExtra")
install.packages("ggpmisc")
install.packages("ggsignif")
install.packages("forcats")


library(plyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggpmisc)
library(ggsignif)
library(forcats)

# Load UK Biobank data

## Input:
## - plotdata: Data of all phenotypes/variables for the scatter and violin plot in Fig.2.

# Extract PDSRll (s-1), PDSRrr (s-1), Age and Sex

unidata<-as.data.frame(cbind(plotdata$`PDSRll (s-1)`, plotdata$`PDSRrr (s-1)`, plotdata$Age, plotdata$Sex))
colnames(unidata)<-c("PDSRll (s-1)", "PDSRrr (s-1)", "Age", "Sex")
unidata<-na.omit(unidata) 
unidata$Sex <- mapvalues(unidata$Sex, from=c("0", "1"), to=c("Female (n=20,097)", "Male (n=18,604)")) 


theme_update(
  
  panel.border = element_blank(),
  axis.line.y = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.text.x  = element_text(size=18),
  axis.text.y  = element_text(size=18),
  axis.title.x  = element_text(size=18, vjust=0.3, face="bold"),
  axis.title.y  = element_text(size=18, face = "bold", vjust=0.9, angle = 90),
  axis.line = element_line(size = 0.8, linetype = "solid"),
  axis.ticks = element_line(size = 1), legend.position="none"
  
) 


# Scatterplot of PDSRll (s-1) with density contours and marginal density plots

p1 <- ggplot(unidata, aes(x=Age, y=`PDSRll (s-1)`)) +
  geom_point(alpha=0.1, colour="steelblue2", pch=19) +
  #scale_shape(solid = FALSE)+
  stat_density2d(geom="density2d", aes(alpha=..level..), colour="magenta3", size=1.5, contour=TRUE) +
  geom_smooth(method="lm", colour="black", se=TRUE, size=1.5, level=0.99)

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  fill = '#9FA4EB'
)

# Scatterplot of PDSRrr (s-1) with density contours and marginal density plots

p2 <- ggplot(unidata, aes(x=Age, y=`PDSRrr (s-1)`)) +
  geom_point(alpha=0.1, colour="steelblue2", pch=19) +
  stat_density2d(geom="density2d", aes(alpha=..level..), colour="magenta3", size=1.5, contour=TRUE) +
  geom_smooth(method="lm", colour="black", se=TRUE, size=1.5, level=0.99)

ggExtra::ggMarginal(
  p = p2,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  fill = '#9FA4EB'
)

# AI - flatten transparency. Print size 10.23 x 6.24


# Violin plot of PDSRll (s-1) by Sex

PDSRll <- ggplot(unidata, aes(Sex, `PDSRll (s-1)`, group=Sex, fill=Sex)) + 
  geom_violin(trim = T, alpha=0.3)+
  geom_boxplot(fill = c("brown4", "#39568CFF"),color = c("#C84B6D","steelblue3"),width=0.1, outlier.size=2, alpha=0.3)+
  scale_fill_manual( values = c('white','white'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x  = element_text(size=18),
        axis.text.y  = element_text(size=18),
        axis.line = element_line(size = 0.8, linetype = "solid"),
        axis.title.x  = element_text(size=18, vjust=0.3, face="bold"),
        axis.title.y  = element_text(size=18, face = "bold", vjust=0.9, angle = 90)
  )+
  labs(x="Sex", y = expression(bold(paste("PDSRll (s"^"-1",")")))) +
  stat_compare_means(aes(label = ..p.signif..), label.x = 1.5, size = 8
  )

PDSRll

# Violin plot of PDSRrr (s-1) by Sex

PDSRrr <- ggplot(unidata, aes(Sex, `PDSRrr (s-1)`, group=Sex, fill=Sex)) + 
  geom_violin(trim = T, alpha=0.3)+
  geom_boxplot(fill = c("brown4", "#39568CFF"),color = c("#C84B6D","steelblue3"),width=0.1, outlier.size=2, alpha=0.3)+
  scale_fill_manual( values = c('white','white'))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.text.x  = element_text(size=18),
        axis.text.y  = element_text(size=18),
        axis.line = element_line(size = 0.8, linetype = "solid"),
        axis.title.x  = element_text(size=18, vjust=0.3, face="bold"),
        axis.title.y  = element_text(size=18, face = "bold", vjust=0.9, angle = 90)
  )+
  labs(x="Sex", y = expression(bold(paste("PDSRrr (s"^"-1",")")))) +
  stat_compare_means(aes(label = ..p.signif..), label.x = 1.5, size = 8
  )

PDSRrr


