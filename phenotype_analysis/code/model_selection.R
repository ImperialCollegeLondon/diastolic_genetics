
install.packages("data.table")
install.packages("BeSS")
install.packages("devtools")
library(devtools)
devtools::install_github("SachaEpskamp/bootnet")
install.packages("stabs")
install.packages("glmnet")
install.packages("car")
install.packages("mctest")

library(data.table)
library(BeSS)
library("bootnet")
library(stabs)
library(glmnet)
library(car)
library(mctest)

### Multivariate LASSO regresion analysis with model  selection for selecting the non-imaging phenotypes 

data_pheno <- read.table("phenotypes_all.txt", header = TRUE)# load all 110 non-imaging and 31 imaging phenotype data

# Method 1 - Best subset selection approach with "beSS"
# using GPDAS algorithm to select the optimal subset selection k and the best model is determined by Extended Bayesian Information Criterion (EBIC)

pheno<-as.matrix(na.omit(data_pheno[,-c(1:2)])) # exclude the IDs

# PDSRll
fit.seqll <- bess(pheno[,-111], pheno[,111], method="sequential", epsilon = 0)
# PDSRrr
fit.seqrr <- bess(pheno[,-112], pheno[,112], method="sequential", epsilon = 0)
# LAVmaxi
fit.seqlav <- bess(pheno[,-127], pheno[,127], method="sequential", epsilon = 0)

K.opt.ebic.ll <- which.min(fit.seqll$EBIC)
K.opt.ebic.rr <- which.min(fit.seqrr$EBIC)
K.opt.ebic.lav <- which.min(fit.seqlav$EBIC)

# PDSRll
fit.one.ll <- bess.one(pheno[,-111], pheno[,111], s = K.opt.ebic.ll, family = "gaussian")
bm.one.ll <- fit.one.ll$bestmodel
#PDSRrr
fit.one.rr <- bess.one(pheno[,-112], pheno[,112], s = K.opt.ebic.rr, family = "gaussian")
bm.one.rr <- fit.one.rr$bestmodel
#LAVmaxi 
fit.one.lav <- bess.one(pheno[,-127], pheno[,127], s = K.opt.ebic.lav, family = "gaussian")
bm.one.lav <- fit.one.lav$bestmodel

pheno_long<-names(bm.one.ll$coefficients)[-1] # PDSRll
pheno_radial<-names(bm.one.rr$coefficients)[-1] #ODSRrr
pheno_lav<-names(bm.one.lav$coefficients)[-1] #LAVmaxi
pheno_all_M1<-c(pheno_long,pheno_radial,pheno_lav) # bind all the variables selected from the beSS
pheno_all_M1<-(unique(pheno_all_M1))
phenonames<-colnames(pheno)
pheno_all_M1<-substring(pheno_all_M1,6) # exclude "xbest"
# pheno_all_M1 all variables selected

pos_pheno_M1<-0
for (iP in 1:length(pheno_all_M1)){pos_pheno_M1[iP]<-grep(pheno_all_M1[iP],phenonames)} # position in the pheno of variables selected

# Method 2 - EBIC on the graphical LASSO model implemented in the R package "bootnet"
# which employs a regularizing penalty that can lead to parameter estimates of exactly zero.
# the graphical model was thresholded, the tuning parameter for EBIC was set to 0.5 and only the non-zero 
# associations between the three diastolic function parameters and all other covariates were selected as covariates.

pheno<-as.matrix(na.omit(data_pheno[,-c(1:2)])) # exclude the IDs
net_thresh <- bootnet_EBICglasso(pheno,
                                 tuning = 0.5, # EBICglasso sets tuning to 0.5
                                 threshold = T,unlock=T)
net<-net_thresh[["results"]][["optnet"]]
nll<-which(net[,111]==0)
nrr<-which(net[,112]==0)
nlav<-which(net[,127]==0)
pheno_long<-names(net[-nll,111]) #PDSRll
pheno_radial<-names(net[-nrr,112]) # PDSRrr
pheno_lav<-names(net[-nlav,127]) # LAVmaxi
pheno_all_M2<-c(pheno_long,pheno_radial,pheno_lav) # bind all the variables selected from the EBICglasso
pheno_all_M2<-(unique(pheno_all_M2))
# pheno_all_M2

phenonames<-colnames(pheno)
pos_pheno_M2<-0
for (iP in 1:length(pheno_all_M2)){pos_pheno_M2[iP]<-grep(pheno_all_M2[iP],phenonames)} # position in the pheno of variables selected

# Method 3 - Stability selection with "stabsel"
# uses resampling to assess the stability of selected imaging phenotypes for a robust  selection of covariates.
# "cutoff" was set to 0.95 allowing more variables to be included in the model, the per-family error rate was 
# set to 1.0 and the "fitfun" parameter was set as "glmnet.lasso"

pheno<-as.matrix(na.omit(data_pheno[,-c(1:2)])) # exclude the IDs

stab.glmnet_ll <- stabsel(x=pheno[,-111], y=pheno[,111] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_rr <- stabsel(x=pheno[,-112], y=pheno[,112] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_lav <- stabsel(x=pheno[,-127], y=pheno[,127] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)


pheno_long<-(stab.glmnet_ll$selected) # for PDSRll
pheno_radial<-(stab.glmnet_rr$selected) # PDSRrr
pheno_lav<-(stab.glmnet_lav$selected) # LAVmaxi
pheno_all_M3<-c(names(pheno_long),names(pheno_radial),names(pheno_lav)) # bind all the variables selected from the stabsel
pheno_all_M<-unique(pheno_all_M)
pheno_all_M3

phenonames<-colnames(pheno)
pos_pheno_M3<-0
for (iP in 1:length(pheno_all_M3)){pos_pheno_M3[iP]<-grep(pheno_all_M3[iP],phenonames)} # position in the pheno of variables selected

# Plot stability selection
# plot(stab.glmnet_rr, main="PDSRrr",type="maxsel")
# plot(stab.glmnet_ll, main="PDSRll",type="maxsel")
# plot(stab.glmnet_lav, main="LAVmaxi",type="maxsel")
# Choose stability selection method
pos_pheno_M3<-pos_pheno_M3[order(pos_pheno_M3)]

multivar_data <- read.table("multivar_datatable.txt", header = TRUE) # load the whole dataset

## Final check for collinearity of the selected variables using Variance Inflation Factor (VIF)
multivar_data<-na.omit(multivar_data[,c(pos_pheno_M3)])
model<-lm(`PDSRll`~., data=as.data.frame(multivar_data))
car::vif(model)
multivar_data_vif<-multivar_data[,-c(4,11,16,17,19,22,24)] # exclude variables with VIF > 5 and include one phenotype of each
# of the four cardiac chambers (LV, LA, RV, RA), one of the relevant strains (Err, Ell) and two aortic sections (AAo, DAo) where possible to avoid collinearity.

model_vif<-lm(`LAVmaxi`~., data=as.data.frame(multivar_data_vif))
car::vif(model_vif)
summary(model_vif)
vif_values<-vif(model_vif)
par(las=2)
par(mar=c(5,10,4,6))
barplot(vif_values,main="VIF values",col="steelblue", horiz = T)
imcdiag(model_vif,method="VIF", vif=5) # 0 if collinearity is not detected by this test

position_final<-0
for (iP in 1:ncol(multivar_data_vif)){position_final[iP]<-grep(colnames(multivar_data_vif)[iP],phenonames)} # position in the pheno of variables selected

write.table(position_final, "postion_final.txt")
# END