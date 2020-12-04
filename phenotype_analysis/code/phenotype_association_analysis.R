
install.packages("data.table")
install.packages("glmnet")
install.packages("BeSS")
install.packages("devtools")
library(devtools)
devtools::install_github("SachaEpskamp/bootnet")
install.packages("stabs")
install.packages("car")
install.packages("mctest")
install.packages("berryFunctions")

library(data.table)
library(BeSS)
library("bootnet")
library(stabs)
library(car)
library(mctest)
library(glmnet)
library(berryFunctions)

# load data for multiple linear regression analysis
multidata <- read.table("20200709_strain_plots_all/data_paper/multiple_datatable.txt", header = TRUE)

beta_ml<-matrix(0,ncol=31, nrow=10)
mat_pv<-matrix(0,ncol=31, nrow=10)
t_BH<-matrix(0,ncol=31,1)
rsq<-matrix(0,ncol=31,1)
conflist<-vector(mode="list",length=31)

iT<-1
for (iN in 11:41){
  data_mul<-as.data.frame(multidata[,c(1:10,iN)])
  cv<-lapply(colnames(data_mul)[11], function(x) lm(formula(paste("`",x,"`","~.", sep="")),data=data_mul))
  pval<-summary(cv[[1]])$coefficients[,"Pr(>|t|)"] # get p-values
  p.cor<-p.adjust(pval,method = "BH") # adjust using Benjamini - Hochberg procedure
  names(p.cor)<-NULL
  beta<-as.vector(t(coef(cv[[1]], complete=TRUE)))[-1] # get beta coefficient
  t_BH[iT]<-get_bh_threshold(pval, 0.05) # get BH threshold
  smcv<-summary(cv[[1]])
  rsq[iT]<-smcv$r.squared # rsquared
  conflist[iT]<-list(confint(cv[[1]],level=0.95)[-1,]) # confidence intervals
  
  beta_ml[,iT]<-beta
  mat_pv[,iT]<-p.cor[-1]
  iT<-iT+1
}
# beta coefficients
colnames(beta_ml)<-colnames(multidata)[11:41]
rownames(beta_ml)<-colnames(multidata)[1:10]
# p-values
colnames(mat_pv)<-colnames(multidata)[11:41]
rownames(mat_pv)<-colnames(multidata)[1:10]


### Multivariate LASSO regresion analysis with model  selection for selecting the non-imaging phenotypes 

# Method 1 - Best subset selection approach with "beSS"
# using GPDAS algorithm to select the optimal subset selection k and the best model is determined by Extended Bayesian Information Criterion (EBIC)

data_pheno <- read.table("20200709_strain_plots_all/data_paper/Phenotypes_40k.txt", header = TRUE)# load only imaging phenotype data
dataf<-(data_pheno[,c(3:114,131,148,155,172:197)])

forscale<-c(1,4:13,85:141)
for (iS in 1:68){
  dataf[,forscale[iS]]<-scale(dataf[,forscale[iS]])
}
pheno<-as.matrix(na.omit(dataf))

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
# [1] "Age"                    "Sex"                    "SBP"                    "Triglycerides"         
# [5] "PDSRrr"                 "Ecc.Global"             "Err.Global"             "Ell.Global"            
# [9] "AAo.distensibility"     "AAo.min.area"           "DAo.min.area"           "LVSVi"                 
# [13] "LAVmaxi"                "LASVi"                  "LVCI"                   "RAEF"                  
# [17] "Pulse_rate"             "Assessment_centre"      "C_reactive_protein_log" "PDSRll"                
# [21] "WT.Global"              "LVEF"                   "RVSVi"                  "RVEF"                  
# [25] "LAVmini"                        
pos_pheno_M1<-0
for (iP in 1:length(pheno_all_M1)){pos_pheno_M1[iP]<-grep(pheno_all_M1[iP],phenonames)} # position in the pheno of variables selected

# Method 2 - EBIC on the graphical LASSO model implemented in the R package "bootnet"
# which employs a regularizing penalty that can lead to parameter estimates of exactly zero.
# the graphical model was thresholded, the tuning parameter for EBIC was set to 0.5 and only the non-zero 
# associations between the three diastolic function parameters and all other covariates were selected as covariates.

pheno<-as.matrix(na.omit(dataf))
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
# [1] "Age"          "DBP"          "PDSRrr"       "Ecc.Global"   "Err.Global"   "Ell.Global"   "AAo.max.area"
# [8] "AAo.min.area" "DAo.max.area" "DAo.min.area" "LVSVi"        "LAVmaxi"      "LVCI"         "LAEF"        
# [15] "Pulse_rate"   "PDSRll"       "LVMi"         "LASVi"        "RVEDVi"       "RVSVi"        "Sex"         
# [22] "BSA"          "Height"       "LVEDVi"       "LVESVi"       "LAVmini"     
phenonames<-colnames(pheno)
pos_pheno_M2<-0
for (iP in 1:length(pheno_all_M2)){pos_pheno_M2[iP]<-grep(pheno_all_M2[iP],phenonames)} # position in the pheno of variables selected

# Method 3 - Stability selection with "stabsel"
# uses resampling to assess the stability of selected imaging phenotypes for a robust  selection of covariates.
# "cutoff" was set to 0.95 allowing more variables to be included in the model, the per-family error rate was 
# set to 1.0 and the "fitfun" parameter was set as "glmnet.lasso"

pheno<-as.matrix(na.omit(dataf))
stab.glmnet_ll <- stabsel(x=pheno[,-111], y=pheno[,111] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_rr <- stabsel(x=pheno[,-112], y=pheno[,112] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_lav <- stabsel(x=pheno[,-127], y=pheno[,127] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)


pheno_long<-(stab.glmnet_ll$selected) # for PDSRll
pheno_radial<-(stab.glmnet_rr$selected) # PDSRrr
pheno_lav<-(stab.glmnet_lav$selected) # LAVmaxi
pheno_all_M3<-c(names(pheno_long),names(pheno_radial),names(pheno_lav)) # bind all the variables selected from the stabsel
pheno_all_M3<-unique(pheno_all_M3)
# pheno_all_M3
# [1] "Age"                "SBP"                "DBP"                "Pulse_rate"         "Triglycerides"     
# [6] "PDSRrr"             "Ecc.Global"         "Err.Global"         "Ell.Global"         "AAo.distensibility"
# [11] "AAo.min.area"       "DAo.distensibility""LVSVi"              "LAVmaxi"            "LAVmini"           
# [16] "LVCI"               "RAVmini"            "Assessment_centre"  "PDSRll"             "LVMi"              
# [21] "LVEF"               "RVEDVi"             "RVSVi"              "LASVi"             

phenonames<-colnames(pheno)
pos_pheno_M3<-0
for (iP in 1:length(pheno_all_M3)){pos_pheno_M3[iP]<-grep(pheno_all_M3[iP],phenonames)} # position in the pheno of variables selected

# Plot stability selection
# plot(stab.glmnet_rr, main="PDSRrr",type="maxsel")
# plot(stab.glmnet_ll, main="PDSRll",type="maxsel")
# plot(stab.glmnet_lav, main="LAVmaxi",type="maxsel")
pos_pheno<-c(pos_pheno_M1,pos_pheno_M2,pos_pheno_M3) # combine all the variables selected 
pos_pheno<-unique(pos_pheno)
pos_pheno<-pos_pheno[order(pos_pheno)]
# colnames(pheno)[pos_pheno]
# [1] "Age"                    "Sex"                    "BSA"                    "SBP"                   
# [5] "DBP"                    "Pulse_rate"             "Height"                 "Assessment_centre"     
# [9] "C_reactive_protein_log" "Triglycerides"          "PDSRll"                 "PDSRrr"                
# [13] "Ecc.Global"             "Err.Global"             "Ell.Global"             "WT.Global"             
# [17] "AAo.distensibility"     "AAo.max.area"           "AAo.min.area"           "DAo.max.area"          
# [21] "DAo.min.area"           "LVEDVi"                 "LVESVi"                 "LVSVi"                 
# [25] "LVMi"                   "LAVmaxi"                "LAVmini"                "LASVi"                 
# [29] "LVEF"                   "LVCI"                   "LAEF"                   "RVEDVi"                
# [33] "RVSVi"                  "RAVmini"                "RVEF"                   "RAEF"                  

multivar_data <- read.table("20200709_strain_plots_all/data_paper/multivar_datatable.txt", header = TRUE) # load the whole dataset
multivar_data_train <- read.table("20200709_strain_plots_all/data_paper/multivar_train_datatable.txt", header = TRUE) # load the training dataset
multivar_data_test <- read.table("20200709_strain_plots_all/data_paper/multivar_test_datatable.txt", header = TRUE) # load the test dataset

## Final check for collinearity of the selected variables using Variance Inflation Factor (VIF)
multivar_data<-na.omit(multivar_data[,pos_pheno])
model<-lm(`LAVmaxi`~., data=as.data.frame(multivar_data))
car::vif(model)
multivar_data_vif<-multivar_data[,-c(3,5,7,13,16,19:24,26,29:30,32:33,36,37)] # exclude variables with VIF > 5 and include one phenotype of each
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


# Apply LASSO regression

# cv.glmnet to train for the lambda parameter
data.train<-as.matrix(multivar_data_train[,position_final])
data.train<-na.omit(data.train)
colnames(data.train)
lambda_min<-matrix(0,ncol = 1, nrow = ncol(data.train))
for (iT in 1:ncol(data.train)){
  if (iT==2){ # for the logistic regression
    cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, family="binomial", alpha=1)$lambda.min
  } else {
    cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, alpha=1)$lambda.min
  }
  lambda_min[iT]<-round(cv,5)
}

# glmnet test
data_test<-as.matrix(multivar_data_test[,position_final])
data_test<-na.omit(data_test)
beta_gl<-matrix(0,ncol = ncol(data_test), nrow = ncol(data_test)-1)
for (iS in 1:ncol(data_test)){
  if (iS==2){ # for the logistic regression
    cv<-glmnet(data_test[,-iS],data_test[,iS], family="binomial", lambda = lambda_min[iS], alpha=1)
    beta<-as.vector(t(coef(cv)))
    beta<-as.matrix(beta[-1])
    
  } else {
    cv<-glmnet(data_test[,-iS],data_test[,iS],lambda = lambda_min[iS], alpha = 1)
    beta<-as.vector(t(coef(cv)))
    beta<-as.matrix(beta[-1])
  }
  beta_gl[,iS]<-beta[,1]
  
}

colnames(beta_gl)<-colnames(data_test)
beta_gl<-as.data.frame(beta_gl)
multivar_lasso<-matrix(0,nrow = ncol(beta_gl),ncol=ncol(beta_gl))
multivar_lasso<-as.data.frame(multivar_lasso)
# insert rows so the rows of correlation matrix will have the same dimensions as the columns
iEx<-1
for (ib in 0:(ncol(multivar_lasso)-1)){
  multivar_lasso[,iEx]<-insertRows(as.data.frame(beta_gl[,iEx]),ib,new=0,rcurrent=TRUE)
  iEx<-iEx+1
}
colnames(multivar_lasso)<-colnames(beta_gl)
rownames(multivar_lasso)<-colnames(multivar_lasso)
# Save
# write.table(multivar_lasso, "multivar_lasso.txt", row.names = T, col.names = T) # save beta coefficients

# Plot 
# col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
# corrplot(as.matrix(multivar_lasso), is.corr = F, method = "color", col = col(200),
#          number.digits = 3, number.cex = .55,
#          tl.cex = .8,
#          addCoef.col = "black",# Add coefficient of correlation
#          tl.col = "black", tl.srt = 90,cl.cex=.8, cl.pos = "r"#, # Text label color and rotation
# )

# END
