
install.packages("data.table")
install.packages("stabs")
install.packages("glmnet")
install.packages("car")
install.packages("mctest")

library(data.table)
library(stabs)
library(glmnet)
library(car)
library(mctest)

### Multivariate LASSO regresion analysis with model  selection for selecting the non-imaging phenotypes 

data_pheno <- read.table("Phenotypes_40k.txt", header = TRUE)# load all 110 non-imaging and 31 imaging phenotype data
pheno<-as.matrix(na.omit(data_pheno))

# Stability selection with "stabsel"
# uses resampling to assess the stability of selected imaging phenotypes for a robust  selection of covariates.
# "cutoff" was set to 0.95 allowing more variables to be included in the model, the per-family error rate was 
# set to 1.0 and the "fitfun" parameter was set as "glmnet.lasso"

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
# write.table(pos_pheno_M3, "position_stabsel.txt")

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

write.table(position_final, "position_final.txt")
# END
