
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

### LASSO regresion analysis with stability selection for selecting the phenotypes 

# We applied feature selection algorithms using stability selection method that uses resampling to assess the stability of selected phenotypes for a 
# robust selection of covariates using the R package "stabs".
# We fitted LASSO models via stability selection procedure setting as predictors the three diastolic function parameters in order to estimate stable 
# statistical associations between the variables. The selection probability "cutoff" was set to 0.95, the per-family error rate was set to 1.0 and the 
# "fitfun" parameter of the "stabsel" function was set as "glmnet.lasso".

# The covariates selected from the model selection method were checked for collinearity by computing the variance inflation factor (VIF). 
# Finally, we excluded the covariates with high VIF values (VIF > 5) in order to avoid collinearity.

## Input:
## - pheno: Data with 110 non-imaging imputed phenotypes/variables and 31 imaging phenotype data.
## Output:
## - position_final: position of the variables selected in the "pheno" data.


stab.glmnet_ll <- stabsel(x=pheno[,-111], y=pheno[,111] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_rr <- stabsel(x=pheno[,-112], y=pheno[,112] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)
stab.glmnet_lav <- stabsel(x=pheno[,-127], y=pheno[,127] ,fitfun = glmnet.lasso, args.fitfun = list(alpha=1), cutoff = 0.95, PFER =1, B=100)


pheno_long<-(stab.glmnet_ll$selected) # for PDSRll
pheno_radial<-(stab.glmnet_rr$selected) # PDSRrr
pheno_lav<-(stab.glmnet_lav$selected) # LAVmaxi
pheno_all_M<-c(names(pheno_long),names(pheno_radial),names(pheno_lav)) # bind all the variables selected from the stabsel
pheno_all_M

phenonames<-colnames(pheno)
pos_pheno_M<-0
for (iP in 1:length(pheno_all_M)){pos_pheno_M[iP]<-grep(pheno_all_M[iP],phenonames)} # position in the pheno of variables selected
# write.table(pos_pheno_M, "position_stabsel.txt")

# Plot stability selection
# plot(stab.glmnet_rr, main="PDSRrr",type="maxsel")
# plot(stab.glmnet_ll, main="PDSRll",type="maxsel")
# plot(stab.glmnet_lav, main="LAVmaxi",type="maxsel")
# Choose stability selection method
pos_pheno_M<-pos_pheno_M[order(pos_pheno_M)]


## Final check for collinearity of the selected variables using Variance Inflation Factor (VIF)
pheno<-na.omit(pheno[,c(pos_pheno_M)])
model<-lm(`PDSRll`~., data=as.data.frame(pheno))
car::vif(model)
pheno_vif<-pheno[,-c(4,11,16,17,19,22,24)] 
# After visual inspection, we exclude variables with VIF > 5 and include one phenotype of each
# of the four cardiac chambers (LV, LA, RV, RA), one of the relevant strains (Err, Ell) and 
# two aortic sections (AAo, DAo) where possible to avoid collinearity.

model_vif<-lm(`LAVmaxi`~., data=as.data.frame(pheno_vif))
car::vif(model_vif)
summary(model_vif)
vif_values<-vif(model_vif)
par(las=2)
par(mar=c(5,10,4,6))
barplot(vif_values,main="VIF values",col="steelblue", horiz = T)
imcdiag(model_vif,method="VIF", vif=5) # 0 if collinearity is not detected by this test

position_final<-0
for (iP in 1:ncol(pheno_vif)){position_final[iP]<-grep(colnames(pheno_vif)[iP],phenonames)} # position in the pheno of variables selected

# Save
# write.table(position_final, "position_final.txt")

# END
