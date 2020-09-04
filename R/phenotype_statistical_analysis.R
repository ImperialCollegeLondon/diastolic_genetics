
install.packages("data.table")

library(data.table)
multidata <- read.table("multiple_datatable.txt", header = TRUE)

# load data for multiple linear regression analysis

lm_l<-lm(long_PDSR~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
co_l<-coef(lm_l, complete=TRUE)
pval<-anova(lm_l)$`Pr(>F)`
slml<-summary(lm_l)

slml$r.squared # r squared
confint(lm_l,level=0.95) # confidence intervals

lm_r<-lm(radial_PDSR~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
co_r<-coef(lm_r, complete=TRUE)
pvalr<-anova(lm_r)$`Pr(>F)`
slmr<-summary(lm_r) 

slmr$r.squared # r squared
confint(lm_r,level=0.95) # confidence intervals

lm_la<-lm(LAVmax~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
co_la<-coef(lm_la, complete=TRUE)
pvalla<-anova(lm_la)$`Pr(>F)`
slmla<-summary(lm_la)

slmla$r.squared # r squared
confint(lm_la,level=0.95) # confidence intervals


### Multivariate LASSO regresion analysis with stability selection for selecting the non-imaging phenotypes 

install.packages("glmnet")
install.packages("stabs")

data_pheno <- read.table("pheno_datatable.txt", header = TRUE)# load only imaging phenotype data

# stabsel

library(stabs)

pheno<-as.matrix(data_pheno)
stab.glmnet <- stabsel(x=pheno[,-1], y=pheno[,1] ,
                       fitfun = glmnet.lasso,
                       args.fitfun = list(alpha=1),
                       cutoff = 0.75, PFER =1, B=100)

pheno_long<-as.data.frame(stab.glmnet$selected) # repeat for PDSRrr
pheno_radial<-as.data.frame(stab.glmnet$selected) # repeat for LAVmaxi
pheno_lav<-as.data.frame(stab.glmnet$selected) # bind all the variables selected from the stabsel
pheno_all<-rbind(pheno_long,pheno_radial,pheno_lav)

pos_pheno<-match(rownames(pheno_all),colnames(pheno)) # order the position of variables selected


library(data.table)
multivariate_data <- read.table("multivar_datatable.txt", header = TRUE) # load the whole dataset
multivariate_data_train <- read.table("multivar_train_datatable.txt", header = TRUE) # load the training dataset
colnames(multivariate_data_train)<-colnames(multivariate_data)
covar <- read.table("cov.txt", header = TRUE) # the position of all non-imaging covariates
pos_pheno <- read.table("fit.txt", header = TRUE) # the position of imaging covariates
## load data for training and for analysis
## covar: position of the covariates to bind with position of phenotypes for analysis
position_stab<-rbind(covar,pos_pheno)

## Final check for collinearity using the selected variables

install.packages("mctest")

model<-lm(`PDSRll (s-1)`~., data=as.data.frame(multivariate_data[,position_stab[,1]]))
library(mctest)
imcdiag(model,method="VIF", vif=5) # 0 if collinearity is not detected by thie test

# Apply LASSO regression

library(glmnet)  

# cv.glmnet to train for the lambda parameter
data.train<-as.matrix(multivariate_data_train[,position_stab[,1]])
data.train<-na.omit(data.train)
colnames(data.train)
lambda_min<-matrix(0,ncol = 1, nrow = ncol(data.train))
for (iT in 1:ncol(data.train)){
  if (iT==2|iT==7){ # for the logistic regression
    cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, family="binomial", alpha=1)$lambda.min
  } else {
    cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, alpha=1)$lambda.min
  }
  lambda_min[iT]<-round(cv,5)
}

# glmnet 
data_selected<-as.matrix(multivariate_data[,position_stab[,1]])
data_selected<-na.omit(data_selected)
beta_gl<-matrix(0,ncol = ncol(data_selected), nrow = ncol(data_selected)-1)
for (iS in 1:ncol(data_selected)){
  if (iS==2|iS==7){ # for the logistic regression
    cv<-glmnet(data_selected[,-iS],data_selected[,iS], family="binomial", lambda = lambda_min[iS], alpha=1)
    beta<-as.vector(t(coef(cv)))
    beta<-as.matrix(beta[-1])
    
  } else {
    cv<-glmnet(data_selected[,-iS],data_selected[,iS],lambda = lambda_min[iS], alpha = 1)
    beta<-as.vector(t(coef(cv)))
    beta<-as.matrix(beta[-1])
  }
  beta_gl[,iS]<-beta[,1]
  
}

colnames(beta_gl)<-colnames(data_selected)
beta_gl<-as.data.frame(beta_gl)
multivariate_beta<-matrix(0,nrow = ncol(beta_gl),ncol=ncol(beta_gl))
multivariate_beta<-as.data.frame(multivariate_beta)
write.table(multivariate_beta, "multivariate_beta.txt", row.names = T, col.names = T) # save beta coefficients

# END