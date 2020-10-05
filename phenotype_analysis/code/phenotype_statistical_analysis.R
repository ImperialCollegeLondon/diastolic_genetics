
install.packages("data.table")

library(data.table)

# load data for multiple linear regression analysis
multidata <- read.table("multiple_datatable.txt", header = TRUE)

beta_ml<-matrix(0,ncol=30, nrow=11)
mat_pv<-matrix(0,ncol=30, nrow=11)
t_BH<-matrix(0,ncol=30,1)
rsq<-matrix(0,ncol=30,1)
conflist<-vector(mode="list",length=30)

iT<-1
for (iN in 12:41){
  cv<-lapply(colnames(data_na)[12], function(x) lm(formula(paste("`",x,"`","~.", sep="")),data=data_na))
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
colnames(beta_ml)<-colnames(multidata)[12:41]
rownames(beta_ml)<-colnames(multidata)[1:11]

### Multivariable LASSO regresion analysis with stability selection for selecting the non-imaging phenotypes 

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
multivar_data <- read.table("multivar_datatable.txt", header = TRUE) # load the whole dataset
multivar_data_train <- read.table("multivar_train_datatable.txt", header = TRUE) # load the training dataset
multivar_data_test <- read.table("multivar_test_datatable.txt", header = TRUE) # load the test dataset
colnames(multivar_data_train)<-colnames(multivar_data)
covar <- read.table("cov.txt", header = TRUE) # the position of all non-imaging covariates
pos_pheno <- read.table("fit.txt", header = TRUE) # the position of imaging covariates
## load data for training and for analysis
## covar: position of the covariates to bind with position of phenotypes for analysis
position_stab<-rbind(covar,pos_pheno)

## Final check for collinearity using the selected variables

install.packages("mctest")

model<-lm(`PDSRll (s-1)`~., data=as.data.frame(multivar_data[,position_stab[,1]]))
library(mctest)
imcdiag(model,method="VIF", vif=5) # 0 if collinearity is not detected by this test

# Apply LASSO regression

library(glmnet)  

# cv.glmnet to train for the lambda parameter
data.train<-as.matrix(multivar_data_train[,position_stab[,1]])
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
data_selected<-as.matrix(multivar_data_test[,position_stab[,1]])
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
multivar_lasso<-matrix(0,nrow = ncol(beta_gl),ncol=ncol(beta_gl))
multivar_lasso<-as.data.frame(multivar_lasso)
write.table(multivar_lasso, "multivar_lasso.txt", row.names = T, col.names = T) # save beta coefficients from LASSO regression analysis

# END
