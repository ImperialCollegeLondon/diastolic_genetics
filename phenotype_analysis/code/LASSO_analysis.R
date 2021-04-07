
install.packages("data.table")
install.packages("glmnet")
install.packages("berryFunctions")

library(data.table)
library(glmnet)
library(berryFunctions)

# We fitted L1-regularized logistic or linear regression (LASSO) that optimises the model coefficient of the linear regression,
# where the λ parameter represents the strength of the regularization.
# We adjusted λ by a ten-fold cross-validation (CV) method on a training set (68%, n=26,893) including all covariates, 
# using the cv.glmnet function from the "glmnet" R package. The lambda.min parameter, which denotes the value that gives 
# minimum mean cross-validated error, was extracted and used for prediction on the test set (32%, n=12,666).

## Input:
## - multivar_data_train: A training set (68%, n=26,893) including all covariates (110 non-imaging imputed phenotypes/variables and 31 imaging phenotype data).
## - multivar_data_test: A test set (32%, n=26,893) including all covariates (110 non-imaging imputed phenotypes/variables and 31 imaging phenotype data).
## - position_final: position of the variables selected in the data. In order to create the plot in Extended Data Fig.5b, "position_final" was not added in the code!
## Output:
## - multivar_lasso: Outcome of the LASSO regression analysis, containing the regression coefficients for all associations used in the Extended Data Fig.4b and circos plot in Fig.3.

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
# save beta coefficients

# END
