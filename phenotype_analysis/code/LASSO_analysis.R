
install.packages("data.table")
install.packages("glmnet")
install.packages("berryFunctions")

library(data.table)
library(glmnet)
library(berryFunctions)

# load data for LASSO regression analysis
multivar_data_train <- read.table("multivar_train_datatable.txt", header = TRUE) # load the training dataset
multivar_data_test <- read.table("multivar_test_datatable.txt", header = TRUE) # load the test dataset
position_final<-as.matrix(read.table("position_final.txt"))
position_final<-as.numeric(position_final) # make numeric
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
