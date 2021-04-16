install.packages("mice")

library(mice)

# An imputation technique for handling missing values with multiple imputations using predictive mean 
# matching as imputation method for continuous variables, with 5 imputations and 5 iterations in the 
# R package "mice".

## Input:
## - data_pheno: Data with all non-imaging phenotypes/variables for imputation.
## Output:
## - data_imp: Data with all non-imaging imputed phenotypes/variables.

mice_imp=mice(data_pheno,maxit=0)
meth<-mice_imp$method
no<-c(3)
meth[no]<-"" # do not impute race

mice_impute=mice(data_pheno,m=5,maxit=5, method = meth, print = T)
impute_data=complete(mice_impute,5) # get the last imputation
data_imp<-as.data.frame(impute_data) # imputed data
