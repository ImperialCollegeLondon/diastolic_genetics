
install.packages("mice")

library(mice)

data_pheno<-read.table("Phenotypes_40k.txt", header=T)
data_p<-as.data.frame(data_pheno[,c(4:11,13,16:17,85:110)]) # impute only certain variables

mice_imp=mice(data_p,maxit=0)
meth<-mice_imp$method
no<-c(3)
meth[no]<-"" # do not impute race

mice_impute=mice(data_p,m=5,maxit=5, method = meth, print = T)
impute_data=complete(mice_impute,5) # get the last imputation
data_imp<-as.data.frame(impute_data) # imputed data

# Save
# write.table(data_imp,"Phenotypes_40k_imp.txt", col.names = T)