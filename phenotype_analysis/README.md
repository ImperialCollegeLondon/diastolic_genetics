# Cardiac Phenotype associations Pipeline

## Multivariable analysis using LASSO models with stability selection for selecting the imaging phenotypes
### Stability selection procedure using 'stabsel'

    # load non-imaging phenotype data
    pheno<-as.matrix(data_pheno)
    library(stabs)
    # for PDSRll (s-1)
    library(stabs)

    stab.glmnet <- stabsel(x=pheno[,-1], y=pheno[,1] ,
                          fitfun = glmnet.lasso,
                          args.fitfun = list(alpha=1),
                          cutoff = 0.75, PFER =1, B=100)

    pheno_long<-as.data.frame(stab.glmnet$selected) # repeat for PDSRrr
    pheno_radial<-as.data.frame(stab.glmnet$selected) # repeat for LAVmaxi
    pheno_lav<-as.data.frame(stab.glmnet$selected) # bind all the variables selected from the stabsel
    pheno_all<-rbind(pheno_long,pheno_radial,pheno_lav)

    pos_pheno<-match(rownames(pheno_all),colnames(pheno)) # order the position of variables selected

 ### LASSO models using 'glmnet'. The parameter of lambda.min was tuned by a 10-fold cross-validation method using 'cv.glmnet' on a training set (~67% of the original dataset).

   
    ## load data for training and for analysis
    ## covar: position of the covariates to bind with position of phenotypes for analysis
    position_stab<-rbind(covar,pos_pheno)
    
    ## Final check for collinearity using the selected variables

    model<-lm(PDSRll~., data=as.data.frame(multivar_data[,position_stab[,1]]))
    library(mctest)
    imcdiag(model,method="VIF", vif=5) # 0 if collinearity is not detected by this test
    
    Call:
    imcdiag(mod = model, method = "VIF", vif = 5)


    Call:
    imcdiag(mod = model1, method = "VIF", vif = 5)


     VIF Multicollinearity Diagnostics

                                 VIF detection
     Age                      1.9780         0
     Sex                      1.6394         0
     SBP                      1.3519         0
     `Pulse rate`             1.6454         0
     Diabetes                 1.1216         0
     Smoking                  1.0122         0
     `Duration of activity`   1.0114         0
     `Medication (n)`         1.1984         0
     `Assessment centre`      1.0688         0
     `C-reactive protein`     1.1589         0
     Cholesterol              1.2081         0
     Triglycerides            1.4019         0
     `eGFR cystatin`          1.4250         0
     PDSRrr                   2.1422         0
     `Err Global`             2.0268         0
     `Ell Global`             1.4646         0
     `AAo distensibility`     2.8909         0
     `DAo distensibility`     2.8495         0
     LVSVi                    4.5815         0
     LVCI                     2.9617         0
     LAVmaxi                  4.4977         0
     LAVmini                  4.0171         0
     RVSVi                    2.7339         0
     RAVmini                  1.5510         0

     NOTE:  VIF Method Failed to detect multicollinearity


     0 --> COLLINEARITY is not detected by the test

    ===================================


    ## Apply LASSO regression

    library(glmnet)  
    
    data.train<-as.matrix(multivar_data_train[,position_stab[,1]])

    lambda_min<-matrix(0,ncol = 1, nrow = ncol(data.train))
    for (iT in 1:ncol(data.train)){
      if (iT==2|iT==7){ # the columns with binary data
        cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, family="binomial", alpha=1)$lambda.min
      } else {
        cv<-cv.glmnet(data.train[,-iT],data.train[,iT],nfolds = 10, alpha=1)$lambda.min
      }
      lambda_min[iT]<-round(cv,5)
    }
    
    # LASSO model using glmnet 
    
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
    multivar_beta<-matrix(0,nrow = ncol(beta_gl),ncol=ncol(beta_gl))
    multivar_beta<-as.data.frame(multivar_beta)


## Circos plot  

<img src="circos_plot.JPG" alt="" class="inline" />
