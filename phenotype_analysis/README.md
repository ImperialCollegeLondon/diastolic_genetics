# Cardiac Phenotype associations Pipeline

## Multivariable analysis using LASSO models with stability selection for selecting the phenotypes
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

   
    ## position_stab: position of phenotypes for analysis
    
    ## Final check for collinearity using the selected variables

    model<-lm(PDSRll~., data=as.data.frame(multivar_data[,position_stab[,1]]))
    library(mctest)
    imcdiag(model,method="VIF", vif=5) # 0 if collinearity is not detected by this test
    
    Call:
    imcdiag(mod = model, method = "VIF", vif = 5)


    VIF Multicollinearity Diagnostics

                                      VIF detection
    Age                            2.4588         0
    Sex                            2.5576         0
    BSA                            3.8769         0
    SBP                            2.4555         0
    DBP                            2.1936         0
    `Pulse rate`                   1.8521         0
    Diabetes                       1.4502         0
    Smoking                        1.0303         0
    `Duration of activity`         1.0144         0
    `Medication (n)`               1.1927         0
    `Assessment centre`            1.0708         0
    HbA1c                          1.4123         0
    `C-reactive protein            1.2030         0
    HDL                            1.6472         0
    Glucose                        1.3148         0
    Triglycerides                  1.4476         0
    `eGFR cystatin`                1.4778         0
     PDSRrr                        2.1762         0
    `Err Global`                   2.6223         0
    `Ell Global`                   1.5321         0
    `AAo distensibility`           2.9105         0
    `DAo distensibility`           2.8771         0
    LVSVi                          3.2973         0
    LAVmaxi                        3.8564         0
    LAVmini                        3.8502         0
    LVEF                           2.8855         0
    LVCO                           3.7458         0
    LVCI                           3.9446         0
    RVSVi                          2.8715         0


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
    multivar_lasso<-matrix(0,nrow = ncol(beta_gl),ncol=ncol(beta_gl))
    multivar_lasso<-as.data.frame(multivar_lasso)


## Circos plot  

<img src="circos_plot.JPG" alt="" class="inline" />
