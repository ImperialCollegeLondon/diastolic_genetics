# Cardiac Phenotype associations Pipeline
## Multiple linear regression analysis

    # load data for multiple linear regression analysis

    lm_l<-lm(long_PDSR~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
    co_l<-coef(lm_l, complete=TRUE)
    pval<-anova(lm_l)$`Pr(>F)`
    slml<-summary(lm_l)
    slml

    ## 
    ## Call:
    ## lm(formula = `PDSRll (s-1)` ~ Age + Sex + BSA + SBP + DBP + Pulse_rate + 
    ##     Diabetes + Smoking + Number_Medication + Moderate_activity + 
    ##     Vigorous_activity + Assessment_centre, data = multidata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5426 -0.5571 -0.0476  0.5093  5.2277 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.1439174  0.0078140  18.418  < 2e-16 ***
    ## Age               -0.3227934  0.0050451 -63.982  < 2e-16 ***
    ## Sex               -0.2896075  0.0116813 -24.793  < 2e-16 ***
    ## BSA               -0.0318173  0.0058983  -5.394 6.92e-08 ***
    ## SBP                0.0259642  0.0063461   4.091 4.30e-05 ***
    ## DBP               -0.1668480  0.0062499 -26.696  < 2e-16 ***
    ## Pulse_rate        -0.3166845  0.0046265 -68.451  < 2e-16 ***
    ## Diabetes          -0.1489125  0.0198740  -7.493 6.88e-14 ***
    ## Smoking           -0.0008847  0.0045455  -0.195   0.8457    
    ## Number_Medication -0.0082938  0.0046728  -1.775   0.0759 .  
    ## Moderate_activity -0.0069714  0.0045313  -1.539   0.1239    
    ## Vigorous_activity -0.0035229  0.0045418  -0.776   0.4380    
    ## Assessment_centre -0.0602252  0.0044104 -13.655  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.853 on 38910 degrees of freedom
    ##   (636 observations deleted due to missingness)
    ## Multiple R-squared:  0.2725, Adjusted R-squared:  0.2723 
    ## F-statistic:  1215 on 12 and 38910 DF,  p-value: < 2.2e-16

    confint(lm_l,level=0.95)

    ##                          2.5 %        97.5 %
    ## (Intercept)        0.128601722  0.1592330253
    ## Age               -0.332681861 -0.3129049608
    ## Sex               -0.312503038 -0.2667119535
    ## BSA               -0.043378081 -0.0202565256
    ## SBP                0.013525683  0.0384028048
    ## DBP               -0.179097966 -0.1545980415
    ## Pulse_rate        -0.325752505 -0.3076165737
    ## Diabetes          -0.187866089 -0.1099589527
    ## Smoking           -0.009794093  0.0080246773
    ## Number_Medication -0.017452639  0.0008649771
    ## Moderate_activity -0.015852824  0.0019099820
    ## Vigorous_activity -0.012424878  0.0053791086
    ## Assessment_centre -0.068869772 -0.0515805981

    lm_r<-lm(`PDSRrr (s-1)`~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
    co_r<-coef(lm_r, complete=TRUE)
    pvalr<-anova(lm_r)$`Pr(>F)`
    slmr<-summary(lm_r)
    slmr

    ## 
    ## Call:
    ## lm(formula = PDSRrr (s-1) ~ Age + Sex + BSA + SBP + DBP + Pulse_rate + 
    ##     Diabetes + Smoking + Number_Medication + Moderate_activity + 
    ##     Vigorous_activity + Assessment_centre, data = multidata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5705 -0.5616 -0.0044  0.5664  3.8891 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.187705   0.007979  23.526  < 2e-16 ***
    ## Age               -0.262127   0.005149 -50.905  < 2e-16 ***
    ## Sex               -0.364371   0.011926 -30.553  < 2e-16 ***
    ## BSA               -0.046981   0.006021  -7.803 6.18e-15 ***
    ## SBP                0.090346   0.006477  13.948  < 2e-16 ***
    ## DBP               -0.161751   0.006375 -25.374  < 2e-16 ***
    ## Pulse_rate        -0.335266   0.004725 -70.959  < 2e-16 ***
    ## Diabetes          -0.141469   0.020313  -6.964 3.35e-12 ***
    ## Smoking           -0.014065   0.004638  -3.033  0.00243 ** 
    ## Number_Medication -0.020682   0.004775  -4.331 1.49e-05 ***
    ## Moderate_activity -0.006512   0.004620  -1.410  0.15867    
    ## Vigorous_activity  0.003042   0.004632   0.657  0.51135    
    ## Assessment_centre -0.069324   0.004504 -15.393  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8681 on 38687 degrees of freedom
    ##   (859 observations deleted due to missingness)
    ## Multiple R-squared:  0.2466, Adjusted R-squared:  0.2464 
    ## F-statistic:  1055 on 12 and 38687 DF,  p-value: < 2.2e-16

    confint(lm_r,level=0.95)

    ##                          2.5 %       97.5 %
    ## (Intercept)        0.172066821  0.203343231
    ## Age               -0.272219472 -0.252033958
    ## Sex               -0.387745717 -0.340996168
    ## BSA               -0.058781282 -0.035180343
    ## SBP                0.077650646  0.103041732
    ## DBP               -0.174246052 -0.149256505
    ## Pulse_rate        -0.344526401 -0.326005005
    ## Diabetes          -0.181282512 -0.101655037
    ## Smoking           -0.023154940 -0.004974810
    ## Number_Medication -0.030041620 -0.011321942
    ## Moderate_activity -0.015567673  0.002542991
    ## Vigorous_activity -0.006037371  0.012122093
    ## Assessment_centre -0.078151341 -0.060496724

    lm_la<-lm(`LAVmaxi (mL/m2)`~Age+Sex+BSA+SBP+DBP+Pulse_rate+Diabetes+Smoking+Number_Medication+Moderate_activity+Vigorous_activity+Assessment_centre, data=multidata)
    co_la<-coef(lm_la, complete=TRUE)
    pvalla<-anova(lm_la)$`Pr(>F)`
    slmla<-summary(lm_la)
    slmla

    ## 
    ## Call:
    ## lm(formula = `LAVmaxi (mL/m2)` ~ Age + Sex + BSA + SBP + DBP + Pulse_rate + 
    ##     Diabetes + Smoking + Number_Medication + Moderate_activity + 
    ##     Vigorous_activity + Assessment_centre, data = multidata)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1673 -0.6416 -0.0689  0.5494  5.9962 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.073534   0.008872   8.288  < 2e-16 ***
    ## Age               -0.078939   0.005739 -13.755  < 2e-16 ***
    ## Sex               -0.141767   0.013274 -10.680  < 2e-16 ***
    ## BSA                0.058865   0.006730   8.746  < 2e-16 ***
    ## SBP                0.194435   0.007198  27.012  < 2e-16 ***
    ## DBP               -0.119673   0.007104 -16.846  < 2e-16 ***
    ## Pulse_rate        -0.228892   0.005269 -43.443  < 2e-16 ***
    ## Diabetes          -0.119099   0.022735  -5.239 1.63e-07 ***
    ## Smoking           -0.002570   0.005169  -0.497   0.6191    
    ## Number_Medication  0.066212   0.005348  12.381  < 2e-16 ***
    ## Moderate_activity  0.011990   0.005171   2.318   0.0204 *  
    ## Vigorous_activity  0.042714   0.005154   8.287  < 2e-16 ***
    ## Assessment_centre -0.078100   0.004997 -15.629  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9577 on 37991 degrees of freedom
    ##   (1555 observations deleted due to missingness)
    ## Multiple R-squared:  0.08303,    Adjusted R-squared:  0.08274 
    ## F-statistic: 286.7 on 12 and 37991 DF,  p-value: < 2.2e-16

    confint(lm_la,level=0.95)

    ##                          2.5 %       97.5 %
    ## (Intercept)        0.056145005  0.090923510
    ## Age               -0.090186864 -0.067690247
    ## Sex               -0.167785074 -0.115748446
    ## BSA                0.045673194  0.072056088
    ## SBP                0.180326738  0.208544178
    ## DBP               -0.133596744 -0.105749409
    ## Pulse_rate        -0.239219048 -0.218564898
    ## Diabetes          -0.163660734 -0.074537417
    ## Smoking           -0.012701331  0.007561499
    ## Number_Medication  0.055730394  0.076693579
    ## Moderate_activity  0.001853694  0.022125399
    ## Vigorous_activity  0.032611690  0.052816534
    ## Assessment_centre -0.087895194 -0.068305563

## Multivariate analysis using LASSO models with stability selection for selecting the imaging phenotypes
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

    model<-lm(`PDSRll (s-1)`~., data=as.data.frame(multivariate_data[,position_stab[,1]]))
    library(mctest)
    imcdiag(model,method="VIF", vif=5) # 0 if collinearity is not detected by this test
    
    Call:
    imcdiag(mod = model2, method = "VIF", vif = 5)


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
    `Cardiac MRI index`            3.9446         0
    PDSRrr                         2.1762         0
    `Err Global`                   2.6223         0
    `Ell Global`                   1.5321         0
    `AAo distensibility`           2.9105         0
    `DAo distensibility`           2.8771         0
    LVSVi                          3.2973         0
    LAVmaxi                        3.8564         0
    LAVmini                        3.8502         0
    LVEF                           2.8855         0
    LVCO                           3.7458         0
    RVSVi                          2.8715         0


    0 --> COLLINEARITY is not detected by the test

    ===================================

    # Apply LASSO regression

    library(glmnet)  
    
    data.train<-as.matrix(multivariate_data_train[,position_stab[,1]])
    colnames(data.train)

    ##  [1] "Age"                              "Sex"                             
    ##  [3] "BSA"                              "SBP"                             
    ##  [5] "DBP"                              "Pulse rate"                      
    ##  [7] "Diabetes"                         "Smoking"                         
    ##  [9] "Duration of activity"             "Medication (n)"                   
    ## [11] "Assessment centre"                "HbA1c"                                          
    ## [13] "C-reactive protein"               "HDL"                                              
    ## [15] "Hb concentration"                 "Triglycerides"                
    ## [17] "eGFR cystatin"                    "Cardiac MRI index"                            
    ## [19] "PDSRll"                           "PDSRrr"                    
    ## [21] "Err Global"                       "Ell Global"                  
    ## [23] "AAo distensibility"               "DAo distensibility"
    ## [25] "LVSVi"                            "LAVmaxi"                 
    ## [27] "LAVmini"                          "LVEF"                        
    ## [29] "LVCO"                             "RVSVi"

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


## Circos plot

     library(circlize)
     library(circlize)
     all_states = rownames(multivariate_beta)
     n_states = nrow(multivariate_beta)
     state_col = c("Age" = "darkgreen",    "Sex" = "darkgreen",
                   "BSA" = "darkgreen",  "SBP" = "darkgreen",
                   "DBP" = "darkgreen",    "Pulse rate" = "darkgreen",
                   "Diabetes" = "darkgreen",     "Smoking" = "darkgreen",
                   "Duration of activity" = "darkgreen",     "Medication (n)" = "darkgreen",
                   "Assessment centre" = "darkgreen", "HbA1c" = "#377EB8",
                   "C-reactive protein" = "#377EB8",  "HDL" = "#377EB8",
                   "Glucose" = "#377EB8", "Triglycerides" = "#377EB8","eGFR cystatin" = "#377EB8",
                   "Cardiac MRI index" = "#E41A1C",    "PDSRll" = "#E41A1C",
                   "PDSRrr" = "#E41A1C",  "Err Global" = "#E41A1C",
                   "Ell Global" = "#E41A1C",    "AAo distensibility" = "#E41A1C",
                   "DAo distensibility" = "#E41A1C",     "LVSVi" = "#E41A1C",
                   "LAVmaxi" = "#E41A1C",     "LAVmini" = "#E41A1C",
                   "LVEF" = "#E41A1C","LVCO" = "#E41A1C", 
                   "RVSVi" = "#E41A1C")

     all_states = names(state_col)
 
     # one for rows and one for columns
     state_col2 = c(state_col, state_col)
     names(state_col2) = c(rownames(multivariate_beta), colnames(multivariate_beta))
     
     colmat = rep(state_col2[rownames(multivariate_beta)], n_states)
     colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
     
     colmat = paste0(colmat, "A0")
     dim(colmat) = dim(multivariate_beta)
     
     circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE) # initialise circos plot

     multivariate_chord = chordDiagram(multivariate_beta, col = colmat, grid.col = state_col2,
                            directional = TRUE, annotationTrack = "grid", 
                            big.gap = 10, small.gap = 1) # plot circos for multivariate_beta
     circos.clear() # clear this cirsos plot

     head(multivariate_chord)
     val<-multivariate_chord$value2
     p<-which(sign(val)==-1) # make all associations absolute
     val[p]<-(-val[p])
     multivariate_chord$value1<-val
     multivariate_chord$value2<-val
     pv<-which(val>=0.4)
     pl<-which(multivariate_chord$rn=='PDSRll')
     pr<-which(multivariate_chord$rn=='PDSRrr')
     pa<-which(multivariate_chord$rn=='LAVmaxi')
     pall<-c(pl,pr,pa)
     pl<-which(multivariate_chord$cn=='PDSRll')
     pr<-which(multivariate_chord$cn=='PDSRrr')
     pa<-which(multivariate_chord$cn=='LAVmaxi')
     pall2<-c(pl,pr,pa)
     pp<-c(pall,pall2)
     pval<-c(pp,pv) # include only the positions with beta coefficient > 0.4 apart from the associations 
                    # between PDSRll, PDSRrr and LAVmaxi and all other phenotypes.
     pval<-unique(pval)
     vpall<-val[pval]
     multivariate_chord$value1[pval]<-vpall
     multivariate_chord$value2[pval]<-vpall
     
     colmat = rep(state_col2[rownames(multivariate_beta)], n_states)
     colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
     dim(colmat) = dim(multivariate_beta)

     multivariate_chord$col[-pval]<-paste0(colmat[-pval], "20") # add faint colour in the links with betas < 0.4
     head(multivariate_chord)

     circos.par(gap.degree=1,canvas.xlim=c(-0.6,0.6), canvas.ylim=c(-1.2,1.2)) # initialise circosplot

     chordDiagram(multivariate_chord, col = multivariate_chord$col, grid.col = state_col2,
                  directional = TRUE, annotationTrack = c("grid"), link.rank = order(multivariate_chord$col),
                  big.gap = 10, small.gap = 1,preAllocateTracks = list(track.height = mm_h(5))) # plot circos for multivariate_chord

     ## Add sector numbers. The numbers in each sector represent the sum of all coefficients for each specific variable. (Optional)
     #
     # for(si in get.all.sector.index()) {
     #   circos.axis(h = "top", labels.cex = 0.4, sector.index = si, track.index = 2)
     # }

     # Add names clockwise
     circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
       xlim = get.cell.meta.data("xlim")
       ylim = get.cell.meta.data("ylim")
       sector.name = get.cell.meta.data("sector.index")
       circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
     }, bg.border = NA)

     # Add small circular rectangles to represent the proportions of different transitions in each variable
     for(i in seq_len(nrow(multivariate_chord))) {
       if(multivariate_chord$value1[i] > 0) {
         circos.rect(multivariate_chord[i, "x1"], -mm_y(1), 
                     multivariate_chord[i, "x1"] - abs(multivariate_chord[i, "value1"]), -mm_y(2), 
                     col = state_col2[multivariate_chord$cn[i]], border = state_col2[multivariate_chord$cn[i]],
                     sector.index = multivariate_chord$rn[i], track.index = 2)
        }
     }

     circos.clear()

     # END    

<img src="plots/circos_final.JPG" alt="" class="inline" />
