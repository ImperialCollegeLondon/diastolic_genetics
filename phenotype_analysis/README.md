# Cardiac Phenotype associations Pipeline

## Multivariable analysis
Multivariable regression analysis was used to explore the phenotype relationship between each diastolic parameter and cardiovascular risk factors. 

## Multivariable analysis using LASSO regression

### Model selection procedure

We employed the model selection approach, setting as predictors the three diastolic function parameters, peak diastolic longitudinal strain rates (PDSRll), peak diastolic radial strain rates (PDSRrr) and left atrial maximum volume indexed to BSA (LAVmaxi), and selecting the variables determined by the stability selection on the least absolute shrinkage and selection operator (LASSO) model, to identify phenotype associated with the trait of interest.

### LASSO regression analysis
We applied LASSO regression that optimises the model coefficient of the linear regression adjusting a 10-fold cross-validation method on a training set using 'cv.glmnet' function tuning the lambda.min parameter, which denotes the value that gives minimum mean cross-validated error, and use it for prediction on the test set.
   
The code has been fully documented and descriptions are available within.