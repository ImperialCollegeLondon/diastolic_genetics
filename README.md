# Genetic and environmental determinants of diastolic heart function

[![DOI](https://zenodo.org/badge/291437734.svg)](https://zenodo.org/badge/latestdoi/291437734)

## Abstract

Diastole is the sequence of physiological events that occur in the heart during ventricular filling and principally depends on myocardial relaxation and chamber stiffness. Abnormal diastolic function is related to many cardiovascular disease processes and is predictive of health outcomes, but its genetic architecture is largely unknown. Here, we use machine learning cardiac motion analysis to measure diastolic functional traits in UK Biobank participants to identify environmental modifiers and associated common genetic variants.  

## Content

* [Cardiac Image analysis (ukbb_cardiac)](https://github.com/baiwenjia/ukbb_cardiac)   
Automated pipeline for image segmentation and motion analysis.

* [Peak detection](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/peak_detection)  
Automated peak detection algorithm to obtain the diastolic peak. This repository includes two approaches for detecting peaks in cardiac strain data:
  - *[Baseline Method](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/peak_detection/BaselineMethod)*: Uses local extrema detection with a predefined look-ahead window for peak identification.  
  - *[Advanced Method](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/peak_detection/AdvancedMethod)*: A new technique using *Expectile Generalized Additive Models (GAMs)* to enhance peak detection accuracy and robustness against noise.
 
* [Phenotype anaysis](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/phenotype_analysis)   
Cardiac phenotype pipeline to investigate the relationships between the myocardial diastolic function parameters and the non-imaging clinical data.

* [Genetic anaysis](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/genetic_analysis)  
Genome-wide association study (GWAS) analysis for the derived diastolic function traits. Contains also the causality assessment and polygenic risk score calculations.

* [Plots](https://github.com/ImperialCollegeLondon/diastolic_genetics/tree/master/plots)  
Scatter and box plots of diastolic function.

## License

Distributed under the GNU GENERAL PUBLIC LICENSE license. See ``LICENSE`` for more information.

## Citation

Thanaj M, Mielke J, McGurk KA, Bai W, Savioli N, de Marvao A, Meyer HV, Zeng L, Sohler F, Lumbers RT, Wilkins MR, Ware JS, Bender C, Rueckert D, MacNamara A, Freitag DF, O'Regan DP. Genetic and environmental determinants of diastolic heart function. _Nat Cardiovasc Res_. 2022;1(4):361-371. doi: [10.1038/s44161-022-00048-2](https://doi.org.10.1038/s44161-022-00048-2)

