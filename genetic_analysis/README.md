# Genetic analysis for diastolic function parameters

In this folder, we collect for replicating the genetic part of the diastolic function paper. Note that the data loading steps, which highly depend on the set-up of UK Biobank data infrastructure are not included since this will be different for each institute. However, we describe the required input for each function.

## Run GWAS

We have run GWAS for quantitative endpoints (e.g., the diastolic function parameters) with BOLT-LMM and for binary endpoints (e.g. heart failure case/control analysis) with PLINK. The command line statements are summarized in run_gwas.txt.

## Mapping of hits to causal genes and annotation
[Aidan to add]

## Causes and consequences of diastolic function traits

We first screen for associations by conduction a PheWAS for PRS of diastolic function trait. Afterwards, we run a MR analyses. Code for both assessments are summarized in PRS_and_MR_analysis_functions.R (required functions) and PRS_and_MR_analysis_run_analysis.R (execute analysis). 
