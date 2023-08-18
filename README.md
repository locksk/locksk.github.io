# Predictors of Absolute Neutrophil Count in the CLOZUK3 sample.

## Mediation and Longitudinal Analysis to interpret the association between clozapine pharmacokinetics, pharmacogenomics, and absolute neutrophil count. 

* CLOZUK3_predictors_ANC.rmd contains the code used for analysis of the CLOZUK3 dataset - see: https://locksk.github.io/clozuk3-anc/
    * .html file contains output
    * .pngs are figures used in the .rmd
    * style.css helps format the .rmd
* PyPGx directory contains code snippets used for calling PGx star alleles.
    * https://github.com/sbslee/pypgx
* CLOZUK3_derive_HLA.R was used to derive HLA genotypes for CLOZUK3 sample, weighting by posterior probabilities to create estimated dosages, and then running linear mixed effect models to test the effect of each HLA genotype in the context of the base model (i.e., ANC ~ dose + clozapine + norclozapine + HLA_genotype + sex + age + age_sqauared + Time between dose and sample + (1|participant ID)).
    * HLA imputation was performed using HIBAG with the InfiniumGlobal-European-HLA4-hg19 prediction model available on their site.
    * https://hibag.s3.amazonaws.com/hlares_index.html

## Data and Packages

The full CLOZUK3 dataset was used for linear mixed-effect regression models and single-mediator analyses, benefitting from multiple FBC/PK measurements per patient over time. 

* lme4 and lmerTest were used to fit the LMMs (Bates et al. 2015; Kuznetsova et al. 2017). The mediation package (Tingley et al. 2014) was used to perform single-mediation analyses in the longitudinal data. 
* The CLOZUK3 dataset was transformed into cross-sectional data through taking the lowest observation of ANC for each participant. This reduced dataset was used for multiple-, and single-mediator analyses using SEM in lavaan (Rosseel 2012), and in replication analyses (i.e., Spearman’s correlations, Linear Models) of previous research (Vaquero-Baez et al. 2019; Willcocks et al. 2021).
* Prior to inclusion in all regression and structural equation models, covariates were standardised (mean centred and scaled) using the datawizard R package (Patil et al. 2022). 


