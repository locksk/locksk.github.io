# Predictors of Absolute Neutrophil Count in the CLOZUK3 sample.

## Details

* CLOZUK3_predictors_ANC.rmd contains the code used for analysis of the CLOZUK3 dataset.
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

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
