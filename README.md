---
title: "MSc Thesis: Multilevel Multivariate Imputation by Chained Equations through Bayesian Additive Regression Trees"
author: 
    - name: Heleen Brüggen
      orcid: 0009-0008-5180-3125
      email: h.bruggen@uu.nl
      affiliation: Utrecht University
      url: https://www.linkedin.com/in/heleen-br%C3%BCggen-98746024b/
format:
  html:
    toc: true
---
 
# Design
This study is a simulation study to compare the performance of different imputation methods in the context of multilevel data. The imputation methods that are compared are: single level predictive mean matching (PMM); multilevel PMM; single level Bayesian Additive Regression Trees (BART); multilevel BART with random intercepts; and multilevel BART with random intercepts, random slopes, and cross-level interactions. Listwise deletion and analysis of the true/complete data are included as additional benchmarks. These methods are compared on their absolute bias, coverage of their 95% confidence intervals and width of their 95% confidence intervals. 

# Data Generation
This study uses simulated data. The data is generated with the script `Scripts\Simulation data.R`. This script is reproducible. The data can also be found in the folder `data` --- not on Github since the files are too large. The following table shows the design factors used in generating the data. 

| Design factors             | Values    | 
|----------------------------|:----------|
| Number of clusters         | 30, 50    | 
| Within-cluster sample size | 15, 50    | 
| Intraclass Correlation     | .5        |
| Missing data mechanism     | MCAR, MAR |
| Amount of missingness      | 0%, 50%   |

The script generates a 1000 datasets for each combination of design factors. However, due to time constraints, the study used only 100 datasets for each combination of design factors for every imputation method, except the imputation method for multilevel BART with random intercepts, random slopes, and cross-level interactions. The next table shows the amount of datasets used for each imputation method. 

|                                                                                     |      | Number of groups: 30, Group size: 15 | Number of groups: 30, Group size: 50 | Number of groups: 50, Group size: 15 | Number of groups: 50, Group size: 50 |
|-------------------------------------------------------------------------------------|------|--------------------------------------|--------------------------------------|--------------------------------------|--------------------------------------|
| Single-level PMM                                                                    | MAR  | 100                                  | 100                                  | 100                                  | 100                                  |
|                                                                                     | MCAR | 100                                  | 100                                  | 100                                  | 100                                  |
| Multilevel PMM                                                                      | MAR  | 100                                  | 100                                  | 100                                  | 100                                  |
|                                                                                     | MCAR | 100                                  | 100                                  | 100                                  | 100                                  |
| Single-level BART                                                                   | MAR  | 100                                  | 100                                  | 100                                  | 100                                  |
|                                                                                     | MCAR | 100                                  | 100                                  | 100                                  | 100                                  |
| Multilevel BART with random intercepts                                              | MAR  | 100                                  | 100                                  | 100                                  | 100                                  |
|                                                                                     | MCAR | 100                                  | 100                                  | 100                                  | 100                                  |
| Multilevel BART with random intercepts, random slopes, and cross-level interactions | MAR  | 20                                   | 20                                   | 20                                   | 20                                   |
|                                                                                     | MCAR | 100                                  | 40                                   | 100                                  | 20                                   |

# Script recipe
the following steps show the order of the scripts that need to be run to reproduce the results of the study. The scripts are located in the folder `Scripts`. Please note that all paths in the scripts are relative to the root folder of the repository. All software and dependencies used can be found in the file `requirements.txt`. A specific version of the `mice`-package should be downloaded with the following code: 
```
devtools::install_github("heleenbrueggen/mice@impute.mbart")
``` 
Furtermore, keep in mind that some scripts are computationally intensive and may take a long time to run. Additionally, some scripts contain parallel processing, so it is recommended to change the code to use the desired number of cores.

1. `Scripts\Simulation data.R`: Generates the simulation data and saves it in appropriate folders.
2. `Scripts\Missing data generation.R`: Generates the missing data and saves it in appropriate folders.
3. `Scripts\Imputation.R`: Imputes the missing data using the different imputation methods and saves the imputed datasets in appropriate folders.
4. `Scripts\Obtaining results.R`: Obtains the results of the imputation methods and saves them in appropriate folders. This script sources multiple other scripts for an easier running procedure. These scripts are:
    - First the substative analysis models are run on the imputed datasets:
        - `Scripts\analyses\2l.pmm analyses.R`': Script for running the substantive model for the 2l.pmm imputation method --- i.e. the multilevel preditive mean matching imputation method.
        - `Scripts\analyses\bart analyses.R`: Script for running the substantive model for the bart imputation method --- i.e. the single-level BART imputation method.
        - `Scripts\analyses\Complete analyses.R`: Script for running the substantive model for the complete/true data --- i.e. the data without missing values.
        - `Scripts\analyses\Listwise deletion analyses.R`: Script for running the substantive model with listwise deletion.
        - `Scripts\analyses\pmm analyses.R`: Script for running the substantive model for the pmm imputation method --- i.e. the single-level preditive mean matching imputation method.
        - `Scripts\analyses\rbart analyses.R`: Script for running the substantive model for the rbart imputation method --- i.e. the multilevel BART imputation method with random intercepts.
        - `Scripts\analyses\stan4bart analyses.R`: Script for running the substantive model for the stan4bart imputation method --- i.e. the multilevel BART imputation method with random intercepts, random slopes, and cross-level interactions.
    - Then, the bias, coverage, and confidence interval width are calculated from the substantive analysis results:
        - `Scripts\analyses\Results 2l.pmm analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the 2l.pmm imputation method.
        - `Scripts\analyses\Results bart analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the BART imputation method.
        - `Scripts\analyses\Results complete analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the complete/true data.
        - `Scripts\analyses\Results listwise deletion analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the listwise deletion imputation method.
        - `Scripts\analyses\Results pmm analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the pmm imputation method.
        - `Scripts\analyses\Results rbart analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the rbart imputation method.
        - `Scripts\analyses\Results stan4bart analyses.R`: Script for calculating the bias, coverage, and confidence interval width for the stan4bart imputation method.
5. *`Scripts\Formatting for visualizations.R`: Formats the results appropriately for the visualizations and saves them in appropriate folders.*[^1]
6. `Scripts\Visualizations.R`: Creates the visualizations and saves them in appropriate folders. This script also sources the script `Scripts\Formatting for visualizations.R` for an easier running procedure.

[^1]: Note that this script is automatically sourced in the script `Scripts\Visualizations.R`. Therefore, it does not need to be run separately.

# Ethics/Privacy/Security

The FETC reference number for this study is: 23-1778. The data is simulated and handled in accordance with faculty protocol and access is limited to the student, supervisor(s) and related UU researchers.

# Permission and access

The archive is accessible on github at the following [link](https://github.com/heleenbrueggen/masterthesis/). I, Heleen Brüggen, am solely responsible for the content of the archive, which is publicly available.

# Manuscript

The manuscript is written in LateX and can be found in the folder `thesis`. This folder also contains the BibteX file for the bibliography. The plots that are used in the manuscript are in the folder `thesis\graphs`.
