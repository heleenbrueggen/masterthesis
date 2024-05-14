# MSc Thesis: Multilevel Multivariate Imputation by Chained Equations through Bayesian Additive Regression Trees
 
## Design
This study is a simulation study to compare the performance of different imputation methods in the context of multilevel data. The imputation methods that are compared are: single level predictive mean matching (PMM); multilevel PMM; single level Bayesian Additive Regression Trees (BART); multilevel BART with random intercepts; and multilevel BART with random intercepts, random slopes, and cross-level interactions.

## Data Generation
This study uses simulated data. The data is generated with the script `Scripts\Simulation data.R`. This script is reproducible. However, the data can also be found in the folder `Data`. The following table shows the design factors used in generating the data. 

| Design factors             | Values    | 
|----------------------------|:----------|
| Number of clusters         | 30, 50    | 
| Within-cluster sample size | 15, 50    | 
| Intraclass Correlation     | .5        |
| Missing data mechanism     | MCAR, MAR |
| Amount of missingness      | 0%, 50%   |

The script generates a 1000 datasets for each combination of design factors. However, due to time constraints, the study used only 100 datasets for each combination of design factors for every imputation method, except the imputation method for multilevel BART with random intercepts, random slopes, and cross-level interactions. The next table shows the amount of datasets used for each imputation method. 

|                                                                                     |      | Number of groups: 50, Group size: 15 | Number of groups: 30, Group size: 50 | Number of groups: 50, Group size: 15 | Number of groups: 50, Group size: 50 |
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
