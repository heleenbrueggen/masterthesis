# MSc Thesis: Multilevel Multivariate Imputation by Chained Equations through Bayesian Additive Regression Trees
 
## Design
This study is a simulation study to compare the performance of different imputation methods in the context of multilevel data. The imputation methods that are compared are: single level predictive mean matching (PMM); multilevel PMM; single level Bayesian Additive Regression Trees (BART); multilevel BART with random intercepts; and multilevel BART with random intercepts, random slopes, and cross-level interactions.

## Data Generation
This study uses simulated data. The data is generated with the script `Scripts\Simulation data.R`. This script is reproducible. However, the data can also be found in the folder `Data`.

| Design factors             | Values    | 
|----------------------------|:----------|
| Number of clusters         | 30, 50    | 
| Within-cluster sample size | 15, 50    | 
| Intraclass Correlation     | .5        |
| Missing data mechanism     | MCAR, MAR |
| Amount of missingness      | 0%, 50%   |

: Simulation design