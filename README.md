# NetworkTargeting
This is an implementation of the methods described in "Policy Targeting under Network Interference," (Viviano, REStud forthcoming).

Given outcomes for each individual, treatments, covariates, and an adjacency matrix with binary entries, this package estimates which individuals should be treated to maximize the expected outcome (welfare).

The package implements an exact optimization procedure and an approximate optimization procedure. The exact procedure is recommended for small $n$, and the approximate procedure for large $n$.

## Installation

```
install.packages("remotes")
library(remotes) 
remotes::install_github("dviviano/NetworkTargeting")
library(NetworkTargeting)
```

## Implementation

This package contains the functions: 

- `NetworkTargeting`
- `NetworkWelfare`
- `PredictPolicies`


See the [vignette](vignettes/vignette_Jan2024.Rmd) for more details on proper usage. Further details on function arguments are available by calling `help` with the name of the function of interest.

## Reference 

Viviano, Davide. "Policy Targeting under Network Interference." The Review of Economic Studies (forthcoming).

## Support 

Davide Viviano: dviviano@fas.harvard.edu

Jacob Carlson: jacob_carlson@g.harvard.edu
