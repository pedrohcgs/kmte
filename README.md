# kmte: An R Package for Treatment Effects with Censored Outcomes

## Description 
This R package includes a variety of policy evaluations when the outcome of interest, typically a duration, is subjected to right censoring. The content includes estimators and tests related to average, quantile and distributional treatment effects under different identifying assumptions including unconfoundedness, local treatment effects, and nonlinear difrences-in-differences. 

In short, the R package implements all estimators proposed in Sant'Anna (2016a) "Program Evaluation with Right-Censored Data", and all tests proposed in Sant'Anna (2016b), "Nonparametric Tests for Treatment Effect Heterogeneity with Duration Outcomes", available at Pedro H.C. Sant'Anna webpage, https://sites.google.com/site/pedrohcsantanna/ .

## Functions
### Exogenous Treatment Allocation (Unconfoundedness setup)
When the treatment is exogenous, i.e. under the unconfoundedness assumption, we estimate the following treatment effect parameters:
* `kmate` - compute the Average Treatment Effect for possibly randomly censored data.
* `kmqte` - compute the Quantile Treatment Effect for possibly randomly censored data.
* `kmdte` - compute the Distributional Treatment Effect for possibly randomly censored data.

In addition, the `kmte` package also implements different nonparametric tests for treatment effects heterogeneity under the unconfoundedness assumption:
* `zcate` - compute different tests for the null hypothesis of Zero Conditional Average Treatment Effetcs.
* `zcdte` - compute different tests for the null hypothesis of Zero Conditional Distributional Treatment Effetcs.
* `hcate` - compute different tests for the null hypothesis of Homogeneous Conditional Average Treatment Effetcs.

### Endogenous Treatment Allocation (Local Treatment Effects setup)    
When the treatment is endogenous and we have a binary instrument, i.e. under the local treatment effects setup, we estimate the 
following treatment effect parameters:
* `kmlate` - compute the Local Average Treatment Effect for possibly randomly censored data.
* `kmlqte` - compute the Local Quantile Treatment Effect for possibly randomly censored data.
* `kmldte` - compute the Local Distributional Treatment Effect for possibly randomly censored data.

Similar to the unconfoundedness case, the following functions implement different nonparametric tests for treatment effects heterogeneity under the local treatment effect setup:
* `zclate` - compute different tests for the null hypothesis of Zero Conditional Local Average Treatment Effetcs.
* `zcldte` - compute different tests for the null hypothesis of Zero Conditional Local Distributional Treatment Effetcs.
* `hclate` - compute different tests for the null hypothesis of Homogeneous Conditional Local Average Treatment Effetcs.

## Installing kmte
This github website hosts the source code. To install the `kmte` package, you simply need to run the following two lines of code in R:

        library(devtools)
        install_github("pedrohcgs/kmte")

A help file with examples will be made available in the near future. Furthermore, the Nonlinear Differences-in-Differences estimator proposed in Sant'Anna (2016a), "Program Evaluation with Right-Censored Data", will also be added.

This package was written by Pedro H. C. Sant'Anna (Vanderbilt University)
