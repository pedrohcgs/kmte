# kmte: An R Package for Treatment Effects with Censored Outcomes

As of now, this n R package implements all the tests proposed in Sant'Anna (2016b), "Nonparametric Tests for Treatment Effect Heterogeneity with Censored data", available at Pedro H.C. Sant'Anna webpage, https://sites.google.com/site/pedrohcsantanna/ .

In short, this package implements different Kolmogorov-Smirnov and Cramér–von Mises type tests. 

When the treatment is exogenous, i.e. under the unconfoundedness assumption, we consider the following null hypotheses:
    
    a. Zero conditional distribution treatment effects;
    
    b. Zero conditional average treatment effects;
    
    c. Homogeneous conditional average treatment effects.
    
When the treatment is endogenous and we have a binary instrument, i.e. under the local treatment effects setup, we consider the 
following null hypotheses:

    a. Zero conditional local distribution treatment effects;
    
    b. Zero conditional local average treatment effects;
    
    c. Homogeneous local conditional average treatment effects.

To install this package, you simply need to run the following two lines of code in R:

        library(devtools)
        install_github("pedrohcgs/kmte")

A help file with examples will be made available in the near future. Furthermore, all the treatment effects measures proposed in Sant'Anna (2016a), "Program Evaluation with Right-Censored Data", will also be added.

This package was written by Pedro H. C. Sant'Anna (Vanderbilt University)


