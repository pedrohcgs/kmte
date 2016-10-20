#' An R Package for Treatment Effects with Censored Outcomes.
#'
#'The kmte R package includes a variety of policy evaluations tools when the outcome of interest, typically a duration, is subjected to right
#'censoring. The content includes estimators and tests related to average, quantile and distributional treatment effects under different
#'identifying assumptions including unconfoundedness, local treatment effects, and nonlinear difrences-in-differences.
#'
#' In short, kmte implements all estimators proposed in Sant'Anna (2016a) "Program Evaluation with Right-Censored Data", and all tests proposed
#' in Sant'Anna (2016b), "Nonparametric Tests for Treatment Effect Heterogeneity with Duration Outcomes". Both articles are available at
#' Pedro H.C. Sant'Anna webpage, \url{http://sites.google.com/site/pedrohcsantanna/} .
#'
#'When the treatment is exogenous, i.e. under the unconfoundedness assumption,
#'one can use the following functions:
#'\itemize{
#'   \item kmate - compute the Average Treatment Effect for possibly randomly censored data.
#'   \item kmqte - compute the Quantile Treatment Effect for possibly randomly censored data.
#'   \item kmdte - compute the Distributional Treatment Effect for possibly randomly censored data.
#' }
#'
#'Alternatively, when the treatment is endogenous and we have a binary instrument, i.e. under the local treatment effect setup,
#'one can use the following functions:
#'#'\itemize{
#'   \item kmlate - compute the Local Average Treatment Effect for possibly randomly censored data.
#'   \item kmlqte - compute the LocalQuantile Treatment Effect for possibly randomly censored data.
#'   \item kmldte - compute the Local Distributional Treatment Effect for possibly randomly censored data.
#' }
#'
#'
#'In addition, the kmte package also implements the following nonparametric tests for treatment effect heterogeneity:
#'
#'\enumerate{
#'\item Under Unconfoundedness:
#'\itemize{
#'   \item zcate - compute different tests for the null hypothesis of Zero Conditional Average Treatment Effetcs.
#'   \item zcdte - compute different tests for the null hypothesis of Zero Conditional Distributional Treatment Effetcs.
#'   \item hcate - compute different tests for the null hypothesis of Homogeneous Conditional Average Treatment Effetcs.
#'}
#'\item Under the Local Treatment Setup:
#'#'\itemize{
#'   \item zclate - compute different tests for the null hypothesis of Zero Conditional Local Average Treatment Effetcs.
#'   \item zcldte - compute different tests for the null hypothesis of Zero Conditional Local Distributional Treatment Effetcs.
#'   \item hclate - compute different tests for the null hypothesis of Homogeneous Conditional Local Average Treatment Effetcs.
#'}
#'}
#'
#' The module was written by Pedro H.C. Sant'Anna (Vanderbilt University).
#'
#' A help file with examples will be made available in the near future.
#'
"_PACKAGE"
