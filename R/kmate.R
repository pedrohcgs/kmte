#' kmate: Two-Step Kaplan-Meier Average Treatment Effect
#'
#' \emph{kmate} computes the Average Treatment Effect for possibly right-censored outcomes.
#' The estimator relies on the unconfoundedness assumption, and on estimating the propensity score.
#' For details of the estimation procedure, see Sant'Anna (2016a), 'Program Evaluation with
#' Right-Censored Data'.
#'
#'
#'@param out vector containing the outcome of interest
#'@param delta vector containing the censoring indicator (1 if observed, 0 if censored)
#'@param treat vector containing the treatment indicator (1 if treated, 0 if control)
#'@param xpscore matrix (or data frame) containing the covariates (and their
#'               transformations) to be included in the propensity score estimation.
#'               Propensity score estimation is based on Logit.
#'@param b 	The number of bootstrap replicates to be performed. Default is 1,000.
#'@param ci A scalar or vector with values in (0,1) containing the confidence level(s)
#'          of the required interval(s). Default is a vector with
#'          0,90, 0.95 and 0.99
#'@param trunc scalar that defined the truncation parameter. Default is NULL, which does not perform any kind of
#'            truncation in the computation of the ATE. When trunc is different than NULL, all outcomes which values greater
#'            than trunc are truncated.
#'@param standardize Default is TRUE, which normalizes propensity score weights to sum to 1 within each treatment group.
#'                    Set to FALSE to return Horvitz-Thompson weights.
#'@param cores number of processesors to be used during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using snow
#'
#'@return a list containing the Average treatment effect estimate, ate,
#'        and the bootstrapped \emph{ci} confidence
#'        confidence interval, ate.lb (lower bound), and ate.ub (upper bound).
#'@export
#'@importFrom stats glm
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#-----------------------------------------------------------------------------
kmate <- function(out, delta, treat, xpscore, b = 1000, ci = c(0.90,0.95,0.99),
                  trunc = NULL, standardize = TRUE, cores = 1) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, xpscore))
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmate <- function(fulldata, i, trunc1 = trunc, standardize1 = standardize){
    #----------------------------------------------------------------------------
    # Select the data for the bootstrap (like the original data)
    df.b=fulldata[i,]
    #----------------------------------------------------------------------------
    # Dimension of data matrix df.b
    dim.b <- dim(df.b)[2]
    # Next, we rename the variable in xpscore to avoid problems
    xpscore1.b <- df.b[, (4:dim.b)]
    datascore.b <- data.frame(y = df.b[, 3], xpscore1.b)
    #-----------------------------------------------------------------------------
    # estimate the propensity score
    pscore.b <- stats::glm(y ~ ., data = datascore.b,
                           family = binomial("logit"))
    df.b$pscore <- pscore.b$fit
    #-----------------------------------------------------------------------------
    # sample size
    n.total.b <- as.numeric(length(df.b[, 1]))
    # subset of treated individuals
    data.treat.b <- subset(df.b, df.b[, 3] == 1)
    # subset of not-treated individuals
    data.control.b <- subset(df.b, df.b[, 3] == 0)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for treated
    data.treat.b <- kmweight(1, 2, data.treat.b)
    n.treat.b <- as.numeric(length(data.treat.b[, 1]))
    data.treat.b$w <- data.treat.b$w * (n.treat.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for control
    data.control.b <- kmweight(1, 2, data.control.b)
    n.control.b <- as.numeric(length(data.control.b[, 1]))
    data.control.b$w <- data.control.b$w * (n.control.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Let's put everything in a single data
    df.b <- data.frame(rbind(data.treat.b, data.control.b))
    #-----------------------------------------------------------------------------
    # Compute weigths for treatment and control groups
    w1km.b <- ((df.b$treat * df.b$w) / df.b$pscore)
    w0km.b <- ((1 - df.b$treat) * df.b$w / (1 - df.b$pscore))

    if (standardize1 == TRUE) {
      w1km.b <- w1km.b / mean(df.b$treat / df.b$pscore)
      w0km.b <- w0km.b / mean((1 - df.b$treat) / (1 - df.b$pscore))
    }
    #-----------------------------------------------------------------------------
    # Compute Counterfactual Average Outcomes, E[Y(1)] and E[Y(0)], and the ATE
    meany1km <- sum(w1km.b * df.b$out)
    meany0km <- sum(w0km.b * df.b$out)

    if (is.null(trunc1) == FALSE){
      meany1km <- sum(w1km.b * df.b$out * (df.b$out <= trunc1))
      meany0km <- sum(w0km.b * df.b$out * (df.b$out <= trunc1))
    }

    ate <- meany1km - meany0km
    #-----------------------------------------------------------------------------
    return(cbind(meany1km, meany0km, ate))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmate <- boot::boot(fulldata, boot1.kmate, R = nboot)
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmate <- boot::boot(fulldata, boot1.kmate, R = nboot, parallel = "snow", ncpus = cores)
    parallel::stopCluster(cl)
  }

  #----------------------------------------------------------------------------
  # Compute Counterfactual Average Outcomes, E[Y(1)] and E[Y(0)], and the ATE
  meany1km <- boot.ci(boot.kmate, type="perc", index=1)$t0
  names(meany1km) <- "E[Y(1)]"
  meany0km <- boot.ci(boot.kmate, type="perc", index=2)$t0
  names(meany0km) <- "E[Y(0)]"
  ate <- boot.ci(boot.kmate, type="perc", index=3)$t0
  names(ate) <- "ATE"
  #----------------------------------------------------------------------------
  #Compute the confidence interval for ate
  if (length(ci) == 1){
    ate.lb <- boot.ci(boot.kmate, type="perc", index=3, conf = ci)$percent[4]
    ate.ub <- boot.ci(boot.kmate, type="perc", index=3, conf = ci)$percent[5]
  }
  if (length(ci) >1){
  ate.lb <- boot.ci(boot.kmate, type="perc", index=3, conf = ci)$percent[,4]
  ate.ub <- boot.ci(boot.kmate, type="perc", index=3, conf = ci)$percent[,5]
  }

  ate.lb <- matrix(ate.lb,length(ci),1)
  ate.ub <- matrix(ate.ub,length(ci),1)
  rownames(ate.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(ate.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  colnames(ate.ub) <- "ATE"
  colnames(ate.lb) <- "ATE"
  #----------------------------------------------------------------------------
  # Return these
  list(ate = ate,
       meany1 = meany1km,
       meany0 = meany0km,
       #boot = boot.kmate,
       ate.lb = ate.lb,
       ate.ub = ate.ub)
}
