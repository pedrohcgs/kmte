#' Kaplan-Meier Distributional Treatment Effect
#'
#' \emph{kmdte} computes the Distributional Treatment Effect for possibly right-censored
#' outcomes. The estimator relies on the unconfoundedness assumption, and on
#' estimating the propensity score. For details of the estimation procedure, see
#' Sant'Anna (2016a), 'Program Evaluation with Right-Censored Data'.
#'
#'
#'@param out vector containing the outcome of interest
#'@param delta vector containing the censoring indicator (1 if observed, 0 if censored)
#'@param treat vector containing the treatment indicator (1 if treated, 0 if control)
#'@param xpscore matrix (or data frame) containing the covariates (and their
#'               transformations) to be included in the propensity score estimation.
#'               Propensity score estimation is based on Logit.
#'@param ysup   scalar or vector of points for which
#'               the distributional treatment effect is computed. If NULL,
#'               all uncensored data points available are used.
#'@param b 	The number of bootstrap replicates to be performed. Default is 1,000.
#'@param ci A scalar or vector with values in (0,1) containing the confidence level(s)
#'          of the required interval(s). Default is a vector with
#'          0,90, 0.95 and 0.99
#'@param standardize Default is TRUE, which normalizes propensity score weights to
#'                   sum to 1 within each treatment group.
#'                    Set to FALSE to return Horvitz-Thompson weights.
#'@param cores number of processesors to be used during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using snow
#'
#'@return a list containing the distributional treatment effect estimate, dte,
#'        and the bootstrapped \emph{ci} confidence
#'        confidence interval, l.dte (lower bound), and u.dte (upper bound).
#'@export
#'@importFrom stats glm quantile approxfun
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#'@importFrom Rearrangement rearrangement
#-----------------------------------------------------------------------------
kmdte <- function(out, delta, treat, ysup = NULL,
                  xpscore, b = 1000, ci = c(0.90,0.95,0.99),
                  standardize = TRUE, cores = 1) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, xpscore))
  #-----------------------------------------------------------------------------
  # set up all unique uncensored points in the data (use it as ysup, if NULL)
  yy <- sort(unique(fulldata[fulldata$delta == 1,]$out))
  if (is.null(ysup) == TRUE){
    ysup <- yy
  }
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmdte <- function(fulldata, i, ysup1 = ysup,
                          standardize1 = standardize){
    #----------------------------------------------------------------------------
    # Select the data for the bootstrap (like the original data)
    df.b=fulldata[i,]
    #----------------------------------------------------------------------------
    # Compute Kaplan-Meier Weigths - data is now sorted!
    df.b <- kmweight(1, 2, df.b)
    # Dimension of data matrix df.b
    dim.b <- dim(df.b)[2]
    # Next, we rename the variable in xpscore to avoid problems
    xpscore1.b <- df.b[, (4:(dim.b - 1))]
    datascore.b <- data.frame(y = df.b[, 3], xpscore1.b)
    #-----------------------------------------------------------------------------
    # estimate the propensity score
    pscore.b <- stats::glm(y ~ ., data = datascore.b,
                           family = binomial("logit"))
    df.b$pscore <- pscore.b$fit
    #-----------------------------------------------------------------------------
    # Create id to help on ordering
    df.b$id <- 1:length(df.b[, 1])
    # Update Dimension of data matrix fulldata
    dim.all <- dim(df.b)[2]
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
    # correct KM weigths
    # First, the datasets
    df.b <- data.frame(rbind(data.treat.b, data.control.b))
    # Sort wrt id
    df.b <- df.b[order(as.numeric(df.b[, dim.all])), ]
    #-----------------------------------------------------------------------------
    # Compute weigths for treatment and control groups
    w1km.b <- ((df.b$treat * df.b$w) / df.b$pscore)
    w0km.b <- ((1 - df.b$treat) * df.b$w / (1 - df.b$pscore))

    if (standardize1 == TRUE) {
      w1km.b <- w1km.b / mean(df.b$treat / df.b$pscore)
      w0km.b <- w0km.b / mean((1 - df.b$treat) / (1 - df.b$pscore))
    }
    #-----------------------------------------------------------------------------
    # Compute Counterfactual distributions, and the DTE
    # First, we KM estimates of the potential outcomes distribution
    kmcdf.y1 <- w.ecdf(df.b$out, w1km.b)
    kmcdf.y0 <- w.ecdf(df.b$out, w0km.b)

    # Next, we rearrange these distributions
    #kmcdf.y1.r <- Rearrangement::rearrangement(data.frame(df.b$out),
    #                                           kmcdf.y1(df.b$out))
    #kmcdf.y1.r[kmcdf.y1.r > 1] <- 1
    #kmcdf.y1.r[kmcdf.y1.r < 0] <- 0
    #kmcdf.y1.r <- r.ecdf(ysup1, kmcdf.y1.r)

    #kmcdf.y0.r <- Rearrangement::rearrangement(data.frame(df.b$out),
    #                                           kmcdf.y0(df.b$out))
    #kmcdf.y0.r[kmcdf.y0.r > 1] <- 1
    #kmcdf.y0.r[kmcdf.y0.r<0] <- 0
    #kmcdf.y0.r <- r.ecdf(ysup1, kmcdf.y0.r)

    #Distribution of y1 and y0, and dte
    #cdfy1 <- kmcdf.y1.r(ysup1)
    #cdfy0 <- kmcdf.y0.r(ysup1)
    cdfy1 <- kmcdf.y1(ysup1)
    cdfy0 <- kmcdf.y0(ysup1)

    dte <- cdfy1 - cdfy0
    #-----------------------------------------------------------------------------
    return(cbind(cdfy1, cdfy0, dte))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmdte <- boot::boot(fulldata, boot1.kmdte, R = nboot,
                             stype = "i", sim = "ordinary")
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmdte <- boot::boot(fulldata, boot1.kmdte, R = nboot, parallel = "snow",
                             ncpus = cores, stype = "i", sim = "ordinary")
    parallel::stopCluster(cl)
  }
  #----------------------------------------------------------------------------
  # Compute Counterfactual distributions and the DTE
  cdfy1 <- matrix(boot.kmdte$t0[,1],1,length(ysup))
  rownames(cdfy1) <- "CDF Y(1)"
  colnames(cdfy1) <- ysup

  cdfy0 <- matrix(boot.kmdte$t0[,2],1,length(ysup))
  rownames(cdfy0) <- "CDF Y(0)"
  colnames(cdfy0) <- ysup

  dte <- matrix(boot.kmdte$t0[,3],1,length(ysup))
  rownames(dte) <- "DTE"
  colnames(dte) <- ysup
  #----------------------------------------------------------------------------
  #Compute the confidence interval for dte
  n.ysup <- length(ysup)
  n.ci <- length(ci)

  if (n.ci == 1 & n.ysup == 1){
    dte.lb <- boot::boot.ci(boot.kmdte, type="perc", index = 3, conf = ci)$percent[4]
    dte.ub <- boot::boot.ci(boot.kmdte, type="perc", index = 3, conf = ci)$percent[5]
  }

  if (n.ci >1 & n.ysup == 1){
    dte.lb <- boot::boot.ci(boot.kmdte, type="perc", index = 3, conf = ci)$percent[,4]
    dte.ub <- boot::boot.ci(boot.kmdte, type="perc", index = 3, conf = ci)$percent[,5]
  }

  if ((n.ci == 1) * (n.ysup > 1) == 1){
    dte.lb <- matrix(NA, n.ci, n.ysup)
    dte.ub <- matrix(NA, n.ci, n.ysup)
    for (i in 1:n.ysup){
      dte.lb[,i] <- boot::boot.ci(boot.kmdte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[4]
      dte.ub[,i] <- boot::boot.ci(boot.kmdte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[5]
    }
  }

  if ((n.ci > 1) * (n.ysup > 1) == 1){
    dte.lb <- matrix(NA, n.ci, n.ysup)
    dte.ub <- matrix(NA, n.ci, n.ysup)
    for (i in 1:n.ysup){
      dte.lb[,i] <- boot::boot.ci(boot.kmdte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[,4]
      dte.ub[,i] <- boot::boot.ci(boot.kmdte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[,5]
    }
  }
  #----------------------------------------------------------------------------
  colnames(dte.lb) <- ysup
  colnames(dte.ub) <- ysup
  rownames(dte.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(dte.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  #----------------------------------------------------------------------------
  # Return these
  list(dte = dte,
       cdfy1 = cdfy1,
       cdfy0 = cdfy0,
       #boot = boot.kmdte,
       dte.lb = dte.lb,
       dte.ub = dte.ub
  )
}
