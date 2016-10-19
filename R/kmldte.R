#' Kaplan-Meier Local Distributional Treatment Effect
#'
#' \emph{kmldte} computes the Local Distributional Treatment Effect for possibly right-censored outcomes.
#' The estimator relies on the availability of an Instrumental variable Z, and on a monotonicity assumption.
#' To implement the estimator, we make use of an instrumental propensity score approach.
#' For details of the estimation procedure, see Sant'Anna (2016a), 'Program Evaluation with
#' Right-Censored Data'.
#'
#'
#'@param out vector containing the outcome of interest
#'@param delta vector containing the censoring indicator (1 if observed, 0 if censored)
#'@param treat vector containing the treatment indicator (1 if treated, 0 if control)
#'@param z  vector containing the binary instrument
#'@param xpscore matrix (or data frame) containing the covariates (and their
#'               transformations) to be included in the instrument propensity score estimation.
#'               Instrument Propensity score estimation is based on Logit.
#'@param ysup   scalar or vector of points for which
#'               the distributional treatment effect is computed. If NULL,
#'               all uncensored data points available are used.
#'@param b 	The number of bootstrap replicates to be performed. Default is 1,000.
#'@param ci A scalar or vector with values in (0,1) containing the confidence level(s)
#'          of the required interval(s). Default is a vector with
#'          0,90, 0.95 and 0.99
#'@param cores number of processesors to be used during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using snow.
#'@param monot Default is TRUE, which impose that the estimated counterfactual distributions are in
#'              proper CDF's, i.e. takes values between [0,1], and are non-decreasing.
#'              Boundedness is imposed by truncantion, and monotonicity is imposed using the
#'              rearrangement procedure proposed by Chernozhukov, Fernandez-Val, and Galichon (2010),
#'              implemented in R through package Rearrangement. If FALSE, no adjustment is made.
#'
#'@return a list containing the local distributional treatment effect estimate, ldte,
#'        and the bootstrapped \emph{ci} confidence
#'        confidence interval, ldte.lb (lower bound), and ldte.ub (upper bound).
#'@export
#'@importFrom stats glm quantile approxfun
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#'@importFrom Rearrangement rearrangement
#-----------------------------------------------------------------------------
kmldte <- function(out, delta, treat, z, xpscore, ysup = NULL, b = 1000,
                   ci = c(0.90,0.95,0.99),
                   cores = 1, monot = TRUE) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, z, xpscore))
  #-----------------------------------------------------------------------------
  # set up all unique uncensored points in the data (use it as ysup, if NULL)
  yy <- sort(unique(fulldata[fulldata$delta == 1,]$out))
  if (is.null(ysup) == TRUE){
    ysup <- yy
  }
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmldte <- function(fulldata, i, ysup1 = ysup, monot1 = monot){
    #----------------------------------------------------------------------------
    # Select the data for the bootstrap (like the original data)
    df.b=fulldata[i,]
    #----------------------------------------------------------------------------
    # Dimension of data matrix df.b
    dim.b <- dim(df.b)[2]
    # Next, we rename the variable in xpscore to avoid problems
    xpscore1.b <- df.b[, (5:dim.b)]
    datascore.b <- data.frame(y = df.b[, 4], xpscore1.b)
    #-----------------------------------------------------------------------------
    # estimate the propensity score
    pscore.b <- stats::glm(y ~ ., data = datascore.b,
                           family = binomial("logit"))
    df.b$pscore <- pscore.b$fit
    #-----------------------------------------------------------------------------
    # sample size
    n.total.b <- as.numeric(length(df.b[, 1]))
    # subset of treated individuals with instrument Z equal to 1
    data.treat.z1.b <- subset(df.b, df.b[, 3] == 1 & df.b[, 4] == 1)
    # subset of treated individuals with instrument Z equal to 0
    data.treat.z0.b <- subset(df.b, df.b[, 3] == 1 & df.b[, 4] == 0)
    # subset of control individuals with instrument Z equal to 1
    data.control.z1.b <- subset(df.b, df.b[, 3] == 0 & df.b[, 4] == 1)
    # subset of control individuals with instrument Z equal to 0
    data.control.z0.b <- subset(df.b, df.b[, 3] == 0 & df.b[, 4] == 0)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for treated with z=1
    data.treat.z1.b <- kmweight(1, 2, data.treat.z1.b)
    n.treat.z1.b <- as.numeric(length(data.treat.z1.b[, 1]))
    data.treat.z1.b$w <- data.treat.z1.b$w * (n.treat.z1.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for treated with z=0
    data.treat.z0.b <- kmweight(1, 2, data.treat.z0.b)
    n.treat.z0.b <- as.numeric(length(data.treat.z0.b[, 1]))
    data.treat.z0.b$w <- data.treat.z0.b$w * (n.treat.z0.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for control with z=1
    data.control.z1.b <- kmweight(1, 2, data.control.z1.b)
    n.control.z1.b <- as.numeric(length(data.control.z1.b[, 1]))
    data.control.z1.b$w <- data.control.z1.b$w * (n.control.z1.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Compute Kaplan-Meier weigth for control with z=0
    data.control.z0.b <- kmweight(1, 2, data.control.z0.b)
    n.control.z0.b <- as.numeric(length(data.control.z0.b[, 1]))
    data.control.z0.b$w <- data.control.z0.b$w * (n.control.z0.b/n.total.b)
    #-----------------------------------------------------------------------------
    # Let's put everything in a single data
    df.b <- data.frame(rbind(data.treat.z1.b, data.treat.z0.b,
                             data.control.z1.b, data.control.z0.b))
    # sort wrt out and if draw, by delta
    df.b <- df.b[order(df.b$out, df.b$delta),]
    #-----------------------------------------------------------------------------
    # Compute weigths for treatment and control groups
    w11km.b <- ((df.b$treat * df.b$z * df.b$w) / df.b$pscore)
    w10km.b <- ((df.b$treat * (1 - df.b$z) * df.b$w) / (1 - df.b$pscore))

    w01km.b <- ((1 - df.b$treat) * df.b$z * df.b$w / df.b$pscore)
    w00km.b <- ((1 - df.b$treat) * (1 - df.b$z) * df.b$w / (1 - df.b$pscore))

    kappa11 <- mean(df.b$treat * df.b$z / df.b$pscore)
    kappa10 <- mean(df.b$treat * (1-df.b$z) / (1 - df.b$pscore))
    kappa01 <- mean((1 - df.b$treat) * df.b$z / df.b$pscore)
    kappa00 <- mean((1 - df.b$treat) * (1-df.b$z) / (1 - df.b$pscore))

    w1km.b <- w11km.b - w10km.b
    w0km.b <- w01km.b - w00km.b
    kappa1 <- kappa11 - kappa10
    kappa0 <- kappa01 - kappa00

    w1km.b <- w1km.b / kappa1
    w0km.b <- w0km.b / kappa0
    #-----------------------------------------------------------------------------
    # Compute Counterfactual quantiles, and the dte
    # First, we KM estimates of the potential outcomes distribution
    kmcdf.y1 <- w.ecdf(df.b$out, w1km.b)
    kmcdf.y0 <- w.ecdf(df.b$out, w0km.b)

    # Warning message if negative weigths
    if (monot1 == FALSE & min(min(w1km.b),min(w0km.b)) < 0){
      warning(" Some of the weights used in computing the Distributions are negative. Consider setting monot == TRUE")
    }

    # Next, we rearrange these distributions
    if (monot1 == TRUE) {
      # get all unique uncensored data points out
      yy <- sort(unique(fulldata[fulldata$delta == 1,]$out))

      #rearrange CDF Y(1) for compliers
      kmcdf.y1 <- Rearrangement::rearrangement(data.frame(yy), kmcdf.y1(yy))
      kmcdf.y1[kmcdf.y1 > 1] <- 1
      kmcdf.y1[kmcdf.y1 < 0] <- 0
      kmcdf.y1 <- r.ecdf(yy, kmcdf.y1)

      #rearrange CDF Y(10) for compliers
      kmcdf.y0 <- Rearrangement::rearrangement(data.frame(yy), kmcdf.y0(yy))
      kmcdf.y0[kmcdf.y0 > 1] <- 1
      kmcdf.y0[kmcdf.y0 < 0] <- 0
      kmcdf.y0 <- r.ecdf(yy, kmcdf.y0)
    }

    #Distribution of y1 and y0, and dte
    cdfy1.c <- kmcdf.y1(ysup1)
    cdfy0.c <- kmcdf.y0(ysup1)

    ldte <- cdfy1.c - cdfy0.c

    #-----------------------------------------------------------------------------
    return(cbind(cdfy1.c, cdfy0.c, ldte))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmldte <- boot::boot(fulldata, boot1.kmldte, R = nboot,
                             stype = "i", sim = "ordinary")
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmldte <- boot::boot(fulldata, boot1.kmldte, R = nboot, parallel = "snow",
                             ncpus = cores, stype = "i", sim = "ordinary")
    parallel::stopCluster(cl)
  }
  #----------------------------------------------------------------------------
  # Compute Counterfactual distributions and the ldte
  cdfy1.c <- matrix(boot.kmldte$t0[,1],1,length(ysup))
  rownames(cdfy1.c) <- "CDF Y(1) for Compliers"
  colnames(cdfy1.c) <- ysup

  cdfy0.c <- matrix(boot.kmldte$t0[,2],1,length(ysup))
  rownames(cdfy0.c) <- "CDF Y(0) for Compliers"
  colnames(cdfy0.c) <- ysup

  ldte <- matrix(boot.kmldte$t0[,3],1,length(ysup))
  rownames(ldte) <- "LDTE"
  colnames(ldte) <- ysup
  #----------------------------------------------------------------------------
  #Compute the confidence interval for ldte
  n.ysup <- length(ysup)
  n.ci <- length(ci)

  if (n.ci == 1 & n.ysup == 1){
    ldte.lb <- boot::boot.ci(boot.kmldte, type="perc", index = 3, conf = ci)$percent[4]
    ldte.ub <- boot::boot.ci(boot.kmldte, type="perc", index = 3, conf = ci)$percent[5]
  }

  if (n.ci >1 & n.ysup == 1){
    ldte.lb <- boot::boot.ci(boot.kmldte, type="perc", index = 3, conf = ci)$percent[,4]
    ldte.ub <- boot::boot.ci(boot.kmldte, type="perc", index = 3, conf = ci)$percent[,5]
  }

  if ((n.ci == 1) * (n.ysup > 1) == 1){
    ldte.lb <- matrix(NA, n.ci, n.ysup)
    ldte.ub <- matrix(NA, n.ci, n.ysup)
    for (i in 1:n.ysup){
      ldte.lb[,i] <- boot::boot.ci(boot.kmldte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[4]
      ldte.ub[,i] <- boot::boot.ci(boot.kmldte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[5]
    }
  }

  if ((n.ci > 1) * (n.ysup > 1) == 1){
    ldte.lb <- matrix(NA, n.ci, n.ysup)
    ldte.ub <- matrix(NA, n.ci, n.ysup)
    for (i in 1:n.ysup){
      ldte.lb[,i] <- boot::boot.ci(boot.kmldte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[,4]
      ldte.ub[,i] <- boot::boot.ci(boot.kmldte, type="perc",
                                  index = (i+ 2* n.ysup), conf = ci)$percent[,5]
    }
  }
  #----------------------------------------------------------------------------
  colnames(ldte.lb) <- ysup
  colnames(ldte.ub) <- ysup
  rownames(ldte.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(ldte.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  #----------------------------------------------------------------------------
  # Return these
  list(ldte = ldte,
       cdfy1.c = cdfy1.c,
       cdfy0.c = cdfy0.c,
       #boot = boot.kmldte,
       ldte.lb = ldte.lb,
       ldte.ub = ldte.ub
  )
}
