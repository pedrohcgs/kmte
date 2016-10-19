#' Kaplan-Meier Local Quantile Treatment Effect
#'
#' \emph{kmlqte} computes the Local Quantile Treatment Effect for possibly right-censored outcomes.
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
#'@param probs   scalar or vector of probabilities with values in (0,1) for which
#'               the quantile treatment effect is computed. Default is 0.5, returning
#'               the median.
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
#'@return a list containing the local quantile treatment effect estimate, lqte,
#'        and the bootstrapped \emph{ci} confidence
#'        confidence interval, lqte.lb (lower bound), and lqte.ub (upper bound).
#'@export
#'@importFrom stats glm quantile approxfun
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#'@importFrom Rearrangement rearrangement
#-----------------------------------------------------------------------------
kmlqte <- function(out, delta, treat, z, xpscore, probs = 0.5, b = 1000,
                   ci = c(0.90,0.95,0.99), cores = 1, monot = TRUE) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, z, xpscore))
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmlqte <- function(fulldata, i, probs1 = probs, monot1 = monot){
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
    # Compute Counterfactual quantiles, and the QTE
    # First, we KM estimates of the potential outcomes distribution
    kmcdf.y1 <- w.ecdf(df.b$out, w1km.b)
    kmcdf.y0 <- w.ecdf(df.b$out, w0km.b)

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

    #quantiles of y1 and y0 for compliers, and lqte
    qy1.c <- stats::quantile(kmcdf.y1, type = 1, probs = probs1)
    qy0.c <- stats::quantile(kmcdf.y0, type = 1, probs = probs1)

    lqte <- qy1.c - qy0.c

    #-----------------------------------------------------------------------------
    return(cbind(qy1.c, qy0.c, lqte))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmlqte <- boot::boot(fulldata, boot1.kmlqte, R = nboot,
                             stype = "i", sim = "ordinary")
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmlqte <- boot::boot(fulldata, boot1.kmlqte, R = nboot, parallel = "snow",
                             ncpus = cores, stype = "i", sim = "ordinary")
    parallel::stopCluster(cl)
  }
  #----------------------------------------------------------------------------
  # Compute Counterfactual quantiles for compliers and the LQTE
  qy1.c <- matrix(boot.kmlqte$t0[,1],1,length(probs))
  rownames(qy1.c) <- "Local Quantile Y(1)"
  colnames(qy1.c) <- names(quantile(1, probs = probs))

  qy0.c <- matrix(boot.kmlqte$t0[,2],1,length(probs))
  rownames(qy0.c) <- "Local Quantile Y(0)"
  colnames(qy0.c) <- names(quantile(1, probs = probs))

  lqte <- matrix(boot.kmlqte$t0[,3],1,length(probs))
  rownames(lqte) <- "LQTE"
  colnames(lqte) <- names(quantile(1, probs = probs))
  #----------------------------------------------------------------------------
  #Compute the confidence interval for lqte
  n.probs <- length(probs)
  n.ci <- length(ci)

  if (n.ci == 1 & n.probs == 1){
    lqte.lb <- boot::boot.ci(boot.kmlqte, type="perc", index = 3, conf = ci)$percent[4]
    lqte.ub <- boot::boot.ci(boot.kmlqte, type="perc", index = 3, conf = ci)$percent[5]
  }

  if (n.ci >1 & n.probs == 1){
    lqte.lb <- boot::boot.ci(boot.kmlqte, type="perc", index = 3, conf = ci)$percent[,4]
    lqte.ub <- boot::boot.ci(boot.kmlqte, type="perc", index = 3, conf = ci)$percent[,5]
  }

  if ((n.ci == 1) * (n.probs > 1) == 1){
    lqte.lb <- matrix(NA, n.ci, n.probs)
    lqte.ub <- matrix(NA, n.ci, n.probs)
    for (i in 1:n.probs){
      lqte.lb[,i] <- boot::boot.ci(boot.kmlqte, type="perc",
                                  index = (i+ 2* n.probs), conf = ci)$percent[4]
      lqte.ub[,i] <- boot::boot.ci(boot.kmlqte, type="perc",
                                  index = (i+ 2* n.probs), conf = ci)$percent[5]
    }
  }

  if ((n.ci > 1) * (n.probs > 1) == 1){
    lqte.lb <- matrix(NA, n.ci, n.probs)
    lqte.ub <- matrix(NA, n.ci, n.probs)
    for (i in 1:n.probs){
      lqte.lb[,i] <- boot::boot.ci(boot.kmlqte, type="perc",
                                  index = (i+ 2* n.probs), conf = ci)$percent[,4]
      lqte.ub[,i] <- boot::boot.ci(boot.kmlqte, type="perc",
                                  index = (i+ 2* n.probs), conf = ci)$percent[,5]
    }
  }
  #----------------------------------------------------------------------------
  colnames(lqte.lb) <- paste(names(quantile(1, probs = probs)), "quantile")
  colnames(lqte.ub) <- paste(names(quantile(1, probs = probs)), "quantile")
  rownames(lqte.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(lqte.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  #----------------------------------------------------------------------------
  # Return these
  list(lqte = lqte,
       qy1.c = qy1.c,
       qy0.c = qy0.c,
       #boot = boot.kmlqte,
       lqte.lb = lqte.lb,
       lqte.ub = lqte.ub
  )
}
