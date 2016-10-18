#' kmqte: Two-Step Kaplan-Meier Quantile Treatment Effect
#'
#' \emph{kmqte} computes the Quantile Treatment Effect for possibly right-censored
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
#'@param probs   scalar or vector of probabilities with values in (0,1) for which
#'               the quantile treatment effect is computed. Default is 0.5, returning
#'               the median.
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
#'@return a list containing the quantile treatment effect estimate, qte,
#'        its associated bootstrapped values, b.ate, and the \emph{ci} confidence
#'        confidence interval, l.ate (lower bound), and u.ate (upper bound).
#'@export
#'@importFrom stats glm quantile approxfun
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#'@importFrom Rearrangement rearrangement
#-----------------------------------------------------------------------------
kmqte <- function(out, delta, treat, probs = 0.5,
                  xpscore, b = 1000, ci = c(0.90,0.95,0.99),
                  standardize = TRUE, cores = 1) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, xpscore))
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmqte <- function(fulldata, i, probs1 = probs,
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
    # Compute Counterfactual quantiles, and the QTE
    # First, we KM estimates of the potential outcomes distribution
    kmcdf.y1 <- w.ecdf(df.b$out, w1km.b)
    kmcdf.y0 <- w.ecdf(df.b$out, w0km.b)

    # Next, we rearrange these distributions
    kmcdf.y1.r <- Rearrangement::rearrangement(data.frame(df.b$out),
                                               kmcdf.y1(df.b$out))
    kmcdf.y1.r[kmcdf.y1.r > 1] <- 1
    kmcdf.y1.r[kmcdf.y1.r < 0] <- 0
    kmcdf.y1.r <- r.ecdf(df.b$out, kmcdf.y1.r)

    kmcdf.y0.r <- Rearrangement::rearrangement(data.frame(df.b$out),
                                               kmcdf.y0(df.b$out))
    kmcdf.y0.r[kmcdf.y0.r > 1] <- 1
    kmcdf.y0.r[kmcdf.y0.r<0] <- 0
    kmcdf.y0.r <- r.ecdf(df.b$out, kmcdf.y0.r)

    #quantiles of y1 and y0, and qte
    qy1 <- stats::quantile(kmcdf.y1.r, type = 1, probs = probs1)
    qy0 <- stats::quantile(kmcdf.y0.r, type = 1, probs = probs1)

    qte <- qy1 - qy0
    #-----------------------------------------------------------------------------
    return(cbind(qy1, qy0, qte))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmqte <- boot::boot(fulldata, boot1.kmqte, R = nboot,
                             stype = "i", sim = "ordinary")
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmqte <- boot::boot(fulldata, boot1.kmqte, R = nboot, parallel = "snow",
                             ncpus = cores, stype = "i", sim = "ordinary")
    parallel::stopCluster(cl)
  }
  #----------------------------------------------------------------------------
  # Compute Counterfactual quantiles and the QTE
  qy1 <- boot.kmqte$t0[,1]
  #names(qy1) <- paste("Counterfactual ", probs, "-quantile for treated", sep="")
  qy0 <- boot.kmqte$t0[,2]
  #names(qy0) <- paste("Counterfactual", probs, "quantile for control", sep="")
  qte <- boot.kmqte$t0[,3]
 # names(qte) <- paste(probs, "-quantile treatment effect", sep="")
  #----------------------------------------------------------------------------
  #Compute the confidence interval for qte
  n.probs <- length(probs)
  n.ci <- length(ci)

  if (n.ci == 1 & n.probs == 1){
    qte.lb <- boot::boot.ci(boot.kmqte, type="perc", index = 3, conf = ci)$percent[4]
    qte.ub <- boot::boot.ci(boot.kmqte, type="perc", index = 3, conf = ci)$percent[5]
  }

  if (n.ci >1 & n.probs == 1){
    qte.lb <- boot::boot.ci(boot.kmqte, type="perc", index = 3, conf = ci)$percent[,4]
    qte.ub <- boot::boot.ci(boot.kmqte, type="perc", index = 3, conf = ci)$percent[,5]
  }

  if ((n.ci == 1) * (n.probs > 1) == 1){
    qte.lb <- matrix(NA, n.ci, n.probs)
    qte.ub <- matrix(NA, n.ci, n.probs)
    for (i in 1:n.probs){
      qte.lb[,i] <- boot::boot.ci(boot.kmqte, type="perc",
                          index = (i+ 2* n.probs), conf = ci)$percent[4]
      qte.ub[,i] <- boot::boot.ci(boot.kmqte, type="perc",
                          index = (i+ 2* n.probs), conf = ci)$percent[5]
    }
  }

  if ((n.ci > 1) * (n.probs > 1) == 1){
    qte.lb <- matrix(NA, n.ci, n.probs)
    qte.ub <- matrix(NA, n.ci, n.probs)
    for (i in 1:n.probs){
      qte.lb[,i] <- boot::boot.ci(boot.kmqte, type="perc",
                                 index = (i+ 2* n.probs), conf = ci)$percent[,4]
      qte.ub[,i] <- boot::boot.ci(boot.kmqte, type="perc",
                                 index = (i+ 2* n.probs), conf = ci)$percent[,5]
    }
  }
  #----------------------------------------------------------------------------
  colnames(qte.lb) <- paste(names(quantile(1, probs = probs)), "quantile")
  colnames(qte.ub) <- paste(names(quantile(1, probs = probs)), "quantile")
  rownames(qte.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(qte.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  #----------------------------------------------------------------------------
  # Return these
  list(qte = qte,
       qy1 = qy1,
       qy0 = qy0,
       #boot = boot.kmqte,
       qte.lb = qte.lb,
       qte.ub = qte.ub
  )
}
