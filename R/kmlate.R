#' Kaplan-Meier Local Average Treatment Effect
#'
#' \emph{kmlate} computes the Local Average Treatment Effect for possibly right-censored outcomes.
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
#'@param b 	The number of bootstrap replicates to be performed. Default is 1,000.
#'@param ci A scalar or vector with values in (0,1) containing the confidence level(s)
#'          of the required interval(s). Default is a vector with
#'          0,90, 0.95 and 0.99
#'@param trunc scalar that defined the truncation parameter. Default is NULL, which does not perform any kind of
#'            truncation in the computation of the ATE. When trunc is different than NULL, all outcomes which values greater
#'            than trunc are truncated.
#'@param standardize Default is TRUE, which normalizes instrument propensity score weights to sum to 1 within each
#'                    treatment/instrument group. Set to FALSE to return Horvitz-Thompson weights.
#'@param cores number of processesors to be used during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using snow
#'
#'@return a list containing the Local Average treatment effect estimate, late,
#'        and the bootstrapped \emph{ci} confidence
#'        confidence interval, late.lb (lower bound), and late.ub (upper bound).
#'@export
#'@importFrom stats glm
#'@importFrom parallel makeCluster stopCluster clusterExport
#'@importFrom boot boot.ci boot
#-----------------------------------------------------------------------------
kmlate <- function(out, delta, treat, z, xpscore, b = 1000, ci = c(0.90,0.95,0.99),
                  trunc = NULL, standardize = TRUE, cores = 1) {
  #-----------------------------------------------------------------------------
  # first, we merge all the data into a single datafile
  fulldata <- data.frame(cbind(out, delta, treat, z, xpscore))
  #-----------------------------------------------------------------------------
  # Next, we set up the bootstrap function
  boot1.kmlate <- function(fulldata, i, trunc1 = trunc, standardize1 = standardize){
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


    if (standardize1 == TRUE) {
      w11km.b <- w11km.b / mean(((df.b$treat * df.b$z) / df.b$pscore))
      w10km.b <- w10km.b / mean(((df.b$treat * (1 - df.b$z)) / (1 - df.b$pscore)))
      w01km.b <- w01km.b / mean(((1 - df.b$treat) * df.b$z / df.b$pscore))
      w00km.b <- w00km.b / mean(((1 - df.b$treat) * (1 - df.b$z) * df.b$w / (1 - df.b$pscore)))
      kappa11 <- kappa11 / mean(df.b$z / df.b$pscore)
      kappa10 <- kappa10 / mean((1-df.b$z) / (1 - df.b$pscore))
      kappa01 <- kappa01 / mean(df.b$z / df.b$pscore)
      kappa00 <- kappa00 / mean((1-df.b$z) / (1 - df.b$pscore))
    }

    w1km.b <- w11km.b - w10km.b
    w0km.b <- w01km.b - w00km.b
    kappa1 <- kappa11 - kappa10
    kappa0 <- kappa01 - kappa00

    w1km.b <- w1km.b / kappa1
    w0km.b <- w0km.b / kappa0
    #-----------------------------------------------------------------------------
    # Compute Counterfactual Local Average Outcomes, E[Y(1)|C] and E[Y(0)|C], and the LATE
    meany1km.c <- sum(w1km.b * df.b$out)
    meany0km.c <- sum(w0km.b * df.b$out)

    if (is.null(trunc1) == FALSE){
      meany1km.c <- sum(w1km.b * df.b$out * (df.b$out <= trunc1))
      meany0km.c <- sum(w0km.b * df.b$out * (df.b$out <= trunc1))
    }

    late <- meany1km.c - meany0km.c
    #-----------------------------------------------------------------------------
    return(cbind(meany1km.c, meany0km.c, late))
  }
  #-----------------------------------------------------------------------------
  # Number of bootstrap draws
  nboot <- b
  #----------------------------------------------------------------------------
  #COmput the bootstrap
  if (cores == 1){
    boot.kmlate <- boot::boot(fulldata, boot1.kmlate, R = nboot)
  }
  if (cores > 1){
    cl <- parallel::makeCluster(cores)
    #clusterExport(cl, "kmweight")
    parallel::clusterSetRNGStream(cl)
    boot.kmlate <- boot::boot(fulldata, boot1.kmlate, R = nboot, parallel = "snow", ncpus = cores)
    parallel::stopCluster(cl)
  }

  #----------------------------------------------------------------------------
  # Compute Counterfactual Average Outcomes, E[Y(1)|C] and E[Y(0)|C], and the LATE
  meany1km.c <- boot.ci(boot.kmlate, type="perc", index=1)$t0
  names(meany1km.c) <- "E[Y(1)|Complier]"
  meany0km.c <- boot.ci(boot.kmlate, type="perc", index=2)$t0
  names(meany0km.c) <- "E[Y(0)|Complier]"
  late <- boot.ci(boot.kmlate, type="perc", index=3)$t0
  names(late) <- "LATE"
  #----------------------------------------------------------------------------
  #Compute the confidence interval for ate
  if (length(ci) == 1){
    late.lb <- boot.ci(boot.kmlate, type="perc", index=3, conf = ci)$percent[4]
    late.ub <- boot.ci(boot.kmlate, type="perc", index=3, conf = ci)$percent[5]
  }
  if (length(ci) >1){
    late.lb <- boot.ci(boot.kmlate, type="perc", index=3, conf = ci)$percent[,4]
    late.ub <- boot.ci(boot.kmlate, type="perc", index=3, conf = ci)$percent[,5]
  }

  late.lb <- matrix(late.lb,length(ci),1)
  late.ub <- matrix(late.ub,length(ci),1)
  rownames(late.ub) <- paste(names(quantile(1, probs = ci)), 'CI: UB')
  rownames(late.lb) <- paste(names(quantile(1, probs = ci)), 'CI: LB')
  colnames(late.ub) <- "LATE"
  colnames(late.lb) <- "LATE"
  #----------------------------------------------------------------------------
  # Return these
  list(late = late,
       meany1.c = meany1km.c,
       meany0.c = meany0km.c,
       #boot = boot.kmlate,
       late.lb = late.lb,
       late.ub = late.ub)
}
