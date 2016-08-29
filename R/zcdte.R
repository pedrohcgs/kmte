
#' zcdte: Testing for Zero Conditional Distribution Treatment Effetcs
#'
#' \emph{zcdte} computes Kolmogorov-Smirnov and Cramer-von Mises type tests
#' for the null hypothesis of zero conditional distribution treatment effects.
#' The test is suitable for both censored and uncensored outcomes, and relies on
#' the unconfoundedness assumption. For details of the testing procedure, see
#' Sant'Anna (2016b),'Nonparametric Tests for Treatment Effect Heterogeneity with
#' Censored data'.
#'
#'@param out vector containing the outcome of interest
#'@param delta vector containing the censoring indicator (1 if observed, 0 if censored)
#'@param treat vector containing the treatment indicator (1 if treated, 0 if control)
#'@param xvector matrix (or data frame) containing the conditioning covariates
#'@param xpscore matrix (or data frame) containing the covariates (and their
#'               transformations) to be included in the propensity score estimation
#'@param b number of bootstrap draws
#'@param cores number of cores to use during the bootstrap (default is 1).
#'        If cores>1, the bootstrap is conducted using parLapply, instead
#'        of lapply type call.
#'
#'@return a list containing the Kolmogorov-Smirnov test statistic (kstest),
#'        the Cramer-von Mises test statistic (cvmtest), and their associated
#'        bootstrapped p-values, pvks and pvcvm, respectively.
#'@export
#'@importFrom stats binomial ecdf glm rbinom
#'@importFrom MASS ginv
#'@importFrom parallel makeCluster parLapply stopCluster


zcdte <- function(out, delta, treat, xvector, xpscore, b, cores=1) {
    # first, we merge all the data into a single datafile
    fulldata <- data.frame(cbind(out, delta, treat, xvector, xpscore))
    # Compute Kaplan-Meier Weigths - data is now sorted!
    fulldata <- kmweight(1, 2, fulldata)
    # Dimension of data matrix fulldata
    dim.all <- dim(fulldata)[2]
    # Next, we rename the variable in xpscore to avoid problems Dimension of Xvector
    dimx <- dim(xvector)[2]
    # Where it ends
    dimxe <- 3 + dimx
    xpscore1 <- fulldata[, ((dimxe + 1):(dim.all - 1))]
    datascore <- data.frame(y = fulldata[, 3], xpscore1)
    # estimate the propensity score
    pscore <- stats::glm(y ~ ., data = datascore, family = binomial("logit"), x = T)
    fulldata$pscore <- pscore$fit

    # Create id to help on ordering
    fulldata$id <- 1:length(fulldata[, 1])

    # Dimension of data matrix fulldata
    dim.all <- dim(fulldata)[2]

    # sample size
    n.total <- as.numeric(length(fulldata[, 1]))

    # subset of treated individuals
    data.treat <- subset(fulldata, fulldata[, 3] == 1)
    # subset of not-treated individuals
    data.control <- subset(fulldata, fulldata[, 3] == 0)

    # Compute Kaplan-Meier weigth for treated
    data.treat <- kmweight(1, 2, data.treat)
    n.treat <- as.numeric(length(data.treat[, 1]))
    data.treat$w <- data.treat$w * (n.treat/n.total)

    # Compute Kaplan-Meier weigth for control
    data.control <- kmweight(1, 2, data.control)
    n.control <- as.numeric(length(data.control[, 1]))
    data.control$w <- data.control$w * (n.control/n.total)

    # Let's put everything in a single data - correct KM weigths First, the datasets
    fulldata <- data.frame(rbind(data.treat, data.control))
    # Sort wrt id
    fulldata <- fulldata[order(as.numeric(fulldata[, dim.all])), ]


    # Compute the indicators needed for the computation of the test
    indx <- outer(fulldata[, 4], fulldata[, 4], "<=")
    indy <- outer(fulldata[, 1], fulldata[, 1], "<=")
    # Dimension of X
    dimx <- dim(xvector)[2]
    dimxe <- 3 + dimx
    if (dimx > 1) {
        for (i in 5:dimxe) {
            indx <- indx * outer(fulldata[, i], fulldata[, i], "<=")
        }
    }
    ind <- indx * indy


    # We compute the test statistic
    testdist <- (((fulldata[, 3]/fulldata[, (dim.all - 1)]) - ((1 - fulldata[, 3])/(1 - fulldata[,
        (dim.all - 1)]))) * fulldata[, (dim.all - 2)])

    testdist <- matrix(rep(testdist, n.total), n.total)
    testdist <- testdist * ind
    testdist <- colSums(testdist)

    kstest <- (n.total^0.5) * max(abs(testdist))
    cvmtest <- sum(testdist^2)

    # Compute linear representations
    linrep.treat <- lr.cdte.treated(fulldata = fulldata, subdata = data.treat, dimx = dimx)
    linrep.control <- lr.cdte.control(fulldata = fulldata, subdata = data.control, dimx = dimx)

    # Let's put the linear represenations in a single matrix
    linrep.tau <- rbind(linrep.treat, -linrep.control)
    # sort linrep.merged wrt id
    linrep.tau <- linrep.tau[order(as.numeric(abs(linrep.tau[, 1]))), ]
    linrep.tau <- linrep.tau[, -1]

    # Remove what we won't need
    rm(linrep.treat, linrep.control)
    # Now, the estimation effect Prepare the matrix
    matest <- pscore$x
    # X'X
    mat1 <- (t(matest)) %*% matest
    mat1 <- mat1/n.total
    # (X'X)^-1
    mat1inv <- MASS::ginv(mat1)

    # Compute the terms in the 'Y'
    wt <- fulldata[, 3]/(fulldata[, (dim.all - 1)]^2)
    wc <- (1 - fulldata[, 3])/((1 - fulldata[, (dim.all - 1)])^2)

    yest <- (wt + wc) * fulldata[, (dim.all - 2)]
    yest <- matrix(rep(yest, n.total), n.total)
    yest <- yest * ind

    # X'Y
    xy <- (t(matest)) %*% yest

    # Beta1=(X'X)^-1 (X'Y)
    beta1 <- mat1inv %*% xy

    # Compute the 'prediction) X*Beta1
    esteff1 <- matest %*% beta1
    # Now, compute the estimation effect of the pscore
    esteff <- fulldata[, 3] - fulldata[, (dim.all - 1)]
    esteff <- matrix(rep(esteff, n.total), n.total)
    esteff <- esteff * esteff1

    # Remove what we won't need
    rm(esteff1, yest)
    # Now, we plug in everything to get the linear represenation that we will bootstrap!
    taudist <- (linrep.tau - esteff)/n.total

    rm(esteff, ind)
    # Now, the bootstrap Number of bootstrap draws
    nboot <- b

    if (cores==1){
      tests <- b.km(n.total = n.total, taudist = taudist,
                    nboot = nboot, kstest = kstest, cvmtest = cvmtest)
    }
    if (cores>1){
      tests <- b.km.mult(n.total = n.total, taudist = taudist,
                         nboot = nboot, kstest = kstest, cvmtest = cvmtest,
                         cores = cores)
    }
    # remove what I do not need
    rm(taudist)

    # Return these
    list(kstest = tests$kstest, cvmtest = tests$cvmtest, pvks = tests$pvks, pvcvm = tests$pvcvm)
}

