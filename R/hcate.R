#' hcate: Testing for Homogeneous Conditional Average Treatment Effetcs
#'
#' \emph{hcate} computes Kolmogorov-Smirnov and Cramer-von Mises type tests
#' for the null hypothesis of homogeneous conditional average treatment effects.
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
#'@importFrom harvestr gather

hcate <- function(out, delta, treat, xvector, xpscore, b, cores = 1) {
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
    pscore <- stats::glm(y ~ ., data = datascore,
                  family = binomial("logit"), x = T)
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

    if (dimx > 1) {
        for (i in 5:dimxe) {
            indx <- indx * outer(fulldata[, i], fulldata[, i], "<=")
        }
    }
    ind <- indx
    rm(indx)

    # Compute unconditional ATE
    wps1 <- fulldata[, 3]/fulldata[, (dim.all - 1)]
    wps0 <- (1 - fulldata[, 3])/(1 - fulldata[, (dim.all - 1)])
    ate <- sum((wps1 - wps0) * fulldata[, 1] * fulldata[, (dim.all - 2)])
    ate <- matrix(rep(ate, n.total), n.total)

    # We now compute the test statistic

    testdist <- ((wps1 - wps0) * fulldata[, 1] - ate) * fulldata[, (dim.all - 2)]
    testdist <- matrix(rep(testdist, n.total), n.total)
    testdist <- testdist * ind
    testdist <- colSums(testdist)

    kstest <- (n.total^0.5) * max(abs(testdist))
    cvmtest <- sum(testdist^2)

    # linear representation for these sub-samples
    linrep.treat <- lr.hom.treated(fulldata = fulldata, subdata = data.treat, dimx = dimx, ate = ate)

    linrep.control <- lr.hom.control(fulldata = fulldata, subdata = data.control, dimx = dimx, ate = ate)
    # Let's put the linear represenations in a single matrix
    linrep.tau <- rbind(linrep.treat, -linrep.control)
    # sort linrep.merged wrt id
    linrep.tau <- linrep.tau[order(as.numeric(abs(linrep.tau[, 1]))), ]
    linrep.tau <- linrep.tau[, -1]

    # Now, the unconditional one
    linrep.unc.treat <- lr.homunc.treated(fulldata = fulldata, subdata = data.treat, dimx = dimx)

    linrep.unc.control <- lr.homunc.control(fulldata = fulldata, subdata = data.control, dimx = dimx)
    # Let's put the linear represenations in a single matrix
    linrep.unc.tau <- rbind(linrep.unc.treat, -linrep.unc.control)
    # sort linrep.merged wrt id
    linrep.unc.tau <- linrep.unc.tau[order(as.numeric(abs(linrep.unc.tau[, 1]))), ]
    linrep.unc.tau <- linrep.unc.tau[, -1]

    # Now, we must match the dimensions for the unconditional ones
    linrep.unc.tau <- matrix(rep(linrep.unc.tau, n.total), n.total, n.total)

    # Remove what we won't need
    rm(linrep.treat, linrep.control, linrep.unc.treat, linrep.unc.control)
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

    yest <- (wt + wc) * fulldata[, (dim.all - 2)] * fulldata[, 1]
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

    # Now the estimation effect of the unconditiona part
    yest.unc <- (wt + wc) * fulldata[, (dim.all - 2)] * fulldata[, 1]

    # X'Y
    xy.unc <- (t(matest)) %*% yest.unc

    # Beta1=(X'X)^-1 (X'Y)
    beta1.unc <- mat1inv %*% xy.unc

    # Compute the 'prediction) X*Beta1
    esteff1.unc <- matest %*% beta1.unc


    # Now, compute the estimation effect of the pscore
    esteff.unc <- fulldata[, 3] - fulldata[, (dim.all - 1)]
    esteff.unc <- esteff.unc * esteff1.unc
    esteff.unc <- matrix(rep(esteff.unc, n.total), n.total)

    # Remove what we won't need
    rm(esteff1.unc, yest.unc, xy.unc, beta1.unc, mat1, matest, mat1inv)

    # Now, we plug in everything to get the linear represenation that we will bootstrap!
    taudist1 <- (linrep.tau - esteff)
    # taudist1=(linrep.tau)
    fx <- colSums(ind * fulldata[, (dim.all - 2)])
    fx <- t(matrix(rep(fx, n.total), n.total, n.total))
    taudist.unc <- (linrep.unc.tau - matrix(ate, n.total, n.total) - esteff.unc) * fx
    # taudist.unc=linrep.unc.tau*fx-as.numeric(ate)*fx
    rm(fx)
    taudist <- (taudist1 - taudist.unc)/n.total
    rm(taudist1, taudist.unc, linrep.tau, linrep.unc.tau)
    # Remove what I won't use
    rm(esteff, ind, fulldata, data.treat, data.control)
    # Now, the bootstrap Number of bootstrap draws
    nboot <- b

    if (cores==1){
      tests <- b.km(n.total = n.total, taudist = taudist,
                    nboot = nboot, kstest = kstest, cvmtest = cvmtest)
    }
    if (cores>1){
      tests <- b.km.mult(n.total = n.total, taudist = taudist,
                         nboot = nboot, kstest = kstest,
                         cvmtest = cvmtest, cores = cores)
    }
    # remove what I do not need
    rm(taudist)

    # Return these
    list(kstest = tests$kstest, cvmtest = tests$cvmtest, pvks = tests$pvks, pvcvm = tests$pvcvm)
}

