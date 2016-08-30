
b.km <- function(n.total, taudist, nboot, kstest, cvmtest, iseed = NULL) {
    # Use the Mammen(1993) binary V's
    k1 <- 0.5 * (1 - 5^0.5)
    k2 <- 0.5 * (1 + 5^0.5)
    pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)

    bootapply <- function(nn) {
        v <- stats::rbinom(n.total, 1, pkappa)
        v <- ifelse(v == 1, k1, k2)
        v <- matrix(rep(v, n.total), n.total)
        # Bootstrap tau
        taudistb <- taudist * v
        # Bootstrap linear represenations
        testdistb <- colSums(taudistb)
        # KS test
        ksb <- (n.total^0.5) * max(abs(testdistb))
        # Cramer-von Mises test
        cvmb <- sum(testdistb^2)
        # Return both tests
        return(cbind(ksb, cvmb))
    }
    if (!is.null(iseed)){
      set.seed(iseed, kind = "L'Ecuyer")
    }
    boottests <- lapply(1:nboot, bootapply)

    # Put the Bootstrap resuls in a matrix
    boottest <- t(matrix(unlist(boottests), 2, nboot))

    # Name the Columns
    colnames(boottest) <- c("ksb", "cvmb")

    # compute the Bootstrap P-value
    pvksb <- sum((boottest[, 1] > kstest))/nboot
    pvcvmb <- sum((boottest[, 2] > cvmtest))/nboot

    list(kstest = kstest, cvmtest = cvmtest, pvks = pvksb,
         pvcvm = pvcvmb)
}
