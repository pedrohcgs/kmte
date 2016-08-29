################################################################ This file is used to compute the multiplicative Bootstrap proposed in the paper 'Nonparametric
################################################################ Tests for Treatment Effect Heterogeneity with Censored data' by Sant'Anna, Pedro H.C.

################################################################ VERSION : 0.1 MODIFIED : August 27, 2016 at 18:39 AUTHOR : Pedro H. C. Sant'Anna WEBPAGE :
################################################################ https://sites.google.com/site/pedrohcsantanna/ AFFILIATION : Vanderbilt University.  EMAIL :
################################################################ pedro.h.santanna@vanderbilt.edu


# n.total = sample size n of the data taudist = asymptotic linear representation to be
# bootstraped nboot = # of bootstrap draws kstest = Kolmogorov-Smirnov test statistic cvmtest =
# Cramer-von Mises test statistic

# Output is a list of ks and cvm test statistics, with their associated bootstrapped p-values,
# pvks and pvcvm, respectively.

b.km <- function(n.total, taudist, nboot, kstest, cvmtest) {
    # Use the Mammen(1993) binary V's
    k1 <- 0.5 * (1 - 5^0.5)
    k2 <- 0.5 * (1 + 5^0.5)
    pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)
    
    bootapply <- function(nn) {
        v <- rbinom(n.total, 1, pkappa)
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
    boottests <- lapply(1:nboot, bootapply)
    
    # Put the Bootstrap resuls in a matrix
    boottest <- t(matrix(unlist(boottests), 2, nboot))
    
    # Name the Columns
    colnames(boottest) <- c("ksb", "cvmb")
    
    # compute the Bootstrap P-value
    pvksb <- sum((boottest[, 1] > kstest))/nboot
    pvcvmb <- sum((boottest[, 2] > cvmtest))/nboot
    
    list(kstest = kstest, cvmtest = cvmtest, pvks = pvksb, pvcvm = pvcvmb)
}
