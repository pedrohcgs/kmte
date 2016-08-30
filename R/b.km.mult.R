b.km.mult <- function(n.total, taudist, nboot, kstest, cvmtest,
                      cores, iseed = NULL) {
  # Use the Mammen(1993) binary V's
  k1 <- 0.5 * (1 - 5^0.5)
  k2 <- 0.5 * (1 + 5^0.5)
  pkappa <- 0.5 * (1 + 5^0.5)/(5^0.5)

  ## Define seeds
  if (!is.null(iseed)){
    seed.temp <- harvestr::gather(nboot, seed = iseed)
  }
  if (is.null(iseed)){
    ss=floor(stats::runif(1)*10000)
    seed.temp <- harvestr::gather(nboot, seed = ss)
  }

  Seed <- matrix(nrow = nboot, ncol = 6)
  for(i in 1:nboot){
    Seed[i,] <- seed.temp[[i]][2:7]
  }
  # Create bootstrap function
  bootapply <- function(nn, n.total, pkappa, k1, k2,
                        taudist,Seed) {
    seed.run <- Seed[nn,]
    set.seed(seed.run, "L'Ecuyer-CMRG") ## to make each run fully reproducible
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

  cl <- parallel::makeCluster(cores)
  if (!is.null(iseed)){
    parallel::clusterSetRNGStream(cl = cl, iseed)
  }
  boottests <- parallel::parLapply(cl, 1:nboot, bootapply,
                         n.total, pkappa, k1, k2,
                         taudist,Seed)
  parallel::stopCluster(cl)

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
