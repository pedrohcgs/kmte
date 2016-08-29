################################################################ This file is used to compute the linear representation of treated individuals for testing the
################################################################ null of zero conditional distribution treatment effects, as proposed in the paper
################################################################ 'Nonparametric Tests for Treatment Effect Heterogeneity with Censored data' by Sant'Anna, Pedro
################################################################ H.C.

################################################################ VERSION : 0.1 MODIFIED : August 27, 2016 at 18:39 AUTHOR : Pedro H. C. Sant'Anna WEBPAGE :
################################################################ https://sites.google.com/site/pedrohcsantanna/ AFFILIATION : Vanderbilt University.  EMAIL :
################################################################ pedro.h.santanna@vanderbilt.edu

# fulldata = the entire dataset subdata = sub-sample of treated individuals dimx = #of covariates
# X

lr.cdte.treated <- function(fulldata, subdata, dimx) {
    
    # Dimension of fulldata
    dimfull <- dim(fulldata)[2]
    
    subdatakm <- subdata
    fulldatakm <- fulldata
    
    n.sub <- dim(subdatakm)[1]
    ################################################################################################### Necessary ingridients for the #################################### linear represenation and
    ################################################################################################### #################################### the bootstrap ####################################
    
    # Compute integrand of gamma0
    h.sub <- ecdf(subdatakm[, 1])
    distq.sub <- h.sub(subdatakm[, 1])
    # Avoid problem of numerator (this is just mechanical, given that the last obs of gamma0 is zero)
    s1.sub <- ifelse(distq.sub == 1, 1e-06, 1 - distq.sub)
    nc.sub <- (1 - subdatakm[, 2])/s1.sub
    
    # Compute Gamma0
    indy1.sub <- outer(subdatakm[, 1], subdatakm[, 1], "<")
    
    # gamma0.sub is a n.sub x 1 matrix
    gamma0.sub <- exp((t(indy1.sub) %*% nc.sub)/n.sub)
    
    # Delete what I won't use to save memory
    rm(nc.sub, h.sub)
    
    ################################################################################################### Compute the indicators needed for the computation
    indx.sub <- outer(subdatakm[, 4], fulldatakm[, 4], "<=")
    indy.sub <- outer(subdatakm[, 1], fulldatakm[, 1], "<=")
    dimxe <- 3 + dimx
    if (dimx > 1) {
        for (i in 5:dimxe) {
            indx.sub <- indx.sub * outer(subdatakm[, i], fulldatakm[, i], "<=")
            
        }
    }
    ind.sub <- indx.sub * indy.sub
    nfull <- length(ind.sub[1, ])
    ################################################################################################### Leading term of the expansion for treated
    tau.sub <- gamma0.sub * subdatakm[, 2]/(subdatakm[, (dimfull - 1)])
    tau.sub <- matrix(rep(tau.sub, nfull), n.sub, nfull)
    tau.sub <- tau.sub * ind.sub
    
    # Drop what we won't use anymore
    rm(indx.sub)
    
    ################################################################################################### Compute gamma1 and the second term of the expansion
    t.sub <- (outer(subdatakm[, 1], subdatakm[, 1], ">"))
    
    gamma1.sub <- ((t(t.sub) %*% tau.sub)/n.sub)
    
    # compute the gamma1
    gamma1.sub <- gamma1.sub/(matrix(rep(s1.sub, nfull), n.sub, nfull))
    
    # Now, compute the second term of the expansion
    term2.sub <- gamma1.sub * matrix(rep((1 - subdatakm[, 2]), nfull), n.sub, nfull)
    
    ################################################################################################### Now, compute gamma2 and the last term of the expansion
    gamma2.sub <- matrix(rep((1 - subdatakm[, 2])/s1.sub, nfull), n.sub, nfull)
    gamma2.sub <- gamma1.sub * gamma2.sub
    gamma2.sub <- (t(indy1.sub) %*% (gamma2.sub))/n.sub
    
    ################################################################################################### This is the term of the asymptotic expansion free of the estimation effect
    terms.sub <- tau.sub + term2.sub - gamma2.sub
    # First term is the id of observation
    terms.sub <- cbind(subdata[, dimfull], terms.sub)
    return(terms.sub)
}
