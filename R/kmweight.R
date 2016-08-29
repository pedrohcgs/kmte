
# dimy = the column of datam in which the outcome of interest
#is located
#dimdelta = the column of datam in which the censoring
#indicator is located
#datam = a data set of interest.

# Function to generate KM weigths
kmweight <- function(dimy, dimdelta, datam) {
    datam <- as.data.frame(datam)
    datam <- datam[order(as.numeric(datam[, dimy]),
                         as.numeric(datam[, dimdelta])), ]

    Y11 <- as.numeric(datam[, dimy])
    delta11 <- as.numeric(datam[, dimdelta])


    srt <- order(Y11)
    sy <- as.double(Y11[srt])
    sdelta <- as.integer(delta11[srt])
    n11 <- length(sdelta)
    kmweights <- numeric(n11)
    kmweights[1] <- 1/n11
    if (n11 > 1) {
        for (i in 2:n11) {
            kmweights[i] <- kmweights[i - 1] *
              (n11 - i + 2)/(n11 - i + 1) * (((n11 - i + 1)/(n11 -
                i + 2))^sdelta[i - 1])
        }
    }
    kmwts <- kmweights * sdelta


    datam$w <- kmwts
    datam <- as.data.frame(datam)
    return(datam)
}

