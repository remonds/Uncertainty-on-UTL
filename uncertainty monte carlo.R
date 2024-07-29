##########################################################################################
# uncertainty on UTL by Monte Carlo using frequentist approach
# 
# iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Uncertainty 
# of measurement Part 3: Guide to the expression 
# of uncertainty in measurement (GUM:1995) 
# Supplement 1: Propagation of distributions using 
# a Monte Carlo method, 5.3.4 and 7.7.2
# 
# 
# Robert Emonds   BSOH           robert.emonds@bsoh.be, robert.emonds@outlook.com
# Theo Scheffers  TSAC           theo.scheffers@tsac.nl
# Peter van Balen PreventPartner peter.van.balen@preventpartner.nl
##########################################################################################
# 28/07/2024
# added adaptive procedure to adapt number of trials 
# to desired number of significant decimal digits
# 
# 23/07/2024
# added p-values plow and phigh to output
# 
# 19/07/2024
# utl.ros.mc(c, detects, CVt)
# for datasets with nondetects
# measurement uncertainty interval on utl by Multiple Linear Regression on detects
# frequentist approach
# see description below of input and output variables
# 
# 26/05/2024
# utl.mc(c, CVt)
# for datasets without nondetects
# measurement uncertainty interval on utl - no nondetects
# frequentist approach
# see description below of input and output variables
##########################################################################################
# utl.mc(c, CVt)
# utl.ros.mc(c, detects, CVt)
# 
# input variables
# c       : numeric vector of measured concentration 
#           values (mg/m³)
# CVt     : expression of uncertainty
#           numeric vector of total coefficient of 
#           variation of measurement 
#           (sampling variation + analytical variation)
# detects : logical vector representing detects (TRUE) 
#           and nondetects (FALSE) 
# ndig    : number of significant decimal digits regarded
#           as meaningful in the numerical value of UTL
#           ndig = 2  computing time is in the order of seconds
#           ndig = 3  computing time is in the order of minutes
# 
# output variables
# yest : estimate of Y, obtained as the average of 
#        the M model values yr from a Monte Carlo run
# ycv  : relative standard uncertainty associated 
#        with ~ yest
# ylow : left-hand endpoint of a coverage interval 
#        for Y
# yhigh: right-hand endpoint of a coverage interval 
#        for Y
# plow : p-value for ylow 
# phigh: p-value for yhigh 
# Mpos : number of Monte Carlo trials with positive
#        values for Y (=UTL)
# ndig : number of significant decimal digits regarded
#        as meaningful in the numerical value of UTL
#        - is same as the input value
##########################################################################################
# Y    : (scalar) output quantity, regarded as a 
#        random variable
#        = 70% one-sided upper tolerance limit of Q95 (UTL95,70)
#        ("confidence limit" is used for the mean)
##########################################################################################
# test values used during development (cottondust)
# c   <- c(0.16, 0.38, 0.20, 0.44, 0.51, 0.60, 0.35, 0.70, 0.18, 0.65)
# CVt <- c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
# CVt <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
# 
# mc <- utl.mc(c, CVt)
######################################################################

utl.mc <- function(c, CVt, ndig = 2) {
  # coverage level for exposure
  pexp <- 0.70
  
  # upper exceedance fraction threshold for exposure
  ueft <- 0.05
  
  # coverage level for measurement uncertainty
  pmu  <- 0.95
  
  # number of Monte Carlo trials
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4 Adaptive procedure
  M    <- max(100/(1-pmu), 1e4)
  
  # number of significant decimal digits regarded
  # as meaningful in the numerical value of UTL
  # ndig <- 3 # computing time is in the order of minutes
  # ndig <- 2 # computing time is in the order of seconds
  
  # number of input quantities
  N    <- length(c)
  
  # coverage factor for quantiles
  UT   <- qt(pexp,
             df = N-1,
             ncp = sqrt(N)*qnorm(1-ueft)
             )/sqrt(N)
  
  # initialize counter for number of iterations
  h    <-1
  
  # Initialize the list to store results of each iteration
  mclist <- list()

  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around c
    X <- data.frame(
      matrix(
        rnorm(M*N, mean = c, sd = CVt*c),
        ncol = N, byrow = TRUE
      )
    )
    names(X) <- paste("x", c(1:N), sep = "")
    
    # filter positive values
    # - because in the following step the log needs to be taken
    # - this is a quick and dirty solution
    # - it could perhaps distort in some way the shape of the 
    #   distribution profile
    # - might need some refinement
    X <- X[apply(X, 1, function(row) all(row > 0)), ]
    
    # Transform to the normally distributed space
    # store log(Xpos) to reduce calculation time a bit
    X <- log(X)
    
    ##########################################################################################
    # Propagation of the input distributions through the model
    # model: calculate by row
    #   1. GM and GSD
    #   2. UTL = GM*GSD^UT
    ##########################################################################################
    # Transform back to the lognormally distributed space and
    # propagate the input distributions through the model
    # The model UTL=GM*GSD^UT is iaw EN 689 Annex F
    # GM <- exp(rowMeans(X))
    GM   <- exp(apply(X,1,mean))
    GSD  <- exp(apply(X,1,sd  ))
    UTL  <- as.double(GM*GSD^UT)
    
    ##########################################################################################
    # Statistics of UTL
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    ##########################################################################################
    yest <- mean(UTL)
    uy   <- sd(UTL)
    ycv  <- uy/yest
    Mpos <- length(UTL)
    
    # index for coverage width in output vector
    q <- as.integer(pmu*Mpos+0.5)
    
    # width of the coverage interval
    UTL   <- as.double(sort(UTL, decreasing = FALSE))
    dy    <- UTL[(q+1):Mpos] - UTL[1:(Mpos-q)]
    
    # minimum width of the coverage interval
    dymin <- min(dy)
    
    # find index of interval with minimum width
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind)+0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow  <- UTL[ind]
    yhigh <- UTL[ind+q]
    # associated p-values
    plow  <- (ind)/Mpos
    phigh <- (ind+q)/Mpos
    
    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y     = UTL,
      yest  = yest,
      uy    = uy,
      ycv   = ycv,
      ylow  = ylow,
      yhigh = yhigh,
      Mpos  = Mpos
    )))
    
    if (h < 2) {
      h <- h + 1
      next
    }
    
    # standard deviation of h*M values of Y
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1
    sd_y <- sd(unlist(lapply(mclist, `[[`, "Y")))

    # numerical tolerance, delta
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.2
    delta <- 0.5 * 10^(floor(log10(sd_y) + 1 - ndig))

    # standard error of the h values of yest, ycv, ylow, yhigh
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4
    yest_se  <- sd(sapply(mclist, `[[`, "yest" ))/sqrt(h)
    # ycv_se   <- sd(sapply(mclist, `[[`, "ycv"  ))/sqrt(h)
    # ycv is not in same order of magnitude as the rest -> use uy
    uy_se    <- sd(sapply(mclist, `[[`, "uy"   ))/sqrt(h)
    ylow_se  <- sd(sapply(mclist, `[[`, "ylow" ))/sqrt(h)
    yhigh_se <- sd(sapply(mclist, `[[`, "yhigh"))/sqrt(h)
    
    # Check if the tolerance condition is met
    if (2 * yest_se  < delta &&
        2 * uy_se    < delta &&
        2 * ylow_se  < delta &&
        2 * yhigh_se < delta ) {
      break
    }
    
    # avoid runaway
    if (h > 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    h <- h + 1
  }

  # endresult based on all the h*M generated samples
  UTL  <- unlist(lapply(mclist, `[[`, "Y"))
  UTL  <- as.double(sort(UTL, decreasing = FALSE))
    yest <- mean(UTL)
  uy   <- sd(UTL)
    ycv  <- uy/yest
    Mpos <- length(UTL)
  q <- as.integer(pmu*Mpos+0.5)
  dy    <- UTL[(q+1):Mpos] - UTL[1:(Mpos-q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind)+0.5)
    ylow  <- UTL[ind]
    yhigh <- UTL[ind+q]
    plow  <- (ind)/Mpos
    phigh <- (ind+q)/Mpos
  
  return(
    list(
      yest  = signif(yest,  digits = ndig),
      ycv   = signif(ycv,   digits = ndig),
      ylow  = signif(ylow,  digits = ndig),
      yhigh = signif(yhigh, digits = ndig),
      plow  = signif(plow,  digits = ndig),
      phigh = signif(phigh, digits = ndig),
      Mpos  = Mpos,
      ndig  = ndig
    )
  )
}

##########################################################################################
# test values used during development (EN 689 Annex H)
# c <-      c( 0.75,  2.00,  0.75,  0.80,  0.92,  1.05,  1.10,  1.45,  0.75,  2.50)
# detects <-c(FALSE,  TRUE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  FALSE, TRUE)
# CVt <-    c( 0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15)
# CVt <-    c( 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00)
# 
# mc <- utl.ros.mc(c, detects, CVt)
##########################################################################################
  
utl.ros.mc <- function(c, detects, CVt, ndig = 2) {
  # coverage level for exposure
  pexp <- 0.70
  
  # upper exceedance fraction threshold for exposure
  ueft <- 0.05
  
  # coverage level for measurement uncertainty
  pmu  <- 0.95
  
  # number of Monte Carlo trials
  M    <- max(100/(1-pmu), 1e4)
  
  # number of significant decimal digits regarded
  # as meaningful in the numerical value of UTL
  # ndig <- 3 # computing time is in the order of minutes
  # ndig <- 2 # computing time is in the order of seconds
  
  # number of input quantities
  N    <- length(c)
  
  # coverage factor for quantiles
  UT   <- qt(pexp,
             df = N-1,
             ncp = sqrt(N)*qnorm(1-ueft)
             )/sqrt(N)

  # create data frame
  Xdf <- data.frame(
    xk = c,
    CVt = CVt,
    detects = detects
  )

  # keep only the detects
  Xdf <- Xdf[detects==TRUE,]
  Ndet <- nrow(Xdf)
  
  # initialize counter for number of iterations
  h    <-1
  
  # Initialize the list to store results of each iteration
  mclist <- list()
  
  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around the detects in c
    X <- data.frame(
      matrix(
        rnorm(M*Ndet, mean = Xdf$xk, sd = Xdf$CVt*Xdf$xk),
        ncol = Ndet, byrow = TRUE
        )
      )
    names(X) <- paste("x", c(1:Ndet), sep = "")
  
    # filter positive values
    # - because in the following step the log needs to be taken
    # - this is a quick and dirty solution
    # - it could perhaps distort in some way the shape of the 
    #   distribution profile
    # - might need some refinement
    X <- X[apply(X, 1, function(row) all(row > 0)), ]
    
    # Transform to the normally distributed space
    # store log(Xpos) to reduce calculation time a bit
    X <- log(X)
    
    ##########################################################################################
    # Propagation of the input distributions through the model
    # model: calculate by row
    #   1. sort rows in ascending order
    #   2. linear regression
    #   3. GM = intercept and GSD = slope
    #   4. UTL = GM*GSD^UT
    ##########################################################################################
    # Sort the data in ascending order to obtain the order statistics
    X <- as.data.frame(t(apply(X, 1, sort)))
    names(X) <- paste("y", c(1:Ndet), sep = "")
  
    # Blom's rankit
    zk <- qnorm((((N-Ndet+1):N) - 3/8) / (N + 1/4))
  
    # Matrix Formulation of the Multiple Linear Regression (MLR) Model
    # Linear regression iaw En 689 Annex H
    # For the regression, only the detects are used
    # the intercept represents the arithmetic mean (AMr)
    # the slope represents the arithmetic standard deviation (ASDr)
    # the index r indicates that the value is obtained by using
    # the regression coefficients
    library(data.table)
    MLR <- function(known_ys, known_xs) {
      X <- cbind(1, known_xs)
      # solve(t(X) %*% X) gives the inverse of (t(X) %*% X)
      beta <- solve(t(X) %*% X) %*% t(X) %*% t(as.matrix(known_ys))
      intercepts <- beta[1, ]
      slopes <- beta[2, ]
      data.frame(AMr = intercepts, ASDr = slopes)
    }
    MLRcoefs <- MLR(X, zk)
    
    # Transform back to the lognormally distributed space and
    # propagate the input distributions through the model
    # The model UTL=GM*GSD^UT is iaw EN 689 Annex H
    # GMr   <- exp(MLRcoefs[,1])
    # GSDr  <- exp(MLRcoefs[,2])
    GMr   <- exp(MLRcoefs$AMr)
    GSDr  <- exp(MLRcoefs$ASDr)
    UTL   <- as.double(GMr*GSDr^UT)
  
    ##########################################################################################
    # Statistics of UTL
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    ##########################################################################################
    yest <- mean(UTL)
    uy   <- sd(UTL)
    ycv  <- uy/yest
    Mpos <- length(UTL)
    
    # index for coverage width in output vector
    q <- as.integer(pmu*Mpos+0.5)
    
    # width of the coverage interval
    UTL   <- as.double(sort(UTL, decreasing = FALSE))
    dy    <- UTL[(q+1):Mpos] - UTL[1:(Mpos-q)]
    
    # minimum width of the coverage interval
    dymin <- min(dy)
    
    # find index of interval with minimum width
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind)+0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow  <- UTL[ind]
    yhigh <- UTL[ind+q]
    # associated p-values
    plow  <- (ind)/Mpos
    phigh <- (ind+q)/Mpos
    
    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y     = UTL,
      yest  = yest,
      uy    = uy,
      ycv   = ycv,
      ylow  = ylow,
      yhigh = yhigh,
      Mpos  = Mpos
    )))
    
    if (h < 2) {
      h <- h + 1
      next
    }
    
    # standard deviation of h*M values of Y
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1
    sd_y <- sd(unlist(lapply(mclist, `[[`, "Y")))
    
    # numerical tolerance, delta
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.2
    delta <- 0.5 * 10^(floor(log10(sd_y)) + 1 - ndig)

    # standard error of the h values of yest, ycv, ylow, yhigh
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4
    yest_se  <- sd(sapply(mclist, `[[`, "yest" ))/sqrt(h)
    # ycv_se   <- sd(sapply(mclist, `[[`, "ycv"  ))/sqrt(h)
    # ycv is not in same order of magnitude as the rest -> use uy
    uy_se    <- sd(sapply(mclist, `[[`, "uy"   ))/sqrt(h)
    ylow_se  <- sd(sapply(mclist, `[[`, "ylow" ))/sqrt(h)
    yhigh_se <- sd(sapply(mclist, `[[`, "yhigh"))/sqrt(h)
    
    # Check if the tolerance condition is met
    if (2 * yest_se  < delta &&
        2 * uy_se    < delta &&
        2 * ylow_se  < delta &&
        2 * yhigh_se < delta ) {
      break
    }
    
    # avoid runaway
    if (h > 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    h <- h + 1
  }
  
  # endresult based on all the h*M generated samples
  UTL  <- unlist(lapply(mclist, `[[`, "Y"))
  UTL  <- as.double(sort(UTL, decreasing = FALSE))
    yest <- mean(UTL)
  uy   <- sd(UTL)
    ycv  <- uy/yest
    Mpos <- length(UTL)
  q <- as.integer(pmu*Mpos+0.5)
  dy    <- UTL[(q+1):Mpos] - UTL[1:(Mpos-q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind)+0.5)
    ylow  <- UTL[ind]
    yhigh <- UTL[ind+q]
    plow  <- (ind)/Mpos
    phigh <- (ind+q)/Mpos
  
  return(
    list(
      yest  = signif(yest,  digits = ndig),
      ycv   = signif(ycv,   digits = ndig),
      ylow  = signif(ylow,  digits = ndig),
      yhigh = signif(yhigh, digits = ndig),
      plow  = signif(plow,  digits = ndig),
      phigh = signif(phigh, digits = ndig),
      Mpos  = Mpos,
      ndig  = ndig
    )
  )
}
