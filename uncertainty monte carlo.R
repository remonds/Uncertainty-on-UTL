##########################################################################################
# Uncertainty on UTL by Monte Carlo using frequentist approach
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
# 15/08/2024
# added the distribution of the output quantity Y to output
# 
# 15/08/2024
# added classic symmetric confidence bounds
# 
# 08/08/2024
# added data frame to output with characteristic quantiles
# of uncertainty distribution of UTL
# for research purposes
# 
# 07/08/2024
# added parameters ueft, pexp and pmu to input, with EN689 defaults
# remains compatible with previous versions
# utl.mc(twa, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95)
# utl.ros.mc(twa, detects, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95)
# 
# 28/07/2024
# added adaptive procedure to adapt number of trials 
# to desired number of significant decimal digits
# utl.mc(twa, CVt, ndig = 2)
# utl.ros.mc(twa, detects, CVt, ndig = 2)
# remains compatible with previous versions
# 
# 23/07/2024
# added p-values plow and phigh to output
# 
# 19/07/2024
# utl.ros.mc(twa, detects, CVt)
# for datasets with nondetects
# measurement uncertainty interval on utl by Multiple Linear Regression on detects
# frequentist approach
# see description below of input and output variables
# 
# 26/05/2024
# utl.mc(twa, CVt)
# for datasets without nondetects
# measurement uncertainty interval on utl - no nondetects
# frequentist approach
# see description below of input and output variables
##########################################################################################
# utl.mc(twa, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95)
# utl.ros.mc(twa, detects, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95)
# 
# input variables
# twa     : numeric vector of 8h time-weighted exposure
#           values (mg/m³ or ppm)
#           These values are obtained from the measured 
#           concentrations from the personal air samplings 
#           and the estimated daily exposure time 
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
# ueft    : upper exceedance fraction threshold for exposure
# pexp    : 1-sided coverage level for exposure
# pmu     : 2-sided coverage level for measurement uncertainty
# 
# output variables
# yest : estimate of Y, obtained as the average of 
#        the M model values yr from a Monte Carlo run
# ycv  : relative standard uncertainty associated 
#        with ~ yest
# ylow : left-hand endpoint of a coverage interval 
#        for Y (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008)
# yhigh: right-hand endpoint of a coverage interval 
#        for Y (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008)
# plow : p-value for ylow  (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008)
# phigh: p-value for yhigh (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008)
# yq   : characteristic quantiles for research purposes
#        p-values:
#        0.000 minimum y-value of Monte Carlo trials
#        plow  ══════╗
#        0.025 ──────╬───┐
#        0.050 ──┐   ║   │  coverage interval
#        0.500  90% 95% 95% of measurement 
#        0.950 ──┘   ║   |  uncertainty
#        phigh ══════╝   |
#        0.975 ──────────┘
#        1.000 maximum y-value of Monte Carlo trials
# Mpos : number of Monte Carlo trials with positive
#        values for Y (= UTL)
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
# twa  <- c(0.16, 0.38, 0.20, 0.44, 0.51, 0.60, 0.35, 0.70, 0.18, 0.65)
# CVt  <- c(0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15)
# CVt  <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
# ueft <- 0.05
# pexp <- 0.70
# pmu  <- 0.95
# ndig <- 2
# 
# mc <- utl.mc(twa, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95)
######################################################################

utl.mc <- function(twa, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95) {
  # classic symmetric probabilities for confidence interval
  plowsym  <- (1 - pmu) / 2
  phighsym <- (1 + pmu) / 2 # = 1 - (1 - pmu) / 2
  
  # number of Monte Carlo trials
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4 Adaptive procedure
  M    <- max(100 / (1 - pmu), 1e4)
  
  # number of input quantities
  N    <- length(twa)
  
  # coverage factor for quantiles
  UT   <- qt(p   = pexp,
             df  = N - 1,
             ncp = sqrt(N) * qnorm(1 - ueft)
             ) / sqrt(N)
  
  # initialize counter for number of iterations
  h    <- 1
  
  # Initialize the list to store results of each iteration
  mclist <- list()
  
  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around twa
    X <- data.frame(
      matrix(
        rnorm(M * N, mean = twa, sd = CVt * twa),
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
    #   2. UTL = GM * GSD^UT
    ##########################################################################################
    # Transform back to the lognormally distributed space and
    # propagate the input distributions through the model
    # The model UTL=GM*GSD^UT is iaw EN 689 Annex F
    # GM <- exp(rowMeans(X))
    GM   <- exp(apply(X, 1, mean))
    GSD  <- exp(apply(X, 1, sd  ))
    UTL  <- as.double(GM * GSD^UT)
    
    ##########################################################################################
    # Statistics of UTL
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    ##########################################################################################
    yest  <- mean(UTL)
    uy    <- sd(UTL)
    ycv   <- uy / yest
    Mpos  <- length(UTL)
    
    # index for coverage width in output vector
    q     <- as.integer(pmu * Mpos + 0.5)
    
    # width of the coverage interval
    UTL   <- as.double(sort(UTL, decreasing = FALSE))
    dy    <- UTL[(q + 1):Mpos] - UTL[1:(Mpos - q)]
    
    # minimum width of the coverage interval
    dymin <- min(dy)
    
    # find index of interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 7.7.2
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind) + 0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow  <- UTL[ind]
    yhigh <- UTL[ind + q]
    # associated p-values
    plow  <- (ind) / Mpos
    phigh <- (ind + q) / Mpos
    
    # classic symmetric confidence bounds
    ylowsym  <- quantile(
      x      = UTL,
      probs  = plowsym,
      names  = FALSE
    )
    yhighsym <- quantile(
      x      = UTL,
      probs  = phighsym,
      names  = FALSE
    )

    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y     = UTL,
      yest  = yest,
      uy    = uy,
      ycv   = ycv,
      ylow  = ylow,
      yhigh = yhigh,
      ylowsym  = ylowsym,
      yhighsym = yhighsym,
      Mpos  = Mpos
    )))
    
    if (h < 2) {
      h <- h + 1
      next
    }
    
    # standard deviation of h*M values of Y
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1
    sd_y  <- sd(unlist(lapply(mclist, `[[`, "Y")))

    # numerical tolerance, delta
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.2
    delta    <- 0.5 * 10^(floor(log10(sd_y)) + 1 - ndig)

    # standard error of the h values of yest, ycv, ylow, yhigh
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4
    yest_se  <- sd(sapply(mclist, `[[`, "yest" )) / sqrt(h)
    # ycv_se   <- sd(sapply(mclist, `[[`, "ycv"  )) / sqrt(h)
    # ycv is not in same order of magnitude as the rest -> use uy
    uy_se    <- sd(sapply(mclist, `[[`, "uy"   )) / sqrt(h)
    ylow_se  <- sd(sapply(mclist, `[[`, "ylow" )) / sqrt(h)
    yhigh_se <- sd(sapply(mclist, `[[`, "yhigh")) / sqrt(h)
    
    # Check if the tolerance condition is met
    if (2 * yest_se  <= delta &&
        2 * uy_se    <= delta &&
        2 * ylow_se  <= delta &&
        2 * yhigh_se <= delta) {
      break
    }
    
    # avoid runaway
    if (h >= 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    h <- h + 1
  }

  # endresult based on all the h*M generated samples
  UTL   <- unlist(lapply(mclist, `[[`, "Y"))
  UTL   <- as.double(sort(UTL, decreasing = FALSE))
  yest  <- mean(UTL)         # output value
  uy    <- sd(UTL)
  ycv   <- uy / yest         # output value
  Mpos  <- length(UTL)       # output value
  q     <- as.integer(pmu * Mpos + 0.5)
  dy    <- UTL[(q + 1):Mpos] - UTL[1:(Mpos - q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind) + 0.5)
  ylow  <- UTL[ind]          # output value
  yhigh <- UTL[ind + q]      # output value
  plow  <- (ind) / Mpos      # output value
  phigh <- (ind + q) / Mpos  # output value
  prob <- unique(sort(
    c(signif(c(plow   , phigh   ), digits = max(ndig, 3)),
      signif(c(plowsym, phighsym), digits = max(ndig, 3)), 
      c(0.000, 0.050, 0.500, 0.950 , 1.000)
    )))
  yq       <- quantile(
    x      = UTL,
    probs  = prob,
    names  = FALSE,
    digits = ndig
  )
  dUTL     <- density(UTL)
  dUTLsub  <- approx(dUTL$x, dUTL$y, xout = yq)$y
  yq       <- data.frame(       # output value
    "p"    = prob,
    "d"    = signif(dUTLsub, digits = max(ndig, 3)),
    "UTL"  = signif(yq     , digits = ndig),
    row.names = format(prob)
  )
  ylowsym  <- yq[yq$p == signif(plowsym,  digits = max(ndig, 3)),]$UTL
  yhighsym <- yq[yq$p == signif(phighsym, digits = max(ndig, 3)),]$UTL
  
  return(
    list(
      yest  = signif(yest,  digits = ndig),
      ycv   = signif(ycv,   digits = ndig),
      ylow  = signif(ylow,  digits = ndig),
      yhigh = signif(yhigh, digits = ndig),
      plow  = signif(plow,  digits = max(ndig, 3)),
      phigh = signif(phigh, digits = max(ndig, 3)),
      ylowsym  = signif(ylowsym,  digits = ndig),
      yhighsym = signif(yhighsym, digits = ndig),
      plowsym  = signif(plowsym,  digits = max(ndig, 3)),
      phighsym = signif(phighsym, digits = max(ndig, 3)),
      yq    = yq,
      Y     = UTL,
      Mpos  = Mpos,
      ndig  = ndig,
      tolerance = delta
    )
  )
}

##########################################################################################
# test values used during development (EN 689 Annex H)
# twa     <- c( 0.75,  2.00,  0.75,  0.80,  0.92,  1.05,  1.10,  1.45,  0.75,  2.50)
# detects <- c(FALSE,  TRUE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  FALSE, TRUE)
# CVt     <- c( 0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15)
# CVt     <- c( 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00)
# ueft    <- 0.05
# pexp    <- 0.70
# pmu     <- 0.95
# ndig <- 2
# 
# mc <- utl.ros.mc(twa, detects, CVt)
##########################################################################################
  
utl.ros.mc <- function(twa, detects, CVt, ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95) {
  # classic symmetric probabilities for confidence interval
  plowsym  <- (1 - pmu) / 2
  phighsym <- (1 + pmu) / 2 # = 1 - (1 - pmu) / 2
  
  # number of Monte Carlo trials
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 Suppl 1, 7.9.4 Adaptive procedure
  M    <- max(100 / (1 - pmu), 1e4)
  
  # number of input quantities
  N    <- length(twa)
  
  # coverage factor for quantiles
  UT   <- qt(p   = pexp,
             df  = N - 1,
             ncp = sqrt(N) * qnorm(1 - ueft)
             ) / sqrt(N)

  # put dataset in data frame for computational purposes
  Xdf  <- data.frame(
    xk = twa,
    CVt = CVt,
    detects = detects
  )
  
  # keep only the detects
  Xdf  <- Xdf[detects == TRUE, ]
  Ndet <- nrow(Xdf)
  
  # initialize counter for number of iterations
  h    <- 1
  
  # Initialize the list to store results of each iteration
  mclist <- list()
  
  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around the detects in twa
    X <- data.frame(
      matrix(
        rnorm(M * Ndet, mean = Xdf$xk, sd = Xdf$CVt * Xdf$xk),
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
    #   4. UTL = GM * GSD^UT
    ##########################################################################################
    # Sort the data in ascending order to obtain the order statistics
    X <- as.data.frame(t(apply(X, 1, sort)))
    names(X) <- paste("y", c(1:Ndet), sep = "")
  
    # Blom's rankit
    # Pk = (k-a)/(n+1-2a) with a=3/8
    # zk = qnorm(Pk)
    zk <- qnorm((((N - Ndet + 1):N) - 3/8) / (N + 1/4))
  
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
    UTL   <- as.double(GMr * GSDr^UT)
  
    ##########################################################################################
    # Statistics of UTL
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    ##########################################################################################
    yest  <- mean(UTL)
    uy    <- sd(UTL)
    ycv   <- uy / yest
    Mpos  <- length(UTL)
    
    # index for coverage width in output vector
    q     <- as.integer(pmu * Mpos + 0.5)
    
    # width of the coverage interval
    UTL   <- as.double(sort(UTL, decreasing = FALSE))
    dy    <- UTL[(q + 1):Mpos] - UTL[1:(Mpos - q)]
    
    # minimum width of the coverage interval
    dymin <- min(dy)
    
    # find index of interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 7.7.2
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind) + 0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow  <- UTL[ind]
    yhigh <- UTL[ind + q]
    # associated p-values
    plow  <- (ind) / Mpos
    phigh <- (ind + q) / Mpos
    
    # classic symmetric confidence bounds
    ylowsym  <- quantile(
      x      = UTL,
      probs  = plowsym,
      names  = FALSE
    )
    yhighsym <- quantile(
      x      = UTL,
      probs  = phighsym,
      names  = FALSE
    )
    
    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y     = UTL,
      yest  = yest,
      uy    = uy,
      ycv   = ycv,
      ylow  = ylow,
      yhigh = yhigh,
      ylowsym  = ylowsym,
      yhighsym = yhighsym,
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
    yest_se  <- sd(sapply(mclist, `[[`, "yest" )) / sqrt(h)
    # ycv_se   <- sd(sapply(mclist, `[[`, "ycv"  )) / sqrt(h)
    # ycv is not in same order of magnitude as the rest -> use uy
    uy_se    <- sd(sapply(mclist, `[[`, "uy"   )) / sqrt(h)
    ylow_se  <- sd(sapply(mclist, `[[`, "ylow" )) / sqrt(h)
    yhigh_se <- sd(sapply(mclist, `[[`, "yhigh")) / sqrt(h)
    
    # Check if the tolerance condition is met
    if (2 * yest_se  <= delta &&
        2 * uy_se    <= delta &&
        2 * ylow_se  <= delta &&
        2 * yhigh_se <= delta) {
      break
    }
    
    # avoid runaway
    if (h >= 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    h <- h + 1
  }
  
  # endresult based on all the h*M generated samples
  UTL   <- unlist(lapply(mclist, `[[`, "Y"))
  UTL   <- as.double(sort(UTL, decreasing = FALSE))
  yest  <- mean(UTL)         # output value
  uy    <- sd(UTL)
  ycv   <- uy / yest         # output value
  Mpos  <- length(UTL)       # output value
  q     <- as.integer(pmu * Mpos + 0.5)
  dy    <- UTL[(q + 1):Mpos] - UTL[1:(Mpos - q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind) + 0.5)
  ylow  <- UTL[ind]          # output value
  yhigh <- UTL[ind + q]      # output value
  plow  <- (ind) / Mpos      # output value
  phigh <- (ind + q) / Mpos  # output value
  prob <- unique(sort(
    c(signif(c(plow   , phigh   ), digits = max(ndig, 3)),
      signif(c(plowsym, phighsym), digits = max(ndig, 3)), 
      c(0.000, 0.050, 0.500, 0.950, 1.000)
    )))
  yq       <- quantile(
    x      = UTL,
    probs  = prob,
    names  = FALSE
  )
  dUTL     <- density(UTL)
  dUTLsub  <- approx(dUTL$x, dUTL$y, xout = yq)$y
  yq <- data.frame(          # output value
    "p"    = prob,
    "d"    = signif(dUTLsub, digits = max(ndig, 3)),
    "UTL"  = signif(yq     , digits = ndig),
    row.names = format(prob)
  )
  ylowsym  <- yq[yq$p == signif(plowsym,  digits = max(ndig, 3)),]$UTL
  yhighsym <- yq[yq$p == signif(phighsym, digits = max(ndig, 3)),]$UTL
  
  return(
    list(
      yest  = signif(yest,  digits = ndig),
      ycv   = signif(ycv,   digits = ndig),
      ylow  = signif(ylow,  digits = ndig),
      yhigh = signif(yhigh, digits = ndig),
      plow  = signif(plow,  digits = max(ndig, 3)),
      phigh = signif(phigh, digits = max(ndig, 3)),
      ylowsym  = signif(ylowsym,  digits = ndig),
      yhighsym = signif(yhighsym, digits = ndig),
      plowsym  = signif(plowsym,  digits = max(ndig, 3)),
      phighsym = signif(phighsym, digits = max(ndig, 3)),
      yq    = yq,
      Y     = UTL,
      Mpos  = Mpos,
      ndig  = ndig,
      tolerance = delta
    )
  )
}
