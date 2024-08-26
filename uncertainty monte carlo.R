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
#           as meaningful in a numerical value
#           ndig = 2  computing time is in the order of seconds
#           ndig = 3  computing time is in the order of minutes
# ueft    : upper exceedance fraction threshold for exposure
# pexp    : 1-sided coverage level for exposure
# pmu     : coverage level for measurement uncertainty
#           both 1-sided and 2-sided coverage intervals
#           will be reported
# 
# output variables
# yest     : estimate of Y, obtained as the average of 
#            the M model values yr from a Monte Carlo run
# cvy      : relative standard uncertainty associated 
#            with yest
# ylow2gum : left-hand endpoint of a 2-sided coverage interval 
#            for Y (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2)
#            with approximately the same density as yhigh2gum
# yhigh2gum: right-hand endpoint of a 2-sided coverage interval 
#            for Y (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2)
#            with approximately the same density as ylow2gum
# plow2gum : p-value for ylow2gum
#            (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2)
#            with approximately the same probability density as phigh2gum
# phigh2gum: p-value for yhigh2gum
#            (iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2)
#            with approximately the same probability density as plow2gum
# 
# ylow2sym : left-hand endpoint of a 2-sided coverage interval 
#            for Y 
# yhigh2sym: right-hand endpoint of a 2-sided coverage interval 
#            for Y 
# plow2sym : p-value for ylow2gum
#            with phigh2sym = 1 - plow2sym
# phigh2sym: p-value for yhigh2gum
#            with plow2sym = 1 - phigh2sym
# 
# ylow1    : left-hand endpoint of a 1-sided coverage interval 
#            for Y (Leidel, Busch and Lynch 1977, 4.1)
# yhigh1   : right-hand endpoint of a 1-sided coverage interval 
#            for Y (Leidel, Busch and Lynch 1977, 4.1)
# plow1    : p-value for ylow1
# phigh1   : p-value for yhigh1
# 
# yq       : characteristic quantiles (for research purposes)
#  p-values:
#  0.000                 │       minimum y-value of Monte Carlo trials
#        plow2gum  ══════╪═══╗
#  0.025 plow2sym  ──────┼───╫───┐
#  0.050 plow1     ──┐   │   ║   │  coverage intervals
#  0.500            95% 95% 95% 95% of measurement 
#  0.950 phigh1    ──┼───┘   ║   |  uncertainty
#        phigh2gum ══╪═══════╝   |
#  0.975 phigh2sym ──┼───────────┘
#  1.000             │           maximum y-value of Monte Carlo trials
# Mpos     : number of Monte Carlo trials with positive
#            values for Y (= UTL)
# ndig     : number of significant decimal digits regarded
#            as meaningful in the numerical value of UTL
#            - is same as the input value
##########################################################################################
# Y    : (scalar) output quantity, regarded as a random variable
#        = UTL
#        = 70% one-sided upper tolerance limit of Q95 (UTL95,70)
#        ("confidence limit" is used for the mean)
##########################################################################################
if (!requireNamespace("wzMisc", quietly = TRUE)) {
  devtools::install_github("slin30/wzMisc")
}
library(wzMisc)
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

utl.mc <- function(twa, CVt = rep(0, length(twa)), ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95) {
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008
  # 7.3   Sampling from probability distributions
  # 7.4   Evaluation of the model
  # 7.5   Discrete representation of the distribution function for the output quantity
  # 7.6   Estimate of the output quantity and the associated standard uncertainty
  # 7.7   Coverage interval for the output quantity
  # 7.9   Adaptive Monte Carlo procedure
  # 7.9.2 Numerical tolerance associated with a numerical value
  # 7.9.4 Adaptive procedure
  
  # number of significant decimal digits
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 a)
  # ndig = appropriate small positive integer (input value)
  
  # number of Monte Carlo trials
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 b)
  M    <- max(100 / (1 - pmu), 1e4)
  
  # number of input quantities
  N    <- length(twa)
  
  # coverage factor for quantiles
  # coefficient of model
  # put outside of loop to reduce calculation time
  UT   <- qt(p   = pexp,
             df  = N - 1,
             ncp = sqrt(N) * qnorm(1 - ueft)
             ) / sqrt(N)
  
  # Initialize the list to store results of each iteration
  mclist <- list()
  
  # initialize counter for number of iterations
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 c)
  h    <- 1
  
  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around twa
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.3 Sampling from probability distributions
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
    # store in X to reduce memory usage
    X <- X[apply(X, 1, function(row) all(row > 0)), ]
    
    # Transform to the normally distributed space
    # store log(Xpos) to reduce calculation time a bit
    # store in X to reduce memory usage
    X <- log(X)
    
    ##########################################################################################
    # Propagation of the input distributions through the model
    # model: calculate by row
    #   1. GM and GSD
    #   2. Y = UTL = GM * GSD^UT
    ##########################################################################################
    # Propagate the input distributions through the model ------- apply()
    # Transform back to the lognormally distributed space and --- exp()
    # Finalize propagation through the model -------------------- Y
    # The model Y=UTL=GM*GSD^UT is iaw EN 689 Annex F
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.4 Evaluation of the model
    # GM <- exp(rowMeans(X))
    GM  <- exp(apply(X, 1, mean))
    GSD <- exp(apply(X, 1, sd  ))
    Y   <- as.double(GM * GSD^UT)
    
    ##########################################################################################
    # Statistics of Y
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 e)
    ##########################################################################################
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.6
    # estimate y of Y and the standard uncertainty u(y) associated with y
    yest  <- mean(Y)
    uy    <- sd(Y)
    cvy   <- uy / yest
    
    Mpos  <- length(Y)
    
    # sort the model values provided by MCM into non-decreasing order
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 e), 7.5.1 a)
    # discrete representation of the distribution function for the output quantity Y
    Y     <- as.double(sort(Y, decreasing = FALSE))

    # index for coverage width in output vector
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    q     <- as.integer(pmu * Mpos + 0.5)
    
    # width of the coverage interval
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    dy    <- Y[(q + 1):Mpos] - Y[1:(Mpos - q)]
    
    # minimum width of the coverage interval
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    dymin <- min(dy)
    
    # find index of 2-sided interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 7.7.2
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind) + 0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow2gum  <- Y[ind]
    yhigh2gum <- Y[ind + q]
    # associated p-values
    plow2gum  <- (ind) / Mpos
    phigh2gum <- (ind + q) / Mpos

    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y         = Y,
      yest      = yest,
      uy        = uy,
      cvy       = cvy,
      ylow2gum  = ylow2gum,
      yhigh2gum = yhigh2gum,
      Mpos      = Mpos
    )))
    
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 f)
    if (h < 2) {
      h <- h + 1
      next
    }
    
    # standard error of the h values of yest, uy, ylow2gum, yhigh2gum
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 g)
    yest_se      <- sd(sapply(mclist, `[[`, "yest"     )) / sqrt(h)

    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 h)
    uy_se        <- sd(sapply(mclist, `[[`, "uy"       )) / sqrt(h)
    ylow2gum_se  <- sd(sapply(mclist, `[[`, "ylow2gum" )) / sqrt(h)
    yhigh2gum_se <- sd(sapply(mclist, `[[`, "yhigh2gum")) / sqrt(h)
    
    # standard deviation of h*M values of Y
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 i)
    sd_y <- sd(unlist(lapply(mclist, `[[`, "Y")))

    # number of decimal digits
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.2 a), 7.9.4 j)
    ndecdig <- -(ceiling(log10(sd_y))-ndig)
    
    # numerical tolerance, delta
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.2 b), 7.9.4 j)
    delta    <- 0.5 * 10^(-ndecdig)

    # Check if the tolerance condition is met
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 k)
    if (2 * yest_se      <= delta &&
        2 * uy_se        <= delta &&
        2 * ylow2gum_se  <= delta &&
        2 * yhigh2gum_se <= delta) {
      break
    }
    
    # avoid runaway
    if (h >= 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 k)
    h <- h + 1
  }

  # endresult based on all the h*M generated samples
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 l)
  Y     <- unlist(lapply(mclist, `[[`, "Y"))
  Y     <- as.double(sort(Y, decreasing = FALSE))
  yest  <- mean(Y)           # output value
  uy    <- sd(Y)
  cvy   <- uy / yest         # output value
  Mpos  <- length(Y)         # output value
  q     <- as.integer(pmu * Mpos + 0.5)
  dy    <- Y[(q + 1):Mpos] - Y[1:(Mpos - q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind) + 0.5)
  ylow2gum  <- Y[ind]            # output value
  yhigh2gum <- Y[ind + q]        # output value
  plow2gum  <- (ind) / Mpos      # output value
  phigh2gum <- (ind + q) / Mpos  # output value
  plow1  <- 1 - pmu              # output value
  phigh1 <- pmu                  # output value

  # in addition to ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7
  # classic symmetric confidence bounds
  plow2sym  <- (1 - pmu) / 2     # output value
  phigh2sym <- (1 + pmu) / 2     # output value = 1 - (1 - pmu) / 2

  prob <- unique(sort(
    c(signif(c(plow2gum, phigh2gum), digits = max(ndig, 3)),
      signif(c(plow1,    phigh1   ), digits = max(ndig, 3)), 
      signif(c(plow2sym, phigh2sym), digits = max(ndig, 3)), 
      c(0.000, 0.025, 0.050, 0.500, 0.950, 0.975, 1.000)
    )))
  yq   <- quantile(
    x      = Y,
    probs  = prob,
    names  = FALSE
  )
  dY       <- density(Y)
  dYsub    <- approx(dY$x, dY$y, xout = yq)$y
  ndecdig_density <- -(ceiling(log10(min(dYsub)))-ndig)
  if(delta == 0) {
    ndecdig <- max(count_sigfigs(twa - floor(twa)))
  }
  yq <- data.frame(              # output value
    "p"    = prob,
    "d"    = round(dYsub, digits = ndecdig_density),
    "UTL"  = round(yq   , digits = ndecdig        ),
    row.names = format(prob)
  )
  
  # 2-sided confidence interval       [ylow1, yhigh1] confidence level pmu
  ylow2sym  <- yq[yq$p == signif(plow2sym,  digits = max(ndig, 3)),]$UTL
  yhigh2sym <- yq[yq$p == signif(phigh2sym, digits = max(ndig, 3)),]$UTL
  
  # Lower 1-sided confidence interval [ylow1, +∞[     confidence level pmu
  ylow1     <- yq[yq$p == signif(plow1,     digits = max(ndig, 3)),]$UTL

  # Upper 1-sided confidence interval [0, yhigh1]     confidence level pmu
  yhigh1    <- yq[yq$p == signif(phigh1,    digits = max(ndig, 3)),]$UTL
  
  return(
    list(
      yest      = round (yest,      digits = ndecdig),
      cvy       = round (cvy,       digits = ndecdig),
      ylow2gum  = round (ylow2gum,  digits = ndecdig),
      yhigh2gum = round (yhigh2gum, digits = ndecdig),
      plow2gum  = signif(plow2gum,  digits = max(ndig, 3)),
      phigh2gum = signif(phigh2gum, digits = max(ndig, 3)),
      ylow2sym  = round (ylow2sym,  digits = ndecdig),
      yhigh2sym = round (yhigh2sym, digits = ndecdig),
      plow2sym  = signif(plow2sym,  digits = max(ndig, 3)),
      phigh2sym = signif(phigh2sym, digits = max(ndig, 3)),
      ylow1     = round (ylow1,     digits = ndecdig),
      yhigh1    = round (yhigh1,    digits = ndecdig),
      plow1     = signif(plow1,     digits = max(ndig, 3)),
      phigh1    = signif(phigh1,    digits = max(ndig, 3)),
      yq        = yq,
      Y         = Y,
      Mpos      = Mpos,
      ndig      = ndig,
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
  
utl.ros.mc <- function(twa, detects, CVt = rep(0, length(twa)), ndig = 2, ueft = 0.05, pexp = 0.70, pmu = 0.95) {
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008
  # 7.3   Sampling from probability distributions
  # 7.4   Evaluation of the model
  # 7.5   Discrete representation of the distribution function for the output quantity
  # 7.6   Estimate of the output quantity and the associated standard uncertainty
  # 7.7   Coverage interval for the output quantity
  # 7.9   Adaptive Monte Carlo procedure
  # 7.9.2 Numerical tolerance associated with a numerical value
  # 7.9.4 Adaptive procedure
  
  # number of significant decimal digits
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 a)
  # ndig = appropriate small positive integer (input value)
  
  # number of Monte Carlo trials
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 b)
  M    <- max(100 / (1 - pmu), 1e4)
  
  # number of input quantities
  N    <- length(twa)
  
  # coverage factor for quantiles
  # coefficient of model
  # put outside of loop to reduce calculation time
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

  # number of detects
  Ndet <- nrow(Xdf)
  
  # Initialize the list to store results of each iteration
  mclist <- list()
  
  # initialize counter for number of iterations
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 c)
  h    <- 1
  
  repeat {
    ##########################################################################################
    # Generation of input distributions
    ##########################################################################################
    # generate M times normal distribution around the detects in twa
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.3 Sampling from probability distributions
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
    # store in X to reduce memory usage
    X <- X[apply(X, 1, function(row) all(row > 0)), ]
    
    # Transform to the normally distributed space
    # store log(Xpos) to reduce calculation time a bit
    # store in X to reduce memory usage
    X <- log(X)
    
    ##########################################################################################
    # Propagation of the input distributions through the model
    # model: calculate by row
    #   1. sort rows in ascending order
    #   2. linear regression of each row
    #   3. GM = intercept and GSD = slope
    #   4. Y = UTL = GM * GSD^UT
    ##########################################################################################
    # Sort the data in ascending order to obtain the order statistics
    X <- as.data.frame(t(apply(X, 1, sort)))
    names(X) <- paste("y", c(1:Ndet), sep = "")
  
    # Blom's rankit
    # Pk = (k-a)/(n+1-2a) with a=3/8
    # zk = qnorm(Pk)
    zk <- qnorm((((N - Ndet + 1):N) - 3/8) / (N + 1/4))
  
    # Matrix Formulation of the Multiple Linear Regression (MLR) Model
    # Linear regression iaw EN 689 Annex H
    # For the regression, only the detects are used
    # as the operation is performed in the normally distributed space,
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
    # Propagate the input distributions through the model ------- MLR()
    MLRcoefs <- MLR(X, zk)
    
    # Transform back to the lognormally distributed space and --- exp()
    # Finalize propagation through the model -------------------- Y
    # The model Y=UTL=GM*GSD^UT is iaw EN 689 Annex F and H
    # GMr   <- exp(MLRcoefs[,1])
    # GSDr  <- exp(MLRcoefs[,2])
    GMr   <- exp(MLRcoefs$AMr)
    GSDr  <- exp(MLRcoefs$ASDr)
    Y     <- as.double(GMr * GSDr^UT)
  
    ##########################################################################################
    # Statistics of Y
    #   1. mean
    #   2. CV
    #   3. lower and upper bounds (endpoints) of coverage interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 e)
    ##########################################################################################
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.6
    # estimate y of Y and the standard uncertainty u(y) associated with y
    yest  <- mean(Y)
    uy    <- sd(Y)
    cvy   <- uy / yest

    Mpos  <- length(Y)
    
    # sort the model values provided by MCM into non-decreasing order
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 e), 7.5.1 a)
    # discrete representation of the distribution function for the output quantity Y
    Y     <- as.double(sort(Y, decreasing = FALSE))

    # index for coverage width in output vector
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    q     <- as.integer(pmu * Mpos + 0.5)
    
    # width of the coverage interval
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    dy    <- Y[(q + 1):Mpos] - Y[1:(Mpos - q)]
    
    # minimum width of the coverage interval
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7.2
    dymin <- min(dy)
    
    # find index of 2-sided interval with minimum width
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008 7.7.2
    ind   <- which(dy == dymin)
    ind   <- as.integer(mean(ind) + 0.5)
    # left-hand and right-hand endpoints of coverage interval
    ylow2gum  <- Y[ind]
    yhigh2gum <- Y[ind + q]
    # associated p-values
    plow2gum  <- (ind) / Mpos
    phigh2gum <- (ind + q) / Mpos
    
    # add results of iteration to list
    mclist <- append(mclist, list(list(
      Y         = Y,
      yest      = yest,
      uy        = uy,
      cvy       = cvy,
      ylow2gum  = ylow2gum,
      yhigh2gum = yhigh2gum,
      Mpos      = Mpos
    )))
    
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 f)
    if (h < 2) {
      h <- h + 1
      next
    }

    # standard error of the h values of yest, uy, ylow2gum, yhigh2gum
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 g)
    yest_se      <- sd(sapply(mclist, `[[`, "yest"     )) / sqrt(h)

    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 h)
    uy_se        <- sd(sapply(mclist, `[[`, "uy"       )) / sqrt(h)
    ylow2gum_se  <- sd(sapply(mclist, `[[`, "ylow2gum" )) / sqrt(h)
    yhigh2gum_se <- sd(sapply(mclist, `[[`, "yhigh2gum")) / sqrt(h)
    
    # standard deviation of h*M values of Y
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 i)
    sd_y <- sd(unlist(lapply(mclist, `[[`, "Y")))
    
    # number of decimal digits
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.2 a), 7.9.4 j)
    ndecdig <- -(ceiling(log10(sd_y))-ndig)
    
    # numerical tolerance, delta
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.2 b), 7.9.4 j)
    delta <- 0.5 * 10^(-ndecdig)
    
    # Check if the tolerance condition is met
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 k)
    if (2 * yest_se      <= delta &&
        2 * uy_se        <= delta &&
        2 * ylow2gum_se  <= delta &&
        2 * yhigh2gum_se <= delta) {
      break
    }
    
    # avoid runaway
    if (h >= 2000) {
      break
    }
    
    # Increment h if tolerance is not met, and repeat
    # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 k)
    h <- h + 1
  }
  
  # endresult based on all the h*M generated samples
  # iaw ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.9.4 l)
  Y     <- unlist(lapply(mclist, `[[`, "Y"))
  Y     <- as.double(sort(Y, decreasing = FALSE))
  yest  <- mean(Y)           # output value
  uy    <- sd(Y)
  cvy   <- uy / yest         # output value
  Mpos  <- length(Y)         # output value
  q     <- as.integer(pmu * Mpos + 0.5)
  dy    <- Y[(q + 1):Mpos] - Y[1:(Mpos - q)]
  dymin <- min(dy)
  ind   <- which(dy == dymin)
  ind   <- as.integer(mean(ind) + 0.5)
  ylow2gum  <- Y[ind]            # output value
  yhigh2gum <- Y[ind + q]        # output value
  plow2gum  <- (ind) / Mpos      # output value
  phigh2gum <- (ind + q) / Mpos  # output value
  plow1  <- 1 - pmu              # output value
  phigh1 <- pmu                  # output value

  # in addition to ISO/IEC GUIDE 98-3/Suppl.1:2008, 7.7
  # classic symmetric confidence bounds
  plow2sym  <- (1 - pmu) / 2     # output value
  phigh2sym <- (1 + pmu) / 2     # output value = 1 - (1 - pmu) / 2

  prob <- unique(sort(
    c(signif(c(plow2gum, phigh2gum), digits = max(ndig, 3)),
      signif(c(plow1,    phigh1   ), digits = max(ndig, 3)), 
      signif(c(plow2sym, phigh2sym), digits = max(ndig, 3)),
      c(0.000, 0.025, 0.050, 0.500, 0.950, 0.975, 1.000)
    )))
  yq   <- quantile(
    x      = Y,
    probs  = prob,
    names  = FALSE
  )
  dY       <- density(Y)
  dYsub    <- approx(dY$x, dY$y, xout = yq)$y
  ndecdig_density <- -(ceiling(log10(min(dYsub)))-ndig)
  if(delta == 0) {
    ndecdig <- max(count_sigfigs(twa - floor(twa)))
  }
  yq <- data.frame(              # output value
    "p"    = prob,
    "d"    = round(dYsub, digits = ndecdig_density),
    "UTL"  = round(yq   , digits = ndecdig        ),
    row.names = format(prob)
  )
  
  # 2-sided confidence interval       [ylow1, yhigh1] confidence level pmu
  ylow2sym  <- yq[yq$p == signif(plow2sym,  digits = max(ndig, 3)),]$UTL
  yhigh2sym <- yq[yq$p == signif(phigh2sym, digits = max(ndig, 3)),]$UTL
  
  # Lower 1-sided confidence interval [ylow1, +∞[     confidence level pmu
  ylow1     <- yq[yq$p == signif(plow1,     digits = max(ndig, 3)),]$UTL
  
  # Upper 1-sided confidence interval [0, yhigh1]     confidence level pmu
  yhigh1    <- yq[yq$p == signif(phigh1,    digits = max(ndig, 3)),]$UTL
  
  return(
    list(
      yest      = round (yest,      digits = ndecdig),
      cvy       = round (cvy,       digits = ndecdig),
      ylow2gum  = round (ylow2gum,  digits = ndecdig),
      yhigh2gum = round (yhigh2gum, digits = ndecdig),
      plow2gum  = signif(plow2gum,  digits = max(ndig, 3)),
      phigh2gum = signif(phigh2gum, digits = max(ndig, 3)),
      ylow2sym  = round (ylow2sym,  digits = ndecdig),
      yhigh2sym = round (yhigh2sym, digits = ndecdig),
      plow2sym  = signif(plow2sym,  digits = max(ndig, 3)),
      phigh2sym = signif(phigh2sym, digits = max(ndig, 3)),
      ylow1     = round (ylow1,     digits = ndecdig),
      yhigh1    = round (yhigh1,    digits = ndecdig),
      plow1     = signif(plow1,     digits = max(ndig, 3)),
      phigh1    = signif(phigh1,    digits = max(ndig, 3)),
      yq        = yq,
      Y         = Y,
      Mpos      = Mpos,
      ndig      = ndig,
      tolerance = delta
    )
  )
}
