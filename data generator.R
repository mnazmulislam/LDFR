
# --- READ ME 
#  title: R code for generating data that have the structure described in: 
#        "Longitudinal dynamic functional regression" (LDFR)
#        by Md Nazmul Islam, Ana-Maria Staicu, and Eric van Heugten
# created by Md Nazmul Islam
# date: 2017/11
# ---


###########################################################################################################
### ------------------------- defining inputs to generate datasets (Gaussian) ------------------------- ###
###########################################################################################################

#            tmin : minimum timepoint
#            tmax : maximum timepoint
#             rng : range of timepoints
#              TT : vector of 41 equidistant timepoints
#               I : Number of subjects
#               J : length of the vector "TT"
#               g : grid points over which functional predictor is observed at each time point 
#              ss : vector of equidistant points
#              SS : length of the vector "ss"
#           minJi : minimum number of repeated observations per subject
#           maxJi : maximum number of repeated observations per subject
#              mm : number of repeated observations per subject in the test set
#           delta : the intensity at which functional predictor evolves over time
#               A : seed number
#           simul : desired number of Monte carlo simulations 
#             nbf : number of knots

#          zetamu : mean of the variables that are used in representing the true functional predictor
#      zetasigma1 : variance of the first gaussian variable used in representing the true functional predictor
#      zetasigma2 : variance of the second gaussian variable used in representing the true functional predictor
#      zetasigma3 : variance of the third gaussian variable used in representing the true functional predictor
#      zetasigma4 : variance of the fourth gaussian variable used in representing the true functional predictor

#         mu.e.x1 : mean of the first noisy term of the observed functional predictor
#         mu.e.x2 : mean of the second noisy term of the observed functional predictor
#         mu.e.x3 : mean of the third noisy term of the observed functional predictor
#        var.e.x1 : variance the first noisy term of the observed functional predictor
#        var.e.x2 : variance the second noisy term of the observed functional predictor
#        var.e.x3 : variance the third noisy term of the observed functional predictor; determined by signal-to-noise ratio (SNR)

#       errYbi0mu : mean of the random intercept
#       errYbi1mu : mean of the random slope
#       errYbijmu : mean of the white noise term
#    errYbi0sigma : variance of the random intercept
#    errYbi1sigma : variance of the random slope
#   errYbi12sigma : covariance between random intercept and slope
#    errYbijsigma : variance of the white noise term
#          method : selection of smoothing parameter; "GCV.Cp", "REML", "ML" or others 
#          Yerror : type of dependence structure for responses; "CS" for compound symmetric type
#                   or "REM" for both subject-specific random intercept and slope. We illustarte
#                   the method using former here.
#       full_traj : logical argument; TRUE if full trajectory is desired; default is TRUE
#       Cov_Error : "High SNR" or "Low SNR"
#             pve : percentage of variance explained; defualt is 0.95





###################################################################################################
### ---------------------- defining functions to generate datasets (Gaussian) ----------------- ###
###################################################################################################



### -------------------------------------------------------------------------------- ###
### ------------------------- mean of functional predictor ------------------------- ###
### -------------------------------------------------------------------------------- ###

mufun <- function (ss, Tij) 
{ 
  mu <- 1 + 2 * ss + 3 * Tij + 4 * ss * Tij
  mu
}


### -------------------------------------------------------------------------------------------------------------------- ###
### ------------------------- true basis functions used in functional predictor representation ------------------------- ###
### -------------------------------------------------------------------------------------------------------------------- ###

phi_fun <- function (ss)
{
  phi1 <- sqrt(2) * cos (2 * pi * ss)          
  phi2 <- sqrt(2) * sin (2 * pi * ss)     
  return = list (phi1 = phi1, phi2 = phi2)
}

### ---------------------------------------------------------------------------------------------------- ###
### ------------------------- complete data (without sparsity) for one subject ------------------------- ###
### ---------------------------------------------------------------------------------------------------- ###

onesubj <- function (A, i, TT, ss, zetamu, zetasigma1, zetasigma2, zetasigma3, 
                     zetasigma4, delta, Cov_Error, mu.e.x1, mu.e.x2,
                     mu.e.x3, var.e.x1, var.e.x2, var.e.x3, 
                     errYbi0mu, errYbijmu, errYbi0sigma, 
                     errYbijsigma, errYbi1mu, errYbi1sigma, errYbi12sigma, J, Yerror)
{
  set.seed (A)
  J <- length(TT) 
  SS <- length(ss)
  Si.i11 <- rnorm (1, mean = zetamu, sd = sqrt(zetasigma1)) 
  Si.i12 <- rnorm (1, mean = zetamu, sd = sqrt(zetasigma2)) 
  Si.i21 <- rnorm(1, mean = zetamu, sd = sqrt(zetasigma3)) 
  Si.i22 <- rnorm(1, mean = zetamu, sd = sqrt(zetasigma4)) 
  zetai1 <- Si.i11 * cos(2 * pi * TT) + Si.i12 *  sin(2 * pi * TT) 
  zetai2 <- Si.i21 * cos(4 * pi * TT) + Si.i22 *  sin(4 * pi * TT) 
  
  MU.Xi <-  do.call ( rbind, lapply(TT, function(cc) as.vector(mufun(ss, Tij = cc ))))
  pf <-  phi_fun(ss)
  true.alpha <- 7 * sin(3 * pi * TT)
  true.beta1  <- exp (- TT * delta)   
  true.beta2 <- delta * TT * sin ( delta * TT)
  
  
  ### --------------------------------------------------------------------------------------- ###
  ### ------------------------- defining true functional predictors ------------------------- ###
  ### --------------------------------------------------------------------------------------- ###
  
  
  true.func.cov <- MU.Xi + do.call(rbind, lapply(zetai2, function (feb) 
    feb * pf$phi2)) + do.call(rbind, lapply(zetai1, function (jan) jan * pf$phi1))  
  
  
  ### ----------------------------------------------------------------------------------------- ###
  ### ------------------------- defining error structure for response ------------------------- ###
  ### ----------------------------------------------------------------------------------------- ###
  
  
 if (Yerror == "IS") {
   
  meanYij3 <- errYbi0mu + errYbi1mu * TT + errYbijmu 
  Z <- cbind(1, TT)
  D <- matrix( c(errYbi0sigma, errYbi12sigma, errYbi12sigma, errYbi1sigma ),2, 2)
  varYij3 <- (Z %*% D) %*% t(Z) + diag(errYbijsigma, J)
  Randi <-  mvrnorm(1, meanYij3, varYij3)
 } 
 else if (Yerror == "CS") {
   meanYij3 <- rep(errYbi0mu + errYbijmu, J)                                           
   varYij3 <- diag(errYbijsigma, J) + matrix(errYbi0sigma, nrow = J, ncol = J) 
   Randi <- mvrnorm(1, meanYij3, varYij3)
 }
  
  rng <- max(TT) -  min(TT) 
  
  ### ---------------------------------------------------------------------------------------------------- ###
  ### ------------------------------------------- Generating responses ----------------------------------- ###
  ### ---------------------------------------------------------------------------------------------------- ###
  
  Yi.r <-  apply(true.func.cov, 1, function(v){mean(v * pf$phi2 ) * rng }) * (true.beta2) + 
    apply(true.func.cov, 1, function(v){mean(v * pf$phi1 ) * rng }) * (true.beta1) + Randi + true.alpha
  
 
  ### ---------------------------------------------------------------------------------------------------- ###
  ### ------------------------- defining error structure for observed predictors ------------------------- ###
  ### ---------------------------------------------------------------------------------------------------- ###
 
  if ( Cov_Error == "High SNR" ) {
    
    SNR <- 2.5
    var.e.x3 <- (0.5 * (zetasigma1 + zetasigma2 + zetasigma3 + zetasigma4) - SNR * (var.e.x1 + var.e.x2)) / SNR  
    
    e1 <- rnorm(J, mean = mu.e.x1, sd = sqrt(var.e.x1))
    e2 <- rnorm(J, mean = mu.e.x2, sd = sqrt(var.e.x2))
    e3 <- matrix(rnorm(J * SS , mean = mu.e.x3, sd = sqrt(var.e.x3)), nrow = J, byrow = T)     
    error1 <-  do.call(rbind, lapply(e1, function (d) d * pf$phi1)) + do.call(rbind, lapply(e2, function (d) d *  pf$phi2)) 
    errorX <- error1 + e3                                                                                      
    obs.func.cov <- true.func.cov +  errorX 
    zetaiw1 <- zetai1 + apply(error1, 1, function(v){mean(v * pf$phi1) * rng })                                              
    zetaiw2 <- zetai2 + apply(error1, 1, function(v){mean(v * pf$phi2) * rng })
  }
  else if ( Cov_Error == "Low SNR" ) {
    
    SNR <- 0.5
    var.e.x3 <- (0.5 * (zetasigma1 + zetasigma2 + zetasigma3 + zetasigma4) - SNR * (var.e.x1 + var.e.x2)) / SNR  
    
    e1 <- rnorm(J, mean = mu.e.x1, sd = sqrt(var.e.x1))
    e2 <- rnorm(J, mean = mu.e.x2, sd = sqrt(var.e.x2))
    e3 <- matrix(rnorm(J * SS , mean = mu.e.x3, sd = sqrt(var.e.x3)), nrow = J, byrow = T)     
    error1 <-  do.call(rbind, lapply(e1, function (d) d * pf$phi1)) + do.call(rbind, lapply(e2, function (d) d *  pf$phi2)) 
    errorX <- error1 + e3                                                                                      
    obs.func.cov <- true.func.cov +  errorX 
    zetaiw1 <- zetai1 + apply(error1, 1, function(v){mean(v * pf$phi1) * rng })                                              
    zetaiw2 <- zetai2 + apply(error1, 1, function(v){mean(v * pf$phi2) * rng })
  }
  
  index <- c((1 + (i - 1) * J) : (i * J))
  
  data <- cbind(i = rep(i, J), j = seq_len( J ), v=seq_len( J ), Tij = TT, 
                Randi = Randi, zetai1 = zetai1, zetai2 = zetai2, true.alpha = true.alpha, 
                true.beta1 = true.beta1, true.beta2 = true.beta2, obs.func.cov = obs.func.cov, 
                index = index, Yij = Yi.r, zetaiw1 = zetaiw1, zetaiw2 = zetaiw2, MU.Xi = MU.Xi, 
                true.func.cov = true.func.cov)  
}


### ------------------------------------------------------------------------------------- ###
### ------------------------- functions creating sparse dataset ------------------------- ###
### ------------------------------------------------------------------------------------- ###


sparsefun <- function ( A, i, minJi, maxJi, TT, data )
{
  set.seed (A )
  J <- length(TT)
  mi <- sample ( minJi : maxJi, size = 1 )
  ith <- sort ( sample ( 1 : J, size = mi), decreasing = FALSE )
  ith.data <- data [ which ( data$i == i ), ]
  k1 <- ith.data [ ith, ]
  k1
}


### ------------------------------------------------------------------------- ###
### -------------------------- data generation ------------------------------ ###
### ------------------------------------------------------------------------- ###

data_func <- function (A , I, TT, ss, zetamu, zetasigma1, zetasigma2, 
                       zetasigma3, zetasigma4, delta, minJi, maxJi, Cov_Error, mu.e.x1, mu.e.x2, 
                       mu.e.x3, var.e.x1, var.e.x2, var.e.x3, errYbi0mu, errYbijmu, errYbi0sigma, 
                       errYbijsigma, errYbi1mu, errYbi1sigma, errYbi12sigma, mm, J, Yerror)
{
  set.seed (A) 
  n <- I
  
  ### ------------------------------------------------------------------------------------- ###
  ### ---------------------------- Full data set for all subjects ------------------------- ###
  ### ------------------------------------------------------------------------------------- ###
  
  dat.full <- data.frame(do.call(rbind, lapply(seq_len(I), function(i) 
    onesubj (A = i + A , i, TT, ss, zetamu, zetasigma1, zetasigma2, zetasigma3, 
             zetasigma4, delta, Cov_Error, mu.e.x1, mu.e.x2, mu.e.x3, var.e.x1, 
             var.e.x2, var.e.x3, errYbi0mu, errYbijmu, errYbi0sigma, errYbijsigma, errYbi1mu, errYbi1sigma, 
             errYbi12sigma, J, Yerror)))) 
  
  ### ------------------------------------------------------------------------------------- ###
  ### ------------------ Sparse data set generated from full data set --------------------- ###
  ### ------------------------------------------------------------------------------------- ###
  
  dat.sparse <- do.call(rbind, lapply(seq_len(I), function(i) 
    sparsefun(A = i + A , i, minJi, maxJi, TT, data = dat.full))) 
  

  ### -------------------------------------------------------------------------------------- ###
  ### ---- creating training and test set used for investigation of prediction accuracy ---- ###
  ### -------------------------------------------------------------------------------------- ###
  
  m1 <- sample(seq_len(I), n, replace = F );                                                     
  m2 <- as.matrix(aggregate(dat.sparse[, 2], by = list(dat.sparse[, 1]), function(a) 
  { 
    sample(a, mm)
  } )[m1, ])
  
  row <- matrix(0, nrow = n, ncol = (dim(m2)[2]-1))
  for (k in 1 : (dim(m2)[2]-1))  {
    row[, k] <- as.vector (do.call(rbind, lapply(1 : dim(m2)[1], 
                                                 function (a) which ( ( dat.sparse[,1] == m2[ a, 1]) * 
                                                                      ( dat.sparse[,2] == m2[ a, (k + 1)]) == 1))))    
  }
  
  row <- as.vector(row)
  dat.test <- dat.sparse[row,];  dat.train <- dat.sparse [- row, ]
  output = list ( dat.full = dat.full, dat.sparse = dat.sparse, 
                  dat.train = dat.train, row = row,
                  dat.test = dat.test , m2 = m2)  
}


