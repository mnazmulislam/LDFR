
# --- READ ME 
#  title: Main function for fitting the model described in: 
#        "Longitudinal dynamic functional regression" (LDFR)
#        by Md Nazmul Islam, Ana-Maria Staicu, and Eric van Heugten
# created by Md Nazmul Islam
# date: 2017/11
# ---
  

############################################################################################
### ------------------------- defining inputs for fitting LDFR ------------------------- ###
############################################################################################


#                I : total number of subjects
#               IW : total number of new + existing subjects; new subject is (IW - I). 
#                     Note that profile information must be available for each IW.  
#     obs.func.cov : observed functional predictor; if IW is NOT null then obs.func.cov is row-concatenation of profiles
#                    for new and existing subjects. If IW is null, then obs.func.cov consists of profiles only for existing subject
#              Tij : sparse/dense time-points corresponding to each response; if IW is NOT null, then Tij is the vector 
#                    of time-points at which profiles are recorded for both existing and new subjects.
#                    Otherwise, Tij contains time-points only for existing subjects
#              Yij : scalar responses at observed Tij's for existing subjects
#               ID : subject id for each response; if IW is NOT null, then ID is the vector of 
#                    subject id for both existing and new subjects
#                V : position of each response within each subject in the domain of time; 
#                    & new subjects. Otherwise, this contains id's only for existing subjects 
#                    e.g. if t_i = (0.20, 0.60) for the i-th subject, where the domain of time is 
#                    {0.0, 0.20, 0.40, 0.60, 0.80, 1.0}, then V_i = c(2, 4) 
#         out_subj : indexing variable that determines the positions of each response across existing subjects;
#                     e.g. if there are 2 subjects which have position V_1 = c(2, 4) and V_2 = c(3, 5). Then
#                     out_subj = c(2, 4, 9, 11)
#              nbf : number of basis functions used to represent temporal non-parametric smooth function
#           method : smoothing parameter selection criteria
#              pve : pre-specified proportion of variance explained
#           Yerror : type of dependence structure for responses; "CS" for subject-specific random intercept
#                    or "IS" for both subject-specific random intercept and slope 
#           family : exponential family distribution
#    pred_interval : TRUE, if prediction interval is needed otherwise FALSE    
#            alpha : level of significance at which prediction interval is desired; e.g. 0.05, 0.10, ... etc. 


####################################################################################
### ---------------------------- defining outputs ------------------------------ ###
####################################################################################



#             fit : an object of bam() or gam() fit
#            pred : predcited responses for observed time-points
#  full_pred_obs  : predcited responses for all time-points; both observed and unobserved time points for observed subject
# full_pred_PA_un : predcited responses for all time-points for observed subject at population-level
#           error : prediction error at observed time points for observed subjects
#           gamma : estimates of functional coefficient
#             Fnp : Number of fPCs; i.e. k = 1, ..., Fnp
#             phi : empirical basis functions
#             xiF : full trajectory of time-varying coefficients for each k
#         evalues : empirical eigenvalues
#             xi0 : time-varying coefficients for each k at observed time-points
#         ID_full : subject number for observed and unobserved time-points
#          ID_obs :  subject number for observed and unobserved time-points
#        SS_inter : prediction of subject-specific random intercept effect
#          SS_slp : prediction of subject-specific random slope effect if Yerror = "IS", otherwise NULL
#      model_adeq : AIC, BIC, Rsq, and deviance explained for model selection
#           fMU.X : mean function of functional predictor for observed and unobserved time-points
#            Mean : mean of functional predictor for observed time points
#       zetaw.hat : raw time-varying loadings (proxy)
#            beta : estimates of dynamic parameters
#           beta0 : estimates of smooth intercept
#        RMPE_trj : Root-mean-prediction-error of response trajectory (If full response is available; otherwise outputs 9999999 as NOT available)
#         boundYF : Prediction band for observed response trajectory
#         boundXF : Prediction band for unobserved response trajectory
#          bounds : Confidence interval for smooth estimates
#       fpca.fit1 : fpca.face fitted object which primarily gives raw scores
#       fpca.fit2 : fpca.sc fitted object which is used on estimating full trajectory of scores

#       NOTE : Model can be fitted with either lme function from nlme package or gam from mgcv R-package.
#              We illustarte the method using the latter approach using mgcv package here. 


LDFR <- function(ID, obs.func.cov, Tij, Yij, V, I, IW, nbf, method, pve, Yerror, family, 
                 pred_interval, alpha, out_subj)  

{  
  
  TT <- sort(unique(Tij)); J <- length(TT)
  ss <- seq_len(dim(obs.func.cov)[2]); 
  SS <- length(ss); ln <- length(Tij)  
  
  if(length(IW) != 0){
  freq <- unlist(lapply(seq_len(IW), function(ii) length(which(ID == ii))))
  }
  
  if(length(IW) == 0){
    freq <- unlist(lapply(seq_len(I), function(ii) length(which(ID == ii))))
  }
  
  
  ### ---------------------------------------------------------------------------------- ###
  ### ------------------------ estimating mean function -------------------------------- ###
  ### ------------ fbps() is used to estimate; and gam() is used to predict ------------ ###
  ### ---------------------------------------------------------------------------------- ###  
  
  
  fbps.fit <-  fbps(data = obs.func.cov, covariates = c(Tij, ss))
  
  muhat0 <- matrix(0, nrow = J, ncol = SS)
  muhat.dat <- do.call(rbind, lapply(seq_len(J), function(de) { 
    t1 <- TT[de]
    r1 <- as.vector(which(Tij == t1))[1]                                                
    muhat0[de,] <- fbps.fit$Yhat[r1, ]}))
  
  newmu.dat <- data.frame (y = as.vector(t(muhat.dat)), x1 = rep(ss, J), x2 = rep(TT, SS)) 
  MU.gam <- gam (y ~ te(x1, x2, k = c(10, 5), bs ="cr"), data = newmu.dat, method = method)                         
  MU.gam.dat <- data.frame(x1 = rep(ss, ln), x2 = rep(Tij, SS) )                                    
  MU.X <- matrix(predict.gam(MU.gam, newdata = MU.gam.dat), nrow = ln, ncol = SS, byrow = T)  
  

  ### ----------------------------------------------------------------------------------------- ###
  ### ------------------------------ demeaning functional predictor --------------------------- ###
  ### ----------------------------------------------------------------------------------------- ###

  
  deM.X <- obs.func.cov - fbps.fit$Yhat        
  
    
  ### ----------------------------------------------------------------------------------------- ###
  ### ----------------  finding the eigen-components of covariance  operator ------------------ ###
  ### ----------------------------------------------------------------------------------------- ###  
  
  
  covsmooth <- fpca.face (Y = deM.X, pve = pve, var =  TRUE)                                                                       
  phi.mat <- covsmooth$efunctions                                 
  phi.hat <- phi.mat * sqrt(SS)                                                                         
  K.hat <- covsmooth$npc;                                         
  
  zetaw.hat <- covsmooth$scores / sqrt (SS)                    
  mu <- covsmooth$mu                                              
  eigenvalues <- covsmooth$evalues / SS                    
  
 
  
  ### ----------------------------------------------------------------------------------------- ###
  ### ------------------------------- Replacing negative eigenvalues -------------------------- ###
  ### ----------------------------------------------------------------------------------------- ###  
  
  
  eigenvalues[which(eigenvalues < 0)] <- 0                      
  
  
  ### ----------------------------------------------------------------------------------------- ###
  ### ----------------------------- checking norm --------------------------------------------- ###
  ### ----------------------------------------------------------------------------------------- ###  
  
  
  k1 <- all.equal ( as.vector ( do.call( rbind, lapply ( seq_len(K.hat), function (feb) 
    sqrt(mean ( length(ss) * phi.mat[,feb] ^ 2 )) ) ) ), rep(1, K.hat) )
  if ( k1 != "TRUE" ) stop ( "l2 norm mismatches")
   
  
  ### ----------------------------------------------------------------------------------------- ###
  ### ----------------- predicting the time-varying loadings (un)observed) -------------------- ###
  ### ----------------------------------------------------------------------------------------- ###  
  
  if(length(IW) == 0){
    IW <- I
  }
    
  xi.hat0 <- list()                                                                                    
    for( k in seq_len (K.hat) ) 
    {  
      xihat0.vec <- zetaw.hat[, k] 
      xi.hat0[[ k ]]<- t( sapply( seq_len(IW), function(i) { 
        xi.subj <- matrix(nrow = 1, ncol = J) 
        xi.subj[V[which(ID == i)]] <- xihat0.vec[which(ID == i)]
        return(xi.subj)
      }))
    }                                                                                                    
  fit2 <- lapply(xi.hat0, function(a) try ( fpca.sc(Y = a, pve = pve, var = TRUE) )) 
  xi.hat<- lapply(fit2, function(a) a$Yhat)    
    
  zetaw.hatF <- do.call(cbind, lapply(xi.hat, function(gg) as.vector(t(gg))))  
  
  
  ### ----------------------------------------------------------------------------------------------------------------- ###
  ### ----------------- generating variables for ID and Time for all subjects at all J time-points -------------------- ###
  ### ----------------------------------------------------------------------------------------------------------------- ###  

  
  IDF <- unlist(lapply(seq_len(IW), function(gh) rep(gh, J)))
  TijF <- rep(TT, IW)
  


  ### ----------------------------------------------------------------------------------------------------------------------- ###
  ### ----------------- If IW == I (implying no new subjects and do analysis on the existing subjects ) -------------------- ###
  ### ---------------------------------------------------------------------------------------------------------------------- ###  
  

  if(IW == I){
    
    Ind <- seq_len(dim(zetaw.hatF)[1]) 
  
    Index <- unlist(lapply(seq_len(IW), function(ui) Ind[which(IDF == ui)][V[which(ID == ui)]]))  # index for time-points at which profiles are available          
    
    zetaw.hatN <- do.call (cbind, lapply(xi.hat, function(gg) as.vector(t(gg))[Index]))           # extracting scores on those time-points             
    zetaw.hatE <- zetaw.hatF
    
    Yij <- Yij
    ID <- ID
    
    TijFN <- TijF
    IDFN <- IDF
  }
  

  ### ---------------------------------------------------------------------------------------------------------------------- ###
  ### ---- If IW != I (implying there exists new subjects; use existing subjects for estimation of model parameters) ---- ###
  ### ---------------------------------------------------------------------------------------------------------------------- ###  
  
  if(IW != I){
    
    index0 <- seq_len(I * J)
    index1 <-  out_subj 
    index2 <- seq_len(IW * J)[- c(index0)]
    
    zetaw.hatE <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg))[index0] ))   # predicted scores for observed subjects at each time-points
    
    zetaw.hatN <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg))[index1] ))   # predicted scores for observed subjects at observed time points
    
    zetaw.hatX <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg))[index2] ))   # redicted scores for unobserved subjects at each time-points
    

    Yij <- Yij
    V <- V[-which(ID > I)] 
    Tij <- Tij[-which(ID > I)]
    ID <- ID[-which(ID > I)]
    
      
    if(length(Yij) != length(Tij)) stop("number of observed responses for observed subjects does not match with number of observed time-points")
  
    TijFN <- TijF[- c(which(IDF > I))]
    TijFX <- TijF[ c(which(IDF > I))]                       

    IDFN <- IDF[- c(which(IDF > I))]          
    IDFX <- IDF[ c(which(IDF > I))]           
}
  
  

  ### -------------------------------------------------------------------------------- ###
  ### --------------- creating design matrix for SS random effects  ------------------ ###
  ### -------------------------------------------------------------------------------- ###  
  
  
  if (Yerror == "CS") { 
    
   # -------------------- Observed subjects  + Observed time-points  --------------------- #      
    Z1 <- matrix(0, nrow = length(Yij), ncol = length(unique(ID)))    
    for (i in 1 : I) {
      Z1[which(ID == unique(ID)[i]), i] = 1                     
    }
    fZ <- cbind(Z1)
    
    
    # -------------------- Observed subjects  + All time-points  -------------------------- #
    Z1F <- matrix(0, nrow = length(TijFN), ncol = length(unique(IDFN)))    
    for (i in 1 : I) {
      Z1F[which(IDFN == sort(unique(IDFN))[i]), i] = 1                     
    }
    fZF <- cbind(Z1F)
    
    
    if(IW != I){
    # ---------------------- Unobserved subjects  + All time-points  --------------------- #
    Z1FX <- matrix(0, nrow = length(TijFX), ncol = length(unique(IDFX)))    
    for (i in 1 : (IW - I)) {
      Z1FX[which(IDFX == sort(unique(IDFX))[i]), i] = 1                     
    }
    fZFX <- cbind(Z1FX)
  }
  }
  
  
  if (Yerror == "IS") {   
    
    # -------------------- Observed subjects  + Observed time-points  --------------------- #    
    Z1 = matrix(0, nrow = length(Yij), ncol = length(unique(ID)))          
    Z2 = matrix(0, nrow = length(Yij), ncol = length(unique(ID)))
    
    for (i in 1 : I) {
      Z1[which(ID == unique(ID)[i]), i] = 1  
      Z2[which(ID == unique(ID)[i]), i] = Tij[which(ID == unique(ID)[i])]
    }
    fZ <- cbind(Z1, Z2)
    
    
    # -------------------- Observed subjects  + All time-points  -------------------------- #
    Z1F <- matrix(0, nrow = length(TijFN), ncol = length(unique(IDFN))) 
    Z2F <- matrix(0, nrow = length(TijFN), ncol = length(unique(IDFN))) 
    
    for (i in 1 : I) {
      Z1F[which(IDFN == unique(IDFN)[i]), i] = 1  
      Z2F[which(IDFN == unique(IDFN)[i]), i] = TijFN[which(IDFN == unique(IDFN)[i])]
    }
    fZF <- cbind(Z1F, Z2F)
  
    
    if(IW != I){
    # ---------------------- Unobserved subjects  + All time-points  --------------------- #
    Z1FX <- matrix(0, nrow = length(TijFX), ncol = length(unique(IDFX)))    
    Z2FX <- matrix(0, nrow = length(TijFX), ncol = length(unique(IDFX)))    
    for (i in 1 : (IW - I)) {
      Z1FX[which(IDFX == sort(unique(IDFX))[i]), i] = 1  
      Z2FX[which(IDFX == sort(unique(IDFX))[i]), i] = TijFX[which(IDFX == sort(unique(IDFX))[i])]
    }
    fZFX <- cbind(Z1FX, Z2FX)
    }
  }


  ### ----------------------------------------------------------------------------------------------------------- ###
  ### ----  creating design matrix smooth intercept using truncated linear splines in mixed model framework  ---- ###
  ### ----------------------------------------------------------------------------------------------------------- ###  
  
  R <- nbf
  qtiles <- seq(0, 1, length = R)[-c(1, R)]
  knots <- quantile(TT, qtiles)
  
  
  # ------------------------ Observed subjects  + Observed time-points  ------------------------------- #    
  Z0 <- matrix(as.vector(cbind(sapply(knots, function(tau) ((Tij > tau) * (Tij - tau) )))), nrow = length(Tij))
  fZ0 <- cbind(1, Tij, Z0)
  bigX <- cbind(fZ, fZ0)
  
  
  # ----------------------------- Observed subjects  + All time-points  ------------------------------ #
  Z0F <- matrix(as.vector(cbind(sapply(knots, function(tau) ((TijFN > tau) * (TijFN - tau) )))), nrow = length(TijFN))
  fZ0F <- cbind(1, TijFN, Z0F)
  bigXF <- cbind(fZF, fZ0F)

  
  if(IW != I){
  # ------------------------------ Unobserved subjects  + All time-points  --------------------------- #
  Z0FX <- matrix(as.vector(cbind(sapply(knots, function(tau) ((TijFX > tau) * (TijFX - tau) )))), nrow = length(TijFX))
  fZ0FX <- cbind(1, TijFX, Z0FX)
  bigXFX <- cbind(fZ0FX)
 }
  
  ### ------------------------------------------------------------------------------------------------------ ###
  ### ----  creating design matrix for K smooth time-varying functions using truncated linear splines   ---- ###
  ### ------------------------------------------------------------------------------------------------------ ###  
  
  
  for(k in 1 : K.hat){
    
    # ------------------------ Observed subjects  + Observed time-points  ------------------------------- #   
    DD <- paste("fZ",k," <- zetaw.hatN[,",k,"] * fZ0", sep = "")
    DD <- eval(parse(text = DD)) 
    bigX <- cbind(bigX, DD)
    rm(DD)
    
    
    # ----------------------------- Observed subjects  + All time-points  ------------------------------ #
    DD <- paste("fZF",k," <- zetaw.hatE[,",k,"] * fZ0F", sep = "")
    DD <- eval(parse(text = DD)) 
    bigXF <- cbind(bigXF, DD)
    rm(DD)
    
    
    if(IW != I){
    # ------------------------------ Unobserved subjects  + All time-points  --------------------------- #
    DD <- paste("fZFX",k," <- zetaw.hatX[,",k,"] * fZ0FX", sep = "")
    DD <- eval(parse(text = DD)) 
    bigXFX <- cbind(bigXFX, DD)
    rm(DD)
    }
  }
  
  if(dim(bigX)[2] > dim(bigX)[1]) stop("WARNING: Number of observations (n) is smaller than number of parameters (p)")
  
  
 
  ### ------------------------------------------------------------------------------------------------------------- ###
  ### --------------------------------  defining the terms to be penalized  --------------------------------------- ###
  ### ------------------------------------------------------------------------------------------------------------- ###  
  
  
  if (Yerror == "CS") { 
    
    E <- list()
    
    E[[1]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[1]])[c(1 : I)] <- 1
    
    E[[2]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[2]])[c( (I + 3) : (I + R) )] <- 1
    
    IN <- list()
    for(kk in 1:K.hat){
      IN[[kk]] <- c( (I + 3 + R + (kk - 1) * R) : (I + R + kk * R))
    }
    
    E[[3]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[3]])[unlist(IN)] <- 1
    
  }
  
  
  if (Yerror == "IS") { 
    
    E <- list()
    
    E[[1]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[1]])[c(1 : I)] <- 1
    
    E[[2]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[2]])[ (I + 1) : (2 * I) ] <- 1
    
    E[[3]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[3]])[c( ((2 * I) + 3) : ( (2 * I) + R) )] <- 1
    
    IN <- list()
    for(kk in 1 : K.hat){
      IN[[kk]] <- c( ( (2 * I) + 3 + R + ((kk - 1) * R) ) : ( (2 * I) + R + (kk * R) ))
    }
    
    E[[4]] <- diag(rep(0, dim(bigX)[2]))
    diag(E[[4]])[unlist(IN)] <- 1
    
  }
  
  
  ### ------------------------------------------------------------------------------------------- ###
  ### --------------------------------------  fitting model  ------------------------------------ ###
  ### ------------------------------------------------------------------------------------------- ###  
  
  
  time <- system.time(applygam <- bam(Yij ~ bigX - 1, paraPen = list(bigX = E), method = method))[3]
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ------ prediction for observed subject corresponding to its observed time points---------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  
  coefs = applygam$coef
  pred <- as.vector(as.matrix(bigX[, 1:length(coefs)]) %*% coefs )
 
  pred_err <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(length(unique(ID))), function(hh)  
    mean((Yij[which(ID == unique(ID)[hh])] - pred[which(ID == unique(ID)[hh])]) ^ 2)))))) 
  
  
  ### ------------------------------------------------------------------------------------------------------- ###
  ### ------ prediction of full (observed  + unobserved time points) trajectory for observed subjects ------- ###
  ### ------------------------------------------------------------------------------------------------------- ###  
  
  
  NY.pred <- as.vector(as.matrix(bigXF[, 1:length(coefs)]) %*% coefs )
  

  ### -------------------------------------------------------=------------------------------------- ###
  ### ----------------------- prediction of full rajectory for unobserved subjects ---------------- ###
  ### --------------------------------------------------------------------------------------------- ###  
  
  if (Yerror == "CS"){
    I1 <- I
  }    
  
  if (Yerror == "IS"){
    I1 <- 2 * I
  } 
  
  if(IW != I){
  NY.pred.X <- as.vector(as.matrix(bigXFX) %*% coefs[-c(1 : I1)] )
  }
  
  if(IW == I){
  NY.pred.X <- 9999999      # not available
  }  
  

  ### ------------------------------------------------------------------------------------------- ###
  ### ---------- prediction of mean of functional predictor at unobserved time points  ---------- ###
  ### ------------------------------------------------------------------------------------------- ###  
  
  
  fMU.gam.dat <- data.frame(x1 = rep(ss, dim(zetaw.hatF)[1]), x2 = rep(TijF, SS))                                    
  fMU.X <- matrix(predict.gam(MU.gam , newdata = fMU.gam.dat), nrow = length(TijF), ncol = SS, byrow = T)
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ------------------------------- Model adequacy (Ad-hoc) ---------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  AIC <- AIC(applygam)
  BIC <- BIC(applygam)
  Rsq <- summary(applygam)$r.sq
  dev <- summary(applygam)$dev.expl
  
  model_adeq <- c(AIC, BIC, Rsq, dev) 
  names(model_adeq) <- c("AIC", "BIC", "R-square", "dev expl")
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ------------------------------ estimation & prediction ----------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  
  qtiles <- seq(0, 1, length = R)[-c(1, R)]
  knots <- quantile(TT, qtiles)
  Zt <- matrix(as.vector(cbind(sapply(knots, function(tau) ((TT > tau) * (TT-tau) ) ))), nrow = length(TT))
  fZt <- cbind(1, TT, Zt)
  
 
  ### ------------------------------------------------------------------------------------------ ###
  ### ---------------------- time-varying smooth intercept ------------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  
  ME <- qnorm(1 - alpha / 2)    
  
  beta0 <- fZt %*% coefs[c( (I1 + 1) : (I1 + R) )]                                                                
  varb <- (as.matrix(fZt)) %*%  applygam$Vp[ c( (I1 + 1) : (I1 + R) ), c( (I1 + 1) : (I1 + R) )]  %*% t(as.matrix(fZt))   
  Bounds0 <-   cbind(beta0 - ME * sqrt(diag(varb)), 
                     beta0 + ME * sqrt(diag(varb)))
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ---------------------- time-varying smooth coefficient ----------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  
  beta <- list(); IN <- list(); varBeta <- list(); varBetaH <- list(); Bounds <- list()
  
  for(jj in 1 : K.hat){
    beta[[jj]] <- rep(0, J)
    IN[[jj]] <- c( (I1 + R + 1 + (jj - 1) * R) : (I1 + R + (jj * R)) )
    varBeta[[jj]] <- applygam$Vp[ unlist(IN[[jj]]), unlist(IN[[jj]])]  
    b0 <- fZt %*% coefs[IN[[jj]]]
    varBetaH[[jj]] <- (as.matrix(fZt)) %*% varBeta[[jj]] %*% t(as.matrix(fZt))
    beta[[jj]] <- b0
    Bounds[[jj]] = cbind(beta[[jj]] - ME * (sqrt(diag(varBetaH[[jj]]))), 
                         beta[[jj]] + ME * (sqrt(diag(varBetaH[[jj]]))))
    rm(b0)
  }
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ---------------------------------- random effects ---------------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  

  if (Yerror == "CS"){
    rand_int <- coefs[c(1 : I)]
    rand_slp <- NULL
  } 
  if (Yerror == "IS"){
    rand_int <- coefs[c(1 : I)]
    rand_slp <- coefs[c( (I + 1) : (2 * I) )]
  }
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### ------------------------------ Functional coefficient ------------------------------------ ###
  ### ------------------------------------------------------------------------------------------ ###  
  

  allb <- as.matrix(do.call(cbind, lapply(seq_len(K.hat), function(ff) beta[[ff]])))
  gammaE <- phi.hat %*% t(allb) 

  
  ### ------------------------------------------------------------------------------------------ ###
  ### ------------------------------------- prediction band ------------------------------------ ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  
  if(pred_interval == TRUE){
  
    
  # ----------------------------- Observed subjects  + All time-points  ------------------------------ #
  varYF <- bigXF %*% applygam$Vp %*% t(bigXF) +   diag(applygam$reml.scale, dim(bigXF)[1])
  boundYF <- cbind(NY.pred - ME * sqrt(diag(varYF)), NY.pred + ME * sqrt(diag(varYF)))
  
  
  if(IW != I){
  # ------------------------------ Unobserved subjects  + All time-points  --------------------------- #
  if (Yerror == "CS"){
    varXF <- bigXFX %*% applygam$Vp[-(1 : I1), -(1 : I1)] %*% t(bigXFX) +  gam.vcomp(applygam)[1]^2 * (fZFX  %*% t( fZFX ))  + 
              diag(applygam$reml.scale, dim(bigXFX)[1])
  }    
  
  if (Yerror == "IS"){
     varXF <- bigXFX %*% applygam$Vp[-(1 : I1), -(1 : I1)] %*% t(bigXFX) +  
     fZFX  %*% bdiag(diag(rep(gam.vcomp(applygam)[1]^2, I)), diag(rep(gam.vcomp(applygam)[2]^2, I)))  %*% t(fZFX)  + 
     diag(applygam$reml.scale, dim(bigXFX)[1])
  } 
  boundXF <- cbind(NY.pred.X - ME * sqrt(diag(varXF)),  NY.pred.X + ME * sqrt(diag(varXF)))
  }
  
  
  if(IW == I){
    boundXF <- 9999999  
  }
  
  if(pred_interval == FALSE){
    boundYF <- boundXF <- 9999999
  }
  }
  
  
  ### ------------------------------------------------------------------------------------------ ###
  ### -------------------------------------- output -------------------------------------------- ###
  ### ------------------------------------------------------------------------------------------ ###  
  
  return = list(fit = applygam, 
                pred = pred, 
                full_pred_obs = NY.pred, 
                full_pred_PA_un = NY.pred.X, 
                error = pred_err, 
                gamma = gammaE, 
                phi = phi.hat, 
                xiF = zetaw.hatF, 
                evalues = eigenvalues,
                xi0 = zetaw.hatN, 
                Fnp = K.hat, 
                ID_full = IDF, 
                ID_obs = ID, 
                SS_inter = rand_int, 
                SS_slp = rand_slp, 
                model_adeq = model_adeq, 
                fMU.X = fMU.X,  
                zetaw_hat = zetaw.hat,  
                MU_X  =  MU.X,
                Mean  = fbps.fit$Yhat, 
                freq = freq, 
                V = V, 
                Y = Yij, 
                Tm = Tij, 
                beta = beta, 
                beta0 = beta0, 
                bounds =  Bounds,
                boundYF = boundYF, 
                boundXF = boundXF, 
                fpr = obs.func.cov,
                fpca.fit1 = covsmooth,
                fpca.fit2 = fit2)
}

  
  
