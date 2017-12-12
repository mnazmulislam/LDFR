# --- READ ME 
#  title: Getting output from the fit of the model described in: 
#        "Longitudinal dynamic functional regression" (LDFR)
#        by Md Nazmul Islam, Ana-Maria Staicu, and Eric van Heugten
# created by Md Nazmul Islam
# date: 2017/11
# ---

############################################################################################
### ------------------------ defining outputs for the Gaussian fit --------------------- ###
############################################################################################


#             IN.MPE.y : IN-sample prediction error
#            OUT.MPE.y : OUT-sample prediction error
#             RMPE_trj : Root-mean-prediction-error for full response trajectory for existing subjects
#            RMPE_trjX : Root-mean-prediction-error for full response trajectory for new subjects; if 999999 then Not available
#                   FT : object of LDFR() fit for full trajectory
#        YijF & Yij.XF : True response trajectories for existing and new subjects  



Yestim <- function (A, I, IW, TT, ss, zetamu, zetasigma1, zetasigma2, zetasigma3, zetasigma4, 
      delta, minJi, maxJi, Cov_Error, mu.e.x1, mu.e.x2, mu.e.x3, var.e.x1, var.e.x2, var.e.x3, 
      errYbi0mu, errYbijmu, errYbi0sigma, errYbijsigma, errYbi1mu, errYbi1sigma, errYbi12sigma, 
      mm, J, Yerror, nbf, method, pve, full_traj, pred_interval, alpha, new_subj)
{ 
  
  set.seed (A)  
  data <- data_func(A, I, TT, ss, zetamu, zetasigma1, zetasigma2, zetasigma3, zetasigma4, 
                     delta, minJi, maxJi, Cov_Error, mu.e.x1, mu.e.x2, mu.e.x3, var.e.x1, 
                    var.e.x2, var.e.x3, errYbi0mu, errYbijmu, errYbi0sigma, 
                     errYbijsigma, errYbi1mu, errYbi1sigma, errYbi12sigma, mm, J, Yerror)
                     

  
  ### -------------------------------------------------------------------------- ###
  ### ------------------------- extracting information ------------------------- ###
  ### -------------------------------------------------------------------------- ###

    dat.full <- data$dat.full; dat.sparse <- data$dat.sparse; dat.train <- data$dat.train; 
    dat.test <- data$dat.test; m2 <- data$m2                                                              
 
    uTij <- sort(unique(dat.full[,4])); lng <- length(uTij); J <- length(TT); SS <- length(ss);                
    if ( lng != J ) stop ( "Error: Differences  in uTij & TT " )                                          

    
  ### ------------------------- observed data ------------------------- ###
  
    IDs <- dat.sparse[,1]; 
    TijS <- dat.sparse[,4];
    Vs <- dat.sparse[,3]
    obs.func.covS <- as.matrix(dat.sparse[,c(11:111)])
    YijS <- dat.sparse[, 113];
    sprse_seq <- dat.sparse[, 112];   
  
  ### ------------------------- train data ------------------------- ###
  
    ID <- dat.train[ , 1]; 
    V <- dat.train[, 3]; 
    Tij <- dat.train[, 4]; 
    obs.func.cov <- as.matrix(dat.train[, c( 11 : 111)]);
    red_seq <- dat.train [ , 112];  
    Yij <- dat.train[ , 113];
  
  ### ------------------------- test data ------------------------- ###
  
    ID.test <- dat.test[, 1];
    Tij.test <- dat.test[, 4]
    obs.func.covT <- as.matrix(dat.test[, c( 11 : 111)]); 
    Yij.test <- dat.test[, 113]; 
    test_seq <- dat.test [, 112]; 
    
  
  ### ------------------------- full data ------------------------- ###
 
    YijF <- dat.full[, 113];
    IDF <- dat.full[, 1];
  
  
  ### --------------------------------------------------------------------------------------------- ###
  ### -------------------------------- estimating mean function ----------------------------------- ###
  ### --- fbps() is used to estimate; and gam() is used to predict the mean for all time points --- ###
  ### --------------------------------------------------------------------------------------------- ###  
    
    
    ln <- length(Tij)
  
    fbps.fit <-  fbps(data = obs.func.cov, covariates = c (Tij, ss))
    muhat0 <- matrix(0, nrow = J, ncol = SS)
    muhat.dat <- do.call(rbind, lapply(seq_len(J), function(de) { 
    t1 <- uTij[de]
    r1 <- as.vector(which(Tij == t1))[1]                                                
    muhat0[de,] <- fbps.fit$Yhat[r1, ] }))
    
    newmu.dat <- data.frame ( y = as.vector(t(muhat.dat)), x1 = rep( ss, J), x2 = rep( uTij, SS) )
    MU.gam <- gam(y ~ te(x1, x2, k = c(10, 5), bs ="cr"), data = newmu.dat, method = method)                         
    MU.gam.dat <- data.frame (x1 = rep(ss, ln ), x2 = rep( Tij, SS) )                                    
    MU.XA <- matrix (predict.gam(MU.gam , newdata = MU.gam.dat ), nrow = ln, ncol = SS, byrow = T )
  
    MU.gam.dat <- data.frame (x1 = rep(ss, length(TT) ), x2 = rep( TT, SS) )                                    
    MU.X <- matrix (predict.gam(MU.gam , newdata = MU.gam.dat ), nrow = length(TT), ncol = SS, byrow = T )
    
    
    
  ### ---------------------------------------------------------------------------------- ###
  ### ---------------------  demeaning functional predictor (train) -------------------- ###
  ### ---------------------------------------------------------------------------------- ###  
    
    
    deM.X <- obs.func.cov - fbps.fit$Yhat  
    
    
  ### ---------------------------------------------------------------------------------------- ###
  ### ----------- finding the eigen-components of covariance operator (train) ---------------- ###
  ### ---------------------------------------------------------------------------------------- ###  
    
    
    covsmooth <- fpca.face (Y = deM.X, pve = pve, var = TRUE)                                                                       
    phi.mat <- covsmooth$efunctions                                 
    phi.hat <- phi.mat * sqrt(SS)                                                                         
    K.hat <- covsmooth$npc;                                         
    
    zetaw.hat <- covsmooth$scores / sqrt (SS)                    
    mu <- covsmooth$mu                                              
    eigenvalues <- covsmooth$evalues / SS                             
    
    
  ### -------------------------------------------------------------------------------- ###
  ### -------------------------- Replacing negative eigenvalues ---------------------- ###
  ### -------------------------------------------------------------------------------- ###  
    
    
    eigenvalues[which(eigenvalues < 0)] <- 0                      
    
    
  ### -------------------------------------------------------------------------------- ###
  ### ----------------------------- checking norm ------------------------------------ ###
  ### -------------------------------------------------------------------------------- ###  
    
    
    k1 <- all.equal ( as.vector ( do.call( rbind, lapply ( seq_len(K.hat), function (feb) 
      sqrt(mean ( length(ss) * phi.mat[,feb] ^ 2 )) ) ) ), rep(1, K.hat) )
    if ( k1 != "TRUE" ) stop ( "l2 norm mismatches with the simulated basis functions")
    
    
  ### ---------------------------------------------------------------------------------------- ###
  ### --------------------- estimating mean of functional predictor (train) ------------------ ###
  ### ---------------------------------------------------------------------------------------- ###  
  
  
    MU.Xi.Hat <- MU.XA + matrix (rep(mu, ln), nrow = ln, ncol = SS, byrow = T  )
    
  
  ### ----------------------------------------------------------------------------------------------- ###
  ### --------------------- predicting the time-varying loadings (train + test) --------------------- ###
  ### ----------------------------------------------------------------------------------------------- ###  
  
  
    xi.hat0 <- list()                                                                                    
    for( k in seq_len (K.hat) ) 
    {
      xihat0.vec <- zetaw.hat[, k]
      xi.hat0[[ k ]]<-t( sapply( seq_len(I), function(i) 
      {xi.subj <- matrix(nrow = 1, ncol = J) 
      xi.subj[V[which( ID == i )]] <- xihat0.vec[ which(ID == i) ]
      return( xi.subj )
      }))
    }                                                                                                    
    fit2 <- lapply(xi.hat0, function(a) try (fpca.sc(Y = a, pve = pve, var = TRUE))) 
    xi.hat<- lapply(fit2, function(a) a$Yhat)    
    
  
  ### --------------------------------------------------------------------------------------------------------- ###
  ### ---------- predicting the time-varying loadings for all (observed + unobserved)  time-points ------------ ###
  ### --------------------------------------------------------------------------------------------------------- ###  
  
  
    zetaw.hatF <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg)) ))
  
  
  ### ---------------------------------------------------------------------------------------------------------- ###
  ### ------------ predicting the time-varying loadings for observed time-points for train set ----------------- ###
  ### ---------------------------------------------------------------------------------------------------------- ###  
  
  
    zetaw.hatN <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg))[red_seq] ))
    
    
  ### ----------------------------------------------------------------------------------------------------------- ###
  ### ------------ predicting the time-varying loadings for observed time-points for test set ------------------- ###
  ### ----------------------------------------------------------------------------------------------------------- ###

  
    zetaw.hatTst <- do.call (cbind, lapply ( xi.hat, function(gg) as.vector(t(gg))[test_seq] ))

    
  ### ------------------------------------------------------------------------------------------------------------ ###
  ### ------------------------ creating design matrix for SS random effects  ------------------------------------- ###
  ### ------------------------------------------------------------------------------------------------------------ ###  
      
    
  if (Yerror == "CS") { 
      
  # ------------------------------------------------ Training ----------------------------------------------------- #  
    Z1 = matrix(0, nrow = length(Yij), ncol = length(unique(ID)))    
    for (i in 1 : I) {
      Z1[which(ID == unique(ID)[i]), i] = 1                     
    }
    fZ <- cbind(Z1)
    
    
    # -------------------------------------------- Test ------------------------------------------------------------ #
    Z1.test = matrix(0, nrow = length(Yij.test), ncol = length(unique(ID.test)))    
    for (i in 1 : length(unique(ID.test))) {
      Z1.test[which(ID.test == sort(unique(ID.test))[i]), i] = 1                     
    }
    fZ.test <- cbind(Z1.test)
    
    }
      
      
  if (Yerror == "IS") {   
    
    # --------------------------------------------- Training --------------------------------------------------------- # 
    Z1 = matrix(0, nrow = length(Yij), ncol = length(unique(ID)))          
    Z2 = matrix(0, nrow = length(Yij), ncol = length(unique(ID)))
  
    for (i in 1 : I) {
    Z1[which(ID == unique(ID)[i]), i] = 1  
    Z2[which(ID == unique(ID)[i]), i] = Tij[which(ID == unique(ID)[i])]
    }
    fZ <- cbind(Z1, Z2)
    
    
    # --------------------------------------------- Test ------------------------------------------------------------- #
    Z1.test = matrix(0, nrow = length(Yij.test), ncol = length(unique(ID.test)))    
    Z2.test = matrix(0, nrow = length(Yij.test), ncol = length(unique(ID.test)))
    
    for (i in 1 : length(unique(ID.test))) {
    Z1.test[which(ID.test == sort(unique(ID.test))[i]), i] = 1                     
    Z2.test[which(ID.test == sort(unique(ID.test))[i]), i] = Tij.test[which(ID.test == sort(unique(ID.test))[i])]
    }
    
    fZ.test <- cbind(Z1.test, Z2.test)
    }
  
  
    ### ---------------------------------------------------------------------------------------------------------- ###
    ### ------------  creating design matrix smooth intercept using truncated linear splines   ------------------- ###
    ### ---------------------------------------------------------------------------------------------------------- ###  
      
      
    R <- nbf
    qtiles <- seq(0, 1, length = R)[- c(1, R)]
    knots <- quantile(TT, qtiles)
    
    # ----------------------------------------------- Training ------------------------------------------------------ # 
    Z0 <- matrix(as.vector(cbind(sapply(knots, function(tau) ((Tij > tau) * (Tij - tau) )))), nrow = length(Tij))
    fZ0 <- cbind(1, Tij, Z0)
    bigX <- cbind(fZ, fZ0)
    
    
    # ------------------------------------------------ Test --------------------------------------------------------- #
    Z0.test <- matrix(as.vector(cbind(sapply(knots, function(tau) ((Tij.test > tau) * (Tij.test - tau) ) ))), nrow = length(Tij.test))
    fZ0.test <- cbind(1, Tij.test, Z0.test)
    bigX.test <- cbind(fZ.test, fZ0.test)
    
    
    ### ----------------------------------------------------------------------------------------------------------- ###
    ### -------  creating design matrix for K smooth time-varying functions using truncated linear splines   ------ ###
    ### ----------------------------------------------------------------------------------------------------------- ###  
    
    
    for(k in 1 : K.hat){
      
      # ----------------------------------------------- Training ---------------------------------------------------- # 
      DD <- paste("fZ",k," <- zetaw.hatN[,",k,"] * fZ0", sep = "")
      DD <- eval(parse(text = DD)) 
      bigX <- cbind(bigX, DD)
      rm(DD)
      
      # ---------------------------------------------- Test ---------------------------------------------------------- #
      DD <- paste("fZ.test",k," <- zetaw.hatTst[,",k,"] * fZ0.test", sep = "")
      DD <- eval(parse(text = DD)) 
      bigX.test <- cbind(bigX.test, DD)
      rm(DD)
    }
    
  
    if(dim(bigX)[2] > dim(bigX)[1]) stop("WARNING: Number of observations (n) is smaller than number of parameters (p)")
    dim(bigX); dim(bigX.test)
    
    
    ### ------------------------------------------------------------------------------------------------------------- ###
    ### -----------------------------  defining the terms to be penalized  ----------------------------------------- ###
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
    
    
    ### ---------------------------------------------------------------------------------------------------------- ###
    ### --------------------------------------------  fitting model  --------------------------------------------- ###
    ### ---------------------------------------------------------------------------------------------------------- ###  
    
    
    time <- system.time(applygam <- bam(Yij ~ bigX - 1, paraPen = list(bigX = E), method = method))[3]

    
    ### ---------------------------------------------------------------------------------------------------------- ###
    ### --------------------------------------  predicting responses in the train set ---------------------------- ###
    ### ---------------------------------------------------------------------------------------------------------- ###  
    
    
    coefs = applygam$coef
    NY.pred <- as.vector(as.matrix(bigX[, 1:length(coefs)]) %*% coefs )
    
    IN.MPE.y <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(I), function(hh)  
      mean((Yij[which(ID == hh)] - NY.pred[which(ID == hh)]) ^ 2))))))
    
    
    ### ----------------------------------------------------------------------------------------------------------- ###
    ### ---------------------------------  predicting responses in the test --------------------------------------- ###
    ### ----------------------------------------------------------------------------------------------------------- ###  
    
    
    NY.pred.test <- as.vector ( as.matrix(bigX.test[, 1:length(coefs)]) %*% coefs )
      
    OUT.MPE.y <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(length(unique(ID.test))), function(hh)  
      mean((Yij.test[which(ID.test == sort(unique(ID.test))[hh])] - NY.pred.test[which(ID.test == sort(unique(ID.test))[hh])]) ^ 2))))))
    
    
    
    ### ------------------------------------------------------------------------------------------------------------ ###
    ### -------------------------------  predicting responses for full trajectory ---------------------------------- ###
    ### ------------------------------------------------------------------------------------------------------------ ###  
    
  
    
    if (full_traj == TRUE){ 
      
      # ----------- Profiles' information of new and existing subjects at sparse time-points -------------------------#
       if(new_subj == TRUE){  ## Generate information for (exisiting + new) subjects 
        
        rm(data); rm(dat.sparse); rm(dat.full)
         
        set.seed (A)  
        data <- data_func(A, I = IW, TT, ss, zetamu, zetasigma1, zetasigma2, zetasigma3, zetasigma4, 
                          delta, minJi, maxJi, Cov_Error, mu.e.x1, mu.e.x2, mu.e.x3, var.e.x1, 
                          var.e.x2, var.e.x3, errYbi0mu, errYbijmu, errYbi0sigma, 
                          errYbijsigma, errYbi1mu, errYbi1sigma, errYbi12sigma, mm, J, Yerror)
        
        
        ### ------------------------------------------------------------------------------------------------------ ###
        ### -------------------------------------- extracting information ---------------------------------------- ###
        ### ------------------------------------------------------------------------------------------------------ ###
        

        dat.full <- data$dat.full; dat.sparse <- data$dat.sparse; 
        dat.obs <- data$dat.sparse;
        
        EX <- which (dat.obs[,1] > I)                                                            # Last IW - I subjects are the new ones
        EXF <- which (dat.full[,1] > I)
        
  
        Yij.ob <- dat.obs[-EX, 113]; red.seq.ob <- dat.obs[-EX, 112]; 
        obs.func.covX <- as.matrix ( dat.obs [, c( 11 : 111)] ); 
        IDX <- dat.obs[, 1]; 
        VX <- dat.obs[, 3]; 
        TijX <- dat.obs[, 4]; 
    

        IDF <- dat.full[-EXF, 1];  YijF <- dat.full[-EXF , 113];                                 # observed subjects full information of response and id
        ID.XF <- dat.full[EXF, 1];  Yij.XF <- dat.full[EXF , 113];                               # Unobserved subjects full information of response and id
        
     
      
      FT <- try(LDFR(ID = IDX, obs.func.cov = obs.func.covX, Tij = TijX, Yij = Yij.ob, V = VX, I, IW, nbf, 
                       method =  method, pve = pve, Yerror = Yerror, family = gaussian, 
                       pred_interval = pred_interval, alpha = alpha, out_subj = red.seq.ob))
        
        
      RMPE_trj <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(length(unique(IDF))), function(hh)  
               mean((YijF[which(IDF == unique(IDF)[hh])] - FT$full_pred_obs[which(IDF == unique(IDF)[hh])]) ^ 2)))))) 
   
      
      RMPE_trjX <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(length(unique(ID.XF))), function(hh)  
            mean((Yij.XF[which(ID.XF == unique(ID.XF)[hh])] - FT$full_pred_PA_un[which(ID.XF == unique(ID.XF)[hh])]) ^ 2)))))) 
        
      }
    
      # -------------- Profiles' information of only existing subjects is available at sparse time-points  ------------# 
      if(new_subj == FALSE){                              
        
        FT <- try(LDFR(ID = IDs, obs.func.cov = obs.func.covS, Tij  = TijS, Yij = YijS,  V = Vs, I, IW = NULL, nbf, 
                       method =  method, pve = pve, Yerror = Yerror, family = gaussian, 
                       pred_interval = pred_interval, alpha = alpha, out_subj = FALSE))
        
        RMPE_trj <- sqrt(mean(as.vector(do.call(rbind, lapply(seq_len(length(unique(IDF))), function(hh)  
          mean((YijF[which(IDF == unique(IDF)[hh])] - FT$full_pred_obs[which(IDF == unique(IDF)[hh])]) ^ 2)))))) 
       
        RMPE_trjX <- 999999                                                                                     # Not available  
      }
    }  
    else if (full_traj == FALSE){
        FT <- NULL
    }
    
    output = list (IN.MPE.y = IN.MPE.y, OUT.MPE.y = OUT.MPE.y, RMPE_trj =  RMPE_trj, 
                   RMPE_trjX = RMPE_trjX, fit = FT, YijF = YijF, Yij.XF = Yij.XF)
  }    
    
    

    
    
    
    