#This program is for Data generation and parameter estimation for
#the Sequential Response Model;
#The details are in the paper:  A Sequential Response Model for 
#Analyzing Process Data on Technology-Based Problem-solving Tasks.

#####################################################################
#                Generate state transition parameters               #
#####################################################################
GenLam <- function(){
  lambda <- c()
  lambda[1] <- rnorm(1,-3,0.2)
  lambda[3:4] <- rnorm(2,-1.5,0.3)
  lambda[6:7] <- rnorm(2,0,0.3)
  lambda[9:11] <- rnorm(3,0,0.5)
  lambda[13:14] <- rnorm(2,0,0.3)
  lambda[16:17] <- rnorm(2,1,0.3)
  
  lambda[2] <- 0 - lambda[1] 
  lambda[5] <- 0 - sum(lambda[3:4])
  lambda[8] <- 0 - sum(lambda[6:7])
  lambda[12]<- 0 - sum(lambda[9:11])
  lambda[15]<- 0 - sum(lambda[13:14])
  lambda[18]<- 0 - sum(lambda[16:17])
  return(lambda)
}

#####################################################################
#                        Calculate probability                      #
#####################################################################
Prob0 <- function(lambda, theta){
  v <- length(lambda)
  n <- length(theta)
  
  a <- rep(-1, v)
  a[c(2, 5, 8, 12, 15, 18)] <- 1

  lam_mat <- matrix(rep(lambda, each=n), nrow = n) # [n,v]
  the_mat <- as.matrix(theta)  %*% t(as.matrix(a)) # [n,v]
  
  Ep <- exp( lam_mat + the_mat) 
  Ea <- rowSums(matrix(Ep[,1:2], ncol = 2)) 
  Eb <- rowSums(matrix(Ep[,3:5], ncol = 3)) 
  Ec <- rowSums(matrix(Ep[,6:8], ncol = 3)) 
  Ed <- rowSums(matrix(Ep[,9:12], ncol = 4))
  Ee <- rowSums(matrix(Ep[,13:15], ncol = 3)) 
  Ef <- rowSums(matrix(Ep[,16:18], ncol = 3))
  p2 = Ep/cbind(Ea, Ea, Eb, Eb, Eb, Ec, Ec, Ec, Ed, Ed, Ed, Ed, Ee, Ee, Ee,Ef,Ef,Ef)
  return(p2 = p2 )}

#####################################################################
#                         Data simulation                           #
#####################################################################
SimData0 <- function(lambda, theta, transSet){
  scv = c()        
  stl = list()     
  n <- length(theta)
  for (i in 1:n){
    s0 = 'A'  
    sc = s0       
    st = c()     
    while( TRUE ){
      p <- Prob0(lambda, theta[i])
      if (s0 == 'A'){  s1 <- sample(c('#','B'), 1, replace=T, p[1:2]) 
      }else if (s0 == 'B'){ s1 <- sample(c('A','#','C'), 1, replace=T, p[3:5]) 
      }else if (s0 == 'C'){ s1 <- sample(c('A','B','D'), 1, replace=T, p[6:8])  
      }else if (s0 == 'D'){ s1 <- sample(c('B','C','#','E'), 1, replace=T, p[9:12]) 
      }else if (s0 == 'E'){ s1 <- sample(c('C','#','F'), 1, replace=T, p[13:15]) 
      }else if (s0 == 'F'){ s1 <- sample(c('C','D','G'), 1, replace=T, p[16:18]) 
      }
      transId  <- grep( paste0(s0,s1) , transSet)  
      st <- c(st, transId )   
      if (s1 == '#' ){break}   
      sc <- paste0(sc, s1)    
      if (s1 == 'G' ){break}  
      s0 <- s1 }             
    scv = c(scv,sc)
    stl[[i]] = st 
  }
  result <- list(scv = scv, stl = stl )
  return(result)
}
#####################################################################
#     Count the frequency of each person on each transition         #                  
#####################################################################
Fre.each <- function(stl, n, u){
  fq.mat <- matrix(0, nrow = n, ncol = u)
  for (i in 1:u){
    fq.mat[,i] <- sapply(stl, function(x) sum(x==i)) 
  }
  return(fq.mat)}

#####################################################################
#                Generate initial value  for lambda                 #                  
#####################################################################
LamStart <- function(stl){
  fq <- rep(1,18)
  fq[as.numeric(names(table(unlist(stl))))] <- as.numeric(table(unlist(stl))) 
  fm <- fq[c(2,2,5,5,5,8,8,8,12,12,12,12,15,15,15,18,18,18)]
  lamS = log(fq/fm)
  lamS[2] = 0- sum(lamS[1])
  lamS[5] = 0- sum(lamS[c(3:4)])
  lamS[8] = 0- sum(lamS[c(6:7)])
  lamS[12]= 0- sum(lamS[c(9:11)])
  lamS[15]= 0- sum(lamS[c(13:14)])
  lamS[18]= 0- sum(lamS[c(16:17)])
  lamS1 <- lamS/3
  
  lamS2 <- c(scale(lamS[1:2]), scale(lamS[3:5]), scale(lamS[6:8]), scale(lamS[9:12]), scale(lamS[13:15]), scale(lamS[16:18]))
  lamS3 <- c(scale(lamS))
  lamS4 <- lamS2/2
  lamS5 <- rep(0,18)
  lamS.Mat <- rbind(lamS1,lamS2,lamS3,lamS4,lamS5)
  
  return(lamS.Mat = lamS.Mat)
}

#####################################################################
#                           Update lambda                           #                  
#####################################################################
UpdateLam <- function(lam0, zq.v, cw.v, sd = rep(0.3,6)){
  z <- length(sd)
  lam1 <- lam0
  for (i in 1:z){
    lam1[cw.v[[i]] ] <- rnorm(length(cw.v[[i]] ),lam0[cw.v[[i]] ],sd[i]) 
    lam1[ zq.v[i] ] <- 0 - sum(  lam1[cw.v[[i]] ]  )  
  }
  return(lam1=lam1)
}

#####################################################################
#                     likelihood for each state                     #                  
#####################################################################
LLABC0 <- function(p, fq.mat, zt.v){
  z <- length(zt.v)
  llABC <- rep(1,z)
  for (i in 1:z){
    zt.id <- zt.v[[i]]  #c(1:2)
    ll.i <- 1
    for (j in zt.id){
      ll.i <- ll.i * p[,j]^fq.mat[,j] 
    }
    llABC[i] <- prod(ll.i)  
  }
  return(llABC)}

#####################################################################
#       Calculate transition probability (Alpha) for lambda         #                  
#####################################################################
AlphaLambda <- function(llABC, zt.v, lam1, lam0){
  z <- length(zt.v)
  pr <- rep(1,z)
  for (i in 1:z){
    pr[i] <-  prod(exp(-0.5*(lam1[ zt.v[[i]] ])^2)/exp(-0.5*(lam0[ zt.v[[i]]  ])^2))     #informative prior: N(0,1)
  # pr[i] <-  prod(exp(-0.005*(lam1[ zt.v[[i]] ])^2)/exp(-0.005*(lam0[ zt.v[[i]]  ])^2)) #non-informative prior: N(0,100) 
    } 
  AlphaLam <- llABC* pr
  return(AlphaLam)
}

#####################################################################
#                     Bayes Parameter Estimation                    #                  
#####################################################################
updatePara <-  function(stl, M, L, B, sd){
  n <- length(stl)
  u <- 18
  z <- 6
  zq.v <- c(2, 5, 8, 12, 15, 18) 
  cw.v <- list(c(1),c(3:4),c(6:7),c(9:11),c(13:14),c(16:17))  
  zt.v <- list(c(1:2),c(3:5),c(6:8),c(9:12),c(13:15),c(16:18))
  
  lamS0 <- LamStart(stl)  # initial value for lambda
  fq.mat <- Fre.each(stl, n, u)  
  theMat <- array(0,c(n,L,M)) #Sampling records
  lamMat <- array(0,c(u,L,M))   
  ifUpdateThe.sum <- matrix(0,n, M) #if update records
  ifUpdateLam <- array(1,c(z,L,M))   
  startTime <- Sys.time() 
  
  for (m in 1:M){ 
    theS <- rnorm(n) 
    lamS <- lamS0[m,] 
    the1 <- theS               
    lam1 <- lamS               
    ifUpdateThe.sum[,m] <- rep(0,n) 
    for (l in 1:L){ 
      if( l%%1000 == 0){print( paste('Chain:',m,'  Iter:',l,'  Time:',Sys.time()-startTime) )}
      #update theta 
      r <- runif(n)
      the0 <- the1              
      the1 <- rnorm(n ,the0, 1.5)  # candidate value ~ proposal distribution 
      p <- Prob0(lam1,the1)/Prob0(lam1,the0)  
      lltheta <- sapply(mapply(function(x,y) x[[1]][y], apply(p,1,list), stl), prod)
      AlphaThe <- lltheta * exp(-0.5*(the1)^2) / exp(-0.5*(the0)^2)  
      the1[AlphaThe < r ] <- the0[AlphaThe < r]    
      ifUpdateThe.sum[,m] <-  ifUpdateThe.sum[,m] + (AlphaThe > r)
      theMat[,l,m] <- the1
      
      #update lambda
      lam0 <- lam1             
      lam1 <- UpdateLam(lam0, zq.v, cw.v, sd) # candidate value ~ proposal distribution 
      p <- Prob0(lam1,the1) / Prob0(lam0,the1) 
      llABC <- LLABC0(p, fq.mat, zt.v)
      AlphaLam <- AlphaLambda(llABC, zt.v,lam1,lam0)
      AlphaLam[is.na(AlphaLam)] <- 0
      r <- runif(z)
      for (i in 1:z){ 
        if (AlphaLam[i] < r[i]){ 
           lam1[ zt.v[[i]] ] = lam0[ zt.v[[i]]] 
           ifUpdateLam[i,l,m] <-0 } 
      }
      lamMat[,l,m] <- lam1
    } 
    print( paste('Chain:',m,'  Time consuming:',Sys.time()-startTime))
    print( 'Update rate:')
    print(colMeans(ifUpdateThe.sum)/L)  
    print(rowMeans(ifUpdateLam[,,m])) 
  }
  return(simpleCoda = list(theMat=theMat, lamMat=lamMat))
  }

#####################################################################
#           Convert sampling matrix into MCMC list                  #                  
#####################################################################
mat.to.mcmc <- function(B,theMat,lamMat){
  the.list <- c()
  lam.list <- c()
  M <- dim(theMat)[[3]]
  L <- dim(theMat)[[2]]
  
  for(m in 1:M){
    the.m <- t(theMat[,B:L,m])
    lam.m <- t(lamMat[,B:L,m])
    
    colnames(the.m) <- 1:n
    colnames(lam.m) <- transSet
    
    the.list[[m]] <- mcmc(the.m)
    lam.list[[m]] <- mcmc(lam.m)
  }
  
  the.mcmc.list <- as.mcmc.list(the.list) # mcmc.list
  lam.mcmc.list <- as.mcmc.list(lam.list)
  return(list(the.mcmc.list=the.mcmc.list, lam.mcmc.list=lam.mcmc.list))
}
#####################################################################
#                        Estimation accuracy                        #                  
#####################################################################
CorMse <- function(theta,lambda,theHat,lamHat,csv){
  patHat <- c()
  pat <- c()
  for (sc in unique(scv)){
    patHat <- c(patHat, mean(theHat[scv == sc])) 
    pat <-  c(pat, mean(theta[scv == sc])) 
  }
  corThe <- cor(patHat,pat)
  corLam <- cor(lamHat,lambda)
  abseThe <- mean(abs(patHat-pat))
  abseLam <- mean(abs(lamHat-lambda))
  mseThe <- mean((patHat-pat)^2) 
  mseLam <- mean((lamHat-lambda)^2) 
  cri <- c(corThe=corThe, corLam=corLam,abseThe=abseThe, abseLam=abseLam, 
           mseThe=mseThe, mseLam= mseLam) 
  cri_lam <- cbind(abseThe = abs(lamHat-lambda), MseLam = (lamHat-lambda)^2)
  return(result = list(cri, cri_lam) )
}

#####################################################################
#                  trace plot for each lambda                       #                  
#####################################################################
plotMCMC <- function(lamMat, lambda, transSet){

  u <- dim(lamMat)[1]
  M <- dim(lamMat)[3]
  for (j in 1:u){
    nam <- paste0('Lambda_',transSet[j])
    plot(lamMat[j,,1] , type="l" , ylim = c(-4,4), ylab = nam, xlab = 'MCMC_Chain')
    if(M>1){ for(m in 2:M){lines(lamMat[j,,m], col = m)} }
    if(sum(lambda)!=0){abline(h=lambda[j] , col="red")}
  }
}

#####################################################################
#                    Model fit : DIC, log.PsML                      #                  
#####################################################################
modFit <- function(theMat,lamMat,theHat,lamHat){
  M <- dim(theMat)[[3]]
  L <- dim(theMat)[[2]]
  n <- dim(theMat)[[1]]
  
  log.like.all <- rep(NA, 5000*M)
  like.theta.all <- matrix(NA,n,5000*M)
  i <- 1
  startTime <- Sys.time()
  for (m in 1:M){
    for (l in 5001:L){ #5001~10000
      if( i%%1000 == 0){print( paste('DIC  iter:',i,'time:',Sys.time()-startTime) )}
      p <- Prob0(lamMat[,l,m],theMat[,l,m])  
      like.theta.l <- sapply(mapply(function(x,y) x[[1]][y], apply(p,1,list), stl), prod)
      log.like.l <- sum(log(like.theta.l))
      log.like.all[i] <-  log.like.l
      like.theta.all[,i] <- like.theta.l 
      i <- i+1
    }
  }
  #DIC
  p <- Prob0(lamHat,theHat)
  lltheta <- sapply(mapply(function(x,y) x[[1]][y], apply(p,1,list), stl), prod) 
  log.like <- sum(log(lltheta))
  D <- -2*log.like 
  D.exp <- -2*mean(log.like.all) 
  DIC <- 2*D.exp - D
  #CPO
  CPO.i <- 1/rowMeans(1/like.theta.all)   
  log.PsML <- sum(log(CPO.i)) #model compare: 2(log.PsML2-log.PsML1)
  print(paste('DIC :',DIC))
  print(paste('log.PsML :',log.PsML))
}

#####################################################################
#                        Model fit : PPMC plot                      #                  
#####################################################################
ppmcPlot <- function(transSet, theMat,lamMat){
  M <- dim(theMat)[[3]]
  L <- dim(theMat)[[2]]
  u <- length(transSet)
  trans.frq.all <- matrix(0, u, 500*M ) 
  rownames(trans.frq.all) <- c(1:u)
  chi.pp <- c()
  i <- 1
  startTime <- Sys.time()
  for (m in 1:M){
    for (l in 9501:L){ #9501~10000
      if( i%%100 == 0){print(paste('PPMC plot  iter:',i,'Time:',Sys.time()-startTime))}
      # PPMC - transition frequency 
      # generate data
      dataList.l <- SimData0(lamMat[,l,m], theMat[,l,m], transSet)
      # transition frequency
      stl.l <- dataList.l$stl
      trans.frq.l <- table(unlist(stl.l))
      trans.frq.all[names(trans.frq.l),i] <- as.numeric(trans.frq.l) 
      chi.pp[i] <- chisq.test(trans.frq.all[,i])$statistic
      i <- i+1
    }
  }
  
  f <- factor(rep(c(1:u),500*M)) 
  trans.frq.all.frame <- data.frame(c(trans.frq.all),f) 
  
  boxplot(c(trans.frq.all) ~ f, trans.frq.all.frame , 
          names = transSet, range = 2, cex=0.5,
          xlab="State Transitions", ylab="Frequency of State Transitions")
  fq <- rep(1,u)
  fq[as.numeric(names(table(unlist(stl))))] <- as.numeric(table(unlist(stl)))
  lines( fq , type = 'b', col = 'red', pch = 16)
  
  chi.obs <- chisq.test(fq)$statistic
  ppp <- mean(chi.pp >= chi.obs)
  return(ppp = ppp)
}




