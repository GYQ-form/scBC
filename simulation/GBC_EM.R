# setwd("/Users/Ziyi/Box Sync/emory/GBC/")
# setwd("/Users/Ziyi/Box Sync/emory/GBC/code_revise_1004")
# library("glmnet")
source("DWL.R")


GBC_EM <- function(X, L, dist.type, param, use.network, edge, init="svd", maxiter=100000, tol=1e-5, nu1=1, nu4=1, cutoff=0.1){
  
  # Function to perform EM algorithm in GBC
  # X:           observed p * n matrix, p is the number of measurements and n is the number of subjects.  
  # L:           desired number of factors
  # dist.type:   a p-element vector, each element represents the distribution type of measurements.  
  #              0 is Gaussian, 1 is Binomial, 2 is Negative Binomial, 3 is Poisson.
  # param:       a p-element vector, each element is the required parameter for the correspondence distribution.
  #              Zeta for Gaussian, n_j for Binomial, r_j for Negative Binomial, N for Poisson.
  # use.network: whether or not use network information.  If specified FALSE, a randomly generated lambda will be used.
  # edge:    provide network information of p measurements.
  # init:        select between "random" and "SVM".  
  #              If "random", randomly generated N(0,1) numbers are used, o.w. use SVD results for initialization.
  # maxiter:     allowed maximum number of iterations
  # tol:         desired total difference of solutions
  
  # start_time <- Sys.time()
  record_convg <- rep(0,maxiter)
  set.seed(1234)
  
  # Initialization
  p <- dim(X)[1]
  n <- dim(X)[2]
  
  ## Specify shrinkage parameters
  kxi <- matrix(runif(L*n, 0.8, 1.5), L, n)
  lambda <- matrix(runif(p*L, 0.8, 1.5), p, L)
  alpha <- log(lambda)
  nu1 <- nu1
  nu2 <- log(2)
  nu3 <- 1
  nu4 <- nu4
  a_omega <- 4
  b_omega <- 1
  
  ## initialize Omega
  if(use.network){
    Omega = matrix(0, p, p)
    for(i in 1:nrow(edge)){
      if(edge[i,1]!=edge[i,2]){
        Omega[edge[i,1], edge[i,2]] = -a_omega/b_omega
        Omega[edge[i,2],edge[i,1]] = Omega[edge[i,1],edge[i,2]]
      }
    }
    for(j in 1:p){
      Omega[j,j] <- 1 - sum(Omega[j,-j])
    }
    Omega_ori <- Omega
  }else{
    Omega <- diag(1, p, p)
    Omega_ori <- Omega
  }
  
  ## process X
  X_mu=init_func(X,dist.type,param)

  ## Initialize m
  m <- apply(X_mu, 1, median)
  
  ## Initialize b  #Pólya-Gamma的时候用到
  b <- matrix(-1, p, n)
  b[dist.type==1, ] <- matrix(rep(param[dist.type==1],n), length(param[dist.type==1]), n)
  b[dist.type==2, ] <- matrix(rep(param[dist.type==2],n), length(param[dist.type==2]), n) + X[dist.type==2, ]
  b[dist.type==3, ] <- matrix(rep(param[dist.type==3],n), length(param[dist.type==3]), n)
  
  ## Initialize k
  k <- matrix(0, p, n)
  k[dist.type==1, ] <- X[dist.type==1, ] - b[dist.type==1, ]/2
  k[dist.type==2, ] <- (X[dist.type==2, ] - matrix(rep(param[dist.type==2],n), length(param[dist.type==2]), n))/2
  k[dist.type==3, ] <- X[dist.type==3, ] - b[dist.type==3, ]/2
  
  k.df = as.data.frame(k)
  k.t.df = as.data.frame(t(k))
  
  ## Initialize psi
  psi <- matrix(0, p, n)
  psi[dist.type==0, ] <- X[dist.type==0, ]
  psi[dist.type==3, ] <- log(matrix(rep(param[dist.type==3],n), length(param[dist.type==3]), n))
  
  psi.df = as.data.frame(psi)
  psi.t.df = as.data.frame(t(psi))
  
  ## Initialize W and Z
  if(init == "random"){
    Wvec <- rnorm(p*L, 0, 0.1)
    W <- matrix(Wvec, p, L)
    Zvec <- rnorm(L*n, 0, 0.1)
    Z <- matrix(Zvec, L, n)
  }else if(init == "svd"){
    svd.res <- svd(X_mu - m %*% matrix(1,1,n),L,L)
    W <- svd.res$u %*% diag(sqrt(svd.res$d),L,L)
    Z <- diag(sqrt(svd.res$d),L,L) %*% t(svd.res$v)
  }
  
  if(L==1){
    W <- matrix(W, p, 1)
    Z <- matrix(Z, 1, n)
  }
  
  ## Initialize miu
  miu <- m %*% matrix(1,1,n) + W %*% Z
  
  # EM iterations
  diff <- 1
  iter <- 0
  likelih <- -1
  # record_convg <- rep(0, maxiter)
  
  rho=rho1=rho2=matrix(0, p, n)
  
  while(diff>tol & iter<maxiter){
    
    likelih_old = likelih
    # W_old = W
    
    # E-step for rho
    rho1 <- b*(exp(miu)-exp(psi))/2/(miu-psi)/(exp(miu)+exp(psi))
    rho2 <- matrix(rep((param+n)/(param+rowSums((X-miu)^2)),n),p,n)
    rho3 <- b*exp(miu)/2/((miu-psi+1)*exp(miu)+exp(psi))
    rho[(dist.type>0),] <- rho1[(dist.type>0),]
    rho[(dist.type==0),] <- rho2[(dist.type==0),]
    rho[is.na(rho)] <- rho3[is.na(rho)]
    
    rho.df = as.data.frame(rho)
    rho.t.df = as.data.frame(t(rho))
    
    # M-step for Z
    Zfunc <- function(rho,k,psi,kxi){
      DWL.Init(sqrt(rho) * W, (1/sqrt(rho)) * (k + rho * (psi-m)))
      DWL(lam=kxi)$coef
    }
    Z <- mapply(Zfunc,rho.df,k.df,psi.df,as.data.frame(kxi))
    
    # M-step for W
    mMatrix <- matrix(rep(m,n),n,p,byrow=T)
    Wfunc <- function(rho,k,psi,lambda,mMatrix){
      DWL.Init(sqrt(rho) * t(Z), (1/sqrt(rho)) * (k + rho * (psi-mMatrix)))
      DWL(lam = lambda)$coef
    }
    W <- t(mapply(Wfunc,rho.t.df,k.t.df,psi.t.df,as.data.frame(t(lambda)),as.data.frame(mMatrix)))
    
    # M-step for m
    mfunc <- function(W,rho,psi,k){
      (-colSums(W*Z)%*%rho + psi%*%rho + sum(k))/sum(rho)
    }
    m <- mapply(mfunc,as.data.frame(t(W)),rho.t.df,psi.t.df,k.t.df)
    
    # update miu
    miu <- m + W %*% Z
    
    # M-step for alpha
    alphafunc <- function(lambda,W,alpha){
      alpha+ chol2inv(chol(diag(exp(alpha)*abs(W),p,p) + Omega/nu2)) %*% (1 - lambda*abs(W) - Omega%*%(alpha-nu1)/nu2)
    }
    alpha <- mapply(alphafunc,as.data.frame(lambda),as.data.frame(W),as.data.frame(alpha))
    lambda <- exp(alpha)
    
    # M-step for kxi
    kxi <- nu3/(abs(Z)+1/nu4)
    
    # E-step for omega
    Omega[Omega_ori==0]=0
    for(j in 1:p){
      for(t in j:p){
        if(j != t & Omega_ori[j,t] != 0){  # do not update the diagonal element now
          Omega[j,t] <- -2*nu2*a_omega/(2*nu2*b_omega + sum((alpha[j,]-alpha[t,])^2))
          Omega[t,j] <- Omega[j,t]
        }
      }
    }
    for(i in 1:p){
      Omega[i,i] <- 1-sum(Omega[i,-i])
    }
    
    likelih <- 0

    if(all(Omega[lower.tri(Omega)] == 0, Omega[upper.tri(Omega)] == 0)){
      like_omega=0
    }else{
      tmp = -(Omega - diag(diag(Omega),p,p))
      lomega = log(tmp)
      lomega[abs(lomega)==Inf] = 0
      like_omega = (a_omega-1)*sum(lomega) - b_omega*sum(tmp)
    }

    like_alpha <- sum(diag(t(alpha - nu1) %*% Omega %*% (alpha - nu1)))
    likelih <- - sum(rho*(miu-psi)^2)/2 + sum(k*miu) + sum(alpha) - sum(lambda*abs(W)) +
      nu3*sum(log(kxi)) - sum(kxi*abs(Z)) - like_alpha/2/nu2 + like_omega -sum(kxi)/nu4
    # 
    # likelih <- - sum(rho*(miu-X)^2)/2 + sum(k*miu) - sum(param*rho[,1]/2) + sum(log(rho[,1])*(p/2+param/2-1)) + 
      # sum(alpha) - sum(lambda*abs(W)) - sum(log(kxi))/2 - sum(kxi*abs(Z)) - like_alpha/2/nu2 + 
      # like_omega + n*L*log(nu3/nu4) - sum((kxi-nu3)^2/2/nu3/kxi)
    
    if(is.na(likelih) | likelih == -Inf | likelih == Inf)
      break
    # diffW <- max(abs(W - W_old))
    # print(diffW)

    iter <- iter + 1
    
    diff = abs((likelih - likelih_old)/likelih)
    print(diff)
    # record_convg[iter] = likelih
  }
  
  S <- list()
  for(l in 1:L){
    S[[l]] <- list(r=which(abs(W[,l])>cutoff),c=which(abs(Z[l,])>cutoff))
  }
  
  # print(Sys.time() - start_time)
  return(list(Z=Z, W=W, mu=miu, m=m, S=S, convg=record_convg, llh=likelih))
}



init_func <- function(X,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  Y = matrix(0,p,n)
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      Y[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      pbar = pmin(pmax(X[j,],1/3),param[j]-1/3)/param[j]
      Y[j,] = log(pbar/(1-pbar))
    }
    if ( type[j] == 2 )
    {
      pbar = pmax(X[j,],1/3)/(param[j]+pmax(X[j,],1/3))
      Y[j,] = log(pbar/(1-pbar))
    }
    if ( type[j] == 3 )
    {
      Y[j,] = log(pmax(X[j,],1))
    }
  }
  Y
}

