
simu_mu <- function(L,p,n){
  W <- matrix(0,p,1)
  Wrow <- list()
  edge <- matrix(c(NA,NA),1)
  for(i in 1:L){
    tmp <- rep(0,p)
    index <- sample(1:p,50)
    Wrow[[i]] <- index
    graph <- sample(index,9)
    tmp_edge <- matrix(c(NA,NA),1)
    for(k in 1:(length(graph)-1)){
      tmp_edge <- rbind(tmp_edge,c(graph[k],graph[k+1]))
    }
    tmp_edge <- tmp_edge[-1,]
    for(j in index){
      tmp[j] <- rnorm(1,1.5,0.1)
      flag <- rbinom(1,1,0.5)
      if(flag)tmp[j] <- -tmp[j]
    }
    W <- cbind(W,tmp)
    edge <- rbind(edge,tmp_edge)
  }
  W <- W[,-1]
  edge <- edge[-1,]
  
  Z <- matrix(0,1,n)
  Zcol <- list()
  colnum <- rpois(L,30)
  for(i in 1:L){
    tmp <- rep(0,n)
    index <- sample(1:n,colnum[i])
    Zcol[[i]] <- index
    for(j in 1:length(index)){
      tmp[index[j]] <- rnorm(1,1.5,0.1)
      flag <- rbinom(1,1,0.5)
      if(flag)tmp[index[j]] <- -tmp[index[j]]
    }
    Z <- rbind(Z,tmp)
  }
  Z <- Z[-1,]
  
  mu <- W%*%Z
  S <- list()
  for(i in 1:L){
    S[[i]] <- list(r=Wrow[[i]],c=Zcol[[i]])
  }
  
  return(list(W=W,Z=Z,mu=mu,S=S,edge=edge))
}

simu_nb_x <- function(mu){
  p <- nrow(mu)
  n <- ncol(mu)
  X <- matrix(0,p,n)
  for(i in 1:p){
    r <- sample(5:20,1)
    for(j in 1:n){
      prob <- 1/(1+exp(mu[i,j]))
      X[i,j] <- rnbinom(1,r,prob)
    }
  }
  return(X)
}

simu_noise_batch <- function(X){
  interval <- ncol(X)/3
  m <- mean(X)/4
  
  for(i in 1:nrow(X)){
    for(j in 1:interval){
      X[i,j] <- X[i,j]+rnorm(1,0,1)
      if(rbinom(1,1,0.4) | X[i,j]<0)X[i,j] <- 0
    }
    
    for(j in (interval+1):(2*interval)){
      X[i,j] <- X[i,j]+rnorm(1,0,2)
      if(rbinom(1,1,0.4) | X[i,j]<0)X[i,j] <- 0
    }
    
    for(j in (2*interval+1):(3*interval)){
      X[i,j] <- X[i,j]+rnorm(1,0,4)
      if(rbinom(1,1,0.4) | X[i,j]<0)X[i,j] <- 0
    }
  }
  
  return(X)
}
