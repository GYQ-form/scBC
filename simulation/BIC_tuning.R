llk <- function(X,mu,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  l = 0
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
      l = l - n*log(mean((X[j,]-mu[j,])^2))/2
    else if ( type[j] == 1 )
      l = l + sum(dbinom(X[j,],param[j],1/(1+exp(-mu[j,])),TRUE))
    else if ( type[j] == 2 )
      l = l + sum(dnbinom(X[j,],param[j],1/(1+exp(mu[j,])),log=TRUE))
    else if ( type[j] == 3 )
      l = l + sum(dpois(X[j,],exp(mu[j,]),TRUE))
  }
  l
}

BIC_tuning <- function(dat,res){
  n = dim(res$mu)[2]
  p = dim(res$mu)[1]
  
  BIC_res = -2*llk(as.matrix(dat$X),as.matrix(res$mu),dat$type,dat$param) + (sum(abs(res$W)>0)+sum(abs(res$Z)>0))*log(n*p)
  
  return(BIC_res)
}




