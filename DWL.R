# DWL Algorithm
# written by Changgee Chang
# ver 20160329
#
# Usage:
# Set the global variables DWL.X and DWL.y
# Call DWL.Init() only once
# Then call DWL(lam) whenever lam changes
#
# lam: the penalty vector
#

DWL.Init <- function(X, y)
{
  DWL.X <<- X
  DWL.y <<- y
  DWL.n <<- nrow(DWL.X)
  DWL.p <<- ncol(DWL.X)
  DWL.m <<- min(DWL.n,DWL.p)
  
  DWL.XxMax <<- min(3*DWL.n,DWL.p)
  DWL.Xx <<- matrix(0,DWL.p,DWL.XxMax)
  DWL.XxIdx <<- rep(0,DWL.p)
  DWL.XxCnt <<- 0
  
  DWL.Xy <<- as.vector(t(DWL.X) %*% DWL.y)
  DWL.yy <<- sum(DWL.y^2)
  
  DWL.lam <<- rep(max(abs(DWL.Xy))*1.2,DWL.p)
  DWL.A <<- c()
  DWL.nA <<- 0
  DWL.B <<- c()
  DWL.S <<- c()
  
  DWL.C <<- DWL.Xy
  DWL.iXXa <<- matrix(0,DWL.m,DWL.m)
  DWL.Idx <<- rep(0,DWL.p)
}




DWL <- function(lam)
{
  # Eat free lunch
  for ( i in 1:DWL.p )
  {
    if ( DWL.Idx[i] == 0 & DWL.lam[i] < lam[i] )
      DWL.lam[i] <<- lam[i]
  }
  
  niter = 0
  repeat
  {
    niter = niter + 1
    dlam = lam - DWL.lam
    
    # calculate dB/dalpha and dC/dalpha
    if ( DWL.nA > 0 )
    {
      dB = -DWL.iXXa[1:DWL.nA,1:DWL.nA] %*% (DWL.S * dlam[DWL.A])
      dC = -DWL.getXXa(DWL.A) %*% dB
    }
    else
      dC = rep(0,DWL.p)
    
    # find breakpoint
    alpha = 1
    
    if ( DWL.nA > 0 )
    {
      pbp0 = -DWL.B/dB
      for ( l in 1:DWL.nA )
        if ( (DWL.B[l]+dB[l])*DWL.S[l] < 0 & pbp0[l] < alpha )
        {
          alpha = pbp0[l]
          type = 0
          idx = DWL.A[l]
        }
    }
    
    pbp1 = (DWL.lam-DWL.C)/(dC-dlam)
    pbp2 = -(DWL.lam+DWL.C)/(dC+dlam)
    for ( k in 1:DWL.p )
      if ( DWL.Idx[k] == 0 )
      {
        if ( DWL.C[k]+dC[k] > DWL.lam[k]+dlam[k] & pbp1[k] < alpha )
        {
          alpha = pbp1[k]
          type = 1
          idx = k
        }
        if ( DWL.C[k]+dC[k] < -DWL.lam[k]-dlam[k] & pbp2[k] < alpha )
        {
          alpha = pbp2[k]
          type = -1
          idx = k
        }
      }
    
    
    # Add or remove var
    if ( alpha < 1 )
      if ( type == 0 )
        DWL.remove(idx)
      else
        DWL.add(idx,type)

        
    # compute B and C at alpha with new A and S
    DWL.lam <<- DWL.lam + dlam*alpha
    if ( DWL.nA > 0 )
    {
      DWL.B <<- DWL.iXXa[1:DWL.nA,1:DWL.nA] %*% ( DWL.Xy[DWL.A] - DWL.S*DWL.lam[DWL.A] )
      DWL.C <<- DWL.Xy - DWL.getXXa(DWL.A) %*% DWL.B
    }
    else
    {
      DWL.B <<- c()
      DWL.C <<- DWL.Xy
    }
    
    if ( alpha ==  1 )
      break
  }
  
  coef = rep(0,DWL.p)
  coef[DWL.A] = DWL.B

  list(coef=coef,niter=niter)
}

DWL.add <- function(k,sgn)
{
  b = DWL.getXXa(k)
  m = DWL.nA
  
  if ( m > 0 )
  {
    a = DWL.iXXa[1:m,1:m]
    del = drop(a %*% b[DWL.A])
    d = drop(b[k] - crossprod(del,b[DWL.A]))
    
    # if ( d < 1e-5 )
      # print("Warning: numerical instability")
  }
  
  # Now add k
  if ( m > 0 )
  {
    DWL.iXXa[1:m,1:m] <<- a + del %*% t(del) / d
    DWL.iXXa[1:m,m+1] <<- -del / d
    DWL.iXXa[m+1,1:m] <<- -del / d
    DWL.iXXa[m+1,m+1] <<- 1/d
  }
  else
  {
    DWL.iXXa[1] <<- 1/b[k]
  }
  DWL.Idx[k] <<- m+1
  DWL.nA <<- m+1
  DWL.A <<- c(DWL.A,k)
  DWL.S <<- c(DWL.S,sgn)
}

DWL.add_iXXa <- function(k)
{
  b = DWL.getXXa(k)
  
  m = DWL.nA
  if ( m > 0 )
  {
    a = DWL.iXXa[1:m,1:m]
    del = drop(a %*% b[DWL.A])
    d = drop(b[k] - crossprod(del,b[DWL.A]))
    
    DWL.iXXa[1:m,1:m] <<- a + del %*% t(del) / d
    DWL.iXXa[1:m,m+1] <<- -del / d
    DWL.iXXa[m+1,1:m] <<- -del / d
    DWL.iXXa[m+1,m+1] <<- 1/d
  }
  else
  {
    DWL.iXXa[1] <<- 1/b[k]
  }
  DWL.Idx[k] <<- m+1
  DWL.nA <<- m+1
  DWL.A <<- c(DWL.A,k)
}


DWL.remove_iXXa <- function(k)
{
  l = DWL.Idx[k]
  m = DWL.nA
  DWL.Idx[k] <<- 0
  if ( l<m )
    DWL.Idx[DWL.A[(l+1):m]] <<- DWL.Idx[DWL.A[(l+1):m]] - 1
  DWL.nA <<- m-1
  DWL.A <<- DWL.A[-l]
  
  if ( m>1 )
  {
    a = DWL.iXXa[1:m,1:m]
    b = a[,l]
    DWL.iXXa[1:(m-1),1:(m-1)] <<- a[-l,-l] - b[-l] %*% t(b[-l]) / b[l]
  }
  
  DWL.iXXa[,m] <<- 0
  DWL.iXXa[m,] <<- 0
}


DWL.remove <- function(k)
{
  l = DWL.Idx[k]
  m = DWL.nA
  DWL.Idx[k] <<- 0
  if ( l<m )
    DWL.Idx[DWL.A[(l+1):m]] <<- DWL.Idx[DWL.A[(l+1):m]] - 1
  DWL.nA <<- m-1
  DWL.A <<- DWL.A[-l]
  DWL.S <<- DWL.S[-l]
  
  if ( m>1 )
  {
    a = DWL.iXXa[1:m,1:m]
    b = a[,l]
    DWL.iXXa[1:(m-1),1:(m-1)] <<- a[-l,-l] - b[-l] %*% t(b[-l]) / b[l]
  }
  
  DWL.iXXa[,m] <<- 0
  DWL.iXXa[m,] <<- 0
}


DWL.getXXa <- function(A)
{
  for ( k in A )
    if ( DWL.XxIdx[k] == 0 )
    {
      DWL.XxCnt <<- DWL.XxCnt + 1
      if ( DWL.XxCnt > DWL.XxMax )
      {
        oldmax = DWL.XxMax
        oldXx = DWL.Xx
        DWL.XxMax <<- min(oldmax*2,p)
        DWL.Xx <<- matrix(0,p,DWL.XxMax)
        DWL.Xx[0,1:oldmax] <<- oldXx
      }
      DWL.XxIdx[k] <<- DWL.XxCnt
      DWL.Xx[,DWL.XxCnt] <<- t(DWL.X) %*% DWL.X[,k]
    }
  
  DWL.Xx[,DWL.XxIdx[A]]
}


DWL.version <- function()
{
  print("0.1")
}