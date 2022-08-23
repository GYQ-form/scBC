#plaid and Xmotif

get_plaid_res <- function(X,L){
  library(biclust)
  plaid_res <- biclust(dat$X,method = BCPlaid())
  res <- list()
  tmp <- list(r=integer(),c=integer())
  for(i in 1:L)res[[i]] <- tmp
  if(plaid_res@Number==0)return(res)
  for(i in 1:plaid_res@Number){
    res[[i]] <- list(r=(which(plaid_res@RowxNumber[,i])),
                     c=(which(plaid_res@NumberxCol[i,])))
  }
  return(res)
}


get_Xmotif_res <- function(X,L){
  library(biclust)
  xmotif_res <- biclust(noiseX,method = BCXmotifs(),number=L)
  res <- list()
  tmp <- list(r=integer(),c=integer())
  for(i in 1:L)res[[i]] <- tmp
  if(xmotif_res@Number==0)return(res)
  for(i in 1:xmotif_res@Number){
    res[[i]] <- list(r=(which(xmotif_res@RowxNumber[,i])),
                     c=(which(xmotif_res@NumberxCol[i,])))
  }
  return(res)
}


