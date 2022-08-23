#permutation


permutation <- function(n){
  judge <<- rep(0,n)
  full_arrange <<- matrix(NA,1,n)
  tmp <<- rep(0,5)
  ni <<- n
  aux_permu(1)
  full_arrange <<- full_arrange[-1,]
}


aux_permu <- function(t){
  if(t>ni){
    full_arrange <<- rbind(full_arrange,tmp)
    return()
  }
  
  for(i in 1:ni){
    if(judge[i]==0){
      tmp[t] <<- i
      judge[i] <<- 1
      aux_permu(t+1)
      judge[i] <<- 0
    }
  }
  
}
