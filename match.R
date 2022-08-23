#compute number of match items
Match <- function(res,dat){
  row_inter <- length(intersect(res[[1]],dat[[1]]))
  col_inter <- length(intersect(res[[2]],dat[[2]]))
  return(row_inter*col_inter)
}


####################
#for evaluation
####################

#CE
###########

#calculate super grid
calcu_supergrid <- function(n,p,L,orig.S,res.S){
  aux_matrix <- matrix(0,n,p)
  for(i in 1:L){
    tmp1 <- orig.S[[i]]
    tmp2 <- res.S[[i]]
    for(j in 1:length(tmp1[[1]])){
      for(k in 1:length(tmp1[[2]]))aux_matrix[tmp1[[1]][j],tmp1[[2]][k]] <- aux_matrix[tmp1[[1]][j],tmp1[[2]][k]]+1
    }
    for(j in 1:length(tmp2[[1]])){
      for(k in 1:length(tmp2[[2]]))aux_matrix[tmp2[[1]][j],tmp2[[2]][k]] <- aux_matrix[tmp2[[1]][j],tmp2[[2]][k]]+1
    }
  }
  super <- sum(aux_matrix!=0)
  super
}

#compute the best match and corresponding CE score
CE_eval <- function(origin.S,res.S,super_grid,full_arrange){
  CEscore <- rep(NA,nrow(full_arrange))
  for(i in 1:nrow(full_arrange)){
    score <- rep(0,5)
    for(j in 1:5){
      score[j] <- Match(res.S[[j]],origin.S[[full_arrange[i,j]]])
    }
    CEscore[i] <- sum(score)/super
  }
  max(CEscore)
  full_arrange[which.max(CEscore),]
  return(list(best_permutation_for_original_data=full_arrange[which.max(CEscore),],
              CEscore=max(CEscore)))
}





###########
#SEN
calcu_sen <- function(orig.S,res.S,L){
  origr <- c()
  origc <- c()
  resr <- c()
  resc <- c()
  for(i in 1:L){
    origr <- union(orig.S[[i]][[1]],origr)
    origc <- union(orig.S[[i]][[2]],origc)
    resr <- union(res.S[[i]][[1]],resr)
    resc <- union(res.S[[i]][[2]],resc)
  }
  length(intersect(origr,resr))*length(intersect(origc,resc))/(length(origr)*length(origc))
}








###########
#SPE
calcu_spe <- function(p,n,L,orig.S,res.S){
  aux_matrix <- matrix(0,p,n)
  for(i in 1:L){
    tmp1 <- orig.S[[i]]
    tmp2 <- res.S[[i]]
    aux_matrix[tmp1[[1]],tmp1[[2]]] <- aux_matrix[tmp1[[1]],tmp1[[2]]]-100
    aux_matrix[tmp2[[1]],tmp2[[2]]] <- aux_matrix[tmp2[[1]],tmp2[[2]]]+1
  }
  sum(aux_matrix==0)/sum(aux_matrix>=0)
}









