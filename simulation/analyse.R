set.seed(1)

source("GBC_EM.R")
source("DWL.R")
source("BIC_tuning.R")
source("plaid_Xmotif.R")
source("permutation.R")
source("match.R")
library(reticulate)
use_condaenv("scvi-tools", required=TRUE)
library(Seurat)

p = 1000
n = 300
type = rep(2,p)  # all variables are assumed to follow nb distribution
L = 5
mu <- simu_mu(5,1000,300)
X <- simu_nb_x(mu$mu)
noiseX <- simu_noise_batch(X)

###############
#GBC and sGBC
###############

dat <- list(X=noiseX,type=type,param=rep(1,1000),S=mu$S,edge=mu$edge)

## select nu1 and nu4 by BIC & apply sGBC

nu1_vec = c(7,9,11,13,15,20,25)
nu4_vec = c(20,40,50,60,70,90,110)

BIC_rec = matrix(999999, length(nu4_vec), length(nu1_vec))
for(i1 in 1:length(nu1_vec)){
  print(paste("i1=",i1))
  for(i2 in 1:length(nu4_vec)){
    print(paste("i2=",i2))
    set.seed(123)
    tmp = try(GBC_EM(X=as.matrix(dat$X), L=L, dist.type=dat$type, param=dat$param, use.network=F, init="svd", nu1=nu1_vec[i1], nu4=nu4_vec[i2], cutoff=0),silent=T)
    if(is.list(tmp)){
      BIC_rec[i2,i1] = BIC_tuning(dat,tmp)
    }
  }
}
minind<-which(BIC_rec == min(BIC_rec), arr.ind=TRUE)[1,]


GBC_res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=dat$type, param=dat$param, use.network=F,init="svd", nu1=nu1_vec[minind[1]], nu4=nu4_vec[minind[2]], cutoff=0)
sGBC_res = GBC_EM(X=as.matrix(dat$X), L=L, dist.type=dat$type, param=dat$param, use.network=T,edge = dat$edge,init="svd", nu1=nu1_vec[minind[1]], nu4=nu4_vec[minind[2]], cutoff=0)

#calculate number of super grid
super_GBC <- calcu_supergrid(1000,300,5,dat$S,GBC_res$S)
super_sGBC <- calcu_supergrid(1000,300,5,dat$S,sGBC_res$S)

#match respective clustering and get the evaluation results
(GBC_CE <- CE_eval(dat$S,GBC_res$S,super_GBC,full_arrange)$CEscore)
(GBC_SEN <- calcu_sen(dat$S,GBC_res$S,5))
(GBC_SPE <- calcu_spe(1000,300,5,dat$S,GBC_res$S))
(sGBC_CE <- CE_eval(dat$S,sGBC_res$S,super_GBC,full_arrange)$CEscore)
(sGBC_SEN <- calcu_sen(dat$S,sGBC_res$S,5))
(sGBC_SPE <- calcu_spe(1000,300,5,dat$S,sGBC_res$S))

###################
#plaid and Xmotif
###################

#plaid
plaid_res <- get_plaid_res(dat$X,L)
super_plaid <- calcu_supergrid(1000,300,5,dat$S,plaid_res)
(plaid_CE <- CE_eval(dat$S,plaid_res,super_plaid,full_arrange)$CEscore)
(plaid_SEN <- calcu_sen(dat$S,plaid_res,5))
(plaid_SPE <- calcu_spe(1000,300,5,dat$S,plaid_res))

#xmotif
Xmotif_res <- get_Xmotif_res(dat$X,L)
super_Xmotif <- calcu_supergrid(1000,300,5,dat$S,Xmotif_res)
(xmotif_CE <- CE_eval(dat$S,Xmotif_res,super_plaid,full_arrange)$CEscore)
(xmotif_SEN <- calcu_sen(dat$S,Xmotif_res,5))
(xmotif_SPE <- calcu_spe(1000,300,5,dat$S,Xmotif_res))


################
#scBC
################

colnames(noiseX) <- paste("cell_",seq(1,300),sep = "")
row.names(noiseX) <- paste("gene_",seq(1,1000),sep = "")
noiseX_seurat <- CreateSeuratObject(counts = round(noiseX), project = "simulation")
batch <- c(rep("1",100),rep("2",100),rep("3",100))
noiseX_seurat[["batch"]] <- batch
adata <- convertFormat(noiseX_seurat, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)

#Setup our AnnData for training
scvi$model$SCVI$setup_anndata(adata,batch_key = "batch")

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

#get denoised expression
denoised <- model$get_normalized_expression(n_samples=as.integer(100),library_size="latent")
reconst_data <- t(py_to_r(denoised))

#GBC
#param = apply(denois_data,1,para_estimate)
scBC_dat <- list(X=reconst_data,type=type,param=rep(1,p),edge=mu$edge)
scBC_res = GBC_EM(X=as.matrix(scBC_dat$X), L=L, dist.type=scBC_dat$type, param=scBC_dat$param, use.network=T,edge = scBC_dat$edge,init="svd", nu1=nu1_vec[minind[1]], nu4=nu4_vec[minind[2]], cutoff=0)

#calculate number of super grid
super_scBC <- calcu_supergrid(1000,300,5,dat$S,scBC_res$S)
#match respective clustering
(scBC_CE <- CE_eval(dat$S,scGBC_res$S,super_scBC,full_arrange)$CEscore)
(scBC_SEN <- calcu_sen(dat$S,scBC_res$S,5))
(scBC_SPE <- calcu_spe(1000,300,5,dat$S,scBC_res$S))


