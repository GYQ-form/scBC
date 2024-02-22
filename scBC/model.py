import numpy as np
from numpy.linalg import cholesky, inv
from .DWL import DWL
import warnings
from scvi.model import SCVI
from anndata import AnnData
from pandas import DataFrame
from typing import List, Optional
from biomart import BiomartServer

class scBC:

    def __init__(self,
                 adata: AnnData,
                 layer: Optional[str] = None,
                 batch_key: Optional[str] = None,
                 labels_key: Optional[str] = None,
                 size_factor_key: Optional[str] = None,
                 categorical_covariate_keys: Optional[List[str]] = None,
                 continuous_covariate_keys: Optional[List[str]] = None) -> None:
        '''
        Parameters
        -----------
        adata:       AnnData object with an observed n * p matrix, p is the number of measurements and n is the number of subjects.  
        layer:       if not None, uses this as the key in adata.layers for raw count data.
        batch_key:    key in adata.obs for batch information. Categories will automatically be converted into integer categories and
                      saved to adata.obs['_scvi_batch']. If None, assigns the same batch to all the data.
        labels_key:   key in adata.obs for label information. Categories will automatically be converted into integer categories and 
                      saved to adata.obs['_scvi_labels']. If None, assigns the same label to all the data.
        size_factor_key:  key in adata.obs for size factor information. Instead of using library size as a size factor, the provided
                          size factor column will be used as offset in the mean of the likelihood. Assumed to be on linear scale.
        categorical_covariate_keys:  keys in adata.obs that correspond to categorical data. These covariates can be added in addition
                                     to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects 
                                     on the latent space). Thus, these should not be used for biologically-relevant factors that you do not want to correct for.
        continuous_covariate_keys:  keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch
                                    covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). 
                                    Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.
        '''   

        self.adata = adata
        self.p = adata.X.shape[1]
        self.n = adata.X.shape[0]

        SCVI.setup_anndata(adata, 
                           layer=layer, 
                           batch_key=batch_key, 
                           labels_key=labels_key, 
                           size_factor_key=size_factor_key, 
                           categorical_covariate_keys=categorical_covariate_keys, 
                           continuous_covariate_keys=continuous_covariate_keys)
        
        self.edge = None
        self.vi_model = None
        self.reconst_data = None
        self.W = None
        self.Z = None
        self.mu = None
        self.convg_trail = None
        self.S = None
        self.niter = None

    
    def get_edge(self, gene_list=None, dataset='hsapiens_gene_ensembl', intensity=5):
        '''
        Parameters
        ----------
        gene_list
            A list containing HGNC gene names used to find prior edge in hsapiens_gene_ensembl dataset. 
            Different actions will be performed based on the data type you provide:
                str --- take self.adata.var[gene_list] as the gene set
                list --- use the provided list as the gene set to be searched
                None(default) --- use the variable names(self.adata.var_names) as the gene set

        dataset
            The dataset to be used for extracting prior information. Default is 'hsapiens_gene_ensembl'.

        intensity
            A integer denoting the intensity of the prior. Generally speaking, the larger the number, the more
            prior information returned (the bigger the edge array is).
        '''
        if isinstance(gene_list, str):
            genes = self.adata.var[gene_list].to_list()
        elif gene_list is not None:
            genes = gene_list
        else:
            genes = self.adata.var_names.to_list()

        server = BiomartServer("http://www.ensembl.org/biomart")

        # Connect to the biomart server and select the dataset
        mart = server.datasets[dataset]
        pathway = []
        gene_idx = 0
        while gene_idx < len(genes):
            filters = {"hgnc_symbol": genes[gene_idx:gene_idx+800]}
            attributes = ["go_id", "hgnc_symbol"]
            response = mart.search({'attributes':attributes, 'filters':filters})
            for line in response.iter_lines():
                line = line.decode('utf-8')
                pathway.append(line.split("\t"))
            gene_idx += 800

        pathway = DataFrame(data=pathway, columns=['go_id','gene'])

        # Convert the response to a pandas dataframe and remove rows with missing values
        pathway = pathway.drop(pathway[pathway['go_id']==''].index)

        # Count the number of genes per pathway
        path_num = pathway["go_id"].value_counts()
        path_num = path_num[(path_num > 10) & (path_num < 100)]

        # Initialize the edge matrix and pool of genes
        edge = []
        pool = np.array([], dtype=int)

        # Select random pathways and add edges between genes in the pathway
        for i in range(intensity):
            path_id = np.random.choice(path_num.index)
            p = pathway[pathway["go_id"] == path_id]
            path = np.where(np.isin(genes, p["gene"]))[0]
            if np.intersect1d(pool, path).size:
                continue
            pool = np.concatenate((pool, path))
            for j in range(len(path)):
                for k in range(len(path)):
                    if k == j:
                        continue
                    pair = [path[j], path[k]] 
                    edge.append(pair)

        self.edge = np.array(edge, dtype=np.int64)


    def train_VI(self, n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1, 
                 dispersion='gene', gene_likelihood='zinb', latent_distribution='normal',
                 max_epochs=None, accelerator='auto', devices='auto', batch_size=128, early_stopping=False):
        '''
        Parameters
        ----------
        n_hidden
            Number of nodes per hidden layer.
        n_latent
            Dimensionality of the latent space.
        n_layers
            Number of hidden layers used for encoder and decoder NNs.
        dropout_rate
            Dropout rate for neural networks.
        dispersion
            One of the following:
            * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
            * ``'gene-batch'`` - dispersion can differ between different batches
            * ``'gene-label'`` - dispersion can differ between different labels
            * ``'gene-cell'`` - dispersion can differ for every gene in every cell
        gene_likelihood
            One of:
            * ``'nb'`` - Negative binomial distribution
            * ``'zinb'`` - Zero-inflated negative binomial distribution
            * ``'poisson'`` - Poisson distribution
        latent_distribution
            One of:
            * ``'normal'`` - Normal distribution
            * ``'ln'`` - Logistic normal distribution (Normal(0, I) transformed by softmax)
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        accelerator
            Supports passing different accelerator types (“cpu”, “gpu”, “tpu”, “ipu”, “hpu”, “mps, “auto”)
            as well as custom accelerator instances.
        devices 
            The devices to use. Can be set to a non-negative index (int or str), a sequence of device indices (list or comma-separated str), 
            the value -1 to indicate all available devices, or “auto” for automatic selection based on the chosen accelerator. 
            If set to “auto” and accelerator is not determined to be “cpu”, then devices will be set to the first available device.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping.
        '''

        self.vi_model = SCVI(adata = self.adata, n_hidden=n_hidden, n_latent=n_latent, n_layers=n_layers, dropout_rate=dropout_rate, 
                 dispersion=dispersion, gene_likelihood=gene_likelihood, latent_distribution=latent_distribution)
        
        self.vi_model.train(max_epochs=max_epochs, accelerator=accelerator, devices=devices, batch_size=batch_size, early_stopping=early_stopping)


    def get_reconst_data(
        self,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
    ):

        '''
        Parameters
        ----------
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model.
        '''

        self.reconst_data = self.vi_model.get_normalized_expression(adata=self.adata, library_size="latent", 
        n_samples=n_samples, batch_size=batch_size, return_mean=True, return_numpy=True)
    

    def __initX(self, X, dist_type, param):
        Y = np.zeros((self.p, self.n))
        for j in range(self.p):
            if dist_type[j] == 0:
                Y[j,:] = X[j,:]
            elif dist_type[j] == 1:
                pbar = np.minimum(np.maximum(X[j,:], 1/3), param[j] - 1/3) / param[j]
                Y[j,:] = np.log(pbar / (1 - pbar))
            elif dist_type[j] == 2:
                pbar = np.maximum(X[j,:], 1/3) / (param[j] + np.maximum(X[j,:], 1/3))
                Y[j,:] = np.log(pbar / (1 - pbar))
            elif dist_type[j] == 3:
                Y[j,:] = np.log(np.maximum(X[j,:], 1))
        return Y

    
    def Biclustering(self, L, mat=None, edge=None, dist_type=None, param=None, initWZ="svd", maxiter=100000, tol=1e-5, nu1=7, nu4=20, cutoff=0.1):
        '''
        Parameters
        ----------
        L:          number of the maximum biclusters
        mat:         a p x n numpy array used for biclustering. We will use the reconstructed data after calling scBC.get_reconst_data(). 
                    However, you can explicitly provide one rather than using the reconstructed data.
        edge:        a 2-colounm matrix providing prior network information of some measurements. 
                    We prefer to use the edge given here, if not provided, we will use scBC.edge as a candidate.
                    Edge can be automaticly generated using scBC.get_edge(). Running Biclustering() without any edge prior is also feasible.
        dist_type:   a p-element vector, each element represents the distribution type of measurements.  
                    0 is Gaussian, 1 is Binomial, 2 is Negative Binomial, 3 is Poisson. 
                    If not given, all measurements are deemed to follow Gaussian distribution.
        param:       a p-element vector, each element is the required parameter for the correspondence distribution.
                    Zeta for Gaussian, n_j for Binomial, r_j for Negative Binomial, N for Poisson. Default is set as a 1 vector.
        initWZ:        select between "random" and "svd".  
                    If "random", randomly generated N(0,1) numbers are used, o.w. use SVD results for initialization.
        maxiter:     allowed maximum number of iterations
        tol:         desired total difference of solutions
        '''
        
        if dist_type is None:
            dist_type = np.zeros(self.p, dtype=np.int64)
        if param is None:
            param = np.ones(self.p)
        if edge:
            use_edge = edge
        else:
            use_edge = self.edge

        if mat:
            X = mat
        else:
            X = self.reconst_data.T
        
        np.random.seed(1234)

        # Specify shrinkage parameters
        kxi = np.random.uniform(0.8, 1.5, size=(L,self.n))
        lambda_ = np.random.uniform(0.8, 1.5, size=(self.p,L))
        alpha = np.log(lambda_)
        nu1 = nu1
        nu2 = np.log(2)
        nu3 = 1
        nu4 = nu4
        a_omega = 4
        b_omega = 1

        ## initialize Omega
        if use_edge is not None:
            Omega = np.zeros((self.p, self.p))
            for i in range(use_edge.shape[0]):
                if use_edge[i,0] != use_edge[i,1]:
                    Omega[use_edge[i,0], use_edge[i,1]] = -a_omega/b_omega
                    Omega[use_edge[i,1], use_edge[i,0]] = Omega[use_edge[i,0], use_edge[i,1]]
            for j in range(self.p):
                Omega[j,j] = 1 - np.sum(Omega[j,:]) + Omega[j,j]
            Omega_ori = Omega
        else:
            Omega = np.eye(self.p, self.p)
            Omega_ori = Omega

        ## process X 
        X_mu = self.__initX(X, dist_type, param)

        ## Initialize m
        m = np.apply_along_axis(np.median, 1, X_mu).reshape(-1,1)

        ## Initialize b
        b = np.full((self.p, self.n), -1)
        b[dist_type==1, :] = np.tile(param[dist_type==1], (self.n, 1)).T
        b[dist_type==2, :] = np.tile(param[dist_type==2], (self.n, 1)).T + X[dist_type==2, :]
        b[dist_type==3, :] = np.tile(param[dist_type==3], (self.n, 1)).T

        ## Initialize k
        k = np.zeros((self.p, self.n))
        k[dist_type==1, :] = X[dist_type==1, :] - b[dist_type==1, :]/2
        k[dist_type==2, :] = (X[dist_type==2, :] - np.tile(param[dist_type==2], (self.n, 1)).T)/2
        k[dist_type==3, :] = X[dist_type==3, :] - b[dist_type==3, :]/2

        ## Initialize psi
        psi = np.zeros((self.p, self.n))
        psi[dist_type==0, :] = X[dist_type==0, :]
        psi[dist_type==3, :] = np.log(np.tile(param[dist_type==3], (self.n, 1)).T)

        ## Initialize W and Z
        if initWZ == "random":
            W = np.random.normal(loc=0, scale=0.1, size=(self.p, L))
            Z = np.random.normal(loc=0, scale=0.1, size=(L, self.n))
        elif initWZ == "svd":
            svd_res = np.linalg.svd(X_mu - np.outer(m.squeeze(), np.ones(self.n)), full_matrices=False)
            W = svd_res[0] @ np.diag(np.sqrt(svd_res[1]))[:, :L]
            Z = np.diag(np.sqrt(svd_res[1]))[:L, :] @ svd_res[2]

        if L == 1:
            W = np.reshape(W, (self.p, 1))
            Z = np.reshape(Z, (1, self.n))

        ## Initialize miu
        miu = np.outer(m.squeeze(), np.ones(self.n)) + W @ Z

        # EM iterations
        diff = 1
        iter = 0
        likelih = -1
        convg = []

        rho = rho1 = rho2 = np.zeros((self.p, self.n))

        while (diff > tol) and (iter < maxiter):
            likelih_old = likelih
            # W_old = W

            # E-step for rho
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rho1 = b*(np.exp(miu)-np.exp(psi))/2/(miu-psi)/(np.exp(miu)+np.exp(psi))
                rho2 = np.tile(((param+self.n)/(param+np.sum((X-miu)**2, axis=1))).reshape(-1,1), (1, self.n))
                rho3 = b*np.exp(miu)/2/((miu-psi+1)*np.exp(miu)+np.exp(psi))
                rho[(dist_type>0), :] = rho1[(dist_type>0), :]
                rho[(dist_type==0), :] = rho2[(dist_type==0), :]
                rho[np.isnan(rho)] = rho3[np.isnan(rho)]

            #similar in R
            def mapply(my_fun,*args):
                cols = args[0].shape[1]
                tmp = []
                for i in range(cols):
                    arg_list = [mat[:,i].reshape((-1,1)) for mat in args]
                    tmp.append(my_fun(*arg_list).reshape(-1,1))
                res = np.concatenate(tmp, axis=1)
                return res

            # M-step for Z
            def Zfunc(rho, k, psi, kxi):
                clf = DWL(np.sqrt(rho) * W, 1/np.sqrt(rho) * (k + rho * (psi-m)))
                clf.fit(kxi)
                return clf.coef_

            Z = mapply(Zfunc, rho, k, psi, kxi)

            # M-step for W
            mMatrix = np.tile(m, (1, self.n)).T
            def Wfunc(rho, k, psi, lambda_, mMat):
                clf = DWL(np.sqrt(rho) * Z.T, (1/np.sqrt(rho)) * (k + rho * (psi-mMat)))
                clf.fit(lambda_)
                return clf.coef_

            W = mapply(Wfunc, rho.T, k.T, psi.T, lambda_.T, mMatrix).T 

            # M-step for m
            def mfunc(W, rho, psi, k):
                return (-W.T @ Z @ rho +  psi.T @ rho + np.sum(k))/np.sum(rho)

            m = mapply(mfunc, W.T, rho.T, psi.T, k.T).reshape(-1,1)

            # update miu
            miu = m + W @ Z

            # M-step for alpha
            def alphafunc(lam, W, alpha):
                omega = Omega/nu2
                L = cholesky(np.diag((np.exp(alpha)*np.abs(W)).squeeze()) + omega)
                L_inv = inv(L)
                alpha_new = alpha + L_inv.T @ L_inv @ (1 - lam*np.abs(W) - omega @ (alpha-nu1))
                return alpha_new

            alpha = mapply(alphafunc, lambda_, W, alpha)
            lambda_ = np.exp(alpha)

            # M-step for kxi
            kxi = nu3/(np.abs(Z)+1/nu4)

            # E-step for omega
            Omega[Omega_ori==0] = 0
            for j in range(self.p):
                for t in range(j, self.p):
                    if j != t and Omega_ori[j, t] != 0:
                        Omega[j, t] = -2*nu2*a_omega/(2*nu2*b_omega + np.sum((alpha[j, :]-alpha[t, :])**2))
                        Omega[t, j] = Omega[j, t]

            for i in range(self.p):
                Omega[i, i] = 1 - np.sum(Omega[i, np.arange(self.p) != i])

            likelih = 0

            if np.all((Omega[np.tril_indices(self.p, k=-1)] == 0) & (Omega[np.triu_indices(self.p, k=1)] == 0)):
                like_omega = 0
            else:
                tmp = -(Omega - np.diag(np.diag(Omega)))
                lomega = np.log(tmp)
                lomega[np.abs(lomega) == np.inf] = 0
                like_omega = (a_omega-1)*np.sum(lomega) - b_omega*np.sum(tmp)

            like_alpha = np.sum(np.diag((alpha - nu1).T @ Omega @ (alpha - nu1)))
            likelih = - np.sum(rho*(miu - psi)**2)/2 + np.sum(k*miu) + np.sum(alpha) - np.sum(lambda_*np.abs(W)) + \
                    nu3*np.sum(np.log(kxi)) - np.sum(kxi*np.abs(Z)) - like_alpha/2/nu2 + like_omega - np.sum(kxi)/nu4

            if np.isnan(likelih) or likelih == -np.inf or likelih == np.inf:
                break

            iter = iter + 1

            diff = np.abs((likelih - likelih_old)/likelih)
            convg.append(likelih)
            S = []
            for l in range(L):
                S.append({'measurements': np.where(np.abs(W[:, l]) > cutoff)[0], 'subjects': np.where(np.abs(Z[l, :]) > cutoff)[0]})

            self.W = W
            self.Z = Z
            self.mu = miu
            self.convg_trail = convg
            self.S = S
            self.niter = iter


    def strong_class(self):
        return np.abs(self.Z).argmax(axis=0)