import anndata as ad
import requests
import os
import numpy as np
from numpy.random import binomial
from scipy.stats import nbinom

def load_data(data, **kwargs):
    '''
    Parameters
    ----------
    data
        the built-in subsampled datasets, including "HEART", "PBMC", "LUAD" and "BC"
    **kwargs
        Keyword args for `anndata.read_h5ad`
    '''


    if data == "HEART":
        if not os.path.exists("heart_sub.h5ad"):
            url = 'https://github.com/GYQ-form/scBC/raw/main/data/heart_sub.h5ad'
            response = requests.get(url)
            with open("heart_sub.h5ad", 'wb') as f:
                f.write(response.content)      
        adata = ad.read_h5ad("heart_sub.h5ad", **kwargs)

    elif data == "PBMC":
        if not os.path.exists("PBMC_sub.h5ad"):
            url = 'https://github.com/GYQ-form/scBC/raw/main/data/PBMC_sub.h5ad'
            response = requests.get(url)
            with open("PBMC_sub.h5ad", 'wb') as f:
                f.write(response.content)            
        adata = ad.read_h5ad("PBMC.h5ad", **kwargs)

    elif data == "LUAD":
        if not os.path.exists("LUAD_sub.h5ad"):
            url = 'https://github.com/GYQ-form/scBC/raw/main/data/LUAD_sub.h5ad'
            response = requests.get(url)
            with open("LUAD_sub.h5ad", 'wb') as f:
                f.write(response.content)            
        adata = ad.read_h5ad("LUAD.h5ad", **kwargs)

    elif data == "BC":
        if not os.path.exists("BC_sub.h5ad"):
            url = 'https://github.com/GYQ-form/scBC/raw/main/data/BC_sub.h5ad'
            response = requests.get(url)
            with open("BC_sub.h5ad", 'wb') as f:
                f.write(response.content)            
        adata = ad.read_h5ad("BC.h5ad", **kwargs)

    else:
        raise ValueError(f"no data named {data}!")

    return adata


def simulate_data(L, p, n, dropout, batch_num):
    '''
    Parameters
    ----------
    L
        number of biclusters
    p
        number of measurements(genes)
    n
        number of objects(cells)
    dropout
        dropout rate (between 0 and 1)
    batch_num
        number of batches you want to simulate
    '''
    W = np.zeros((p, 1))
    Wrow = []
    edge = np.array([[np.nan, np.nan]])
    indexall = np.random.choice(range(1, p+1), int(p/20*L), replace=False)
    for i in range(1, L+1):
        tmp = np.zeros(p)
        index = indexall[int(p/20*(i-1)):int(p/20*i)]
        
        # add prior
        for j in range(len(index)):
            for k in range(len(index)):
                if k == j:
                    continue
                edge = np.vstack((edge, [index[j], index[k]]))
        
        Wrow.append(list(index))
        for j in index:
            tmp[j-1] = np.random.normal(1.5, 0.1, size=1)
            flag = np.random.binomial(1, 0.5, size=1)
            if flag:
                tmp[j-1] = -tmp[j-1]
        W = np.hstack((W, tmp.reshape(-1, 1)))
    W = W[:, 1:]
    edge = edge[1:, :]
    
    Z = np.zeros((1, n))
    Zcol = []
    colnum = np.random.poisson(n/10, size=L)
    for i in range(L):
        tmp = np.zeros(n)
        index = np.random.choice(range(n), colnum[i], replace=False)
        Zcol.append(list(index))
        for j in range(len(index)):
            tmp[index[j]] = np.random.normal(1.5, 0.1, size=1)
            flag = np.random.binomial(1, 0.5, size=1)
            if flag:
                tmp[index[j]] = -tmp[index[j]]
        Z = np.vstack((Z, tmp.reshape(1, -1)))
    Z = Z[1:, :]
    
    mu = W @ Z
    S = []
    for i in range(L):
        S.append({'r': Wrow[i], 'c': Zcol[i]})
    
    #simulate X
    X = np.zeros((p, n))
    for i in range(p):
        r = np.random.randint(5, 21, size=1)
        for j in range(n):
            prob = 1 / (1 + np.exp(mu[i, j]))
            X[i, j] = nbinom.rvs(r, prob, size=1)

    #simulate noise and dropout
    interval = n // batch_num
    m = np.mean(X)*2 / batch_num
    for i in range(batch_num):
        for j in range(interval*i, interval*(i+1)):
            X[:, j] += i * m
            X[np.intersect1d(np.where(binomial(1, dropout, p))[0],np.where(X[:, j] < 0)[0]), j] = 0

    dat = ad.AnnData(X.T)
    return {'W': W, 'Z': Z, 'dat': dat, 'S': S, 'edge': edge}
