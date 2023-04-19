# DWL Algorithm
# proposed by Changgee Chang

import numpy as np
import warnings

class DWL:

    def __init__(self, X, y):
        self.DWL_X = X
        self.DWL_y = y
        self.DWL_n = X.shape[0]
        self.DWL_p = X.shape[1]
        self.DWL_m = min(self.DWL_n, self.DWL_p)

        self.DWL_XxMax = min(3*self.DWL_n, self.DWL_p)
        self.DWL_Xx = np.zeros((self.DWL_p, self.DWL_XxMax))
        self.DWL_XxIdx = np.full(self.DWL_p, -1, dtype=np.int64)
        self.DWL_XxCnt = 0

        self.DWL_Xy = np.dot(X.T, y)
        self.DWL_yy = np.sum(y**2)

        self.DWL_lam = np.repeat(max(np.abs(self.DWL_Xy))*1.2, self.DWL_p).reshape(-1,1)
        self.DWL_A = np.array([],dtype=np.int32)
        self.DWL_nA = 0
        self.DWL_B = np.array([])
        self.DWL_S = np.array([])

        self.DWL_C = self.DWL_Xy
        self.DWL_iXXa = np.zeros((self.DWL_m, self.DWL_m))
        self.DWL_Idx = np.full(self.DWL_p, -1, dtype=np.int64)

        self.coef_ = None
        self.niter_ = None


    def fit(self, lam):
        for i in range(self.DWL_p):
            if self.DWL_Idx[i] == -1 and self.DWL_lam[i] < lam[i]:
                self.DWL_lam[i] = lam[i]

        niter = 0
        while True:
            niter += 1
            dlam = lam - self.DWL_lam

            if self.DWL_nA > 0:
                dB = -np.dot(self.DWL_iXXa[:self.DWL_nA,:self.DWL_nA], self.DWL_S * dlam[self.DWL_A])
                dC = -np.dot(self.__DWL_getXXa(self.DWL_A), dB)
            else:
                dC = np.zeros((self.DWL_p,1))

            alpha = 1

            if self.DWL_nA > 0:
                pbp0 = -self.DWL_B / dB
                for l in range(self.DWL_nA):
                    if (self.DWL_B[l] + dB[l]) * self.DWL_S[l] < 0 and pbp0[l] < alpha:
                        alpha = pbp0[l]
                        tp = 0
                        idx = self.DWL_A[l]
                        
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pbp1 = (self.DWL_lam - self.DWL_C) / (dC - dlam)
                pbp2 = -(self.DWL_lam + self.DWL_C) / (dC + dlam)

            for k in range(self.DWL_p):
                if self.DWL_Idx[k] == -1:
                    if self.DWL_C[k] + dC[k] > self.DWL_lam[k] + dlam[k] and pbp1[k] < alpha:
                        alpha = pbp1[k]
                        tp = 1
                        idx = k
                    if self.DWL_C[k] + dC[k] < -self.DWL_lam[k] - dlam[k] and pbp2[k] < alpha:
                        alpha = pbp2[k]
                        tp = -1
                        idx = k
            
            # Add or remove var
            if alpha < 1:
                if tp == 0:
                    self.__DWL_remove(idx)
                else:
                    self.__DWL_add(idx, tp)


            # compute B and C at alpha with new A and S
            self.DWL_lam += dlam * alpha
            if self.DWL_nA > 0:
                self.DWL_B = np.dot(self.DWL_iXXa[:self.DWL_nA,:self.DWL_nA], self.DWL_Xy[self.DWL_A] - self.DWL_S * self.DWL_lam[self.DWL_A])
                self.DWL_C = self.DWL_Xy - np.dot(self.__DWL_getXXa(self.DWL_A), self.DWL_B)
            else:
                self.DWL_B = np.array([])
                self.DWL_C = self.DWL_Xy

            if alpha == 1:
                break

        coef = np.zeros(self.DWL_p)
        coef[self.DWL_A] = self.DWL_B.squeeze()

        self.coef_ = coef
        self.niter_ = niter


    def __DWL_add(self, k, sgn):
        b = self.__DWL_getXXa(k)
        m = self.DWL_nA

        if m > 0:
            a = self.DWL_iXXa[:m, :m]
            del_val = np.dot(a, b[self.DWL_A]).squeeze()
            d = b[k] - np.dot(del_val, b[self.DWL_A]).squeeze()
            #add k
            self.DWL_iXXa[:m, :m] += np.outer(del_val, del_val) / d
            self.DWL_iXXa[:m, m] = -del_val / d
            self.DWL_iXXa[m, :m] = -del_val / d
            self.DWL_iXXa[m, m] = 1.0 / d
        else:
            self.DWL_iXXa[0] = 1.0 / b[k]

        self.DWL_Idx[k] = m
        self.DWL_nA = m + 1
        self.DWL_A = np.append(self.DWL_A, k)
        self.DWL_S = np.append(self.DWL_S, sgn).reshape(-1,1)


    def __DWL_remove(self, k):
        l = self.DWL_Idx[k]
        m = self.DWL_nA
        self.DWL_Idx[k] = -1
        if l+1 < m:
            self.DWL_Idx[self.DWL_A[(l+1):m]] -= 1

        self.DWL_nA = m-1
        self.DWL_A = np.delete(self.DWL_A, l)
        self.DWL_S = np.delete(self.DWL_S, l).reshape(-1,1)

        if m > 1:
            a = self.DWL_iXXa[:m, :m]
            b = a[:, l]
            idx = np.delete(np.arange(m),l)
            self.DWL_iXXa[:m-1, :m-1] = a[np.ix_(idx, idx)] - np.outer(b[idx], b[idx]) / b[l]

        self.DWL_iXXa[:, m-1] = 0
        self.DWL_iXXa[m-1, :] = 0


    def __DWL_getXXa(self, A):
        
        if isinstance(A, int):
            if self.DWL_XxIdx[A] == -1:
                self.DWL_XxCnt += 1
                if self.DWL_XxCnt > self.DWL_XxMax:
                    oldmax = self.DWL_XxMax
                    oldXx = self.DWL_Xx
                    self.DWL_XxMax = min(oldmax*2, self.DWL_p)
                    self.DWL_Xx = np.zeros((self.DWL_p, self.DWL_XxMax))
                    self.DWL_Xx[:,:oldmax] = oldXx
                self.DWL_XxIdx[A] = self.DWL_XxCnt - 1
                self.DWL_Xx[:, self.DWL_XxCnt-1] = np.dot(self.DWL_X.T, self.DWL_X[:, A])
        
        else:
            for k in A:
                if self.DWL_XxIdx[k] == -1:
                    self.DWL_XxCnt += 1
                    if self.DWL_XxCnt > self.DWL_XxMax:
                        oldmax = self.DWL_XxMax
                        oldXx = self.DWL_Xx
                        self.DWL_XxMax = min(oldmax*2, self.DWL_p)
                        self.DWL_Xx = np.zeros((self.DWL_p, self.DWL_XxMax))
                        self.DWL_Xx[:,:oldmax] = oldXx
                    self.DWL_XxIdx[k] = self.DWL_XxCnt - 1
                    self.DWL_Xx[:, self.DWL_XxCnt-1] = np.dot(self.DWL_X.T, self.DWL_X[:, k])

        return self.DWL_Xx[:, self.DWL_XxIdx[A]]