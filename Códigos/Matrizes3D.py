
import numpy as np

class matriz3D():  

    def __init__(self, X, Y, Z, vi, vj, vk, vl):

        self.bi = (Y[vj] - Y[vl])*(Z[vk] - Z[vl]) - (Y[vk] - Y[vl])*(Z[vj] - Z[vl])
        self.bj = (Y[vk] - Y[vl])*(Z[vi] - Z[vl]) - (Y[vi] - Y[vl])*(Z[vk] - Z[vl])
        self.bk = (Y[vi] - Y[vl])*(Z[vj] - Z[vl]) - (Y[vj] - Y[vl])*(Z[vi] - Z[vl])
        self.bl = -(self.bi + self.bj + self.bk)

        self.ci = (X[vk] - X[vl])*(Z[vj] - Z[vl]) - (X[vj] - X[vl])*(Z[vk] - Z[vl])
        self.cj = (X[vi] - X[vl])*(Z[vk] - Z[vl]) - (X[vk] - X[vl])*(Z[vi] - Z[vl])
        self.ck = (X[vj] - X[vl])*(Z[vi] - Z[vl]) - (X[vi] - X[vl])*(Z[vj] - Z[vl])
        self.cl = -(self.ci + self.cj + self.ck)

        self.di = (X[vj] - X[vl])*(Y[vk] - Y[vl]) - (X[vk] - X[vl])*(Y[vj] - Y[vl])
        self.dj = (X[vk] - X[vl])*(Y[vi] - Y[vl]) - (X[vi] - X[vl])*(Y[vk] - Y[vl])
        self.dk = (X[vi] - X[vl])*(Y[vj] - Y[vl]) - (X[vj] - X[vl])*(Y[vi] - Y[vl])
        self.dl = -(self.di + self.dj + self.dk)


        self.vol = (1.0/6.0)*np.linalg.det(np.array([[1.0, X[vi], Y[vi], Z[vi]],
                                                     [1.0, X[vj], Y[vj], Z[vj]],
                                                     [1.0, X[vk], Y[vk], Z[vk]],
                                                     [1.0, X[vl], Y[vl], Z[vl]]]))
    def volCalc(self):
        return self.vol
    
    def matrizm(self):
        melem = (self.vol/20.0)*np.array([[2.0, 1.0, 1.0, 1.0],
                                          [1.0, 2.0, 1.0, 1.0],
                                          [1.0, 1.0, 2.0, 1.0],
                                          [1.0, 1.0, 1.0, 2.0]])
        return melem

    def matrizk(self):
        # Para material isotropico (kx=ky=kz=k)
        k = 1.0  # W/m.K - Aco SAE 1020
        
        B = (1.0/(6.0*self.vol))*np.array([[self.bi, self.bj, self.bk, self.bl],
                                           [self.ci, self.cj, self.ck, self.cl],
                                           [self.di, self.dj, self.dk, self.dl]])
        BT = np.transpose(B)

        kelem = self.vol * k * np.dot(BT, B)
        return kelem
