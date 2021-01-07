
class matriz3D():  
    def __init__(self, X, Y, Z, v1, v2, v3, v4):
        import numpy as np

        self.bi = (Y[v2] - Y[v4])*(Z[v3] - Z[v4]) - (Y[v3] - Y[v4])*(Z[v2] - Z[v4])
        self.bj = (Y[v3] - Y[v4])*(Z[v1] - Z[v4]) - (Y[v1] - Y[v4])*(Z[v3] - Z[v4])
        self.bk = (Y[v1] - Y[v4])*(Z[v2] - Z[v4]) - (Y[v2] - Y[v4])*(Z[v1] - Z[v4])
        self.bl = -(self.bi + self.bj + self.bk)

        self.ci = (X[v3] - X[v4])*(Z[v2] - Z[v4]) - (X[v2] - X[v4])*(Z[v3] - Z[v4])
        self.cj = (X[v1] - X[v4])*(Z[v3] - Z[v4]) - (X[v3] - X[v4])*(Z[v1] - Z[v4])
        self.ck = (X[v2] - X[v4])*(Z[v1] - Z[v4]) - (X[v1] - X[v4])*(Z[v2] - Z[v4])
        self.cl = -(self.ci + self.cj + self.ck)

        self.di = (X[v2] - X[v4])*(Y[v3] - Y[v4]) - (X[v3] - X[v4])*(Y[v2] - Y[v4])
        self.dj = (X[v3] - X[v4])*(Y[v1] - Y[v4]) - (X[v1] - X[v4])*(Y[v3] - Y[v4])
        self.dk = (X[v1] - X[v4])*(Y[v2] - Y[v4]) - (X[v2] - X[v4])*(Y[v1] - Y[v4])
        self.dl = -(self.di + self.dj + self.dk)


        self.vol = (1.0/6.0)*np.linalg.det(np.array([[1.0, X[v1], Y[v1], Z[v1]],
                                                     [1.0, X[v2], Y[v2], Z[v2]],
                                                     [1.0, X[v3], Y[v3], Z[v3]],
                                                     [1.0, X[v4], Y[v4], Z[v4]]]))
    def volCalc(self):
        return self.vol
    
    def matrizm(self):
        import numpy as np
        melem = (self.vol/20.0)*np.array([[2.0, 1.0, 1.0, 1.0],
                                          [1.0, 2.0, 1.0, 1.0],
                                          [1.0, 1.0, 2.0, 1.0],
                                          [1.0, 1.0, 1.0, 2.0]])
        return melem

    def matrizk(self):
        import numpy as np
        # Para material isotropico (kx=ky=kz=k)
        k = 1.0
        
        B = (1.0/(6.0*self.vol))*np.array([[self.bi, self.bj, self.bk, self.bl],
                                           [self.ci, self.cj, self.ck, self.cl],
                                           [self.di, self.dj, self.dk, self.dl]])
        BT = np.transpose(B)

        kelem = self.vol * k * np.dot(BT, B)
        return kelem
