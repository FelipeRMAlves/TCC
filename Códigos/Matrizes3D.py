
class matriz3D():  
    def __init__(self, X, Y, v1, v2, v3): # Z e v4
        import numpy as np

        self.bi = (Y[v2] - Y[v4])*(Z[v3] - Z[v4]) - (Y[v3] - Y[v4])*(Z[v2] - Z[v4])
        self.bj = (Y[v3] - Y[v4])*(Z[v1] - Z[v4]) - (Y[v1] - Y[v4])*(Z[v3] - Z[v4])
        self.bk = (Y[v1] - Y[v4])*(Z[v2] - Z[v4]) - (Y[v2] - Y[v4])*(Z[v1] - Z[v4])
        self.bl = -(self.bi + self.bj + self.bk)

        self.ci = (X[v3] - X[v4])*(Z[v2] - Z[v4]) - (X[v2] - X[v4])*(Z[v3] - Z[v4])
        self.cj = (X[v1] - X[v4])*(Z[v3] - Z[v4]) - (X[v3] - X[v4])*(Z[v1] - Z[v4])
        self.ck = (X[v2] - X[v4])*(Z[v1] - Z[v4]) - (X[v1] - X[v4])*(Z[v2] - Z[v4])
        self.cl = -(self.Ci + self.cj + self.ck)

        self.di = (Y[v2] - Y[v4])*(Z[v3] - Z[v4]) - (Y[v3] - Y[v4])*(Z[v2] - Z[v4])
        self.dj = (Y[v3] - Y[v4])*(Z[v1] - Z[v4]) - (Y[v1] - Y[v4])*(Z[v3] - Z[v4])
        self.dk = (Y[v1] - Y[v4])*(Z[v2] - Z[v4]) - (Y[v2] - Y[v4])*(Z[v1] - Z[v4])
        self.dl = -(self.di + self.dj + self.dk)


        self.vol = (1.0/6.0)*np.linalg.det(np.array([[1.0, X[v1], Y[v1]],
                                                     [1.0, X[v2], Y[v2]],
                                                     [1.0, X[v3], Y[v3]],
                                                     [1.0, X[v4], Y[v4]]]))
    def volCalc(self):
        return self.vol
    
    def matrizm(self):
        import numpy as np
        melem = (self.area/12.0)*np.array([[2.0, 1.0, 1.0],
                                           [1.0, 2.0, 1.0],
                                           [1.0, 1.0, 2.0]])
        return melem

    def matrizk(self):
        import numpy as np
        B = (1.0/(2.0*self.area))*np.array([[self.bi, self.bj, self.bk],
                                            [self.ci, self.cj, self.ck]])
        BT = np.transpose(B)

        kelem = self.area * np.dot(BT, B)
        return kelem



