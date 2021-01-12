
import numpy as np 

class matriz2D():
    
    def __init__(self, X, Y, v1, v2, v3):

        self.bi = Y[v2] - Y[v3]
        self.bj = Y[v3] - Y[v1]
        self.bk = Y[v1] - Y[v2]
        self.ci = X[v3] - X[v2]
        self.cj = X[v1] - X[v3]
        self.ck = X[v2] - X[v1]

        self.area = (1.0/2.0)*np.linalg.det(np.array([[1.0, X[v1], Y[v1]],
                                                      [1.0, X[v2], Y[v2]],
                                                      [1.0, X[v3], Y[v3]]]))
    def areaCalc(self):
        return self.area
    
    def matrizm(self):
        melem = (self.area/12.0)*np.array([[2.0, 1.0, 1.0],
                                           [1.0, 2.0, 1.0],
                                           [1.0, 1.0, 2.0]])
        return melem

    def matrizk(self):
        B = (1.0/(2.0*self.area))*np.array([[self.bi, self.bj, self.bk],
                                            [self.ci, self.cj, self.ck]])
        BT = np.transpose(B)

        kelem = self.area * np.dot(BT, B)
        return kelem



