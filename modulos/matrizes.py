
import numpy as np
from numpy.linalg import det
from numpy import array


class matriz3D():
    def __init__(self, X, Y, Z, v1, v2, v3, v4):

        self.b1 = (Y[v2] - Y[v4])*(Z[v3] - Z[v4]) - \
            (Y[v3] - Y[v4])*(Z[v2] - Z[v4])
        self.b2 = (Y[v3] - Y[v4])*(Z[v1] - Z[v4]) - \
            (Y[v1] - Y[v4])*(Z[v3] - Z[v4])
        self.b3 = (Y[v1] - Y[v4])*(Z[v2] - Z[v4]) - \
            (Y[v2] - Y[v4])*(Z[v1] - Z[v4])
        self.b4 = -(self.b1 + self.b2 + self.b3)

        self.c1 = (X[v3] - X[v4])*(Z[v2] - Z[v4]) - \
            (X[v2] - X[v4])*(Z[v3] - Z[v4])
        self.c2 = (X[v1] - X[v4])*(Z[v3] - Z[v4]) - \
            (X[v3] - X[v4])*(Z[v1] - Z[v4])
        self.c3 = (X[v2] - X[v4])*(Z[v1] - Z[v4]) - \
            (X[v1] - X[v4])*(Z[v2] - Z[v4])
        self.c4 = -(self.c1 + self.c2 + self.c3)

        self.d1 = (X[v2] - X[v4])*(Y[v3] - Y[v4]) - \
            (X[v3] - X[v4])*(Y[v2] - Y[v4])
        self.d2 = (X[v3] - X[v4])*(Y[v1] - Y[v4]) - \
            (X[v1] - X[v4])*(Y[v3] - Y[v4])
        self.d3 = (X[v1] - X[v4])*(Y[v2] - Y[v4]) - \
            (X[v2] - X[v4])*(Y[v1] - Y[v4])
        self.d4 = -(self.d1 + self.d2 + self.d3)

        self.vol = (1.0/6.0)*det(array([[1.0, X[v1], Y[v1], Z[v1]],
                                        [1.0, X[v2], Y[v2], Z[v2]],
                                        [1.0, X[v3], Y[v3], Z[v3]],
                                        [1.0, X[v4], Y[v4], Z[v4]]])
                                 )

    def volCalc(self):
        return self.vol

    def matrizm(self):
        melem = (self.vol/20.0)*array([[2.0, 1.0, 1.0, 1.0],
                                       [1.0, 2.0, 1.0, 1.0],
                                       [1.0, 1.0, 2.0, 1.0],
                                       [1.0, 1.0, 1.0, 2.0]])
        return melem

    def matrizk(self, k):
        # condutividade termica k
        # Para material isotropico (kx=ky=kz=k)

        B = (1.0/(6.0*self.vol))*array([[self.b1, self.b2, self.b3, self.b4],
                                        [self.c1, self.c2, self.c3, self.c4],
                                        [self.d1, self.d2, self.d3, self.d4]])
        BT = np.transpose(B)

        kelem = self.vol * k * np.dot(BT, B)
        return kelem


class matriz2D():

    def __init__(self, X, Y, v1, v2, v3):

        self.b1 = Y[v2] - Y[v3]
        self.b2 = Y[v3] - Y[v1]
        self.b3 = Y[v1] - Y[v2]
        self.c1 = X[v3] - X[v2]
        self.c2 = X[v1] - X[v3]
        self.c3 = X[v2] - X[v1]

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
        B = (1.0/(2.0*self.area))*np.array([[self.b1, self.b2, self.b3],
                                            [self.c1, self.c2, self.c3]])
        BT = np.transpose(B)

        kelem = self.area * np.dot(BT, B)
        return kelem
