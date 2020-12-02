# -*- coding: utf-8 -*-

class mesh2d():

    def __init__(self, Lx, Ly, nx, ny, lados):
        import numpy as np

        # Variáveis
        if lados == 4:
            f = 1
        if lados == 3:
            f = 2

        self.ne = f * (nx - 1) * (ny - 1)  # Número de elementos
        self.npoints = nx * ny  # Número de pontos
        self.dx = Lx / (nx - 1)  # Comprimento de cada elemento em x
        self.dy = Ly / (ny - 1)  # Comprimento de cada elemento em y

        # Criação das Listas
        # Matriz de coordenada x de cada ponto da malha
        self.X = np.zeros((self.npoints), dtype='float')
        # Matriz de coordenada y de cada ponto da malha
        self.Y = np.zeros((self.npoints), dtype='float')

        # 4. Preenchimento das listas X e Y
        for i in range(0, self.npoints):
            self.X[i] = Lx * i / (nx - 1)
            if i > (nx - 1):
                self.X[i] = self.X[i - nx]

        for i in range(nx, self.npoints):
            #  Malha com perturbação:
            A = 0.0  # "A = 0.0" é a malha sem perturbacao
            fi = 2 * np.pi / 4.0
            self.Y[i] = (self.Y[i - nx] + self.dy) + A * \
                np.sin((2 * np.pi / 5.0) * self.X[i] - fi)

        # Inner e Bound
        P = np.ones((self.npoints), dtype='int')
        for i in range(self.npoints):
            P[i] = i

        self.bound = list(P)
        inner = np.zeros((nx-2)*(ny-2), dtype='int')
        s = 0
        for j in range(1, ny-1):
            for i in range(self.npoints):
                if j * nx < i < (((j+1) * nx) - 1):
                    self.bound.remove(i)
                    inner[s] = i
                    s = s + 1

        self.inner = list(inner)
        self.tipo = lados
        self.nx = nx
        self.ny = ny

    def matriz_IEN(self):
        import numpy as np
        self.IEN = np.zeros((self.ne, self.tipo), dtype='int')
        # 5. Matriz IEN
        # 5.1. Malha de quadriláteros
        if self.tipo == 4:
            s = -1
            for e in range(0, self.ne):
                if e % (self.nx - 1) == 0:
                    s = s + 1
                self.IEN[e] = [e + s, e + 1 + s, e +
                               self.nx + 1 + s, e + self.nx + s]

        # 5.2. Malha de triângulos
        if self.tipo == 3:
            i = 1
            for a in range(self.ny - 1):
                for b in range(self.nx - 1):
                    self.IEN[i] = [self.nx * a + b, self.nx +
                                   1 + b + self.nx * a, self.nx * a + b + 1]
                    self.IEN[i - 1] = [self.nx * a + b,
                                       self.nx + self.nx * a + b, 
                                       self.nx + self.nx * a + b + 1]
                    i += 2
        return self.IEN

    def plotMalha(self):
        import numpy as np
        import matplotlib.pyplot as plt
        plt.plot(self.X, self.Y, 'ro')
        if self.tipo == 3:
            plt.triplot(self.X, self.Y, self.IEN)

        if self.tipo == 4:  # OBS: Tambem funciona para tipo == 3
            for i in range(0, len(self.IEN)):
                cx = np.zeros((self.tipo + 1), dtype='float')
                cy = np.zeros((self.tipo + 1), dtype='float')
                for l in range(0, self.tipo):
                    cx[l] = self.X[self.IEN[i][l]]
                    cy[l] = self.Y[self.IEN[i][l]]
                cx[self.tipo] = cx[0]
                cy[self.tipo] = cy[0]
                plt.plot(cx, cy, 'b-')

        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()
        # plt.savefig("2dmesh")

    def plotSol(self, var):
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        if self.tipo == 3:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            plot = ax.plot_trisurf(self.X, self.Y, var, cmap='jet')
            ax.set_title("Temperatura em uma placa plana utilizando MEF (ºC)")

        elif self.tipo == 4:
            levels = 20
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            plot = ax.tricontourf(self.X, self.Y, var, levels, cmap='jet')
            fig.colorbar(plot)
            ax.set_title("Temperatura em uma placa plana utilizando MEF (ºC)")

        return plot
