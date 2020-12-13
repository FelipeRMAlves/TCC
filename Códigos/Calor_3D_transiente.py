

from Mesh_2D import mesh2d
from pyevtk.hl import gridToVTK
import meshio
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg
from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from mpl_toolkits.mplot3d import Axes3D


# 1) Definicoes da simulacao
alpha = 1
time = 0.0
dt = 0.1
nIter = 50


# 2) Importacao da malha (GMSH)
filename = 'minha_malha'
formato = '.msh'
msh = meshio.read(filename + formato)
X = msh.points[:, 0]
Y = msh.points[:, 1]
Z = msh.points[:, 2]
npoints = len(X)

IEN = msh.cells[-1][1]   # msh.cells['tetra']
ne = len(IEN)
IENbound = []
for elem in msh.cells:
    if elem[0] == 'line':
        IENbound.append(elem[1])
print(IENbound)
IENboundTypeElem = list(msh.cell_data['gmsh:physical'])
boundNames = list(msh.field_data.keys())

IENboundElem = []
for elem in IENboundTypeElem:
    IENboundElem.append(boundNames[elem-1])

# cria lista de nos do contorno
cc = np.unique(IENbound.reshape(IENbound.size))
ccName = [[] for i in range(len(X))]
for elem in range(0, len(IENbound)):
    ccName[IENbound[elem][0]] = IENboundElem[elem]
    ccName[IENbound[elem][1]] = IENboundElem[elem]

# Consertando bug para malha do aerofolio
if filename == 'airfoil_naca':
    for i in range(6, 60):
        ccName[i] = 'cylinder'

#  2.1) Plot malha e pontos de contorno
plt.triplot(X, Y, IEN, 'k-')
plt.plot(X[cc], Y[cc], 'bo')
plt.show()



# Paraview
temp = np.random.rand(npoints).reshape((nx + 1, ny + 1, nx + 1)) 
gridToVTK("./structured", X, X, Y, pointData = {"temp" : temp})



# 2) Condicao de contorno
bval = np.zeros((npoints), dtype='double')
for b in range(nx):
    bval[b] = X[b]
for b in range(nx, (nx*(ny-1)+1), nx):
    bval[b] = Y[b]
for b in range((2*nx)-1, npoints, nx):
    bval[b] = ((Y[b]) ** 2 + 1)
for b in range(nx * (ny-1), npoints):
    bval[b] = ((X[b]) ** 2 + 1)
# print('bval=',bval)

# 3) Assembling
# LIL is a convenient format for constructing sparse matrices
K = lil_matrix((npoints, npoints), dtype='double')
M = lil_matrix((npoints, npoints), dtype='double')
for e in range(0, ne):
    # construir as matrizes do elemento
    v1 = IEN[e, 0]
    v2 = IEN[e, 1]
    v3 = IEN[e, 2]
    v4 = IEN[e, 2] # fix this

    from Matrizes2D import matriz2D
    sq = matriz2D(v1=v1, v2=v2, v3=v3, v4=v4, X=malha.X, Y=malha.Y)
    area = sq.areaCalc()
    melem = sq.matrizm()
    kelem = sq.matrizk()

    for ilocal in range(0, 3):
        iglobal = IEN[e, ilocal]
        for jlocal in range(0, 3):
            jglobal = IEN[e, jlocal]

            K[iglobal, jglobal] = K[iglobal, jglobal] + kelem[ilocal, jlocal]
            M[iglobal, jglobal] = M[iglobal, jglobal] + melem[ilocal, jlocal]

    print(f'{round(100*e/ne, 0)} % - calculando as matrizes...')

# change to csr: efficient arithmetic operations CSR + CSR, CSR * CSR, etc.
M = M.tocsr()
K = K.tocsr()

# lado direito do sistema linear eh fixo
H = M + dt*alpha*K

print(type(H))
print(type(M))
print(type(K))

# imposicao das condicoes de contorno de Dirichlet
f = np.zeros((npoints), dtype='double')
T = np.zeros((npoints), dtype='double')

# save H into H2
H2 = H.copy()
H2 = H2.todense()

# imposicao das condicoes de contorno de Dirichlet
# deixando a matriz k simetrica (passa os valores para o outro lado)
for i in bound:
    H[i, :] = 0.0  # zera a linha toda
    H[:, i] = 0.0
    H[i, i] = 1.0
    T[i] = bval[i]

# check if H is symmetric
# def check_symmetric(a, tol=1e-8):
#     return np.all(np.abs(a-a.T) < tol)
# print(check_symmetric(H, tol=1e-8))

# visualizacao da condicao inicial
plotado = malha.plotSol(var=T)
plt.show()

for n in range(0, nIter):
    f = M * T

    # aplicar c.c. de Dirichlet no vetor f a cada iteração
    for i in bound:
        for j in range(npoints):
            # passa os valores para o outro lado da equacao
            f[j] = f[j] - H2[j, i]*bval[i]

        # f[i] = bval[i]  ## PQ NAO FUNCIONA ???

    for i in bound:
        f[i] = bval[i]

    # solucao do sistema linear
    T = cg(H,f)[0]
    # print('T=', T)

    if lados == 3:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        plot = ax.plot_trisurf(X, Y, T,cmap='jet')
        ax.set_title(
            "Temperatura em uma placa plana utilizando MEF (ºC)")
        plt.pause(1)

    elif lados == 4:
        levels = 20
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        plot = ax.tricontourf(X, Y, T, levels, cmap='jet')
        fig.colorbar(plot)
        ax.set_title(
            "Temperatura em uma placa plana utilizando MEF (ºC)") 
        plt.pause(1)

    # plt.savefig(f'iteracao {n} triang.png')
    
