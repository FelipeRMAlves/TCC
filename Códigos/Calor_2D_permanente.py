# Solucao do problema termico 2D permanente usando MEF
#
#   Lap(T) = 0
#
# Com c.c.s de Dirichlet


from Mesh_2D import mesh2d
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg
from scipy.sparse import lil_matrix
# from scipy.sparse import csr_matrix

# entrada da malha
nx = 10
ny = 10
Lx = 1
Ly = 1
lados = 4  # num de lados do elemento da malha

malha = mesh2d(Lx=Lx, Ly=Ly, nx=nx, ny=ny, lados=lados)

# lil_matrix
K = lil_matrix((malha.nx * malha.ny, malha.nx * malha.ny), dtype='double')
M = lil_matrix((malha.nx * malha.ny, malha.nx * malha.ny), dtype='double')

for e in range(0, malha.ne):
    # construir as matrizes do elemento
    v1 = malha.matriz_IEN()[e, 0]
    v2 = malha.matriz_IEN()[e, 1]
    v3 = malha.matriz_IEN()[e, 2]

    from Matrizes2D import matriz2D
    sq = matriz2D(v1=v1, v2=v2, v3=v3, X=malha.X, Y=malha.Y)
    area = sq.areaCalc()
    kelem = sq.matrizk()

    for ilocal in range(0, 3):
        iglobal = malha.matriz_IEN()[e, ilocal]
        for jlocal in range(0, 3):
            jglobal = malha.matriz_IEN()[e, jlocal]

            K[iglobal, jglobal] = K[iglobal, jglobal] + kelem[ilocal, jlocal]

    print(f'{round(100*e/malha.ne, 0)} % - calculando  matrizes...')

# vetores de identificacao da malha e condicao de contorno
bound = malha.bound
inner = malha.inner
bval = np.zeros((malha.npoints), dtype='double')
for b in range(nx):
    bval[b] = malha.X[b]
for b in range(nx, (nx*(ny-1) + 1), nx):
    bval[b] = malha.Y[b]
for b in range((2*nx)-1, malha.npoints, nx):
    bval[b] = ((malha.Y[b]) ** 2 + 1)
for b in range(nx * (ny-1), malha.npoints):
    bval[b] = ((malha.X[b]) ** 2 + 1)

# imposicao das condicoes de contorno de Dirichlet
# deixando a matriz k simetrica (passa os valores para o outro lado)
b = np.zeros((malha.npoints), dtype='double')
for i in bound:
    K[i, :] = 0.0  # zera a linha toda
    b[i] = bval[i]
    for j in range(malha.npoints):
        # passa os valores para o outro lado da equacao
        b[j] = b[j] - K[j, i]*bval[i]
    K[:, i] = 0.0
    K[i, i] = 1.0


# check if H is symmetric
def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)
print(check_symmetric(K, tol=1e-8))


# solucao do sistema linear
T = cg(K, b)[0]
print('T=',T)
plotado = malha.plotSol(var=T)
plt.show()
