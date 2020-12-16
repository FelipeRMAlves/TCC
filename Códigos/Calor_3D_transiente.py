

from Mesh_3D import mesh3d
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
# 2.1) Definicao da malha pelo usuario
Lx = 1
Ly = 1
Lz = 0.5
le = 0.1    # tamanho medio do elemento
filename = 'minha_malha.msh'
malha = mesh3d(Lx, Ly, Lz, le, filename)


msh = meshio.read(filename)
X = msh.points[:, 0]
Y = msh.points[:, 1]
Z = msh.points[:, 2]
npoints = len(X)         # numero de nos

IENbound = []            # nos do contorno
for elem in msh.cells:
    if elem[0] == 'triangle':
        IENbound.append(elem[1])
    elif elem[0] == 'tetra':
        IEN = elem[1]
print('IEN: \n',IEN)
ne = len(IEN)            # numero de elementos tetraedricos
    
print(IENbound)
bound1 = []  # lista com os nos da primeira superficie das ccs
for elem in IENbound[0]:
    for no in elem:
        bound1.append(no)

IENboundTypeElem = list(msh.cell_data['gmsh:physical'])
boundNames = list(msh.field_data.keys())

IENboundElem = []
for elem in IENboundTypeElem:
    try:
        IENboundElem.append(boundNames[elem[0]])
    except IndexError:
        pass


# LINKAR AS TAGS DOS NOS COM OS NOMES DAS CCS A QUAL ELE PERTENCE
# cria lista de nos do contorno
# cc = np.unique(IENbound.reshape(IENbound.size))
# ccName = [[] for i in range(len(X))]
# for elem in range(0, len(IENbound)):
#     ccName[IENbound[elem][0]] = IENboundElem[elem]
#     ccName[IENbound[elem][1]] = IENboundElem[elem]

# # Paraview
# temp = np.random.rand(npoints).reshape((nx + 1, ny + 1, nx + 1)) 
# gridToVTK("./structured", X, X, Y, pointData = {"temp" : temp})


# 2) Condicao de contorno
bval = np.zeros((npoints), dtype='int')
for b in range(len(bval)):
    if b in bound1:
        bval[b] = 10
print('bval=',bval)

# ok ok ok ok ok ok ok ok ok ok next step





# 3) Assembling
# LIL is a convenient format for constructing sparse matrices
K = lil_matrix((npoints, npoints), dtype='double')
M = lil_matrix((npoints, npoints), dtype='double')
for e in range(0, ne):
    # construir as matrizes do elemento
    v1 = IEN[e, 0]     # vertice 1
    v2 = IEN[e, 1]     # vertice 2
    v3 = IEN[e, 2]     # vertice 3
    v4 = IEN[e, 3]     # vertice 4

    from Matrizes3D import matriz3D
    m = matriz3D(v1=v1, v2=v2, v3=v3, v4=v4, X=X, Y=Y, Z=Z)
    volume = m.volCalc()
    melem = m.matrizm()
    kelem = m.matrizk()

    for ilocal in range(0, 4):
        iglobal = IEN[e, ilocal]
        for jlocal in range(0, 4):
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

f = np.zeros((npoints), dtype='double')
T = np.zeros((npoints), dtype='double')

# save H into H2
H2 = H.copy()
H2 = H2.todense()

# imposicao das condicoes de contorno de Dirichlet
# deixando a matriz k simetrica (passa os valores para o outro lado)
for i in bound1:
    H[i, :] = 0.0  # zera a linha toda
    H[:, i] = 0.0
    H[i, i] = 1.0
    T[i] = bval[i]

# check if H is symmetric
# def check_symmetric(a, tol=1e-8):
#     return np.all(np.abs(a-a.T) < tol)
# print(check_symmetric(H, tol=1e-8))

for n in range(0, nIter):
    f = M * T

    # aplicar c.c. de Dirichlet no vetor f a cada iteração
    for i in bound1:
        for j in range(npoints):
            # passa os valores para o outro lado da equacao
            f[j] = f[j] - H2[j, i]*bval[i]

        # f[i] = bval[i]  ## PQ NAO FUNCIONA ???

    for i in bound1:
        f[i] = bval[i]

    # solucao do sistema linear
    T = cg(H,f)[0]
    # print('T=', T)

    # plt.savefig(f'iteracao {n} triang.png')
    
