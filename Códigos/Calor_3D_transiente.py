
import meshio
import numpy as np
import pandas as pd
from scipy.sparse.linalg import cg, spsolve
from scipy.sparse import lil_matrix, csr_matrix, issparse
from Mesh_3D import mesh3d
from Matrizes3D import matriz3D


'''
##############################################################################
# 1) Input - Definicoes da simulacao
##############################################################################
'''
Q = 0.0           # geracao de calor
dt = 0.0001       # time step
nIter = 100000    # numero de iteracoes
teta = 1.0        # metodo dif. finitas - implicito      = 1.0;
#                                       - explicito      = 0.0;
#                                       - crank nicolson = 0.5.

'''
##############################################################################
# 2) Input Malha
##############################################################################
'''
Lx = 1
Ly = 1
Lz = 1
le = 0.1        # tamanho medio do elemento
nome_arquivo = 'minha_malha'
formato = '.msh'


##############################################################################
# 2.1) Malha gerada no API do GMSH
##############################################################################
arquivo = nome_arquivo + formato
malha = mesh3d(Lx, Ly, Lz, le, arquivo)


##############################################################################
# 2.2) Leitura das matrizes da malha
##############################################################################
msh = meshio.read(arquivo)
X = msh.points[:, 0]                # coordenada x dos nos
Y = msh.points[:, 1]                # coordenada y dos nos
Z = msh.points[:, 2]                # coordenada z dos nos
npoints = len(X)                    # numero de nos


##############################################################################
# 2.3) Matriz de conectividade (IEN) e de contorno
##############################################################################
IENbound = []                       # nos do contorno
for elem in msh.cells:
    if elem[0] == 'triangle':
        IENbound.append(elem[1])
    elif elem[0] == 'tetra':        # elementos tetraedricos
        IEN = elem[1]               # matriz de conectivide IEN
ne = len(IEN)                       # numero de elementos

# print('Para verificacao da IEN, somar 1 nas tags dos nos')
# print('IEN \n', IEN)
# print('IENbound: \n', IENbound)

bound1 = []  # lista com os nos da primeira superficie das ccs
for elem in IENbound[0]:
    for no in elem:
        bound1.append(no)

bound2 = []  # lista com os nos da segunda superficie das ccs
# for elem in IENbound[5]:
#     for no in elem:
#         bound2.append(no)

bound = bound1 + bound2

##############################################################################
# 3) Condicao de contorno
##############################################################################
bval = np.zeros((npoints), dtype='float')
for b in range(len(bval)):
    if b in bound1:
        bval[b] = 10.0
    # elif b in bound2:
    #     bval[b] = 0.0
# print('bval=',bval)


##############################################################################
# 4) Assembling (matrizes K e M)
##############################################################################
# LIL is a convenient format for constructing sparse matrices
K = lil_matrix((npoints, npoints), dtype='double')
M = lil_matrix((npoints, npoints), dtype='double')
for e in range(0, ne):
    # construir as matrizes do elemento
    v1 = IEN[e, 0]          # vertice 1
    v2 = IEN[e, 1]          # vertice 2
    v3 = IEN[e, 2]          # vertice 3
    v4 = IEN[e, 3]          # vertice 4

    # importando matrizes do modulo Matrizes3D
    m = matriz3D(vi=v1, vj=v2, vk=v3, vl=v4, X=X, Y=Y, Z=Z)
    melem = m.matrizm()
    kelem = m.matrizk()

    for ilocal in range(0, 4):
        iglobal = IEN[e, ilocal]
        for jlocal in range(0, 4):
            jglobal = IEN[e, jlocal]
            K[iglobal, jglobal] = K[iglobal, jglobal] + kelem[ilocal, jlocal]
            M[iglobal, jglobal] = M[iglobal, jglobal] + melem[ilocal, jlocal]

    print(f'Assembling - {round(100*e/ne, 1)} % ...')


##############################################################################
# 5) Montagem do sistema linear
##############################################################################
# change to csr: efficient arithmetic operations as CSR + CSR, CSR * CSR, etc.
rho = 7870  # kg/m^3
cv = 486  # J/kg.K
M = rho*cv*M
M = M.tocsr()
K = K.tocsr()

# lado esquerdo do sistema
H = M + (teta)*dt*K
H = H.tolil()

# salvar H em H2
H2 = H.copy()
H2 = H2.todense()


# criacao das listas das variaveis
f = np.zeros((npoints), dtype='double')     # lado direito da eq
T = np.zeros((npoints), dtype='double')     # Temperaturas


##############################################################################
# 5.1) Imposicao das condicoes de contorno de Dirichlet
##############################################################################
# deixar a matriz H simetrica (passa os valores para o outro lado da equacao)
for i in bound:
    H[i, :] = 0.0                           # zera a linha toda
    H[:, i] = 0.0                           # zera a coluna toda
    H[i, i] = 1.0                           # 1 na diagonal
    T[i] = bval[i]                          # Temperatura de contorno
print('H eh esparsa?', issparse(H))


##############################################################################
# 6) Iteracoes no tempo
##############################################################################
T_in = T                                    # Temperatura inicial
T_time = [T_in]                             # Lista de temperaturas por
#                                             iteracao de tempo

for n in range(0, nIter):
    print(f'{n} - {round(100*n/nIter,2)}%')
    # lado direito da equacao
    f = M.dot(T) - dt*(1-teta)*K.dot(T)   # + M*dt*Q

    # # aplicar c.c. de Dirichlet no vetor f a cada iteração
    # for i in bound:
    #     for j in range(npoints):
    #         # passa os valores para o outro lado da equacao
    #         f[j] = f[j] - H2[j, i]*bval[i]

    for i in bound:
        f[i] = bval[i]  # mantem os nos de contorno com a temperatura inicial

    # solucao do sistema linear
    T = spsolve(H.tocsc(), f) # cg(H, f)[0]
    T_time.append(T)


##############################################################################
# 7) Salva resultados para visualizacao no Paraview
##############################################################################
header = ['X', 'Y', 'Z', 'T1', 'T2', 'T3', 'Tfinal']
df = pd.DataFrame([X, Y, Z, T_time[0], T_time[100], T_time[200], T_time[-1]]).T

df.to_excel(f'Temperaturas.xlsx',
            header=header,
            float_format="%.2f",
            index=False
            )
df.to_csv('Temperaturas_csv.csv', encoding='utf-8', index=False)

print('done')
