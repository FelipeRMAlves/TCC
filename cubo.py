
import meshio
import numpy as np
from datetime import datetime
from tqdm import tqdm
from scipy.sparse.linalg import cg, spsolve
from scipy.sparse import lil_matrix, csr_matrix, issparse
from modulos.malha import cube
from modulos.montagem import assembling3D, assembling2D
from modulos.matrizes import matriz3D, matriz2D

startTime = datetime.now()

'''
##############################################################################
# 1) Input - Parametros da simulacao
##############################################################################
'''
dt = 0.1                        # Time step
nIter = 50                      # Numero de iteracoes
theta = 1.0                     # Metodo dif. finitas - implicito      = 1.0;
#                                                     - explicito      = 0.0;
#                                                     - crank nicolson = 0.5.
condCont = 1                    # 1-dirichlet, 2-neumann ou neumann/robin

T2 = 100
k = 49.8
rho = 7850
cv = 486
Tin = 35



'''
##############################################################################
# 2) Input - gerar malha no API do GMSH
##############################################################################
'''
# Inputs de um cubo
Lx = 0.05             # Tamanho da aresta no eixo x  [m]
Ly = 0.05             # Tamanho da aresta no eixo y  [m]
Lz = 0.1              # Tamanho da aresta no eixo z  [m]
le = 0.025            # Tamanho medio do elemento    [m]
filename = 'cube.msh'
cube(Lx, Ly, Lz, le, filename=filename)

# Salvando o tempo levado para construcao da malha
meshTime = datetime.now() - startTime


##############################################################################
# 2.1) Leitura das matrizes da malha
##############################################################################
msh = meshio.read(filename)
X = msh.points[:, 0]                # Coordenada x dos nohs
Y = msh.points[:, 1]                # Coordenada y dos nohs
Z = msh.points[:, 2]                # Coordenada z dos nohs
npoints = len(X)                    # Numero de nos


##############################################################################
# 2.3) Matriz de conectividade (IEN)
##############################################################################
IENbound = []                       # Nohs do contorno
for elem in msh.cells:
    if elem[0] == 'triangle':
        IENbound.append(elem[1])
    elif elem[0] == 'tetra':        # Elementos tetraedricos
        IEN = elem[1]               # Matriz de conectivide IEN
ne = len(IEN)                       # Numero de elementos


##############################################################################
# 2.4) Nohs de contorno
##############################################################################
bound1 = []          # Lista com os nohs da primeira superficie das ccs
bound2 = []          # Lista com os nohs da segunda superficie das ccs

for elem in IENbound[0]:
    for noh in elem:
        bound1.append(noh)

for elem in IENbound[5]:
    for noh in elem:
        bound2.append(noh)

bound = bound1 + bound2


##############################################################################
# 3) Condicao de contorno
##############################################################################
print('\nPreparando simulacao:')
if condCont == 1:  # Dirichlet
    bval = np.zeros((npoints), dtype='float')
    for b in tqdm(range(len(bval))):
        if b in bound1:
            bval[b] = T2
        elif b in bound2:
            bval[b] = T2

else:
    qi = np.zeros((npoints), dtype='float')
    Tc = np.zeros((npoints), dtype='float')
    for b in tqdm(range(len(qi))):
        if b in bound1:
            qi[b] = 1.0     # Neumann
            Tc[b] = -Tinf   # Robin
        elif b in bound2:
            qi[b] = -1.0    # Neumann
            Tc[b] = Tinf    # Robin


##############################################################################
# 4) Assembling (matrizes K, M, MG e MC)
##############################################################################
print('\nMain Assembling:')
main = assembling3D(IEN,npoints,ne,X=X,Y=Y,Z=Z,k=k)
K = main[0]
M = main[1]

# print('\nBoundary Assembling:')
# boundAss = assembling2D(IENboundG,npoints,X=X,Y=Y,Z=Z)
# MC = boundAss


##############################################################################
# 5) Montagem do sistema linear
##############################################################################
# change to csr: efficient arithmetic operations as CSR + CSR, CSR * CSR, etc.
M = M.tocsr()
K = K.tocsr()

# Lado esquerdo do sistema
A = rho*cv*M + (theta)*dt*(K)  # + MC nas iteracoes

# # Salvar A em A2
# A2 = A.copy()
# A2 = A2.todense()

A = A.tolil()

# Criacao das listas das variaveis
b = np.zeros((npoints), dtype='double')       # Lado direito da eq
T = Tin*np.ones((npoints), dtype='double')    # Temperaturas iniciais


##############################################################################
# 5.1) Imposicao das condicoes de contorno
##############################################################################
# Deixar a matriz H simetrica (passa os valores para o outro lado da equacao)
if condCont == 1:
    for i in bound:
        A[i, :] = 0.0                             # zera a linha toda
        # for j in range(npoints):
        #     # passa os valores para o outro lado da equacao
        #     b[j] = b[j] - A2[j, i]*bval[i]
        # A[:, i] = 0.0                           # zera a coluna toda
        A[i, i] = 1.0                             # 1 na diagonal
        T[i] = bval[i]                            # Temperatura de contorno
    # print('A eh esparsa?', issparse(H))


# Salva resultados iniciais para visualizacao no Paraview
point_data = {'temp' : T}
meshio.write_points_cells(f'sol-0.vtk',msh.points,
                        msh.cells,point_data=point_data,)


##############################################################################
# 6) Iteracoes no tempo
##############################################################################
print('\nEquacao transiente:')

if condCont == 1:  # Dirichlet
    for n in tqdm(range(nIter)):
        # Lado direito da equacao
        b = rho*cv*M.dot(T) - dt*(1-theta)*K.dot(T) # + M*dt*Q
        for i in bound:
            b[i] = bval[i]  # mantem os nohs de contorno c/ a temp de contorno

    # Solucao do sistema linear (AT=b)
    T = spsolve(A.tocsc(), b) # cg(H, f)[0] -> apenas zerando coluna:
                                            # simetrica positiva definida.

    # Salva resultados para visualizacao no Paraview
    point_data = {'temp' : T}
    meshio.write_points_cells(f'sol-{n+1}.vtk',msh.points,
                            msh.cells,point_data=point_data,)

# else:
#     i = 0
#     for n in tqdm(range(nIter)):
#         # qm = q0*(1-(n*dt/t_stop))    # Fluxo de calor momentaneo    [W/m^2]
#         # q = qm*qi                    # Vetor de fluxo de energia
#         # h = 100

#         if n >= (i+1)*100:
#             i = i+1

#         h = hList[i]
#         qm = qList[i]/(2*Ad)

#         q = qm*qi                      # Vetor de fluxo de energia
#         MC = h*MG
#         # Lado direito da equacao
#         f = rho*cv*M.dot(T) + dt*MG.dot(q) + dt*MC.dot(Tc)  # Implicito e Q=0


#         # Solucao do sistema linear
#         T = spsolve(H.tocsc(), f)

#         # Salva resultados para visualizacao no Paraview
#         point_data = {'temp' : T}
#         meshio.write_points_cells(f'sol-{n+1}.vtk',msh.points,
#                                 msh.cells,point_data=point_data,)

print('\nSimulacao finalizada')


##############################################################################
# 7) Tempos de simulacao
##############################################################################
totalTime = datetime.now() - startTime
print(f'\n\nTempos da simulacao:')
print(f'Construcao da malha -> {meshTime}')
print(f'Simulacao -> {totalTime - meshTime}')
print(f'Tempo total -> {totalTime}')
