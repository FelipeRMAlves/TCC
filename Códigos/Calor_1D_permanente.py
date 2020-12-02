import numpy as np
import matplotlib.pyplot as plt

#  1) parametros da simulacao
Ti = 0.0
Tf = 1.0
alpha = 1.0
Q = 5.0
grau = 3  #  define o tipo de elemento 
# (linear=1, quadratico=2, cubico=3)

#  2) malha
npoints = 20
nelem = int(npoints-1)
L = 1.0
npoints = npoints + (grau-1) * nelem  # atualiza numero de pontos

#  3) funcao que retorna as matrizes necessarias 
def matrizes(grau_elem):
    IEN = np.zeros((nelem, grau + 1), dtype='int')
    if grau_elem == 1:
        for e in range(0, nelem):
            IEN[e] = [e, e+1]
    elif grau_elem == 2:
        for e in range(0, nelem):
            IEN[e] = [2*e, 2*e+1, 2*e+2]
    elif grau_elem == 3:
        for e in range(0, nelem):
            IEN[e] = [3*e, 3*e+1, 3*e+2, 3*e+3]

    X = np.linspace(0, L, npoints)
    # para elementos não uniformes, eh necessario
    # sempre calcular tudo pra cada elemento
    for e in range(0, nelem):  
        v1 = IEN[0, 0]
        v2 = IEN[0, grau]
        length = X[v2] - X[v1]  
        # Modulo com as matrizes de MEF
        from Matrizes import matriz
        matrix = matriz(grau, nelem, L, length)
        kelem = matrix.matrizk()
        melem = matrix.matrizm()
        gelem = matrix.matrizg()
        felem = matrix.matrizf(fonte=Q)
    return kelem, felem, IEN, X

#  4) chamada das matrizes
matrizes = matrizes(grau_elem=grau)  
kelem = matrizes[0]
felem = matrizes[1]
IEN = matrizes[2]
X = matrizes[3]

#  5) assembling
K = np.zeros((npoints, npoints), dtype='float')
F = np.zeros((npoints, 1), dtype='float')
for e in range(0, nelem):
    for i in range(0, grau + 1):
        iglobal = IEN[e, i]
        F[iglobal] += felem[i]
        for j in range(0, grau + 1):
            jglobal = IEN[e, j]
            K[iglobal, jglobal] += kelem[i, j]

#  6) imposicao das c.cs:
#  6.1) Dirichlet:
# noh 0
K[0, :] = 0.0
K[0, 0] = 1.0
F[0] = Ti
K[-1, :] = 0.0
K[-1, -1] = 1.0
F[-1] = Tf

# #  6.2) Neumann homogeneo:
# # noh 0
# K[0, :] = 0.0
# K[0, 0] = 1.0
# F[0] = Ti

#  6.3) Neumann nao homogeneo: 
# # noh 0
# K[0, :] = 0.0
# K[0, 0] = 1.0
# F[0] = Ti
# # ultimo noh
# F[-1] = F[-1] - (4*L-3)  # (para elemento linear)


#  7) solucao do sistema linear KT=F
T = np.linalg.solve(K,F)

#  8) plot
plt.plot(X, T, 'r.-')
plt.title(f'Temperatura com elemento de grau {grau}', fontsize=16)
plt.show()
