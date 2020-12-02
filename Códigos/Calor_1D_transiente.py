

import numpy as np
import matplotlib.pyplot as plt

#  1) parametros da simulacao
Ti = 0.0
Tf = 1.0
alpha = 1.0
Q = 5.0
dt = 0.1
time = 0.0
grau = 3  # tipo de elemento(linear=1, quadratico=2, cubico=3)
teta = 1  # metodo (implicito=1, explicito=0 ou crank nicolson=0.5)

#  2) malha
npoints = 20
nelem = int(npoints-1)
L = 1.0
npoints = npoints + (grau-1) * nelem  # atualizacao do numero de pontos


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
    for e in range(0, nelem):  
    # para elementos não uniformes, 
    # eh necessario sempre calcular tudo pra cada elemento
        v1 = IEN[0, 0]
        v2 = IEN[0, grau]
        length = X[v2] - X[v1]  # tamanho do elemento

        # Modulo com as matrizes de MEF
        from Matrizes import matriz
        matrix = matriz(grauElem=grau, nelem=nelem, L=L, length=length)
        kelem = matrix.matrizk()
        melem = matrix.matrizm()
        felem = matrix.matrizf(fonte=Q)
    return kelem, felem, melem, IEN, X

#  4) chamada das matrizes
matrizes = matrizes(grau_elem=grau)  
kelem = matrizes[0]
felem = matrizes[1]
melem = matrizes[2]
IEN = matrizes[3]
X = matrizes[4]

#  5) assembling
K = np.zeros((npoints, npoints), dtype='float')
F = np.zeros((npoints, 1), dtype='float')
M = np.zeros((npoints, npoints), dtype='float')
for e in range(0, nelem):
    for i in range(0, grau + 1):
        iglobal = IEN[e, i]
        F[iglobal] += felem[i]
        for j in range(0, grau + 1):
            jglobal = IEN[e, j]
            K[iglobal, jglobal] += kelem[i, j]
            M[iglobal, jglobal] += melem[i, j]

#  6) contrucao do sistema linear e aplicacao das c.cs:
T = np.zeros((npoints, 1), dtype='float')
T[0] = Ti
T[-1] = Tf

# matriz do lado esquerdo
H = M + (teta) * dt*K   

# noh 0
H[0,:] = 0.0
H[0,0] = 1.0
# ultimo noh
H[-1,:] = 0.0
H[-1,-1] = 1.0

# inversa da matriz do lado esquerdo
invH = np.linalg.inv(H)
Tan = Q*L/(2*alpha) * (X - ((X**2)/L)) + (Tf/L)*X + Ti
for n in range(0,10):
    b = np.dot(M,T) - (1 - teta) * dt*np.dot(K,T) + dt*F 
    # b = np.dot(M, T) + dt*F

    # noh 0
    b[0] = Ti
    # ultimo noh
    b[-1] = Tf

    T = np.dot(invH,b)
    time = time + dt
    aux = round(n * dt, 3)
    plt.title(f'Temperatura em {aux} segundos', fontsize=16)
    plt.plot(X, Tan, 'b-')
    plt.plot(X, T, 'r.-')
    plt.pause(dt)
    plt.savefig(f'step {n} cubico')
    plt.clf()
    # print(n * dt)
 

