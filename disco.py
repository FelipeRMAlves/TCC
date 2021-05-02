
import time
import meshio
import numpy as np
from tqdm import tqdm
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix, issparse
from modulos.malha import disc, polar, boundf, contornoDisco
from modulos.montagem import assembling3D, assembling2D
from modulos.input import inputInfo

startTime = time.time()
open("Results.txt", "w").close()

'''
##############################################################################
# 1) Leitura do input
##############################################################################
'''
fileIpt = 'Input.xlsx'          # Nome do arquivo de input
ipt     = inputInfo(fileIpt)    # Chamada da funcao
param   = ipt.getParam()        # Parametros da simulacao
geom    = ipt.getGeom()         # Parametros da geometria do disco
prop    = ipt.getProperties()   # Propriedades termicas do disco
ambt    = ipt.getAmb()          # Temperatura inicial e do ar
qcoefs  = ipt.getCoefs('q')     # Fluxo de calor
hcoefs  = ipt.getCoefs('h')     # Coeficientes de conveccao

##############################################################################
# 1.1) Parametros da simulacao
##############################################################################
dt    = param[0]                # Time step
nIter = int(param[1])           # Numero de iteracoes
theta = param[2]                # Metodo dif. finitas - implicito      = 1.0;
#                                                     - explicito      = 0.0;
#                                                     - crank nicolson = 0.5.

##############################################################################
# 1.2) Parametros fisicos
##############################################################################
Tin  = ambt[0]                  # Temp inicial do solido             [ºC]
Tinf = ambt[1]                  # Temp do fluido na conveccao        [ºC]

##############################################################################
# 1.3) Dados do Disco de freio
##############################################################################
e   = geom[0]/1000              # Espessura do disco                 [m]
re  = geom[1]/1000              # Raio externo do disco              [m]
ri  = geom[2]/1000              # Raio interno do disco              [m]
rt  = geom[3]/1000              # Raio interno da pastilha           [m]
r   = geom[4]/1000              # Raio dos furinhos                  [m]
k   = prop[0]                   # Condutividade termica disco        [W/m.K]
cv  = prop[1]                   # Calor especifico do disco          [J/kg.K]
rho = prop[2]                   # Massa especifica do disco          [kg/m^3]
Ad  = np.pi*((re**2)-(rt**2))   # Area de troca de calor do disco    [m^2]
Afur = 32*np.pi*(r**2)          # Area dos furinhos
Ad = Ad - Afur                  # Area disco descontando os furinhos


'''
##############################################################################
# 2) Malha
##############################################################################
'''
# Gerar malha no API do GMSH
le_min = param[3]/1000          # Tamanho minimo dos elementos [m]
le_max = param[3]/1000          # Tamanho maximo dos elementos [m]
filename = 'disc.msh'
disc(re, rt, ri, e, le_min, le_max, filename, furos=[r,(11.25/180)*np.pi])

# Salvando o tempo levado para construcao da malha
meshTime = time.time() - startTime


##############################################################################
# 2.1) Leitura das matrizes da malha
##############################################################################
msh = meshio.read(filename)
X = msh.points[:, 0]                # Coordenada x dos nohs
Y = msh.points[:, 1]                # Coordenada y dos nohs
Z = msh.points[:, 2]                # Coordenada z dos nohs
npoints = len(X)                    # Numero de nos

pol  = polar(X,Y,Z)
ang  = pol[0]
raio = pol[1]


##############################################################################
# 2.2) Matriz de conectividade (IEN)
##############################################################################
boundParam = boundf(filename)
IENbound = boundParam[0]
IEN = boundParam[1]
ne = len(IEN)


##############################################################################
# 2.3) Nohs de contorno
##############################################################################
boundDisco = contornoDisco(IENbound,raio,Z,rt,e)
bound1 = boundDisco[0]
bound2 = boundDisco[1]
IENboundG1 = boundDisco[2]
IENboundG2 = boundDisco[3]
IENboundG = IENboundG1 + IENboundG2
IENconv1 = boundDisco[4]
IENconv2 = boundDisco[5]


'''
##############################################################################
# 3) Vetores para condicao de contorno
##############################################################################
'''
print('\nPreparando simulacao:')
Tc = np.zeros((npoints), dtype='float')
qi = np.zeros((npoints), dtype='float')
for b in tqdm(bound1):
    qi[b] = 1.0     # Neumann - vetor para fluxo de calor
for b in bound2:
    qi[b] = -1.0    # Neumann - vetor para fluxo de calor
for b in IENconv1:
    Tc[b] = -Tinf   # Robin
for b in IENconv2:
    Tc[b] = Tinf    # Robin


'''
##############################################################################
# 4) Assembling (matrizes K, M, MC e Mh)
##############################################################################
'''
print('\nMain Assembling:')
main = assembling3D(IEN,npoints,ne,X=X,Y=Y,Z=Z,k=k)
K = main[0]
M = main[1]

print('\nBoundary Assembling:')
MC = assembling2D(IENboundG,npoints,X=X,Y=Y,Z=Z,c=1)

# Matrizes de conveccao precisam ter sinal
Mh1 = assembling2D(IENconv1,npoints,X=X,Y=Y,Z=Z,c= 1)
Mh2 = assembling2D(IENconv2,npoints,X=X,Y=Y,Z=Z,c=-1)
Mhin = Mh1 + Mh2


'''
##############################################################################
# 5) Sistema linear
##############################################################################
'''
# change to csr: efficient arithmetic operations as CSR + CSR, CSR * CSR, etc.
M = M.tocsr()
K = K.tocsr()
MC = MC.tocsr()
Mhin = Mhin.tocsr()

# Lado esquerdo do sistema
A = rho*cv*M + (theta)*dt*(K)  # + (theta)*dt*Mh nas iteracoes
A = A.tolil()
A2 = A.copy()

# Criacao das listas das variaveis
b = np.zeros((npoints), dtype='double')       # Lado direito da eq
T = Tin*np.ones((npoints), dtype='double')    # Temperaturas iniciais

# Salva resultados iniciais para visualizacao no Paraview
point_data = {'temp' : T}
meshio.write_points_cells(f'sol-0.vtk',msh.points,
                        msh.cells,point_data=point_data,)


##############################################################################
# 5.1) Iteracoes no tempo
##############################################################################
print('\nEquacao transiente:')

i = 0
maxT = []
rel = nIter/(len(qcoefs))
for n in tqdm(range(nIter)):
    if n >= (i+1)*rel:
        i = i+1
    h = hcoefs[i][0]
    qm = qcoefs[i][0]/(2*Ad)
    q = qm*qi                      # Vetor de fluxo de energia

    Mh = h*MC
    Mh2 = h*Mhin

    # Lado esquerdo da equacao (incluindo a conveccao)
    A = A2 + (theta)*dt*Mh2

    # Lado direito da equacao
    f = rho*cv*M - (1-theta)*dt*(K+Mh2)
    b = f.dot(T) + dt*Mh.dot(Tc) + dt*MC.dot(q)  # Q=0

    # Solucao do sistema linear (AT=b)
    T = spsolve(A.tocsc(), b)

    maxT.append(max(T))
    # Salva resultados para visualizacao no Paraview
    point_data = {'temp' : T}
    meshio.write_points_cells(f'sol-{n+1}.vtk',msh.points,
                            msh.cells,point_data=point_data,)

print('\nSimulacao finalizada')


'''
##############################################################################
# 6) Tempos de simulacao
##############################################################################
'''
totalTime = time.time()/60 - startTime/60
print(f'\n\nTempos da simulacao:')
print(f'Tempo total -> {round(totalTime, 2)} min')
print(f'Construcao da malha -> {round(meshTime, 2)} seg')
print(f'Simulacao -> {round(totalTime - meshTime/60, 2)} min')

with open("Results.txt","a") as f:
    f.write(f'Simulacao realizada por: {ipt.header[0]} em {ipt.header[2]}\n\n')
    f.write(f'Temperatura maxima = {max(maxT)} \n')
    f.write(f'Temperatura max final = {max(T)} \n\n')
    f.write(f'Tempo total -> {round(totalTime, 2)} min \n')
    f.write(f'Construcao da malha -> {round(meshTime, 2)} seg \n')
    f.write(f'Simulacao -> {round(totalTime - meshTime/60, 2)} min \n\n')
    f.write('Para melhor visualizacao dos resultados, utilize o Paraview')