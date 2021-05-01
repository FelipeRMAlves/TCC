
from tqdm import tqdm
from modulos.matrizes import matriz3D, matriz2D
from scipy.sparse import lil_matrix, csr_matrix, issparse


def assembling3D(IEN,npoints,ne,X,Y,Z,k):
    # LIL is a convenient format for constructing sparse matrices
    K = lil_matrix((npoints, npoints), dtype='double')
    M = lil_matrix((npoints, npoints), dtype='double')

    for e in tqdm(range(ne)):
        # construir as matrizes do elemento
        v1 = IEN[e, 0]          # vertice 1
        v2 = IEN[e, 1]          # vertice 2
        v3 = IEN[e, 2]          # vertice 3
        v4 = IEN[e, 3]          # vertice 4

        # importando matrizes do modulo matrizes
        m = matriz3D(v1=v1, v2=v2, v3=v3, v4=v4, X=X, Y=Y, Z=Z)
        melem = m.matrizm()
        kelem = m.matrizk(k)

        for ilocal in range(0, 4):
            iglobal = IEN[e, ilocal]
            for jlocal in range(0, 4):
                jglobal = IEN[e, jlocal]
                K[iglobal, jglobal] = K[iglobal, jglobal] + kelem[ilocal, jlocal]
                M[iglobal, jglobal] = M[iglobal, jglobal] + melem[ilocal, jlocal]

    return K, M


def assembling2D(IENboundG,npoints,X,Y,Z,c):
    # Matriz de massa do contorno (de Neumann e Robin)
    MC = lil_matrix((npoints, npoints), dtype='double')
    for e in tqdm(IENboundG):
        # construir as matrizes do elemento
        v1 = e[0]
        v2 = e[1]
        v3 = e[2]
        v = [v1,v2,v3]

        # importando matrizes do modulo matrizesMEF
        # obs: o codigo pode ser usado pois os elementos estao sempre na mesma
        # coordenada Z. Caso nao estivesse, teria que fazer diferente.
        sq = matriz2D(v1=v1, v2=v2, v3=v3, X=X, Y=Y)
        melemg = sq.matrizm()

        for ilocal in range(0, 3):
            iglobal = v[ilocal]
            for jlocal in range(0, 3):
                jglobal = v[jlocal]
                MC[iglobal, jglobal] = MC[iglobal, jglobal] + c*melemg[ilocal, jlocal]
                
    return MC
