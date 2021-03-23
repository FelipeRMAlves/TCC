import meshio
import numpy as np
import pandas as pd
from Disc_mesh import disc


'''
##############################################################################
# 2) Input Malha
##############################################################################
'''
Lx = 0.1
Ly = 0.1
Lz = 0.1
le = 0.05        # tamanho medio do elemento
nome_arquivo = 'disc'
formato = '.msh'


##############################################################################
# 2.1) Malha gerada no API do GMSH
##############################################################################
arquivo = nome_arquivo + formato
malha = disc(Lx, Ly, Lz, le, arquivo)