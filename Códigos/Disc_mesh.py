

# class disc():

#     def __init__(self, Lx, Ly, Lz, le, filename):
#         # le = tamanho medio do elemento
#         import numpy as np
#         import gmsh


import gmsh
import numpy as np



'''
##############################################################################
# 1) Input - Definicoes da simulacao
##############################################################################
'''
re = 0.015        # raio externo [m]
ri = 0.005        # raio interno [m]
e = 0.005          # espessura [m]


##############################################################################
# Codigo
##############################################################################
gmsh.initialize()
gmsh.model.add("DFG 3D")

# addCylinder(x,y,z)
cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, e, 0, re, tag = 1)
cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, e, 0, ri, tag = 2)
gmsh.model.occ.cut([(3, 3)], [(3, 7)], 8)

gmsh.model.occ.synchronize()
volumes = gmsh.model.getEntities(dim=3)





# sincronizar o modelo
gmsh.model.geo.synchronize()

# Gerar malha
gmsh.model.mesh.generate(3)  # 3D

# Salvar arquivo .msh
gmsh.write("mesh3D.msh")
# gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.finalize()
