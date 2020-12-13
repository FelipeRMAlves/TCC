
import gmsh

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# in order to output messages on the terminal, just set the 
# "General.Terminal" option to 1:
gmsh.option.setNumber("General.Terminal", 1)

# iniciando a malha
gmsh.model.add("minha_malha")

# pontos:
# gmsh.model.geo.addPoint(x, y, z, target mesh size, tag)
lc = 0.5
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 0, 0, lc, 2)
gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
gmsh.model.geo.addPoint(0, 1, 0, lc, 4)

# linhas:
# gmsh.model.geo.addLine(ponto inicial, ponto final, tag)
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# superficies:
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)


# physical group
ps = gmsh.model.addPhysicalGroup(2, [1], 6)
gmsh.model.setPhysicalName(2, ps, "My surface")


# extrudar:
# e = gmsh.model.geo.extrude([(dim, tag)], 0, 0, h)
h = 0.1  # geometry height in the z-direction
e = gmsh.model.geo.extrude([(2, 1)], 0, 0, h)


# Physical groups are collections of model entities and are identified 
# by their dimension and by a tag.
# Whole domain:
domain_tag = e[1][1]
domain_physical_tag = 1001  # the volume
gmsh.model.addPhysicalGroup(dim=3, tags=[domain_tag], tag=domain_physical_tag)
gmsh.model.setPhysicalName(dim=3, tag=domain_physical_tag, name="Whole domain")

# Find domain boundary tags
boundary_dimtags = gmsh.model.getBoundary(dimTags=[(3, domain_tag)],
                                          oriented=False)
boundary_tags = []
n = -1
for tag in boundary_dimtags:
    n = n + 1
    boundary_tags.append(tag[1])
    gmsh.model.addPhysicalGroup(dim=2, tags=[tag[1]], tag=n)
    gmsh.model.setPhysicalName(dim=2,tag=n,name=f'cc{n}')

# entities = gmsh.model.getEntities()
# print(entities)


# sincronizar o modelo
gmsh.model.geo.synchronize()

# Gerar malha
gmsh.model.mesh.generate(3)  # 3D

# Salvar arquivo .msh
gmsh.write("minha_malha.msh")
# gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.finalize()