

class disc():

    def __init__(self, re, rt, ri, e, le_min, le_max, filename):
        # le = tamanho medio do elemento
        import numpy as np
        import gmsh


        # Before using any functions in the Python API, Gmsh must be initialized:
        gmsh.initialize()
        gmsh.model.add("minha_malha")



        includeBoundary = True


        

        # in order to output messages on the terminal, just set the 
        # "General.Terminal" option to 1:
        gmsh.option.setNumber("General.Terminal", 1)

        # Adicionar cilindros para contruir o disco
        # addCylinder(x,y,z, dx, dy, dz)

        # construcao da superficie do trilho da pastilha
        cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, re, tag = 1)  # exterior trilho
        cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, ri, tag = 2)  # interior trilho
        gmsh.model.occ.cut([(3, 1)], [(3, 2)], 3)  # disco apenas com o trilho
        # gmsh.model.setPhysicalName(3, 1, "The volume")


        # # construcao do interior do disco alem do trilho
        # cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, rt, tag = 4)
        # cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, ri, tag = 5)
        # gmsh.model.occ.cut([(3, 4)], [(3, 5)], 6)


        # sincronizar o modelo
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(dim=3)


        # Definir os tamanhos dos elementos da malha
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", le_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", le_max)


        # Gerar malha
        gmsh.model.mesh.generate(3)  # 3D


        # Salvar arquivo .msh
        gmsh.write("disc.msh")
        # gmsh.option.setNumber("Mesh.SaveAll", 1)
        gmsh.finalize()
