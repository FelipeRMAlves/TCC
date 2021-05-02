
import numpy as np
import meshio
import gmsh


def disc(re, rt, ri, e, le_min, le_max, filename, furos):
    # Before using any functions in the Python API, Gmsh must be initialized:
    gmsh.initialize()
    gmsh.model.add("minha_malha")

    # in order to output messages on the terminal, just set the 
    # "General.Terminal" option to 1:
    gmsh.option.setNumber("General.Terminal", 1)

    # Adicionar cilindros para contruir o disco
    # addCylinder(x, y, z, dx, dy, dz, raio, tag)
    cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, re, tag=1)
    cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, e, ri, tag=2)

    if furos != 'no':
        # furos eh um lista com [raio do furinho, angulo entre furos]
        r   = furos[0]
        ang = furos[1]
        d = re-(13.6/1000)
        p = ri+(13.6/1000)
        coordIn = [d,p]
        
        coord = []
        for i in range(1,9):
            coord.append([d*np.cos(i*ang), d*np.sin(i*ang)])
        for i in range(1,9):
            coord.append([p*np.cos(i*ang), p*np.sin(i*ang)])

        t = 4
        for c in coordIn:
            gmsh.model.occ.addCylinder( c, 0, 0, 0, 0, e, r, tag=t)
            gmsh.model.occ.addCylinder(-c, 0, 0, 0, 0, e, r, tag=t+1)
            gmsh.model.occ.addCylinder(0, c, 0, 0, 0, e, r, tag=t+2)
            gmsh.model.occ.addCylinder(0,-c, 0, 0, 0, e, r, tag=t+3)
            t = t+4

        for c in coord:
            gmsh.model.occ.addCylinder( c[0], c[1],0,0,0,e,r,tag=t)
            gmsh.model.occ.addCylinder(-c[0], c[1],0,0,0,e,r,tag=t+1)
            gmsh.model.occ.addCylinder(-c[0],-c[1],0,0,0,e,r,tag=t+2)
            gmsh.model.occ.addCylinder( c[0],-c[1],0,0,0,e,r,tag=t+3)
            t = t+4

        gmsh.model.occ.cut([(3, 1)], [(3, i) for i in range(4,t)], 3)  # furinhos
        domain = gmsh.model.occ.cut([(3, 3)], [(3, 2)], 4)  # miolo do disco

    else:
        gmsh.model.occ.cut([(3, 1)], [(3, 2)], 4)  # apenas miolo do disco


    # sincronizar o modelo
    gmsh.model.occ.synchronize()
    volumes = gmsh.model.getEntities(dim=3)

    # Definir os tamanhos dos elementos da malha
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", le_min)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", le_max)

    # Gerar malha
    gmsh.model.mesh.generate(3)  # 3D

    # Salvar arquivo .msh
    gmsh.write(filename)
    # gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.finalize()


def polar(X,Y,Z):
    npoints = len(X)                    # Numero de nos

    # listas com angulo e raio de cada noh
    ang = np.zeros((npoints), dtype='float')
    raio = np.zeros((npoints), dtype='float')
    for p in range(len(X)):
        raio[p] = np.sqrt(X[p]**2 + Y[p]**2)
        a = np.arctan2(Y[p], X[p])
        if a < 0:
            a = a + 2*np.pi
        ang[p] = a

    return ang, raio


def boundf(filename):
    msh = meshio.read(filename)
    IENbound = []                       # IEN de nohs do contorno
    for elem in msh.cells:
        if elem[0] == 'triangle':       # Elementos triangulares
            IENbound.append(elem[1])
        elif elem[0] == 'tetra':        # Elementos tetraedricos
            IEN = elem[1]               # Matriz de conectivide IEN

    return IENbound, IEN


def contornoDisco(IENbound, raio, Z, rt, e):
    bound1 = []          # Lista com os nohs da 1a superficie das ccs
    bound2 = []          # Lista com os nohs da 2a superficie das ccs
    IENboundG1 = []      # IEN de fluxo de calor (contorno de neumann)
    IENboundG2 = []      # IEN de fluxo de calor (contorno de neumann)
    IENconv1= []         # IEN de conveccao (contorno de robin)
    IENconv2= []         # IEN de conveccao (contorno de robin)

    for elem in IENbound[1]:
        aux1 = []
        aux2 = []
        for noh in elem:
            if Z[noh] == e:
                aux2.append(noh)
                # apenas nohs no contato com pastilha
                if raio[noh] > rt:
                    bound1.append(noh)
                    aux1.append(noh)

        if len(aux1) == 3:
            IENboundG1.append(aux1)
        if len(aux2) == 3:
            IENconv1.append(aux2)

    for elem in IENbound[2]:
        aux1 = []
        aux2 = []
        for noh in elem:
            if Z[noh] == 0:
                aux2.append(noh)
                # apenas nohs no contato com pastilha
                if raio[noh] > rt:
                    bound2.append(noh)
                    aux1.append(noh)
        if len(aux1) == 3:
            IENboundG2.append(aux1)
        if len(aux2) == 3:
            IENconv2.append(aux2)

    return bound1, bound2, IENboundG1, IENboundG2, IENconv1, IENconv2


def cube(Lx, Ly, Lz, le, filename):

    # Before using any functions in the Python API, Gmsh must be initialized:
    gmsh.initialize()

    # in order to output messages on the terminal, just set the 
    # "General.Terminal" option to 1:
    gmsh.option.setNumber("General.Terminal", 1)

    # iniciando a malha
    gmsh.model.add("minha_malha")

    # pontos:
    # gmsh.model.geo.addPoint(x, y, z, target mesh size, tag)
    gmsh.model.geo.addPoint(0, 0, 0, le, 1)
    gmsh.model.geo.addPoint(Lx, 0, 0, le, 2)
    gmsh.model.geo.addPoint(Lx, Ly, 0, le, 3)
    gmsh.model.geo.addPoint(0, Ly, 0, le, 4)

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
    h = Lz  # geometry height in the z-direction
    e = gmsh.model.geo.extrude([(2, 1)], 0, 0, h)

    # Physical groups are collections of model entities and are identified 
    # by their dimension and by a tag.
    # Whole domain:
    domain_tag = e[1][1]
    domain_physical_tag = 1001  # the volume
    gmsh.model.addPhysicalGroup(dim=3, tags=[domain_tag], tag=domain_physical_tag)
    gmsh.model.setPhysicalName(dim=3, tag=domain_physical_tag, name="Whole domain")

    # O volume tem planos de contorno. Os planos tem retas. As retas tem pontos
    planos = gmsh.model.getBoundary(dimTags=[(3, domain_tag)],
                                            oriented=False)
    retas = []
    pontos = []
    geometria = {}
    n = -1
    for plano in planos:
        n = n + 1
        gmsh.model.addPhysicalGroup(dim=2, tags=[plano[1]], tag=n)
        gmsh.model.setPhysicalName(dim=2,tag=n,name=f'cc{n}')
        retas_lim = gmsh.model.getBoundary(dimTags=[plano], oriented=False)
        retas.append(retas_lim)
        geometria[plano[1]] = {}
        for r in retas:
            noh = gmsh.model.getBoundary(dimTags=[r], oriented=False)
            pontos.append(noh)
            nohs = []
            for p in noh:
                nohs.append(p[1])
            geometria[plano[1]][r[1]] = nohs

    # sincronizar o modelo
    gmsh.model.geo.synchronize()

    # Gerar malha
    gmsh.model.mesh.generate(3)  # 3D

    # Salvar arquivo .msh
    gmsh.write(filename)
    # gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.finalize()
