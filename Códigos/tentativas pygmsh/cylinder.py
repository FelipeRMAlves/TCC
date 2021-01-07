import meshio
import pygmsh

def createCylinder(centre, radius, height, lcar, verbose=False, gmshFile=None, vtkFile=None):
    """
    Creates an unstructured mesh of tetrahedra inside a cylinder.
    The height of the cylinder is along the z axis.

    Parameters
    ----------
        centre: 1D array
            The two coordinates of the centre of the base disk (yx).

        radius: float
            The radius of the base disk.

        height: float
            The height of the cylinder of the z direction.

        lcar: float
            length of the mesh.
            The average distance between two nodes.

        gmshFile: string, optional
            If not None, save the clear text gmsh file with name ``gmshFile`` and suffix ``.msh``
            Default = None

        vtkFile: string, optional
            If not None, save the clear text vtk file with name ``vtkFile`` and suffix ``.vtk``
            Defaut = None

        verbose: bool, optionsl
            Verbose mode.
            Default = False

    Returns
    -------
        points: 2D numpy array
            The coordinates of the mesh nodes (zyx)
            Each line is [zPos, yPos, xPos]

        connectivity: 2D numpy array
            The connectivity matrix of the tetrahedra elements
            Each line is [node1, node2, node3, node4]

    Example
    -------
        >>> points, connectivity =  spam.mesh.createCylinder((0.0,0.0), 0.5, 2.0, 0.5)
        create a mesh in a cylinder of centre 0,0,0 radius, 0.5 and height 2.0 with a characteristic length of 0.5
    """

    import pygmsh
    import meshio

    # unpack
    cy, cx = centre
    r = radius

    # raw code
    code = []
    code.append("Point(1) = {{ {x}, {y},  0, {lcar} }};".format(x=cx, y=cy, lcar=lcar))
    code.append("Point(2) = {{ {x}, {y},  0, {lcar} }};".format(x=cx + r, y=cy, lcar=lcar))
    code.append("Point(3) = {{ {x}, {y},  0, {lcar} }};".format(x=cx, y=cy + r, lcar=lcar))
    code.append("Point(4) = {{ {x}, {y},  0, {lcar} }};".format(x=cx - r, y=cy, lcar=lcar))
    code.append("Point(5) = {{ {x}, {y},  0, {lcar} }};".format(x=cx, y=cy - r, lcar=lcar))
    code.append("Circle(1) = { 2, 1, 3 };")
    code.append("Circle(2) = { 3, 1, 4 };")
    code.append("Circle(3) = { 4, 1, 5 };")
    code.append("Circle(4) = { 5, 1, 2 };")
    code.append("Line Loop(5) = { 1, 2, 3, 4 };")
    code.append("Plane Surface(6) = { 5 };")
    code.append("Extrude {{ 0, 0, {h} }} {{ Surface{{ 6 }}; }}".format(h=height))

    # geom = pygmsh.opencascade.Geometry(characteristic_length_min=lcar, characteristic_length_max=lcar,)

    # # add raw code to geometry
    # geom = pygmsh.built_in.Geometry()
    # geom.add_raw_code(code)

    # mesh
    # points, cells, _, _, _ = pygmsh.generate_mesh(geom, verbose=verbose, extra_gmsh_arguments=["-optimize_netgen"])
    # points, cells, _, _, _ = pygmsh.generate_mesh(geom, verbose=verbose, extra_gmsh_arguments=["-optimize","-algo","del3d","-clmin",str(lcar),"-clmax",str(lcar)])
    mesh = pygmsh.generate_mesh(geom, verbose=verbose, extra_gmsh_arguments=["-optimize", "-algo", "del3d", "-clmin", str(lcar), "-clmax", str(lcar)])
    points = mesh.points
    cells = mesh.cells
    connectivity = cells['tetra']

    # write gmsh/vtk file
    if gmshFile is not None:
        meshio.write_points_cells("{}.msh".format(gmshFile), points, cells, file_format='gmsh2-ascii')
    if vtkFile is not None:
        meshio.write_points_cells("{}.vtk".format(vtkFile), points, {'tetra': connectivity}, file_format='vtk-ascii')

    # NOTE: pygmsh returns coordinates in xyz. This means that:
    # 1. the coordinates array must be switched back to zyx
    # 2. the connectivity array must change so that the nodes of each tetrahedron are numbered
    #    in such a way that the Jacobian is positive. Which means a perturbation of the node numbering.

    points = points[:, ::-1]
    tmp = connectivity.copy()
    connectivity[:, 1] = tmp[:, 3]
    connectivity[:, 3] = tmp[:, 1]

    # return coordinates and connectivity matrix
    return points, connectivity

createCylinder([0,0], 1.0, 0.1, 0.1, verbose=False, gmshFile='yes', vtkFile=None)