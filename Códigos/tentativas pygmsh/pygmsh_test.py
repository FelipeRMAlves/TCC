import pygmsh
import meshio

with pygmsh.geo.Geometry() as geom:
    poly = geom.add_polygon(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ],
        mesh_size=0.5,
    )
    geom.extrude(poly, [0.0, 0.0, 0.1], num_layers=2)
    # geom.add_boundary_layer
    mesh = geom.generate_mesh()

    points = mesh.points
    cells = mesh.cells
    connectivity = cells[2][1]

    filename = 'test.msh'
    meshio.write(filename, mesh, file_format="gmsh22")
    mesh.write(f"{filename}.vtk")