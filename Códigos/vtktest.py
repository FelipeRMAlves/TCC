from pyevtk.hl import gridToVTK 
import numpy as np 
import random as rnd 
# Dimensions 
nx, ny, nz = 6, 6, 2 
lx, ly, lz = 1.0, 1.0, 1.0 
dx, dy, dz = lx/nx, ly/ny, lz/nz 
ncells = nx * ny * nz 
npoints = (nx + 1) * (ny + 1) * (nz + 1) 
# Coordinates 
X = np.arange(0, lx + 0.1*dx, dx, dtype='float64') 
Y = np.arange(0, ly + 0.1*dy, dy, dtype='float64') 
Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64') 

# Variables 
pressure = np.random.rand(ncells).reshape((nx, ny, nz)) 
temp = np.random.rand(npoints).reshape((nx + 1, ny + 1, nz + 1)) 
gridToVTK("./structured", X, Y, Z, cellData = {"pressure" : pressure}, pointData = {"temp" : temp})
