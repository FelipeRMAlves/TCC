
import meshio
import numpy as np
import pandas as pd
from vtk import *
import glob


'''
##############################################################################
# 1) Arquivos a serem lidos
##############################################################################
'''
filenames = ['sol-0.vtk', 'sol-99.vtk', 'sol-199.vtk', 'sol-299.vtk',
             'sol-399.vtk', 'sol-499.vtk', 'sol-749.vtk', 'sol-999.vtk',
             'sol-1499.vtk', 'sol-1999.vtk', 'sol-2499.vtk', 'sol-2999.vtk']


##############################################################################
# 6) Salvando dados em uma lista
##############################################################################
T_time = []
for arq in filenames:
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(arq)
    reader.Update()
    vtkdata = reader.GetOutput()
    pointData = vtkdata.GetPointData()
    vectorData = pointData.GetArray(3)
    T_time.append(pointData)

print(T_time)


##############################################################################
# 7) Salva resultados para visualizacao no Paraview
##############################################################################
header = ['X', 'Y', 'Z', 'T_in', 'T10s', 'T20s', 'T30s', 'T40s', 'T50s',
          'T75s', 'T100s', 'T150s', 'T200s', 'T250s', 'T300s']
df = pd.DataFrame([X, Y, Z, T_time[0], T_time[99], T_time[199], T_time[299], 
                    T_time[399], T_time[499], T_time[749], 
                    T_time[999], T_time[1499], T_time[1999], T_time[2499], 
                    T_time[-1]]).T

df.to_excel(f'Temperaturas.xlsx',
            header=header,
            float_format="%.2f",
            index=False
            )
df.to_csv('Temperaturas_csv.csv', encoding='utf-8', index=False)

print('done')


