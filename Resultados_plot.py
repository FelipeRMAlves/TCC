import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


###############################################################################
# Lista de arquivos a serem lidos
###############################################################################
filenames = ['Resultados_ansys.xlsx', 'Resultados_tcc.xlsx']


###############################################################################
# Lista de colunas de temperaturas a serem lidas
###############################################################################
cols = ['10s', '20s', '30s', '40s', '50s', '75s', '100s', '150s', '200s',
        '250s', '300s']


###############################################################################
# Transformando resultados do Ansys em listas
###############################################################################
ansys_df = pd.read_excel('Resultados_ansys.xlsx', sheet_name='py_tab')
ansys_Z = ansys_df['Z'].values.tolist()

ansys_temp = []
for col in cols:
    ansys_temp.append(ansys_df[col].values.tolist())


###############################################################################
# Transformando resultados do TCC em listas
###############################################################################
tcc_df = pd.read_excel('Resultados_tcc.xlsx', sheet_name='py_tab')
tcc_Z = tcc_df['Z'].values.tolist()

tcc_temp = []
for col in cols:
    tcc_temp.append(tcc_df[col].values.tolist())


###############################################################################
# Plot
###############################################################################
plt.plot(ansys_Z[::1], ansys_temp[5][::1], 'bx', label='Ansys', linewidth=0.5, 
         markersize=2.8)
plt.plot(tcc_Z[::1], tcc_temp[5][::1], 'ro', label='codigo', linewidth=0.5,
         markersize=2.8)

ax = plt.axes()
ax.set_xlabel('Coordenada Z [mm]', fontsize=14)
ax.set_ylabel('Temperatura [ºC]', fontsize=13)

plt.legend()
plt.show()






# ###############################################################################
# # Plot dados |Z| vs f:
# ###############################################################################
# plt.plot(frequencias[::1],impedancias[::1],'b*',linewidth=0.5,markersize=2.8)
# ax = plt.axes()
# ax.set_xlabel('frequencia [kHz]', fontsize=14)
# ax.set_ylabel('Amplitude da Impedancia [Ω]', fontsize=13)
# ax.annotate("Maior valor", 
#             xy=(10, np.amax(impedancias)),
#             xytext=(10.5, 2.3),
#             arrowprops=dict(arrowstyle='simple',connectionstyle="arc3"))
# plt.show()


