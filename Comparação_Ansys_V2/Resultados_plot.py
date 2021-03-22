import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


###############################################################################
# Lista de colunas de temperaturas a serem lidas
###############################################################################
cols = ['10s', '20s', '30s', '40s', '50s', '100s']
ind = 5  # indice de interesse


###############################################################################
# Transformando resultados do Ansys em listas
###############################################################################
ansys_df = pd.read_excel('Resultados_ansys.xlsx', sheet_name='py_tab')
ansys_X = ansys_df['X'].values.tolist()
ansys_Y = ansys_df['Y'].values.tolist()
ansys_Z = ansys_df['Z'].values.tolist()

ansys_temp = []
for col in cols:
    ansys_temp.append(ansys_df[col].values.tolist())


###############################################################################
# Transformando resultados do TCC em listas
###############################################################################
tcc_df = pd.read_excel('Resultados_tcc.xlsx', sheet_name='py_tab')
tcc_X = tcc_df['X'].values.tolist()
tcc_Y = tcc_df['Y'].values.tolist()
tcc_Z = tcc_df['Z'].values.tolist()

tcc_temp = []
for col in cols:
    tcc_temp.append(tcc_df[col].values.tolist())


###############################################################################
# erro
###############################################################################
erro = []
z_erro = []
for a in range(len(ansys_Z)):
    for t in range(len(tcc_Z)):
        # print(f'AnsysXYZ - {ansys_X[a]} / {ansys_Y[a]} / {ansys_Z[a]}')
        # print(f'TCC_XYZ - {tcc_X[t]} / {tcc_Y[t]} / {tcc_Z[t]} \n')
        zdif = abs(ansys_Z[a] - tcc_Z[t])
        xdif = abs(ansys_X[a] - tcc_X[t])
        ydif = abs(ansys_Y[a] - tcc_Y[t])
        if zdif < 0.01 and xdif < 0.5 and ydif < 0.5:
            dif = abs(ansys_temp[ind][a]-tcc_temp[ind][t])
            erro.append(dif)
            z_erro.append(tcc_Z[t])
# print(erro)


###############################################################################
# Plot
###############################################################################
fig, ax = plt.subplots()

color = 'black'
ax.plot(ansys_Z[::1], ansys_temp[ind][::1], 'bo', label=f'Ansys em {cols[ind]}',
         linewidth=0.05, markersize=4.0)
ax.plot(tcc_Z[::1], tcc_temp[ind][::1], 'r.', label=f'Código em {cols[ind]}',
         linewidth=0.7, markersize=3.0)
ax.set_xlabel('Coordenada Z [mm]', fontsize=13)
ax.set_ylabel('Temperatura [ºC]', fontsize=13, color=color)         
ax.tick_params(axis='y', labelcolor=color)
ax.legend()


color = 'slategrey'
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('erro absoluto', fontsize=13, color=color)
ax2.plot(z_erro[::1], erro[::1], 'mx', label='erro',
         linewidth=2.0, markersize=3.4)
ax2.tick_params(axis='y', labelcolor=color)
ax2.legend()
fig.tight_layout()  # otherwise the right y-label is slightly clipped

# ax2.annotate("Maior erro absoluto", 
#             xy=(83.5, np.amax(erro)),
#             xytext=(63.5, np.amax(erro)),
#             arrowprops=dict(arrowstyle='simple',connectionstyle="arc3"))

ttl = ax.set_title(f'Temperaturas ao longo do eixo Z no instante {cols[ind]}',
             fontsize=15)
ttl.set_weight('bold')
# ax.spines['right'].set_color((.8,.8,.8))
# ax.spines['top'].set_color((.8,.8,.8))
ax.grid()

plt.show()
