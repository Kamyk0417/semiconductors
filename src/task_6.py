import numpy as np
import matplotlib.pyplot as plt

from utils import linear_int, binary_comps, strain

X_DATA = np.linspace(0,1,100)

vbo_data = np.zeros_like(X_DATA)
Eg_data = np.zeros_like(X_DATA)

for i in range(len(X_DATA)):
    vbo_data[i] = linear_int("VBO", "Ga", "As", "P", X_DATA[i], params=binary_comps)
    Eg_data[i] = linear_int("Eg", "Ga", "As", "P", X_DATA[i], bow=0.19, params=binary_comps)


eps_xy = (binary_comps["GaAs"]["a"]-binary_comps["GaP"]["a"])/binary_comps["GaAs"]["a"]
eps_z = -2*binary_comps["GaAs"]["c12"]/binary_comps["GaAs"]["c11"]*eps_xy

strain_effects_gaas = []
delta_Ec_list_gaas = []
delta_Ev_list_gaas = []

strain_effects_gap = []
delta_Ec_list_gap = []
delta_Ev_list_gap = []

for x in X_DATA:
    eps_xy_gaas, eps_z_gaas, delta_Ec_gaas, delta_Ev_gaas, eps_xy_gap, eps_z_gap, delta_Ec_gap, delta_Ev_gap = strain(x)
    strain_effects_gaas.append(100*2*eps_xy_gaas+eps_z_gaas)
    c12 = linear_int("c12", "Ga", "As", "P", x)
    c11 = linear_int("c11", "Ga", "As", "P", x)
    delta_Ec_list_gaas.append(delta_Ec_gaas)
    delta_Ev_list_gaas.append(delta_Ev_gaas)

    strain_effects_gap.append(100*2*eps_xy_gap+eps_z_gap)
    delta_Ec_list_gap.append(delta_Ec_gap)  
    delta_Ev_list_gap.append(delta_Ev_gap)

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
fig.subplots_adjust(hspace=0)
strain_zero_line = np.zeros_like(X_DATA)

axs[0].plot(X_DATA, vbo_data, label="Valence band")
axs[0].plot(X_DATA, vbo_data+Eg_data, label="Conduction band")
axs[0].set_ylabel("Energy [eV]")
axs[0].legend()
axs[0].grid()

axs[1].plot(X_DATA, strain_zero_line, 'k--', color='gray', label="Zero strain")
axs[1].plot(X_DATA, strain_effects_gaas, label="GaAs substrate")
axs[1].plot(X_DATA, strain_effects_gap, label="GaP substrate")
axs[1].set_ylabel("Strain effect")
axs[1].legend()
axs[1].grid()


Eg_new_data_gaas = Eg_data + np.array(delta_Ec_list_gaas) - np.array(delta_Ev_list_gaas)
Eg_new_data_gap = Eg_data + np.array(delta_Ec_list_gap) - np.array(delta_Ev_list_gap)

axs[2].plot(X_DATA, Eg_data, 'k--', color='gray', label="Unstrained Eg")
axs[2].plot(X_DATA, Eg_new_data_gaas, label="GaAs substrate")
axs[2].plot(X_DATA, Eg_new_data_gap, label="GaP substrate")
axs[2].set_ylabel("Eg [eV]")
axs[2].legend()
axs[2].grid()
axs[2].set_xlabel("x (GaAs$_{x}$P$_{1-x}$)")
plt.savefig("./results/task_6_strain.png", dpi=300)
plt.close()
