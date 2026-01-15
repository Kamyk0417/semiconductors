import numpy as np
import matplotlib.pyplot as plt
from utils import binary_comps, linear_int

###############################################################################################
# Task 2: Trends in properties of III-V semiconductors
###############################################################################################

Z_data = []
data = {
    "a_data": [],
    "Eg_data": [],
    "dSO_data": []
}

for comp in binary_comps:
    Z_data.append(binary_comps[comp]["Z"])
    data["a_data"].append(binary_comps[comp]["a"])
    data["Eg_data"].append(binary_comps[comp]["Eg"])
    data["dSO_data"].append(binary_comps[comp]["dSO"])

for d in data:
    plt.scatter(Z_data, data[d])
    plt.legend(list(data.keys()))

plt.xlabel("Atomic Number Z")
plt.grid()
plt.savefig("results/task2_trends.png", dpi=300)
plt.close()

#################################################################################################
# Task 3: Linear interpolation of alloy properties
#################################################################################################

X_DATA = np.linspace(0,1,100)
fig, ax = plt.subplots(1,3, figsize=(15,6))

Eg_lin = []
a_lin = []
dSO_lin = []

for x in X_DATA:
    Eg_lin.append(linear_int("Eg", "Ga", "As", "P", x, bow=0.19))
    a_lin.append(linear_int("a", "Ga", "As", "P", x))
    dSO_lin.append(linear_int("dSO", "Ga", "As", "P", x))


ax[0].plot(X_DATA, Eg_lin)
ax[0].set_xlabel("x")
ax[0].set_ylabel("Eg")
ax[0].grid()

ax[1].plot(X_DATA, a_lin)
ax[1].set_xlabel("x")
ax[1].set_ylabel("a")
ax[1].grid()

ax[2].plot(X_DATA, dSO_lin)
ax[2].set_xlabel("x")
ax[2].set_ylabel("dSO")
ax[2].grid()

fig.suptitle("Linear interpolation of GaAsP alloy properties", fontsize=16)
plt.savefig("results/task3_interpolation.png", dpi=300)
plt.close()

