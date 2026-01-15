import numpy as np
import matplotlib.pyplot as plt

from utils import linear_int, binary_comps

#K = 8.617*10**(-5)
#m0 = 9.11*10**(-31)
#e = 2.7
K=8.617*10**-5
m0=9.11
e=1.6

X = [0.25, 0.5, 0.75]
T_task4 = np.linspace(0, 300, 1000)
Tinv_task5 = 1/np.linspace(50, 300, 1000)

#GaAsP

def Eg_f(Eg0, alpha, beta, temp):
    return Eg0 - alpha*temp**2/(temp+beta)

def sigma_f(t, x):
    mee = linear_int("m*e", "Ga", "As", "P", x)

    gamma_eqs = []
    gamma_eqs.append(binary_comps["GaAs"]["gamma1"] - 2*binary_comps["GaAs"]["gamma2"])
    gamma_eqs.append(1/2*(2*binary_comps["GaAs"]["gamma1"] - binary_comps["GaAs"]["gamma2"] - 3*binary_comps["GaAs"]["gamma3"]))
    gamma_eqs.append(binary_comps["GaAs"]["gamma1"] - 2*binary_comps["GaAs"]["gamma3"])
    meh1 = m0/np.mean(gamma_eqs)

    gamma_eqs = []
    gamma_eqs.append(binary_comps["GaP"]["gamma1"] - 2*binary_comps["GaP"]["gamma2"])
    gamma_eqs.append(1/2*(2*binary_comps["GaP"]["gamma1"] - 2*binary_comps["GaP"]["gamma2"] - 3*binary_comps["GaP"]["gamma3"]))
    gamma_eqs.append(binary_comps["GaP"]["gamma1"] - 2*binary_comps["GaP"]["gamma3"])
    meh2 = m0/np.mean(gamma_eqs)

    meh=x*meh1 + (1-x)*meh2

    Eg0 = linear_int("Eg", "Ga", "As", "P", x, bow=0.19)
    alpha = linear_int("alpha", "Ga", "As", "P", x)
    beta = linear_int("beta", "Ga", "As", "P", x)

    #Eg = Eg_f(Eg0, alpha, beta, t)
    Eg = Eg0

    Nc = (mee*K*t)**(3/2)
    Nv = (meh*K*t)**(3/2)
    #sigma = e*t**(-3/2)*np.sqrt(Nc*Nv)*np.exp(-Eg/(2*K*t))
    sigma = np.exp(-Eg/(2*K*t))
    return sigma, Eg

eps_xy = (binary_comps["GaAs"]["a"]-binary_comps["GaP"]["a"])/binary_comps["GaAs"]["a"]
eps_z = -2*binary_comps["GaAs"]["c12"]/binary_comps["GaAs"]["c11"]*eps_xy

def strain(x):
    eps_xy = (binary_comps["GaAs"]["a"]-binary_comps["GaP"]["a"])/binary_comps["GaAs"]["a"]
    eps_z_gaas = -2*binary_comps["GaAs"]["c12"]/binary_comps["GaAs"]["c11"]*eps_xy
    eps_z_gap = -2*binary_comps["GaP"]["c12"]/binary_comps["GaP"]["c11"]*eps_xy

    eps_z = x*eps_z_gaas + (1-x)*eps_z_gap

    ac = x*binary_comps["GaAs"]["ac"] + (1-x)*binary_comps["GaP"]["ac"]
    av = x*binary_comps["GaAs"]["av"] + (1-x)*binary_comps["GaP"]["av"]

    delta_Ec = ac*(2*eps_xy + eps_z)
    delta_Ev = av*(2*eps_xy + eps_z)
    return eps_xy, eps_z, delta_Ec, delta_Ev

strain_effects = []
delta_list = []
for x in X:
    eps_xy, eps_z, delta_Ec, delta_Ev = strain(x)
    strain_effects.append(2*eps_xy+eps_z)
    c12 = linear_int("c12", "Ga", "As", "P", x)
    c11 = linear_int("c11", "Ga", "As", "P", x)
    delta_3 = c12*2*eps_xy+ c11*eps_z
    delta_list.append(delta_3)

plt.plot(X, strain_effects)
plt.xlabel("x") 
plt.ylabel("Strain effect")
plt.savefig("strain.png", dpi=300)
plt.close()

fig,ax = plt.subplots(1,3, figsize=(10,6))

for x in X:
    data_sigma = []
    data_eg = []
    for t in T_task4:
        sigma, eg = sigma_f(t,x)
        data_sigma.append(sigma)
        data_eg.append(eg)
    log_sig_inv = [np.log(1/s) for s in data_sigma]
    s_inv = [np.log(s) for s in data_sigma]
    print(data_sigma)
    ax[0].plot(Tinv_task5, s_inv)
    ax[1].plot(T_task4, data_eg)
    ax[2].plot(Tinv_task5, log_sig_inv)


ax[0].set_title("sigma")
#ax[0].set_yscale('log')
ax[1].set_title("Eg")
ax[2].set_yscale("log")
legend = [f"x={x:.2f}" for x in X]
fig.legend(legend)
plt.savefig("sigma-eg.png", dpi=300)
plt.close()




