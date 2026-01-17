import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from utils import linear_int, binary_comps

K=8.617*10**-5
m0=9.11*10**-31
e=2.7

X = [0.25, 0.5, 0.75]
T_task4 = np.linspace(0, 300, 1000)
T_task5 = 1/np.linspace(50, 300, 1000)

def sigma_f(t, x):
    mee = linear_int("m*e", "Ga", "As", "P", x) #masa efektywna elektronu

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

    Eg = linear_int("Eg", "Ga", "As", "P", x, bow=0.19)

    Nc = (mee*K*t)**(3/2)
    Nv = (meh*K*t)**(3/2)
    sigma = e*t**(-3/2)*np.sqrt(Nc*Nv)*np.exp(-Eg/(2*K*t))

    return sigma, Eg


fig,ax = plt.subplots(1,2, figsize=(10,6))

for x in X:
    data_sigma = []
    data_eg = []

    for t in T_task4:
        sigma, eg = sigma_f(t,x)
        data_eg.append(eg)

    for t in T_task5:
        sigma, eg = sigma_f(1/t,x)
        data_sigma.append(sigma)

    log_sig_inv = [np.log(1/s) for s in data_sigma]

    ax[0].plot(T_task4, data_eg)
    ax[1].plot(T_task5, log_sig_inv)

ax[0].set_title("Eg vs Temperature")
ax[0].set_xlabel("Temperature [K]")
ax[0].set_ylabel("Eg [eV]")
ax[0].grid()

ax[1].set_title("ln(1/sigma) vs 1/T")
ax[1].set_xlabel("1/Temperature [1/K]")
ax[1].set_ylabel("ln(1/sigma) [-]")
ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
ax[1].set_yscale("log")
ax[1].grid()

legend = [f"x={x:.2f}" for x in X]
fig.legend(legend)

plt.tight_layout()
plt.savefig("results/task_4_5_eg_res.png", dpi=300)
plt.close()




