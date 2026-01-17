from numpy.linalg import eigvalsh
import matplotlib.pyplot as plt
import numpy as np

from utils import linear_int, binary_comps, strain

n=1000    #liczba punktow

d=2       #szerokosc studni w nm   
ratio=1/10. #jaki stosunek zajmuje studnia w stosunku do calej struktury
d=d/0.52917721*10
#bariera GaP
Veb=binary_comps["GaP"]["VBO"]
print(Veb)
Veb=Veb/27.21
Vdb=binary_comps["GaP"]["VBO"] + binary_comps["GaP"]["Eg"]
print(Vdb)
Vdb=Vdb/27.21

eps_xy_gaas, eps_z_gaas, delta_Ec_gaas, delta_Ev_gaas, eps_xy_gasb, eps_z_gasb, delta_Ec_gasb, delta_Ev_gasb = strain(0.5)
#studnia GaAs_0.5P_0.5
Ves=linear_int("VBO", "Ga", "P", "Sb", 0.5, params=binary_comps) + delta_Ev_gasb
print(Ves)
Ves=Ves/27.21
Vds=linear_int("VBO", "Ga", "P", "Sb", 0.5, params=binary_comps) + linear_int("Eg", "Ga", "P", "Sb", 0.5, bow=0.19, params=binary_comps) + delta_Ec_gasb
print(Vds)
Vds=Vds/27.21
me=0.05     #masa elektronu
md=0.5      #masa dziury

a=0.
b=d/ratio
h=(b-a)/(n-1)
con=1/2/h**2 #hbar**2/2m/h**2

X=[]
Ve=[]
Vd=[]

#generuje wektor polozenia i potencjalu CB i VB
for i in range(0,n):
    if a+i*h<a+(b-a)*(1.-ratio)/2. or a+i*h>b-(b-a)*(1.-ratio)/2:
        Ve.append(Veb)
        Vd.append(Vdb)
    else:
        Ve.append(Ves)
        Vd.append(Vds)
    X.append(a+i*h)

#generuje macierz
def M(V,m):
    M=con/m*2*np.identity(n)
    for i in range(0,n-1):
        M[i][i+1]=-con/m
        M[i+1][i]=-con/m
        M[i][i]=M[i][i]+V[i]
    M[n-1][n-1]=M[n-1][n-1]+V[n-1]
    return M

#diagonalizuje
el=eigvalsh(M(Ve,me))
dz=eigvalsh(M(np.array(Vd),-1.*md))
print(el[:3])
print(dz[:-4:-1])

ax = plt.axes()
fig = plt.gcf()
fig.set_size_inches(9, 6)
plt.plot(np.array(X)*0.52917721/10,np.array(Ve)*27.21)
plt.plot(np.array(X)*0.52917721/10,np.array(Vd)*27.21)
plt.xlabel(r'$X[nm]$', fontsize=18)
plt.ylabel(r'$V[eV]$', fontsize=18) 
plt.plot(np.array([a,b])*0.52917721/10, 2*[(el[:3])*27.21] , np.array([a,b])*0.52917721/10, 2*[(dz[:-4:-1])*27.21])
plt.savefig("./results/task_7_8_QW.png", dpi=300)
plt.close()

print("Eprzejscia="+str(((el[0])*27.21)-(dz[-1])*27.21))
print("Ebariery="+str((Veb - Vdb)*27.21))
