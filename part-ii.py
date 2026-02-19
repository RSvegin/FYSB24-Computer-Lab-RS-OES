import Radial_functions as Rf 
import matplotlib.pyplot as plt
import numpy as np
import radial as rad
import radiallog as radlog
import pandas as pd
import scipy.optimize as opt
#numbers for Na I
Z = 11
N = 11
zeta = Z - N + 1    # = 1 for Na I
a = 0.2683          # lab value

orbitals=[

    (1, 0, "1s"),
    (2, 0, "2s"),
    (2, 1, "2p"),
    (3, 0, "3s"),
    (3, 1, "3p"),
    (3, 2, "3d"),
    (4, 0, "4s"),
    (4, 1, "4p"),
    (4, 2, "4d"),
    (4, 3, "4f"),
    (5, 0, "5s"),
] #which orbitals to compute with n and l and their name (i keep forgetting)



energies=[]
for (n,l, lab) in orbitals:
    E=radlog.radiallog(l=l,n=n,Z=Z,zeta=zeta,a=a, plot=False)
    energies.append(E[2])

data={"orbital":lab,
    "Computed Energies":energies,}
df=pd.DataFrame(data)
print(df)


#plotting radial functions for 1s, 2s, 3p, 3d,
r = np.linspace(0,20,1000) 
Z = 3

P1s=Rf.P1s(r,Z)
P2s=Rf.P2s(r,Z)
P3p=Rf.P3p(r,Z)
P3d=Rf.P3d(r,Z)


plt.plot(r,P2s)
plt.plot(r,P3p)
plt.plot(r,P3d)
plt.plot(r,P1s)
plt.xlabel("r (a.u.)")
plt.ylabel("P(r)")
plt.title("Radial Density Functions for 1s, 2s, 3p, 3d. Z=3")
plt.legend(["P1s","P2s","P3p","P3d"])
plt.show()
plt.savefig('part2_radial_functions.png')
#---------task sieben----------


nist_li_centroids_au = {
    '2p':  0.067907105,
    '3s':  0.123960204,
    '3p':  0.140906400,
    '3d':  0.142536310,
    '4s':  0.159526684,
    '4p':  0.166167497,
    '4d':  0.166868453,
    '4f':  0.166899473,
}

Z = 3
N = 3
zeta = Z - N + 1

l_map = {'s':0, 'p':1, 'd':2, 'f':3}


def get_energy(n, l, a):
    # radiallog returns (r, P, E, grid_points)
    out = radlog.radiallog(l=l, n=n, Z=Z, zeta=zeta, a=a, plot=False)

    E = out[2]     # <-- ENERGY IS INDEX 2

    return float(E)


def residuals(a_array):

    a = float(a_array[0])

    # ground state (2s)
    E2s = get_energy(2, 0, a)

    res = []

    for key in nist_li_centroids_au:

        Enist = float(nist_li_centroids_au[key])

        n = int(key[0])
        l = l_map[key[1]]

        Ecalc = get_energy(n, l, a)

        deltaE = Ecalc - E2s

        res.append(deltaE - Enist)

    return np.array(res)


a0 = np.array([0.3])

a_opt, _ = opt.leastsq(residuals, a0)

a_best = float(a_opt[0])

print("\nOptimized a =", a_best)

print("\nEffect of changing a (example: 3p level):") #shows how E varys with a 

for a_test in [0.15, 0.25, 0.35, 0.45]:

    E2s = get_energy(2, 0, a_test)
    E3p = get_energy(3, 1, a_test)

    delta = E3p - E2s

    print(f"a = {a_test:.2f}   ΔE(3p) = {delta:.6f}")

# ------------------
# Comparison table
# ------------------

E2s = get_energy(2, 0, a_best)

print("\nState   Computed (a.u.)   NIST (a.u.)")

for key in nist_li_centroids_au:

    Enist = nist_li_centroids_au[key]

    n = int(key[0])
    l = l_map[key[1]]

    Ecalc = get_energy(n, l, a_best)
    deltaE = Ecalc - E2s

    print(f"{key:4s}   {deltaE:12.6f}   {Enist:10.6f}")



#calc energies for all levels with radlog
#subtract against nist values
#minmise this difference by varying a
#plot difference v a
