import Radial_functions as Rf 
import matplotlib.pyplot as plt
import numpy as np
import radial as rad
import radiallog as radlog


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


for (n, l, lab) in orbitals: #calcing the energies
    r, P, E = radlog.radiallog(l=l, n=n, Z=Z, zeta=zeta, a=a, plot=False)
    energies[lab] = E
    radials[lab] = (r, P)


#Print orbital energies in a.u
print("\nNa I orbital energies (parametric V_a, a=0.2683) in a.u.:")
for lab in orbitals:
    key = lab[2]
    print(f"{key:>3s}: {energies[key]: .6f}") #look at me using keys

# configuration energies and excitation energies 
def E_config(occ):
    # occ is dict like {"1s":2,"2s":2,"2p":6,"3s":1}
    return sum(occ[k]*energies[k] for k in occ)

# Ground configuration and some excited ones from NIST
cfg_ground = {"1s":2,"2s":2,"2p":6,"3s":1}
cfg_3p     = {"1s":2,"2s":2,"2p":6,"3p":1}
cfg_3d     = {"1s":2,"2s":2,"2p":6,"3d":1}
cfg_4s     = {"1s":2,"2s":2,"2p":6,"4s":1}
cfg_4p     = {"1s":2,"2s":2,"2p":6,"4p":1}
cfg_5s     = {"1s":2,"2s":2,"2p":6,"5s":1}
cfg_4d     = {"1s":2,"2s":2,"2p":6,"4d":1}
cfg_4f     = {"1s":2,"2s":2,"2p":6,"4f":1}

configs = [
    ("2p^6 3s", cfg_ground),
    ("2p^6 3p", cfg_3p),
    ("2p^6 3d", cfg_3d),
    ("2p^6 4s", cfg_4s),
    ("2p^6 4p", cfg_4p),
    ("2p^6 5s", cfg_5s),
    ("2p^6 4d", cfg_4d),
    ("2p^6 4f", cfg_4f),
]

E0 = E_config(cfg_ground)
AU2CM = 219474.63   # cm^-1 per a.u.best name for variable trust me

print("\nConfiguration energies and excitation energies:")
print(f"{'Config':<12} {'E (a.u.)':>12} {'ΔE to g.s. (a.u.)':>16} {'ΔE (cm^-1)':>14}")
for name, occ in configs:
    Ec = E_config(occ)
    dE = Ec - E0
    print(f"{name:<12} {Ec:12.6f} {dE:16.6f} {dE*AU2CM:14.0f}")

# --- Plot multiple radial functions on one figure ---
to_plot = ["1s","2s","2p","3s","3p"]
plt.figure(figsize=(7,5))
for lab in to_plot:
    r, P = radials[lab]
    plt.plot(r, np.sqrt(r)*P, label=lab)  # same quantity plotted in radiallog
plt.xlabel("r (a.u.)")
plt.ylabel("radial function  $\\sqrt{r}\\,P$")
plt.title("Na I radial functions (parametric $V_a$)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()