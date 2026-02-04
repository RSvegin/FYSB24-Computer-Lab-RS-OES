import Radial_functions as Rf 
import matplotlib.pyplot as plt
import numpy as np
import radial as rad
import radiallog as radlog


r = np.linspace(0,20,1000) 
Z = 1


# Commenting out previous radial density plots
# P = Rf.P1s(r,Z)
# P3s=Rf.P3s(r,Z)
# P3p=Rf.P3p(r,Z)
# P3d=Rf.P3d(r,Z)

# radial density plots

# plt.plot(r,P)
# plt.plot(r,P3s)
# plt.plot(r,P3p)
# plt.plot(r,P3d)
# plt.legend(["p","p3s","p3p","p3d"])
# plt.show()
# plt.savefig('P1s.png')

# task 2
# effective potential
r=np.linspace(0.0001,20,1000)
def V(r,l):
    return -(Z/r)+l*(l+1)/(2*r**2)
Vs=V(r,0)
Vp=V(r,1)
Vd=V(r,2)

Vs=np.clip(Vs,-2,1)
Vp=np.clip(Vp,-2,1)
Vd=np.clip(Vd,-2,1)
# plt.plot(r,Vs)
# plt.plot(r,Vp)
# plt.plot(r,Vd)
# plt.legend(["Vs","Vp","Vd"])
# plt.show,plt.savefig("effective_potential")

#testing radial function
l=0
n=1
Z=1
r,P,E=rad.radial(l,n,Z)
#plt.savefig("radial_function.png")

#testing radiallog function
r,P,E,grid_points=radlog.radiallog(l,n,Z)
#plt.savefig("radiallog_function.png")

#task 3
import pandas as pd
l=0
n=[1,2,3,6,9]
ana_energies=[]
rad_energies=[]
radlog_energies=[]
rad_gridpoints=[]
radlog_gridpoints=[]
rs = []
Ps = []
labels = []

for i in n:
    # --- radial.py ---
    r_lin, P_lin, E = rad.radial(l, i, Z)
    rad_energies.append(E)
    rad_gridpoints.append(10000)

    # --- radiallog.py ---
    r_log, P_log_trans, E_log, grid_points = radlog.radiallog(l, i, Z)
    radlog_energies.append(E_log)
    radlog_gridpoints.append(grid_points)

    # analytical
    ana_energies.append(-Z**2/(2*i**2))

    #store log solution for plotting
    P_true = np.sqrt(r_log) * P_log_trans   # convert transformed P(ρ) to real P(r)
    rs.append(r_log)
    Ps.append(P_true)
    labels.append(f"{i}s")


data={"Analytical Energies":ana_energies,
      "Radial Energies":rad_energies,
      "RadialLog Energies":radlog_energies,
      "Radial Grid Points":rad_gridpoints, #making table
      "RadialLog Grid Points":radlog_gridpoints}
df=pd.DataFrame(data) #creating dataframe
print(df)


plt.figure(figsize=(8,6))

for r, P, label in zip(rs, Ps, labels):
    plt.plot(r, P, label=label)

plt.xlabel("r (a.u.)")
plt.ylabel("P(r)")
plt.title("Radial Wave Functions for 1s, 2s, 3s, 6s, 9s")
plt.legend(title="Orbitals")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig("radial_wavefunctions_1s_2s_3s_6s_9s.png")


# Create subplots once before the loop
plt.subplots(2,2, figsize=(12, 10))

n = [1,2,3,4]
for i in n:
    r,P,E=rad.radial(l,i,Z, plot=False) #getting radial function for each n
    dr = r[1]-r[0]
    r_expectation_value = sum(r*P**2)*dr
    
    plt.subplot(2,2,n.index(i)+1)
    plt.plot(r,P**2)
    plt.vlines(r_expectation_value, ymin=0, ymax=max(P**2), colors='r', linestyles='dashed', label='<r>')
    plt.title(f"Radial Probability Density for n={i}, l={l}\n<r>={r_expectation_value:.2f}")
    plt.xlabel("r")
    plt.ylabel("Radial Probability Density")

plt.tight_layout()
plt.savefig("radial_probability_n1-4_l0.png")
plt.close()
