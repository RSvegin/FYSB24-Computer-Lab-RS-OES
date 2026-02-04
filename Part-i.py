import Radial_functions as Rf 
import matplotlib.pyplot as plt
import numpy as np
import radial as rad
import radiallog as radlog


r = np.linspace(0,20,1000) 
Z = 1
P = Rf.P1s(r,Z)
P3s=Rf.P3s(r,Z)
P3p=Rf.P3p(r,Z)
P3d=Rf.P3d(r,Z)

#radial density plots

#plt.plot(r,P)
#plt.plot(r,P3s)
#plt.plot(r,P3p)
#plt.plot(r,P3d)
#plt.legend(["p","p3s","p3p","p3d"])
#plt.show
#plt.savefig('P1s.png')

#task 2
#effective potential
r=np.linspace(0.0001,20,1000)
def V(r,l):
    return -(Z/r)+l*(l+1)/(2*r**2)
Vs=V(r,0)
Vp=V(r,1)
Vd=V(r,2)

Vs=np.clip(Vs,-2,1)
Vp=np.clip(Vp,-2,1)
Vd=np.clip(Vd,-2,1)
#plt.plot(r,Vs)
#plt.plot(r,Vp)
#plt.plot(r,Vd)
plt.legend(["Vs","Vp","Vd"])
plt.show,plt.savefig("effective_potential")

#testing radial function
l=0
n=1
Z=1
r,P,E=rad.radial(l,n,Z)
plt.savefig("radial_function.png")

#testing radiallog function
r,P,E,grid_points=radlog.radiallog(l,n,Z)
plt.savefig("radiallog_function.png")

#task 3
import pandas as pd
l=0
n=[1,2,3,6,9]
ana_energies=[]
rad_energies=[]
radlog_energies=[]
rad_gridpoints=[]
radlog_gridpoints=[]
for i in n:
    r,P,E=rad.radial(l,i,Z)
    rad_energies.append(E)
    r,p,E,grid_points=radlog.radiallog(l,i,Z) #adding all the data for the table
    radlog_energies.append(E)
    rad_gridpoints.append(10000)
    radlog_gridpoints.append(grid_points)
    ana_energies.append(-Z**2/(2*i**2))
    

data={"Analytical Energies":ana_energies,
      "Radial Energies":rad_energies,
      "RadialLog Energies":radlog_energies,
      "Radial Grid Points":rad_gridpoints, #making table
      "RadialLog Grid Points":radlog_gridpoints}
df=pd.DataFrame(data) #creating dataframe
print(df)




 


#task 4


