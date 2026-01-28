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
r,P,E=radlog.radiallog(l,n,Z)
plt.savefig("radiallog_function.png")

#task 3
l=0
n=[1,2,3,6,9]
rad_energies=[]
radlog_energies=[]
for i in n:
    r,P,E=rad.radial(l,i,Z)
    rad_energies.append(E)
    r,p,E=radlog.radiallog(l,i,Z)
    radlog_energies.append(E)

gridpoints=np.linspace(1,10000,1)

