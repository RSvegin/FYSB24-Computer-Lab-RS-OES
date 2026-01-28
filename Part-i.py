import Radial_functions as Rf 
import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(0,20,1000) 
Z = 1
P = Rf.P1s(r,Z)
plt.plot(r,P)
plt.show
