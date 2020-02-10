import numpy as np
import matplotlib.pyplot as plt


data1 = np.loadtxt('tovqmc.dat')
data2 = np.loadtxt('tovtw.dat')

ener1 = []
press1 = []

ener2 = []
press2 = []

#ener1, press1 = np.loadtxt(filename, unpack =True)

for k in range((len(data1))):
    ener1.append(data1[k,1])  
    press1.append(data1[k,2])  
        
for k in range((len(data2))):
    ener2.append(data2[k,1])  
    press2.append(data2[k,2])  



plt.plot(ener1,press1,label="QMC")
plt.plot(ener2,press2,label="TW")
#plt.axis([-10, 60, 0, 11])
plt.title("EoS")
plt.grid(True)
plt.xlabel("$\u03B5$ (fm$^{-4}$)")#energy density
plt.ylabel("P (fm$^{-4}$)")#pressure
plt.legend(loc='upper left')
plt.savefig("eos.eps")
plt.show()
