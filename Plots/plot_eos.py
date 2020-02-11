import numpy as np
import matplotlib.pyplot as plt


rho1, ener1, press1 = np.loadtxt('tovqmc.dat', unpack=True)
data2 = np.loadtxt('tovtw.dat')

ener2 = []
press2 = []

for k in range((len(data2))):
    ener2.append(data2[k, 1])
    press2.append(data2[k, 2])

plt.plot(ener1, press1, label="QMC")
plt.plot(ener2, press2, label="TW")
#plt.axis([-10, 60, 0, 11])
plt.title("EoS")
plt.grid(True)
plt.xlabel("$\u03B5$ (fm$^{-4}$)")#energy density
plt.ylabel("P (fm$^{-4}$)")#pressure
plt.legend(loc='upper left')
plt.savefig("eos.eps")
plt.show()
