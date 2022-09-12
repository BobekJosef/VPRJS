from typing import List

import numpy as np
import matplotlib.pyplot as plt

data_file = np.loadtxt('/home/josef/School/smash/build/data/1/output.dat')
y = data_file[0, :]
pt_pi = data_file[1, :]
pt_K = data_file[2, :]
pt_p = data_file[3, :]
v2_pi = data_file[4, :]
v2_K = data_file[5, :]
v2_p = data_file[6, :]

pt_X = []
for k in range(0, 20):
    pt_X.append(k * 0.1 + 0.05)


y_X = []
for k in range(0, 20):
    y_X.append(k * 0.4 + 0.1 - 4)


plt.plot(pt_X, pt_pi, 'b--', label=r'$\pi^{\pm}$')
plt.plot(pt_X, pt_K, 'g--', label=r'$K^{\pm}$')
plt.plot(pt_X, pt_p, 'r--', label=r'$p^{\pm}$')
plt.legend(loc='best')
plt.xlabel(r'$p_{T}$ [GeV]')
plt.ylabel(r'$\frac{1}{N_{ev}}\frac{1}{p_{T}}\frac{dN}{dp_{T}\ dy\ d\phi}\ [A.U.]$', fontsize=12)
plt.yscale("log")
plt.show()

plt.plot(y_X,  y, 'r--', label=r'$N_{ch}$')
plt.xlabel("y")
plt.legend(loc='best')
plt.ylabel(r'$\frac{1}{N_{ev}}\frac{dN}{dy\ d\phi}\ [A.U.]$', fontsize=12)
plt.show()

plt.plot(pt_X, v2_pi,  'b--', label=r'$\pi^{\pm}$')
plt.plot(pt_X, v2_K,  'g--', label=r'$K^{\pm}$')
plt.plot(pt_X, v2_p,  'r--', label=r'$p^{\pm}$')
plt.legend(loc='best')
plt.xlabel(r'$p_{T}$ [GeV]')
plt.ylabel(r'$v_{2}$')
plt.show()