import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['text.usetex'] = True
#plt.style.use('bmh')

N, T, R1, R2, I1, I2 = np.loadtxt(dir_path+"/../data/output.txt", unpack=True)
P_num = R1**2 + I1**2

m12 = 1.0
m22 = 2.0
dm2 = m22 - m12
th = np.pi/6

def P_e(t):
    return 1 - (np.sin(2*th) * np.sin(dm2*t/4))**2

plt.plot(T, P_e(T), label=r'$P_{exato}$')
plt.plot(T, P_num,  label=r'$P_{num}$')
plt.xlabel(r'$t$', fontsize=20)
plt.ylabel(r'$P_e(t)$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'Neutrino Oscillations')
plt.savefig(dir_path+"/../fig/vacuum.png", dpi=300, format='png', bbox_inches="tight")
