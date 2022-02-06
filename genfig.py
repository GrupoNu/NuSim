#!/usr/bin/env python3

import sys

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['text.usetex'] = True
#plt.style.use('bmh')

N, T, R1, R2, R3, I1, I2, I3 = np.loadtxt(sys.stdin, unpack=True)
P_num = R1**2 + I1**2

# vacuum constants
m12 = 1.0
m22 = 2.0
dm2 = m22 - m12
th = np.pi/6

# matter constants
E = 1
Gf = 3
Ne = np.exp(1)
D = 2 * np.sqrt(2) * Gf * E * Ne

sin2THm = (dm2**2 * np.sin(2*th))**2 / ( (D - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )
Delta = np.sqrt( (D - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )

def P_e(t):
    return 1 - sin2THm * np.sin(Delta * t / (4*E))**2

#np.set_printoptions(threshold=sys.maxsize)
#print(P_num.T)

plt.plot(T, P_e(T), label=r'$P_{exact}$')
plt.plot(T, P_num,  label=r'$P_{num}$')
plt.xlabel(r'$t$', fontsize=20)
plt.ylabel(r'$P_e(t)$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'Neutrino Oscillations - Interpolated Density')
plt.savefig("./fig/nusim.png", dpi=300, format='png', bbox_inches="tight")
