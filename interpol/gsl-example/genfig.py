#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['text.usetex'] = True
#plt.style.use('bmh')

BRUTE_X, BRUTE_Y = np.loadtxt('./data-brute.txt', unpack=True)
X, Y_CUBIC, Y_AKIMA, Y_STEFF = np.loadtxt('./data-interpol.txt', unpack=True)

plt.plot(BRUTE_X, BRUTE_Y, color='#000000', marker='o',
         label=r'dados')
plt.plot(X, Y_CUBIC, label=r'$y_{cubic}$')
plt.plot(X, Y_AKIMA, label=r'$y_{akima}$')
plt.plot(X, Y_STEFF, label=r'$y_{steffen}$')
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$y$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'GSL Interpolation Example')
plt.savefig('examp_interpol.png', dpi=300, format='png', bbox_inches="tight")
