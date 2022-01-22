import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['text.usetex'] = True
#plt.style.use('bmh')

BRUTE_X, BRUTE_Y = np.loadtxt(str(dir_path)+"/../data/ln-elecDens-clean.txt", unpack=True)
X, Y_STEFF = np.loadtxt(str(dir_path)+"/../data/ln-elecDens-clean.txt", unpack=True)

plt.plot(BRUTE_X, BRUTE_Y, color='#000000', marker='o',
         label=r'dados')
plt.plot(X, Y_STEFF, label=r'$y_{stef}$')
plt.xlabel(r'$R/R_{\odot}$', fontsize=20)
plt.ylabel(r'$\log(N_e/N_A)$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'Interpolation - Solar Electron Density')
plt.savefig(str(dir_path)+"/../fig/2-steff.png", dpi=300, format='png', bbox_inches="tight")
