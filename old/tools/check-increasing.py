#!/usr/bin/env python3

import sys
import numpy as np

x, y = np.loadtxt(sys.stdin, unpack=True)

for i in range(len(x)-1):
    oldx = x[i]
    newx = x[i+1]
    if newx <= oldx:
        print("Failed! Not increasing!")
        quit()

print("Success! It is increasing!")
