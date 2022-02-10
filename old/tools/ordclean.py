#!/usr/bin/env python3

import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))

def importa():
    L = []
    for linha in sys.stdin:
        e = linha[:].split()
        f1 = float(e[0])
        f2 = float(e[1])
        F = [f1,f2]
        L.append(F)
    return L

def mean(lista):
    s = sum(lista)
    media = s/len(lista)
    return media

def limpa_lista(L):
    limpa = []
    i = 0
    while i < len(L):
        var = True
        e_aux = []
        e_aux.append(L[i][1])
        j = 1

        while var and (i+j < len(L)):
            if L[i][0] == L[i+j][0]:
                e_aux.append(L[i+j][1])
                j += 1
            else:
                var = False

        limpa.append([L[i][0], mean(e_aux)])

        i += j

    return limpa

def escreve(limpa):
    for n in range(len(limpa)):
        x = float(limpa[n][0])
        y = float(limpa[n][1])
        if y < 0:
            print(("{:.7f}".format(x))+("  {:.6f}".format(y)))
        else:
            print(("{:.7f}".format(x))+("   {:.6f}".format(y)))

def main():
    Ldesord = importa()
    L = sorted(Ldesord)
    limpa = limpa_lista(L)
    escreve(limpa)

if __name__ == '__main__':
    main()
