#!/bin/sh

PROGDIR="$(dirname "$(realpath $0)")"

gcc -I/usr/include -g -Wall -O3 interpol.c -L/usr/lib -lpthread -lgsl -lgslcblas -lm -o "$PROGDIR"/bin-interpol
"$PROGDIR"/bin-interpol > data-interpol.txt
python3 "$PROGDIR"/genfig.py
