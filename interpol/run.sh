#!/bin/sh

PROGDIR="$(dirname "$(realpath $0)")"

touch "$PROGDIR"/data/ln-elecDens-clean.txt
python3 "$PROGDIR"/aux/ordclean.py
gcc -I/usr/include -g -Wall -O3 interpol.c -L/usr/lib -lpthread -lgsl -lgslcblas -lm -o "$PROGDIR"/bin-interpol
cat "$PROGDIR"/data/ln-elecDens-clean.txt | "$PROGDIR"/bin-interpol > "$PROGDIR"/data/ln-elecDens-interpol.txt
python3 "$PROGDIR"/aux/genfig.py
