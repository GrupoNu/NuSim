#!/bin/sh

PROG="main"
PROGDIR="$(dirname "$(realpath $0)")"

touch "$PROGDIR"/data/ln-elecDens-clean.txt
python3 "$PROGDIR"/aux/ordclean.py
gcc -I/usr/include -g -Wall -O3 "$PROG".c -L/usr/lib -lpthread -lgsl -lgslcblas -lm -o "$PROGDIR"/"$PROG".out
cat "$PROGDIR"/data/ln-elecDens-clean.txt | "$PROGDIR"/"$PROG".out > "$PROGDIR"/data/output.txt
python3 "$PROGDIR"/aux/genfig.py
