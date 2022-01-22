#!/bin/sh

PROG="vacuum"
PROGDIR="$(dirname "$(realpath "$0")")"

gcc -I/usr/include -g -Wall -O3 "$PROGDIR"/"$PROG".c -L/usr/lib -lgsl -lgslcblas -lm -o "$PROGDIR"/"$PROG".out
"$PROGDIR"/"$PROG".out > "$PROGDIR"/data/output.txt
python3 "$PROGDIR"/aux/genfig.py
