#!/bin/sh

PROG="gsl-interpol"
PROGDIR="$(dirname "$(realpath $0)")"

cd "$PROGDIR"
make
"$PROGDIR"/"$PROG" | "$PROGDIR"/genfig.py
