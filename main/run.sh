#!/bin/sh

PROG="main"
PROGDIR="$(dirname "$(realpath $0)")"

cd "$PROGDIR"
make
cat "$PROGDIR"/ln-elecDens.txt | "$PROGDIR"/"$PROG" | "$PROGDIR"/genfig.py
