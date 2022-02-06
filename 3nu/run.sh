#!/bin/sh

PROG="main"
PROGDIR="$(dirname "$(realpath $0)")"

make -s clean
cd "$PROGDIR"
make -s
cat "$PROGDIR"/ln-elecDens.txt | "$PROGDIR"/"$PROG" | "$PROGDIR"/genfig.py
make -s clean
