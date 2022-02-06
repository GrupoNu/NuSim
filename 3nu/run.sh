#!/bin/sh

PROG="mass"
PROGDIR="$(dirname "$(realpath $0)")"

make -s clean
cd "$PROGDIR"
make -s
cat "$PROGDIR"/const.txt | "$PROGDIR"/"$PROG" | "$PROGDIR"/genfig.py
make -s clean
