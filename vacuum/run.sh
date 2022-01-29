#!/bin/sh

PROG="vacuum"
PROGDIR="$(dirname "$(realpath "$0")")"

cd "$PROGDIR"
make
"$PROGDIR"/"$PROG" | "$PROGDIR"/genfig.py
