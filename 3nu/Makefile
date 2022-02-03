PROG=interpol
OBJECTS=
CFLAGS=-std=gnu11 -g -Wall -O2 `pkg-config --cflags gsl`
LDLIBS=`pkg-config --libs gsl`
CC=gcc

$(PROG): $(OBJECTS)

clean:
	rm -rf $(PROG)
