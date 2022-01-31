PROG=main
OBJECTS=
CFLAGS=`pkg-config --cflags gsl` -g -Wall -std=c99 -O3
LDLIBS=`pkg-config --libs gsl`
CC=gcc

$(PROG): $(OBJECTS)

clean:
	rm -rf $(PROG)
