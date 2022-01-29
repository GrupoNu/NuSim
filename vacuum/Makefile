PROG=vacuum
OBJECTS=
CFLAGS=`pkg-config --cflags gsl` -g -Wall -O3
LDLIBS=`pkg-config --libs gsl`
CC=gcc

$(PROG): $(OBJECTS)

clean:
	rm -rf $(PROG)
