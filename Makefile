PREFIX = /usr/local
FC = gfortran
FFLAGS ?= -g -O2 -Wall -mieee-fp
OBJECTS = geo.o kernel.o globals.o model.o
VPATH = src

all: libslezanbear.so

include deps.inc

%.o: %.f90
	$(FC) $(FFLAGS) -fpic -c $< -o $@

libslezanbear.so: $(OBJECTS)
	$(FC) -shared $^ -o $@

clean:
	$(RM) *.o *.mod *.so
