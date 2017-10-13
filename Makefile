PREFIX = /usr/local
FC = gfortran
FFLAGS = -g -O2 -Wall -mieee-fp -fopenmp -Warray-temporaries -Wrealloc-lhs-all
OBJECTS = geo.o kernel.o globals.o model.o
VPATH = src

all: libslezanbear.so

include deps.inc

%.o: %.f90
	$(FC) $(FFLAGS) -fpic -c $< -o $@

libslezanbear.so: $(OBJECTS)
	$(FC) $(FFLAGS) -shared $^ $(LDLIBS) -o $@

clean:
	$(RM) *.o *.mod *.so
