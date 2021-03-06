
FC = gfortran
ARCH_FLAGS ?= -march=native
FFLAGS = -O3 -funsafe-math-optimizations $(ARCH_FLAGS) -fopenmp -flto
# FFLAGS = -g -Og
FFLAGS += -fimplicit-none -ffree-line-length-none \
				  -pedantic -Wall -Wno-unused-dummy-argument -Warray-temporaries

#---------------------------------------------------------------------------

VPATH = ../src ../src/gdal ../src/test
OBJECTS = globals.o geo.o kernel.o model.o mapio.o gdal.o simplex.o \
	datahandling.o
PROGRAMS = sbmap sbcalib sbcalibrand sbverify
TESTS = testw testr
LDLIBS = -lgdal

#---------------------------------------------------------------------------

all: libslezanbear.so $(PROGRAMS) $(TESTS)

#---------------------------------------------------------------------------

include deps.inc

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

#---------------------------------------------------------------------------

$(PROGRAMS): %: %.o $(OBJECTS)
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

libslezanbear.so: globals.f90 geo.f90 kernel.f90 model.f90
	$(FC) $(FFLAGS) -shared -fPIC $^ -o $@

testw: gdal_test_common.o gdal_write_test.o $(OBJECTS)
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@
testr: gdal_test_common.o gdal_read_test.o $(OBJECTS)
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

#---------------------------------------------------------------------------

install: $(PROGRAMS) $(TESTS)
	install -d ../bin
	install $(PROGRAMS) ../bin
	install -d ../testbin
	install $(TESTS) ../testbin

#---------------------------------------------------------------------------

clean:
	$(RM) *.o *.mod *.so $(PROGRAMS) $(TESTS)
	$(RM) -r ../bin ../testbin
