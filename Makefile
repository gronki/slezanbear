PREFIX = /usr/local
FC = gfortran
LDLIBS = -lgdal

#---------------------------------------------------------------------------

ifeq ($(FC),gfortran)
FFLAGS = -g -O3 -march=native -flto -fopenmp
FFLAGS += -Wall -Warray-temporaries -Wrealloc-lhs-all
endif

ifeq ($(FC),ifort)
FFLAGS = -g -O3 -xHost -ipo -qopenmp -fp-model precise
FFLAGS += -Winline
endif

FFLAGS ?= -g -O2

#---------------------------------------------------------------------------

VPATH = src:src/gdal:src/test
OBJECTS = globals.o geo.o kernel.o model.o mapio.o gdal.o

#---------------------------------------------------------------------------

all: libslezanbear.so bin/sbmap bin/sbcalib bin/sbcalibrand bin/sbverify \
	 testbin/testr testbin/testw

#---------------------------------------------------------------------------

include deps.inc

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

#---------------------------------------------------------------------------

bin/%: %.f90 $(OBJECTS)
	mkdir -p bin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

libslezanbear.so: globals.f90 geo.f90 kernel.f90 model.f90
	$(FC) $(FFLAGS) -shared -fPIC $^ -o $@

#---------------------------------------------------------------------------

testbin/testw: gdal_test_common.o gdal_write_test.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@
testbin/testr: gdal_test_common.o gdal_read_test.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

#---------------------------------------------------------------------------

clean:
	$(RM) *.o *.mod *.so
	$(RM) -r bin testbin
