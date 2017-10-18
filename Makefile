PREFIX = /usr/local
FC = gfortran
FFLAGS = -g -O3 -march=native -finline-functions -fopenmp
LDLIBS = -lgdal

ifeq ($(FC),gfortran)
FFLAGS += -flto -funsafe-math-optimizations
FFLAGS += -Wall -Warray-temporaries -Wrealloc-lhs-all
endif
ifeq ($(FC),ifort)
FFLAGS += -ipo
FFLAGS += -warn all -Winline
endif

VPATH = src:src/gdal:src/test
OBJECTS = globals.o geo.o kernel.o model.o maputils.o gdal.o

all: libslezanbear.so bin/sbmap testbin/testr testbin/testw

include deps.inc

testbin/testw: gdal_test_common.o gdal_write_test.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@
testbin/testr: gdal_test_common.o gdal_read_test.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

bin/%: %.o $(OBJECTS)
	mkdir -p bin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

libslezanbear.so: globals.f90 geo.f90 kernel.f90 model.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -shared -fPIC $^ -o $@

clean:
	$(RM) *.o *.mod *.so
	$(RM) -r bin testbin
