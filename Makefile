PREFIX = /usr/local
FC = ifort
FFLAGS = -g -O2

ifeq ($(FC),gfortran)
FFLAGS += -mieee-fp
FFLAGS += -Wall -Warray-temporaries -Wrealloc-lhs-all
FFLAGS += -fopenmp
endif
ifeq ($(FC),ifort)
FFLAGS += -fp-model precise
FFLAGS += -warn all -Winline
FFLAGS += -qopenmp
endif

VPATH = src
SOURCES = globals.f90 geo.f90 kernel.f90 model.f90
OBJECTS = $(SOURCES:.f90=.o)

all: libslezanbear.so testbin/testr testbin/testw

include deps.inc

testbin/testw: LDLIBS += -lgdal
testbin/testw: test/gdal_test_common.f90 test/gdal_write_test.f90 gdal.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@
testbin/testr: LDLIBS += -lgdal
testbin/testr: test/gdal_test_common.f90 test/gdal_read_test.f90 gdal.o $(OBJECTS)
	mkdir -p testbin
	$(FC) $(FFLAGS) $^ $(LDLIBS) -o $@

gdal.o: gdal/gdal.F90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

libslezanbear.so: $(SOURCES)
	$(FC) $(FFLAGS) -shared -fPIC $^ $(LDLIBS) -o $@

clean:
	$(RM) *.o *.mod *.so
	$(RM) -r bin testbin
