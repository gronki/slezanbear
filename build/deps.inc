sbcalibrand.o: mapio.o geo.o gdal.o globals.o model.o
sbverify.o: mapio.o geo.o gdal.o globals.o model.o
model.o: globals.o kernel.o geo.o
kernel.o: globals.o
gdal.o: gdalproto_type.f90 gdalproto_proc.f90 gdalproto_interf.f90
mapio.o: geo.o gdal.o globals.o
gdal_read_test.o: gdal_test_common.o geo.o gdal.o
geo.o: globals.o
gdal_write_test.o: gdal_test_common.o geo.o gdal.o
gdal_test_common.o: geo.o
sbmap.o: mapio.o geo.o gdal.o globals.o model.o
sbcalib.o: model.o mapio.o geo.o gdal.o globals.o datahandling.o simplex.o
simplex.o: globals.o
datahandling.o: globals.o
