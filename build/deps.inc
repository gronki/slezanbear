sbcalibrand.o: globals.o gdal.o model.o geo.o mapio.o
geo.o: globals.o
sbmap.o: globals.o gdal.o model.o geo.o mapio.o
model.o: geo.o globals.o kernel.o
gdal.o: gdalproto_proc.f90 gdalproto_type.f90 gdalproto_interf.f90
sbverify.o: globals.o gdal.o model.o geo.o mapio.o
gdal_read_test.o: geo.o gdal_test_common.o gdal.o
mapio.o: geo.o globals.o gdal.o
sbcalib.o: globals.o gdal.o model.o simplex.o geo.o mapio.o
gdal_test_common.o: geo.o
gdal_write_test.o: geo.o gdal_test_common.o gdal.o
kernel.o: globals.o
simplex.o: globals.o
