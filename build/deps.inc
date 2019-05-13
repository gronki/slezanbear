sbcalibrand.o: model.o globals.o mapio.o gdal.o geo.o
geo.o: globals.o
sbmap.o: model.o globals.o mapio.o gdal.o geo.o
model.o: kernel.o globals.o geo.o
gdal_read_test.o: gdal.o gdal_test_common.o geo.o
gdal_test_common.o: geo.o
gdal.o: gdalproto_interf.f90 gdalproto_type.f90 gdalproto_proc.f90
kernel.o: globals.o
sbverify.o: model.o globals.o mapio.o gdal.o geo.o
gdal_write_test.o: gdal.o gdal_test_common.o geo.o
mapio.o: gdal.o globals.o geo.o
sbcalib.o: model.o globals.o mapio.o gdal.o geo.o
