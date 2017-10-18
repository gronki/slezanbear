gdal_test_common.o: geo.o 
geo.o           : globals.o 
kernel.o        : globals.o 
maputils.o      : globals.o geo.o gdal.o 
gdal_read_test.o: gdal_test_common.o geo.o gdal.o 
gdal_write_test.o: gdal_test_common.o geo.o gdal.o 
sbmap.o         : maputils.o globals.o geo.o gdal.o model.o 
model.o         : kernel.o globals.o geo.o model.o 
