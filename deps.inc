gdal_test_common.o: geo.o 
geo.o           : globals.o 
kernel.o        : globals.o 
mapio.o         : globals.o geo.o gdal.o 
gdal_read_test.o: gdal_test_common.o geo.o gdal.o 
gdal_write_test.o: gdal_test_common.o geo.o gdal.o 
sbcalib.o       : model.o globals.o geo.o gdal.o mapio.o 
sbcalibrand.o   : model.o globals.o geo.o gdal.o mapio.o 
sbmap.o         : model.o globals.o geo.o gdal.o mapio.o 
model.o         : kernel.o globals.o geo.o model.o 
