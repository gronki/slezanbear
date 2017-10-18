import gdal
from sys import argv, stdout

f = open('demlist.txt','w')
for fn in argv[1:]:
    gd = gdal.Open(fn)
    gt = gd.GetGeoTransform()
    # fortran format: A32,1X,I6,1X,I6,4(1X,F9.4)
    f.write("{fn:32} {nx:6} {ny:6} {x0:9.4f} {x1:9.4f} {y0:9.4f} {y1:9.4f}\n".format(
        fn = fn,
        nx = gd.RasterXSize,
        ny = gd.RasterYSize,
        x0 = gt[0],
        y0 = gt[3],
        x1 = gt[0] + gt[1] * (gd.RasterXSize),
        y1 = gt[3] + gt[5] * (gd.RasterYSize),
    ))
f.close()
