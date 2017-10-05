import numpy as np
import matplotlib.pyplot as plt
from gdal import Open as OpenGdal

#-------------------------------------------------------------------------------

fn0 = 'SVDNB_npp_20150101-20151231_75N060W_{}_v10_c201701311200.avg_rade9.tif'

# vcm - viirs cloud mask
# -orm = outlier removed
# -ntl = background (non-lights) removed
gd = OpenGdal(fn0.format('vcm-orm'))

print gd.RasterXSize, gd.RasterYSize
gt = gd.GetGeoTransform()

#-------------------------------------------------------------------------------

xy2geo = lambda g,x,y: (g[3] + g[5] * y, g[0] + g[1] * x)
geo2xy = lambda g,f,l: ((l - g[0]) / g[1], (f - g[3]) / g[5])

#-------------------------------------------------------------------------------

# geo0 = 51.11, 17.03 # wroclaw
# geo0 = 52.22, 21.01 # warszawa
# geo0 = 50.816, 15.383 # Orle
geo0 = 50.846015, 16.698650 # tapadla
geor = 0.4
f0, l0 = geo0[0] + geor, geo0[1] - geor
f1, l1 = geo0[0] - geor, geo0[1] + geor

#-------------------------------------------------------------------------------

x0,y0 = geo2xy(gt, f0, l0)
x1,y1 = geo2xy(gt, f1, l1)
x0,y0,x1,y1 = [ int(round(x)) for x in (x0,y0,x1,y1) ]

print x0, y0
print x1, y1
print x1-x0, y1-y0

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize = (12,6))

#-------------------------------------------------------------------------------

axes[0].imshow(np.sqrt(gd.ReadAsArray(x0,y0,x1-x0,y1-y0)),
    extent = [l0,l1,f1,f0],
    interpolation = 'none', cmap = 'hot')

#-------------------------------------------------------------------------------

gd1 = OpenGdal('eudem_dem_5deg_n50e015.tif')
gt1 = gd1.GetGeoTransform()

x0,y0 = geo2xy(gt1, f0, l0)
x1,y1 = geo2xy(gt1, f1, l1)
x0,y0,x1,y1 = [ int(round(x)) for x in (x0,y0,x1,y1) ]

#-------------------------------------------------------------------------------

axes[1].imshow(gd1.ReadAsArray(x0, y0, x1-x0, y1-y0),
    extent = [l0,l1,f1,f0],
    interpolation = 'none', cmap = 'BrBG_r')

#-------------------------------------------------------------------------------

plt.show()
