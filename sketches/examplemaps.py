import numpy as np
import matplotlib.pyplot as plt
from gdal import Open as OpenGdal

#-------------------------------------------------------------------------------

fn0 = 'SVDNB_npp_20150101-20151231_75N060W_{}_v10_c201701311200.avg_rade9.tif'

# vcm - viirs cloud mask
# vcm-orm = outlier removed
# vcm-ntl = background (non-lights) removed
# vcm-orm-ntl = both
gd = OpenGdal(fn0.format('vcm'))

print gd.RasterXSize, gd.RasterYSize
gt = gd.GetGeoTransform()

#-------------------------------------------------------------------------------

xy2geo = lambda g,x,y: (g[3] + g[5] * y, g[0] + g[1] * x)
geo2xy = lambda g,f,l: ((l - g[0]) / g[1], (f - g[3]) / g[5])

#-------------------------------------------------------------------------------

# geo0 = 51.11, 17.03 # wroclaw
# geo0 = 52.22, 21.01 # warszawa
# geo0 = 50.816, 15.383 # Orle
# geo0 = 50.846015, 16.698650 # tapadla
geo0 = 50.995681, 16.901729
geor = 0.6
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

map_i = np.ndarray((x1-x0, y1-y0), order = 'F', dtype = np.float32)
print map_i.flags.f_contiguous
gd.ReadAsArray(x0, y0, x1-x0, y1-y0, buf_obj = map_i)
axes[0].imshow(np.sqrt(map_i),
    extent = [l0,l1,f1,f0],
    interpolation = 'none', cmap = 'hot')

#-------------------------------------------------------------------------------

gd1 = OpenGdal('eudem_dem_5deg_n50e015.tif')
gt1 = gd1.GetGeoTransform()

x0,y0 = geo2xy(gt1, f0, l0)
x1,y1 = geo2xy(gt1, f1, l1)
x0,y0,x1,y1 = [ int(round(x)) for x in (x0,y0,x1,y1) ]

#-------------------------------------------------------------------------------

map_h = np.ndarray((x1-x0, y1-y0), order = 'F', dtype = np.float32)
gd1.ReadAsArray(x0, y0, x1-x0, y1-y0, buf_obj = map_h)
print map_h.flags.f_contiguous
axes[1].imshow(map_h,
    extent = [l0,l1,f1,f0],
    interpolation = 'none', cmap = 'BrBG_r')

#-------------------------------------------------------------------------------

plt.show()
