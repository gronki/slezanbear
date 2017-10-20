# coding: utf-8
import gdal
import numpy as np
import matplotlib.pyplot as plt

# fn0 = 'model_izer_hires'
fn0 = 'model_wroc_hires'

gd = gdal.Open(fn0 + '.tif')
gt = gd.GetGeoTransform()
nx = gd.RasterXSize
ny = gd.RasterYSize
arr = gd.ReadAsArray()

fs = ( 8 * np.sqrt(nx/float(ny)), 6 / np.sqrt(nx/float(ny)) )

ex = [ gt[0], gt[0] + nx * gt[1], 
       gt[3] + ny * gt[5], gt[3] ]

i2m = lambda i: 7.5 - 2.5 * np.log10(i)
m2i = lambda m: 10 ** ((7.5 - m) / 2.5)

fig,ax = plt.subplots(figsize = fs)
m = ax.imshow(arr[2,:,:], cmap = 'BrBG_r')
plt.colorbar(m, ax = ax)
plt.savefig(fn0 + '.1.png', dpi = 144)

fig,ax = plt.subplots(figsize = fs)
m = ax.imshow(i2m(arr[1,:,:] + m2i(22)), vmin = 18, vmax = 22, cmap = 'inferno_r')
plt.colorbar(m, ax = ax)
m = ax.contour(i2m(arr[1,:,:] + m2i(22)), levels = [18, 19, 20, 20.5, 20.75, 21, 21.25, 21.5, 21.75, 22 ], colors = 'white')
plt.clabel(m, fontsize = 8)
plt.savefig(fn0 + '.2.png', dpi = 144)

fig,ax = plt.subplots(figsize = fs)
m = ax.imshow(i2m(arr[0,:,:] + m2i(22)), vmin = 18, vmax = 22, cmap = 'inferno_r')
plt.colorbar(m, ax = ax)
m = ax.contour(i2m(arr[0,:,:] + m2i(22)), levels = [18, 19, 20, 20.5, 20.75, 21, 21.25, 21.5, 21.75, 22 ], colors = 'white')
plt.clabel(m, fontsize = 8)
plt.savefig(fn0 + '.2b.png', dpi = 144)

fig,ax = plt.subplots(figsize = fs)
m = ax.imshow(-2.5 * np.log10((arr[1,:,:] + m2i(22)) / (arr[0,:,:] + m2i(22))), cmap = 'inferno_r', vmin = 0)
plt.colorbar(m, ax = ax)
plt.savefig(fn0 + '.3.png', dpi = 144)
