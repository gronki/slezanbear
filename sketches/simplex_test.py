#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

from numpy import sqrt, exp, log, log10, sin, cos
from numpy.random import uniform, normal
a0 = uniform(-5, 5)
b0 = uniform(-5, 5)
xd = uniform(-10, 10, 120)
yd = a0 * xd + b0 + normal(0, 0.1, xd.shape[0])
def f(a,b):
    s = a * 0
    for i in range(xd.shape[0]):
        s += (a * xd[i] + b - yd[i])**2
    return sqrt(s / xd.shape[0])
    # return (x - 1)**2 + (y - 2)**2 + (x + y)**3 / 60

#-------------------------------------------------------------------------------

from numpy import linspace, meshgrid
x = linspace(-8,8,128)
y = linspace(-8,8,128)

xx, yy = meshgrid(x, y)

import matplotlib.pyplot as plt
plt.contourf(x, y, f(xx,yy), 16, cmap = 'PuBu_r')

#-------------------------------------------------------------------------------

from matplotlib.cm import YlOrRd_r
from numpy import average, array, argsort, argmin
x = uniform(-5, 5, (2,3))

plt.plot(
    [ x[0,0], x[0,1], x[0,2], x[0,0] ],
    [ x[1,0], x[1,1], x[1,2], x[1,0] ],
    color = YlOrRd_r(0.0),
)

a = [1.0, 2.0, 0.5, 0.5]

niter = 32
for i in range(niter):
    y = f(x[0,:], x[1,:])
    s = argsort(y)

    x_0 = array([average(x[0,s[0:2]]), average(x[1,s[0:2]])])

    # reflection
    x_r = x_0 + a[0] * (x_0 - x[:,s[2]])
    y_r = f(x_r[0], x_r[1])
    if y_r < y[s[1]]:
        if y_r > y[s[0]]:
            x[:,s[2]] = x_r
        else:
            # expansion
            x_e = x_0 + a[1] * (x_r - x_0)
            y_e = f(x_e[0], x_e[1])
            if y_e < y_r:
                x[:,s[2]] = x_e
            else:
                x[:,s[2]] = x_r
    else:
        # contraction
        x_c = x_0 + a[2] * (x[:,s[2]] - x_0)
        y_c = f(x_c[0], x_c[1])
        if y_c < y[s[2]]:
            x[:,s[2]] = x_c
        else:
            # shrink
            x[:,s[1]] = x[:,s[0]] + a[3] * (x[:,s[1]] - x[:,s[0]])
            x[:,s[2]] = x[:,s[0]] + a[3] * (x[:,s[2]] - x[:,s[0]])

    print f(x[0,:], x[1,:])
    plt.plot(
        [ x[0,0], x[0,1], x[0,2], x[0,0] ],
        [ x[1,0], x[1,1], x[1,2], x[1,0] ],
        color = YlOrRd_r((i + 1.0) / (niter)),
    )

print a0, b0
ibest = argmin(f(x[0,:], x[1,:]))
print x[0,ibest], x[1,ibest]
#-------------------------------------------------------------------------------
plt.show()
#-------------------------------------------------------------------------------
