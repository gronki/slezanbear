#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------------------------

from numpy import exp, sqrt, linspace, logspace, pi, log
from numpy.random import uniform, normal
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#-------------------------------------------------------------------------------

B = lambda x,mu,sig: exp(-(log(x) - mu)**2 / (2*sig**2)) / (sig * x * sqrt(2*pi))
B = lambda x,mu,sig: exp(-(x - mu)**2 / (2*sig**2)) / (sig * sqrt(2*pi))
B1 = lambda x,A1,A2: 2*x**3/(exp(x/A2) - 1)
B2 = lambda x,A1,A2: 2*A1*x**4*exp(x/A2)/(A2**2*(exp(x/A2) - 1)**2)
#-------------------------------------------------------------------------------

x = logspace(-2,2,200)
plt.plot(x,B(x,log(0.5),1.0))
plt.plot(x,B(x,log(3.0),1.0))
plt.plot(x,B(x,log(3.0),2.0))
# plt.yscale('log')

#-------------------------------------------------------------------------------
plt.show()
