from numpy import logspace
from numpy import cos, sin, sqrt, arccos, log, log10, pi
import matplotlib.pyplot as plt
# import latexrc

#-------------------------------------------------------------------------------

R = 1
A = pi * R**2
L = R * logspace(-2,2,2**10)

W1 = A / L**2
W2 = 2 * pi * A / (2 * pi * L**2 + A)
W0 = 2 * pi * (1 - L / sqrt(L**2 +R**2) )

fig,ax = plt.subplots()
ax.set_title('symmetrical cap')
ax.set_xlabel('distance to the center ($L$)')
ax.set_ylabel('solid angle ($\\Omega$)')
ax.plot(L, W1, color = '#DBA80F', label = 'approx. $A / L^2$')
ax.plot(L, W2, color = '#D90D0D', label = 'approx. $2\\pi A / (2\\pi L^2 +A)$')
ax.plot(L, W0, '--', color = 'black', label = 'exact')

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend()

#-------------------------------------------------------------------------------

from numpy.random import uniform, normal, seed
from numpy import array, abs, ndarray, linspace

seed()

n = 2**9
W = ndarray((4,n))

def katbryl(A,B,C):
    from numpy import abs, arctan, dot, pi
    from numpy.linalg import norm, det
    X0 = det([A,B,C])
    X1 = norm(A) * norm(B) * norm(C)
    X2 = dot(A,B) * norm(C)
    X3 = dot(A,C) * norm(B)
    X4 = dot(B,C) * norm(A)
    Y = 2 * arctan(abs(X0) / (X1 + X2 + X3 + X4))
    return Y if Y >= 0 else Y + 2 * pi

t = linspace(0,1,n)
LL = ndarray(n)
COSTH = ndarray(n)

rcirc = 1 / sqrt(pi)
nvert = 24
area_poly = nvert * rcirc**2 * sin(2 * pi / nvert) / 2
rcirc /= sqrt(area_poly)

getarea = lambda VA,VB,VC: abs(VA[0] * (VB[1] - VC[1]) \
      + VB[0] * (VC[1] - VA[1]) \
      + VC[0] * (VA[1] - VB[1])) / 2

for i in range(n):
    xc,yc,zc = 0.0, 2*t[i], 0.2

    L = sqrt(xc**2 + yc**2 + zc**2)
    costh = zc / L

    wangl = 0.0
    area = 0.0

    VC = array([xc,yc,zc])
    for j in range(nvert):
        VA = VC + rcirc * array([cos(2*pi*j / nvert),sin(2*pi*j / nvert),0])
        VB = VC + rcirc * array([cos(2*pi*(j+1) / nvert),sin(2*pi*(j+1) / nvert),0])
        area += getarea(VA,VB,VC)
        wangl += katbryl(VA,VB,VC)

    LL[i] = L
    COSTH[i] = costh

    W[0,i] = wangl
    W[1,i] = costh * area / L**2
    W[2,i] = costh * 2 * pi * area / (2 * pi * L**2 + area)
    W[3,i] = costh * 2 * pi * area / (2 * pi * L**2 + area * costh)


fig,axes = plt.subplots(2,2,figsize=(12,12))

xax = LL

axes[0,0].set_title('exact')
axes[0,0].plot(xax, W[0,:] / pi, color = 'black')

axes[0,1].set_title('$A \\cos\\theta / L^2$')
axes[0,1].plot(xax, W[0,:] / pi, color = '#CFCFCF')
axes[0,1].plot(xax, W[1,:] / pi, color = 'black')

axes[1,0].set_title('$ 2 \\pi A \\cos\\theta / (2 \\pi L^2 + A)$')
axes[1,0].plot(xax, W[0,:] / pi, color = '#CFCFCF')
axes[1,0].plot(xax, W[2,:] / pi, color = 'black')

axes[1,1].set_title('$ 2 \\pi A \\cos\\theta / (2 \\pi L^2 + A \\cos \\theta)$')
axes[1,1].plot(xax, W[0,:] / pi, color = '#CFCFCF')
axes[1,1].plot(xax, W[3,:] / pi, color = 'black')

from numpy import mean
print mean((W[0,:]-W[1,:])**2)
print mean((W[0,:]-W[2,:])**2)
print mean((W[0,:]-W[3,:])**2)

for axe in axes:
    for ax in axe:
        # ax.set_xlim(0,1)
        ax.set_ylim(0.0001,2)
        # ax.set_yscale('log')
        ax.set_xlabel('$L$')
        ax.set_ylabel('$\\Omega$')

# plt.savefig('/tmp/a.png', dpi = 144)

#-------------------------------------------------------------------------------

plt.show()
