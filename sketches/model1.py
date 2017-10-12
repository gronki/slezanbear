from numpy import linspace, logspace, meshgrid, where, arange
from numpy import sin, cos, log, sqrt, arccos, exp, abs, pi, arcsinh, log10
import matplotlib.pyplot as plt
# import latexrc

#-------------------------------------------------------------------------------

dist_max = 40.0
x0 = linspace(-dist_max, dist_max, 2 * dist_max / 0.2)
h0 = linspace(0.01, min(dist_max / 1.5, 40), 500)
x,H = meshgrid(x0, h0, indexing = 'xy')
A0 = 1.0 #(x0[1] - x0[0])**2

#-------------------------------------------------------------------------------


Re = 6371.0
Hscale = 5.0

albedo = 0.2

from threshold import th
cutoff = lambda cosgam,cut: th(cosgam / cos(pi * (90 - cut) / 180) - 1, 0.5)

E = lambda cosgam: 1.0
cosgam_max = cos( pi / 180 * (90 - 15) )
E = lambda cosgam: cosgam + th(cosgam_max - cosgam, 0.05) / albedo

g = lambda costh: 1.0
g = lambda costh: 0.75 * (1 + costh**2)


def FF(x):
    from numpy import where, abs, exp
    F1 = lambda x: (exp(x) - 1) / x
    F2 = lambda x: 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
    return where(abs(x) < 1e-2, F2(x), F1(x))

#-------------------------------------------------------------------------------

def kern(alpha, rabs, Hscale, I0, A0, R, D, H, Hobs):
    L = sqrt(H**2 + (1 + H / R) * D**2)
    Q = D**2 / (2*R*H)
    costh = H / L * (1 + Q)
    cosgam = H / L * (1 - Q * (1 + H / R))
    sct0 = alpha / Hscale
    chi0 = alpha / Hscale * rabs
    tau1 = chi0 * L * FF(-H / Hscale)
    tau2 = chi0 * Hscale * exp(-H / Hscale)  * (exp((H-Hobs) / Hscale) - 1)
    J0 = 0.5 * I0 * E(cosgam) * cutoff(cosgam, 0.1) \
        * A0 / (2 * pi * L**2 + A0) * exp(-tau1)
    J1 = sct0 * exp(-H/Hscale) * g(costh)
    J2 = exp(-tau2)
    return J0 * J1 * J2


def tau(alpha, rabs, Hscale, R, D, H, Hobs):
    L = sqrt(H**2 + (1 + H/R) * D**2)
    chi0 = alpha / Hscale * rabs
    tau1 = chi0 * L * FF(-H / Hscale)
    tau2 = chi0 * Hscale * (exp(-Hobs / Hscale) - exp(-H / Hscale))
    return tau1, tau2


#-------------------------------------------------------------------------------

def berry(s,b,k,H,D):
    L = sqrt(H**2 + D**2)
    return s * (H + b*L) / L**2 * exp(-k*L)

#-------------------------------------------------------------------------------

alpha = 0.2
rabs = 0.75
I = kern(alpha, rabs, Hscale, 1.0 / A0, A0, Re, abs(x), H, 0.0)
tau1, tau2 = tau(alpha, rabs, Hscale, Re, abs(x), H, 0.0)
Ifl = kern(alpha, rabs, Hscale, 1.0 / A0, A0, Re * 1e8, abs(x), H, 0.0)
Ifl_noa = kern(alpha, 0.0, Hscale, 1.0 / A0, A0, Re * 1e8, abs(x), H, 0.0)

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(2, 1)

ax = axes[0]
cn = arange(-24,-1)
ax.pcolor(  x0, h0[::2], log10(I[::2,:]),
            vmin = log10(alpha) - 6,
            vmax = log10(alpha),
            cmap = 'jet')
cs = ax.contour( x0, h0, tau1,
    levels = logspace(-2,1,13),
    linewidths = 0.6, colors = 'black')
ax.clabel(cs, fontsize = 6)
ax.set_aspect(1)

integrate = lambda y,x: ((y[1:] + y[:-1])) * (x[1:] - x[:-1]) / 2
dh = h0[1] - h0[0]

ax = axes[1]
ax.plot(x0, sum(Ifl,0)*dh, color = '#03A890')
ax.plot(x0, sum(Ifl_noa,0)*dh, color = '#9ECBC4')

# ax.plot(x0, berry(alpha, alpha * chi / Hscale, Hscale / 3, abs(x0)), color = '#FF6E11')

ax.plot(x0, sum(I,0)*dh, linestyle = '--', color = '#141414')
# ax.axhline(alpha, color = '#BEBEBE', linewidth = 0.6)
ax.set_yscale('log')
ax.set_xlim(x0.min(), x0.max())

#-------------------------------------------------------------------------------

plt.savefig('/tmp/model1.png', dpi = 192)
# plt.show()

#-------------------------------------------------------------------------------
