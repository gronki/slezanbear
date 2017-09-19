from numpy import linspace, logspace, meshgrid, where, arange
from numpy import sin, cos, log, sqrt, arccos, exp, abs, pi, arcsinh, log10
import matplotlib.pyplot as plt
# import latexrc

#-------------------------------------------------------------------------------

dist_max = 20.0
x0 = linspace(-dist_max * 1.5, dist_max * 1.5, 401)
h0 = linspace(0.01, 30, 500)
x,H = meshgrid(x0, h0, indexing = 'xy')
A0 = (x0[1] - x0[0])**2

#-------------------------------------------------------------------------------

def FF(x):
    from numpy import where, abs, exp
    F1 = lambda x: (exp(x) - 1) / x
    F2 = lambda x: 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
    return where(abs(x) < 1e-2, F2(x), F1(x))

#-------------------------------------------------------------------------------

Re = 6371.0
Hscale = 5.0

albedo = 0.3

from threshold import th
cutoff = lambda cosgam,cut: th(cosgam / cos(pi * (90 - cut) / 180) - 1, 0.5)

E = lambda cosgam: 1.0
cosgam_max = cos( pi / 180 * (90 - 20) )
E = lambda cosgam: cosgam + th(cosgam_max - cosgam, 0.05) / albedo

g = lambda costh: 1.0
# g = lambda costh: 0.75 * (1 + costh**2)


#-------------------------------------------------------------------------------

def kern(sig0, chi0, I0, A0, R, D, H, Hobs):
    L = sqrt(H**2 + (1 + H / R) * D**2)
    Q = D**2 / (2*R*H)
    costh = H / L * (1 + Q)
    cosgam = H / L * (1 - Q * (1 + H / R))
    tau = chi0 * L * FF(-H / Hscale) \
        + chi0 * Hscale * exp(-H / Hscale)  * (exp((H-Hobs) / Hscale) - 1)
    return 0.5 * I0 * exp(-tau) * sig0 * exp(-H/Hscale) \
        * g(costh) * E(cosgam) * cutoff(cosgam, 2) \
        * A0 / (2 * pi * L**2 + A0)

def tau(chi0, R, D, H, Hobs):
    L = sqrt(H**2 + (1 + H/R) * D**2)
    tau1 = chi0 * L * FF(-H / Hscale)
    tau2 = chi0 * Hscale * (exp(-Hobs / Hscale) - exp(-H / Hscale))
    return tau1, tau2

#-------------------------------------------------------------------------------

sig = 0.15 / Hscale
chi = sig / 2
I = kern(sig, chi, 1.0 / A0, A0, Re, abs(x), H, 0.0)
tau1, tau2 = tau(chi, Re, abs(x), H, 0.0)
Ifl = kern(sig, chi, 1.0 / A0, A0, Re * 1e8, abs(x), H, 0.0)
Ifl_noa = kern(sig, 0.0, 1.0 / A0, A0, Re * 1e8, abs(x), H, 0.0)

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(2, 1)

ax = axes[0]
cn = arange(-24,-1)
ax.pcolor(  x0, h0[::2], log10(I[::2,:]),
            vmin = -8,
            vmax = -2,
            cmap = 'jet')
cs = ax.contour( x0, h0, tau1,
    levels = logspace(-2,1,13),
    linewidths = 0.6, colors = 'black')
ax.clabel(cs, fontsize = 6)
ax.set_aspect(1)

ax = axes[1]
ax.plot(x0, sum(Ifl,0), color = '#9ECBC4')
# ax.plot(x0, sum(Ifl,0), color = '#03A890')

ax.plot(x0, sum(I,0), linestyle = '--', color = '#141414')
ax.axhline(sig*Hscale/2.7, color = '#BEBEBE', linewidth = 0.6)
ax.set_yscale('log')
ax.set_xlim(x0.min(), x0.max())

#-------------------------------------------------------------------------------

plt.savefig('/tmp/model1.png', dpi = 192)
# plt.show()

#-------------------------------------------------------------------------------
