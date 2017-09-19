from numpy import linspace, logspace
from numpy import sin, cos, log, sqrt, arccos, exp, abs, diff
import matplotlib.pyplot as plt
import latexrc

#-------------------------------------------------------------------------------

x = linspace(-2,2,2**10)

#-------------------------------------------------------------------------------


F1 = lambda x: (exp(x) - 1) / x
F1x = lambda x: exp(x)/x - (exp(x) - 1)/x**2

F2 = lambda x: 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
F2x = lambda x: x**3/30 + x**2/8 + x/3 + 0.5

def FF(x):
    from numpy import where, abs, exp
    F1 = lambda x: (exp(x) - 1) / x
    F2 = lambda x: 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
    return where(abs(x) < 1e-2, F2(x), F1(x))

#-------------------------------------------------------------------------------

fig, axes = plt.subplots(1,3, figsize = (3.5*3,3.5))

#-------------------------------------------------------------------------------

ax = axes[0]
ax.axvline(0, linewidth = 0.3, color = '#C3C3C3')
ax.axhline(1, linewidth = 0.3, color = '#C3C3C3')
ax.plot(x, F2(-x), color = '#2AA0FF', label = 'approx')
ax.plot(x, F1(-x), color = '#002542', linestyle = '--', label = 'orig')
ax.legend()
ax.set_title('$F(x)$')
ax.set_xlim(x.min(), x.max())

#-------------------------------------------------------------------------------

ax = axes[1]
ax.axvline(0, linewidth = 0.3, color = '#C3C3C3')
ax.axhline(0, linewidth = 0.3, color = '#C3C3C3')
ax.plot(x, F2(-x) - F1(-x), color = 'black')
ax.set_title('$F_{\\rm approx}(x) - F(x)$')
ax.set_xlim(x.min(), x.max())

#-------------------------------------------------------------------------------

ax = axes[2]
ax.axvline(0, linewidth = 0.3, color = '#C3C3C3')
ax.axhline(0.5, linewidth = 0.3, color = '#C3C3C3')
ax.plot(x, F2x(-x), color = '#2AA0FF', label = 'approx')
ax.plot(x, F1x(-x), color = '#002542', linestyle = '--', label = 'orig')
ax.set_title('$\\frac{d}{dx}F(x)$')
ax.legend()

ax.set_xlim(x.min(), x.max())

#-------------------------------------------------------------------------------

plt.tight_layout()
plt.savefig('/tmp/fig.png')
# plt.show()

#-------------------------------------------------------------------------------
