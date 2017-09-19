from numpy import linspace, logspace
from numpy import sin, cos, log, sqrt, arccos, exp, abs
import matplotlib.pyplot as plt
import latexrc


def main(R, H, ax):

    #---------------------------------------------------------------------------

    Dmax = R * arccos(R / (R + H))
    D = linspace(0, 1.5 * Dmax, 2**10)
    L = sqrt(H**2 + (1 + H/R) * D**2)

    # exact length from cosine theorem
    Lx = sqrt(H**2 + 2*R*(R+H)*(1 - cos(D/R)))

    costh = (Lx**2 + H**2 + 2*R*H) / (2 * Lx * (R + H))
    cosgam = cos(arccos(costh) + D / R)

    #---------------------------------------------------------------------------

    Q = D**2 / (2 * R * H)
    A = L**2 / (2 * R * H)

    #---------------------------------------------------------------------------

    ax.axvline(Dmax / R, color = '#8DD1CF', linestyle = '-', linewidth = 0.6)
    ax.axhline(0, color = '#8DD1CF', linestyle = '-', linewidth = 0.6)

    ax.plot(D / R, H / L, ':', color = '#656565', label = '$H / L$')

    ax.plot(D / R, H / L * (1 - Q), color = '#40FFAF', label = '$\\frac{H}{L}(1 - Q)$')
    ax.plot(D / R, H / L * (1 - Q*(1+H/R)), color = '#068500', label = '$\\frac{H}{L}\\left[1 - Q \\left(1 + \\frac{H}{R} \\right)\\right]$')
    ax.plot(D / R, cosgam, color = 'black', label = 'exact $\\cos \\gamma$')

    ax.plot(D / R, H / L * (1 + Q), color = '#E9720D', label = '$\\frac{H}{L}(1 + Q)$')
    ax.plot(D / R, costh, '--', color = 'black', label = 'exact $\\cos \\theta$')

    ax.set_xlim(0, 1.5 * Dmax / R)
    ax.set_ylim(None, 1)

    ax.set_xlabel('$D/R$')
    ax.set_ylabel('cosine')
    ax.set_title('$H / R = 1 / {}$'.format(int(round(R/H))))

fig, axes = plt.subplots(1, 2, figsize = (6,3.2))
R = 6371.0
main(R, R / 2, axes[0])
main(R, R / 20, axes[1])
# axes[1].legend()
fig.subplots_adjust(0.10,0.15,0.97,0.91,0.24,0.20)
# plt.tight_layout()
plt.savefig('/tmp/angles.pdf')

#------------------------------------------------------------------------

# fig,ax = plt.subplots()
# ax.plot(D, Q)
# ax.set_xlabel('$D$')
# ax.set_ylabel('$D^2 / (2RH)$')

#------------------------------------------------------------------------

# plt.show()

#------------------------------------------------------------------------
