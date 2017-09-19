from numpy import linspace, logspace
from numpy import sin, cos, log, sqrt, arccos, exp, abs
import matplotlib.pyplot as plt
import latexrc

def main(R, D, H1, H2, ax):

    #---------------------------------------------------------------------------

    H = linspace(H1, H2, 2**10)
    # Dmax = R * arccos(R / (R + H))
    L = sqrt(H**2 + (1 + H/R) * D**2)

    # exact length from cosine theorem
    Lx = sqrt(H**2 + 2*R*(R+H)*(1 - cos(D/R)))

    costh = (Lx**2 + H**2 + 2*R*H) / (2 * Lx * (R + H))
    cosgam = cos(arccos(costh) + D / R)

    #---------------------------------------------------------------------------

    Q = D**2 / (2 * R * H)
    A = L**2 / (2 * R * H)

    #---------------------------------------------------------------------------

    ax.axvline(0, color = '#C3DEDD', linewidth = 0.6)
    ax.axhline(0, color = '#C3DEDD', linewidth = 0.6)

    ax.plot(H / R, H / L, ':', color = '#656565', label = '$H / L$')

    ax.plot(H / R, H / L * (1 - Q), color = '#2DD9FF', label = '$\\frac{H}{L}(1 - Q)$')
    ax.plot(H / R, H / L * (1 - Q*(1+H/R)), color = '#068500', label = '$\\frac{H}{L}\\left[1 - Q \\left(1 + \\frac{H}{R} \\right)\\right]$')
    ax.plot(H / R, cosgam, color = 'black', label = 'exact $\\cos \\gamma$')

    ax.plot(H / R, H / L * (1 + Q), color = '#E9720D', label = '$\\frac{H}{L}(1 + Q)$')
    ax.plot(H / R, costh, '--', color = 'black', label = 'exact $\\cos \\theta$')

    ax.set_xlim(H1/R, H2/R)

    ax.set_xlabel('$H / R$')
    ax.set_ylabel('cosine')
    ax.set_title('$D / R = {:.2f}$'.format(D/R))

fig, axes = plt.subplots(1, 2, figsize = (6,3.2))
R = 6371.0
main(R, R / 3, -R/3, R/3, axes[0])
main(R, R / 10, -R/3, R/3, axes[1])
# axes[1].legend()
fig.subplots_adjust(0.10,0.15,0.97,0.91,0.24,0.20)
# plt.tight_layout()
plt.savefig('/tmp/angles-2.png')

#------------------------------------------------------------------------

# fig,ax = plt.subplots()
# ax.plot(D, Q)
# ax.set_xlabel('$D$')
# ax.set_ylabel('$D^2 / (2RH)$')

#------------------------------------------------------------------------

# plt.show()

#------------------------------------------------------------------------
