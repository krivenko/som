from h5 import HDFArchive
from triqs.gf import *
from matplotlib import pyplot as plt
import numpy as np

# Read data from archive
ar = HDFArchive('results.h5', 'r')

w_mesh = ar["0"]["w_mesh"]
w_points = np.array([float(w) for w in w_mesh])

#
# Plot averaged spectra with error bars derived from estimated dispersions
#

fig, axes = plt.subplots(3, 2, figsize=(8, 10))
for i in (0, 1):
    gr = ar[str(i)]

    # Rectangular resolution function
    axes[0, i].errorbar(w_points,
                        gr["avg_rect"],
                        xerr=w_mesh.delta,
                        yerr=np.sqrt(gr["disp_rect"]),
                        lw=0.1)
    axes[0, i].set_title(r"$A_%d(\omega)$, Rectangular" % i)

    # Lorentzian resolution function
    axes[1, i].errorbar(w_points,
                        gr["avg_lorentz"],
                        xerr=w_mesh.delta,
                        yerr=np.sqrt(gr["disp_lorentz"]),
                        lw=0.1)
    axes[1, i].set_title(r"$A_%d(\omega)$, Lorentz" % i)

    # Gaussian resolution function
    axes[2, i].errorbar(w_points,
                        gr["avg_gauss"],
                        xerr=w_mesh.delta,
                        yerr=np.sqrt(gr["disp_gauss"]),
                        lw=0.1)
    axes[2, i].set_title(r"$A_%d(\omega)$, Gaussian" % i)

plt.setp(axes, xlim=(w_points[0], w_points[-1]), xlabel=r"$\omega$")
plt.tight_layout()
plt.show()

#
# Plot two-point correlators
#

fig, axes = plt.subplots(3, 2, figsize=(8, 10))
for i in (0, 1):
    gr = ar[str(i)]

    extent = [w_points[0], w_points[-1], w_points[-1], w_points[0]]

    # Rectangular resolution function
    corr = axes[0, i].imshow(gr["corr_rect"], extent=extent, cmap="hot")
    fig.colorbar(corr, ax=axes[0, i])
    axes[0, i].set_title(r"$\sigma_{\omega\omega'}^{(%d)}$, Rectangular" % i)

    # Lorentzian resolution function
    corr = axes[1, i].imshow(gr["corr_lorentz"], extent=extent, cmap="hot")
    fig.colorbar(corr, ax=axes[1, i])
    axes[1, i].set_title(r"$\sigma_{\omega\omega'}^{(%d)}$, Lorentz" % i)

    # Gaussian resolution function
    corr = axes[2, i].imshow(gr["corr_gauss"], extent=extent, cmap="hot")
    fig.colorbar(corr, ax=axes[2, i])
    axes[2, i].set_title(r"$\sigma_{\omega\omega'}^{(%d)}$, Gaussian" % i)

plt.tight_layout()
plt.show()
