# Figure 2: relative change in the linear matter power spectrum P(k, z) when varying the dark energy sound speed,
# for PPF dark energy in GR (left column) and MG with alpha_K propto Omega_DE without the quasi-static approximation
# (right column). Replaces the P(k) plots of cs2_impact_fluid.ipynb and cs2_impact_mg.ipynb.
# Run from this folder: python figure2_pk_noqsa.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import camb

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

# Cosmological parameters
H0 = 67
omegam = 0.319
omegab = 0.049
As = 2.1e-9
ns = 0.96
tau = 0.07
h = H0/100
ombh2 = omegab*h*h
omch2 = (omegam - omegab)*h*h

# DESI DR2 BAO + CMB + Pantheon+ mean
w0 = -0.838
wa = -0.62

# MG: alpha_K propto Omega_DE with lambda_MG = 1; this case has a pole in the scale-dependent mu for cs2 < 0.73
alpha_K_parametrization = 1
alpha_K_0 = 1

# k, z grid for power spectrum
ks = np.logspace(-4, -0.9, 256) # h/Mpc
redshifts = np.linspace(0, 3, 256)

def get_pk(cs2_0, use_cs2):
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, As=As, ns=ns, tau=tau,
        mnu=0.06, nnu=3.044, WantTransfer=True, WantCls=False, w=w0, wa=wa, dark_energy_model="ppf", cs2_0=cs2_0,
        use_cs2=use_cs2, use_qsa=False, alpha_K_parametrization=alpha_K_parametrization, alpha_K_0=alpha_K_0
    )
    cosmo.set_matter_power(redshifts=redshifts, kmax=20, silent=True, nonlinear=None)
    results = camb.get_results(cosmo)
    return results.get_matter_power_interpolator(nonlinear=None).P(redshifts, ks)

# Both columns are compared to GR with cs2 = 1
pk_base = get_pk(1, use_cs2=False)
columns = [
    {
        "title": "GR",
        "cbar_label": "$100 \\times \\Delta P(k)/P(k)_{c_s^2 = 1}$",
        "cases": [("$c_s^2 = 10^{-1}$", get_pk(1e-1, use_cs2=False)),
                  ("$c_s^2 = 10^{-3}$", get_pk(1e-3, use_cs2=False)),
                  ("$c_s^2 = 10^{-5}$", get_pk(1e-5, use_cs2=False))],
    },
    {
        "title": "MG",
        "cbar_label": "$100 \\times \\Delta P(k)/P(k)_\\mathrm{GR}$",
        "cases": [(f"$c_s^2 = {cs2_0}$", get_pk(cs2_0, use_cs2=True)) for cs2_0 in [0.75, 1, 1.25]],
    },
]

fig, axs = plt.subplots(3, 2, figsize=(11, 7), sharex=True, sharey=True, gridspec_kw={"hspace": 0.11, "wspace": 0.5})

for col, column in zip(axs.T, columns):
    # Same color scale for all panels in a column, symmetric around zero
    relative_errors = [100*(pk/pk_base - 1) for _, pk in column["cases"]]
    scale = max(np.amax(np.abs(rel)) for rel in relative_errors)
    norm = Normalize(vmin=-scale, vmax=scale)
    for ax, (label, _), rel in zip(col, column["cases"], relative_errors):
        mesh = ax.pcolormesh(ks, redshifts, rel, cmap="seismic", norm=norm)
        ax.set_xscale("log")
        ax.text(0.97, 0.77, label, transform=ax.transAxes, bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5}, fontsize=20, ha="right")
        ax.set_yticks([0, 1, 2, 3])
        ax.tick_params(labelsize=17, labelleft=True)
        ax.set_ylabel("$z$", fontsize=20)
    col[0].set_title(column["title"], fontsize=30)
    col[-1].set_xlabel("$k \\; (h/\\mathrm{Mpc})$", fontsize=20)
    cbar = fig.colorbar(mesh, ax=col, orientation="vertical", fraction=0.05, pad=0.02)
    cbar.set_label(column["cbar_label"], fontsize=20)
    cbar.ax.tick_params(labelsize=18)

plt.savefig("plots/impact_pk_noqsa.png", bbox_inches="tight", dpi=200)
