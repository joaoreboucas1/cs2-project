# Figure 3: relative change in the CMB TT and lensing potential power spectra when varying the dark energy sound speed,
# for PPF dark energy in GR (left column) and MG with alpha_K propto Omega_DE without the quasi-static approximation
# (right column). Replaces the CMB plots of cs2_impact_fluid.ipynb and cs2_impact_mg.ipynb.
# Run from this folder: python figure3_cmb_noqsa.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

lmax = 2000
ells = np.arange(2, lmax + 1)

def get_cls(cs2_0, use_cs2):
    # Returns the lensed TT and the lensing potential spectra for ells
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, As=As, ns=ns, tau=tau,
        mnu=0.06, nnu=3.044, w=w0, wa=wa, dark_energy_model="ppf", cs2_0=cs2_0,
        use_cs2=use_cs2, use_qsa=False, alpha_K_parametrization=alpha_K_parametrization, alpha_K_0=alpha_K_0
    )
    cosmo.set_for_lmax(lmax, lens_potential_accuracy=1)
    results = camb.get_results(cosmo)
    cl_tt = results.get_lensed_scalar_cls(CMB_unit="muK", lmax=lmax)[ells, 0]
    cl_pp = results.get_lens_potential_cls(lmax=lmax)[ells, 0]
    return cl_tt, cl_pp

# Both columns are compared to GR with cs2 = 1
cl_tt_base, cl_pp_base = get_cls(1, use_cs2=False)
columns = [
    {
        "title": "GR",
        "cases": [("$c_s^2 = 10^{-1}$", get_cls(1e-1, use_cs2=False)),
                  ("$c_s^2 = 10^{-3}$", get_cls(1e-3, use_cs2=False)),
                  ("$c_s^2 = 10^{-5}$", get_cls(1e-5, use_cs2=False))],
    },
    {
        "title": "MG",
        "cases": [(f"$c_s^2 = {cs2_0}$", get_cls(cs2_0, use_cs2=True)) for cs2_0 in [0.75, 1, 1.25]],
    },
]

# Generated from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=4
styles = [
    {"color": "#fdcc8a", "lw": 2, "ls": "-"},
    {"color": "#fc8d59", "lw": 2, "ls": "--"},
    {"color": "#d7301f", "lw": 2, "ls": "-."},
]

fig, axs = plt.subplots(2, 2, figsize=(15, 8), sharex=True, gridspec_kw={"hspace": 0.05, "wspace": 0.35})

for col, column in zip(axs.T, columns):
    for style, (label, (cl_tt, cl_pp)) in zip(styles, column["cases"]):
        col[0].semilogx(ells, 100*(cl_tt/cl_tt_base - 1), label=label, **style)
        col[1].semilogx(ells, 100*(cl_pp/cl_pp_base - 1), label=label, **style)
    col[0].set_title(column["title"], fontsize=30)
    col[0].set_ylabel("$100 \\times \\Delta C_\\ell^{TT}/C_{\\ell, \\mathrm{GR}}^{TT}$", fontsize=27)
    col[1].set_ylabel("$100 \\times \\Delta C_\\ell^{\\phi\\phi}/C_{\\ell, \\mathrm{GR}}^{\\phi\\phi}$", fontsize=27)
    col[1].set_xlabel("$\\ell$", fontsize=27)
    col[0].legend(fontsize=22, frameon=True, framealpha=1, edgecolor="black")

for ax in axs.flatten():
    ax.tick_params(axis="both", which="major", length=6, labelsize=25)
    ax.tick_params(axis="both", which="minor", length=3)
    ax.grid()
    ax.set_xlim([2, lmax])

plt.savefig("plots/impact_cmb_noqsa.pdf", bbox_inches="tight")
