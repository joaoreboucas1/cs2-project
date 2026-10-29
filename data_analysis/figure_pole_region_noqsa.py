# Region of the (w0, wa) plane excluded by the pole of the scale-dependent mu(k, a) without the quasi-static
# approximation (Cataneo+2024 eq. 31a): CAMB-cs2 rejects models with mu_p*c_sN2 < 0 at any time.
# One column per alpha_K parametrization (k-essence-like left out); top row varies cs2 at lambda_MG = 1,
# bottom row varies lambda_MG at cs2 = 1. Gray contours: w0waCDM (PPF) CMB+BAO+SN posterior (MCMC48).
# Run from this folder: python figure_pole_region_noqsa.py (writes plots/pole_region_noqsa.pdf)
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import camb
import utils

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

# Cosmological parameters (same as observables_impact/figure1_mu_noqsa.py)
H0 = 67
omegam = 0.319
omegab = 0.049
h = H0/100
ombh2 = omegab*h*h
omch2 = (omegam - omegab)*h*h

# alpha_K_parametrization: 1 = Omega_DE, 3 = cubic Galileon-like, 4 = proportional to alpha_B
parametrizations = [
    (1, "$\\alpha_K \\propto \\Omega_\\mathrm{DE}$"),
    (4, "$\\alpha_K \\propto \\alpha_B$"),
    (3, "Cubic Galileon-like"),
]
# Each row: list of (cs2_0, alpha_K_0, label)
rows = [
    [(cs2_0, 1, f"$c_s^2 = {cs2_0}$") for cs2_0 in [0.25, 0.5, 1]],
    [(1, alpha_K_0, f"$\\lambda_\\mathrm{{MG}} = {alpha_K_0}$") for alpha_K_0 in [0.5, 1, 2]],
]
row_labels = ["$\\lambda_\\mathrm{MG} = 1$", "$c_s^2 = 1$"]

w0s = np.linspace(-1.4, -0.4, 81)
was = np.linspace(-2.0, 1.0, 81)

# Generated from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=4
styles = [
    {"color": "#fdcc8a", "ls": "-"},
    {"color": "#fc8d59", "ls": "--"},
    {"color": "#d7301f", "ls": "-."},
]

def has_pole(alpha_K_parametrization, cs2_0, alpha_K_0, w0, wa):
    # Returns True if CAMB-cs2 rejects the model because of the pole, NaN if the model is invalid (w0 + wa >= 0)
    if w0 + wa >= 0: return np.nan
    try:
        cosmo = camb.set_params(
            H0=H0, ombh2=ombh2, omch2=omch2, nnu=3.044, mnu=0.06, w=w0, wa=wa, dark_energy_model="ppf",
            use_cs2=True, use_qsa=False, alpha_K_parametrization=alpha_K_parametrization, cs2_0=cs2_0, alpha_K_0=alpha_K_0
        )
        camb.get_background(cosmo)
        return False
    except camb.CAMBError as e:
        if "pole" not in str(e): raise
        return True

# GR posterior contours (68% and 95%), same burn-in as triangle_ds1_subluminal_noqsa.ipynb
gr_chain = utils.load_chain(48, burn_in=0.2, smooth_2D=0.3)
density = gr_chain.get2DDensityGridData("w", "wa", num_plot_contours=2)

fig, axs = plt.subplots(len(rows), len(parametrizations), figsize=(15, 10), sharex=True, sharey=True, gridspec_kw={"wspace": 0.05, "hspace": 0.05})

props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
for row_axs, row, row_label in zip(axs, rows, row_labels):
    for ax, (alpha_K_parametrization, title) in zip(row_axs, parametrizations):
        handles = []
        for (cs2_0, alpha_K_0, label), style in zip(row, styles):
            pole = np.array([[has_pole(alpha_K_parametrization, cs2_0, alpha_K_0, w0, wa) for w0 in w0s] for wa in was], dtype=float)
            print(f"alpha_K_parametrization = {alpha_K_parametrization}, {label}: {np.nanmean(pole)*100:.0f}% of the valid grid excluded")
            ax.contourf(w0s, was, pole, levels=[0.5, 1.5], colors=[style["color"]], alpha=0.25)
            ax.contour(w0s, was, pole, levels=[0.5], colors=[style["color"]], linestyles=[style["ls"]], linewidths=3)
            # With no allowed model there is no boundary to draw, so say it in the legend
            if np.nanmin(pole) == 1: label += " (fully excluded)"
            handles.append((Patch(facecolor=style["color"], alpha=0.5), Line2D([], [], color=style["color"], ls=style["ls"], lw=3), label))
        ax.legend(handles=[(patch, line) for patch, line, _ in handles], labels=[label for _, _, label in handles], fontsize=14, loc="lower left")
        ax.contourf(density.x, density.y, density.P, levels=sorted(density.contours) + [np.inf], colors=["#b0b0b0", "#606060"], alpha=0.8, zorder=3)
        ax.axvline(-1, color="black", lw=1, ls=":")
        ax.plot(w0s, -1 - w0s, color="black", lw=1, ls=":")
        ax.fill_between(w0s, -w0s, was[-1], color="white", hatch="//", edgecolor="gray", lw=0)
        ax.set_xlim(w0s[0], w0s[-1])
        ax.set_ylim(was[0], was[-1])
        ax.grid()
        ax.tick_params(axis="both", which="major", labelsize=20)
    row_axs[0].set_ylabel("$w_a$", fontsize=30)
    row_axs[-1].text(1.05, 0.5, row_label, transform=row_axs[-1].transAxes, fontsize=25, ha="left", va="center", rotation=90, bbox=props)

for ax, (_, title) in zip(axs[0], parametrizations):
    ax.set_title(title, fontsize=25)
for ax in axs[-1]:
    ax.set_xlabel("$w_0$", fontsize=30)
fig.legend(handles=[Patch(color="#606060", alpha=0.8, label="GR, CMB+BAO+SN (68%, 95%)"),
                    Patch(facecolor="white", hatch="//", edgecolor="gray", label="$w_0 + w_a > 0$ (not allowed)")],
           loc="lower center", bbox_to_anchor=(0.5, 0.9), ncol=2, fontsize=20, frameon=False)

fig.savefig("plots/pole_region_noqsa.pdf", bbox_inches="tight")
