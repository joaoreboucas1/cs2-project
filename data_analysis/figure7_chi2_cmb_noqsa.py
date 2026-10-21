# Figure 7: distribution of chi2_CMB in the CMB+BAO+SN chains for GR (PPF) and each MG parametrization of alpha_K
# without the quasi-static approximation. Replaces the plot of chi2_cmb_comparison.ipynb (which used the QSA chain).
# The paper shows alpha_K propto Omega_DE (plots/chi2_dist_comparison_noqsa.pdf); the other models are saved with a suffix.
# Run from this folder: python figure7_chi2_cmb_noqsa.py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import utils

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"

# (MG chain, GR chain, label, color, output suffix); k-essence-like is compared to w0waCDM with w >= -1, as in Table II
models = [
    (113, 48, "$\\alpha_K \\propto \\Omega_\\mathrm{DE}$", utils.colors[0], ""),
    (115, 48, "$\\alpha_K \\propto \\alpha_B$",            utils.colors[2], "_prop"),
    (116, 48, "Cubic Galileon-like",                       utils.colors[3], "_cugal"),
    (114, 83, "K-essence-like",                            utils.colors[1], "_kessence"),
]

# Same burn-in as triangle_ds1_subluminal_noqsa.ipynb and table2_noqsa.py
chains = {i: utils.load_chain(i, burn_in=0.2) for i in {i for model in models for i in model[:2]}}

def chi2_cmb(chain):
    return chain["chi2__planck_2020_hillipop.TTTEEE"] + chain["chi2__planck_2018_lowl.TT"] + chain["chi2__planck_2020_lollipop.lowlEB"]

chi2_range = [30623, 30662]

for i_mg, i_gr, label, color, suffix in models:
    cases = [
        (chains[i_gr], {"facecolor": "gray", "alpha": 0.8, "hatch": "//",   "label": "GR"}),
        (chains[i_mg], {"facecolor": color,  "alpha": 0.6, "hatch": "\\\\", "label": label}),
    ]
    fig, ax = plt.subplots(figsize=(10, 5))

    for chain, style in cases:
        chi2 = chi2_cmb(chain)
        mask = np.isfinite(chi2)
        counts, bins = np.histogram(chi2[mask], bins=50, range=chi2_range, weights=chain.weights[mask], density=True)
        # Step-wise coordinates for a single filled polygon
        x = np.repeat(bins, 2)[1:-1]
        y = np.repeat(counts, 2)
        ax.fill_between(x, y, edgecolor="black", linewidth=0.8, **style)
        ax.plot(x, y, color="black", linewidth=0.8)
        # Weighted median (the mean is dominated by a few samples with very large chi2)
        order = np.argsort(chi2[mask])
        cdf = np.cumsum(chain.weights[mask][order])/np.sum(chain.weights[mask])
        median = chi2[mask][order][np.searchsorted(cdf, 0.5)]
        print(f"MCMC{i_mg} plot, {style['label']}: weighted median chi2_CMB = {median:.2f}, minimum = {np.min(chi2[mask]):.2f}")

    patches = [Patch(facecolor=style["facecolor"], edgecolor="black", hatch=style["hatch"], label=style["label"]) for _, style in cases]
    ax.set_xlabel("$\\chi^2_\\mathrm{CMB}$", fontsize=25)
    ax.set_ylabel("Sample density", fontsize=25)
    ax.tick_params(axis="both", which="major", labelsize=22)
    ax.set_xlim(chi2_range)
    ax.set_ylim(bottom=0)
    ax.legend(handles=patches, fontsize=22, loc="upper left")
    fig.tight_layout()
    fig.savefig(f"plots/chi2_dist_comparison_noqsa{suffix}.pdf", bbox_inches="tight")
