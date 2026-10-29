# (w0, wa) posteriors for alpha_K propto Omega_DE with and without the quasi-static approximation, compared to PPF,
# for CMB+BAO+SN (left) and CMBL+BAO+SN+CS (right). 2D version of triangle_ds{1,2}_subluminal_qsa_vs_noqsa_omega.pdf.
# Run from this folder: python figure_w0wa_qsa_vs_noqsa.py (writes plots/w0wa_qsa_vs_noqsa_omega.pdf)
import numpy as np
from getdist import plots
import utils

# (no QSA, QSA, PPF) chains and loading settings, same as triangle_ds{1,2}_subluminal_noqsa.ipynb
datasets = [
    {"name": "CMB+BAO+SN",     "indices": [113, 94, 48], "load_args": {"burn_in": 0.2, "smooth_2D": 0.3}},
    {"name": "CMBL+BAO+SN+CS", "indices": [117, 98, 76], "load_args": {"burn_in": 0.6, "smooth_2D": 0.35, "smooth_1D": 0.35}},
]
label = "$\\alpha_K \\propto \\Omega_\\mathrm{DE}$"
color = utils.colors[0]

p = plots.get_subplot_plotter(subplot_size=5)
p.settings.axes_fontsize = 18
p.settings.axes_labelsize = 25
p.settings.legend_fontsize = 18
p.make_figure(nplot=2, nx=2)

for i, dataset in enumerate(datasets):
    chains = []
    for index in dataset["indices"]:
        chain = utils.load_chain(index, **dataset["load_args"])
        chain.filter(np.isfinite(chain.samples).all(axis=1))
        chains.append(chain)
    p.plot_2d(chains, "w", "wa", ax=(0, i), filled=[False, False, True],
              colors=[color, color, "gray"], lws=3, ls=["-", "--", "-"], lims=[-1.05, -0.55, -1.6, 0.2])
    p.get_axes((0, i)).text(0.95, 0.95, dataset["name"], transform=p.get_axes((0, i)).transAxes, fontsize=20, ha="right", va="top")

p.finish_plot(legend_labels=[f"{label}, no QSA", f"{label}, QSA", "PPF"], legend_ncol=3, no_extra_legend_space=False)
p.fig.savefig("plots/w0wa_qsa_vs_noqsa_omega.pdf", bbox_inches="tight")
