# Table II: constraints (mean and 68% intervals) on Omega_m, H0, S8, w0, wa and mu0 = mu(z = 0, k -> inf),
# from the no-QSA chains, for CMB+BAO+SN and CMBL+BAO+SN+CS.
# Run from this folder: python table2_noqsa.py (writes table2_noqsa.tex)
import numpy as np
import utils

# Columns: w0waCDM, alpha_K propto Omega_DE, alpha_K propto alpha_B, Cubic Galileon-like, K-essence-like, w0waCDM with w >= -1
header = ["$w_0w_a$CDM", "$\\alpha_\\mathrm{K} \\propto \\Omega_\\mathrm{DE}$", "$\\alpha_\\mathrm{K} \\propto \\alpha_\\mathrm{B}$",
          "Cubic Galileon-like", "K-essence-like", "$w_0w_a$CDM $w \\geq -1$"]
# alpha_K_parametrization of each column (None for GR)
aktypes = [None, 1, 4, 3, 2, None]
# Same chains and loading settings as triangle_ds{1,2}_subluminal_noqsa.ipynb
datasets = [
    {"name": "CMB+BAO+SN",     "indices": [48, 113, 115, 116, 114, 83], "load_args": {"burn_in": 0.2, "smooth_2D": 0.3}},
    {"name": "CMBL+BAO+SN+CS", "indices": [76, 117, 119, 120, 118, 84], "load_args": {"burn_in": 0.6, "smooth_2D": 0.35, "smooth_1D": 0.35}},
]
params = [("omegam", "$\\Omega_\\mathrm{m}$"), ("H0", "$H_0$"), ("S8", "$S_8$"), ("w", "$w_0$"), ("wa", "$w_a$")]

def format_value(mean, upper, lower, decimals=None):
    # Round to the precision of the smallest error: 2 significant figures if its leading digits are 10-35, else 1
    if decimals is None:
        err = min(upper, lower)
        exponent = int(np.floor(np.log10(err)))
        leading = int(round(err/10**(exponent - 1)))
        decimals = max(-exponent + (1 if leading <= 35 else 0), 0)
    up, low = f"{upper:.{decimals}f}", f"{lower:.{decimals}f}"
    if up == low:
        return f"${mean:.{decimals}f} \\pm {up}$"
    return f"${mean:.{decimals}f}^{{+{up}}}_{{-{low}}}$"

def param_value(chain, param):
    # Mean and 68% limits of the 1D marginalized posterior (same numbers as getInlineLatex)
    par = chain.getMargeStats().parWithName(param)
    return format_value(par.mean, par.limits[0].upper - par.mean, par.mean - par.limits[0].lower)

def mu0_value(chain, aktype):
    # mu(z = 0) in the k -> inf limit, from samples drawn according to the chain weights
    _, mu_infs, _ = utils.get_mu_k_from_chain(chain, aktype=aktype, k_values=[], num_samples=1000)
    mu0 = mu_infs[:, -1]
    mean = np.mean(mu0)
    return format_value(mean, np.percentile(mu0, 84) - mean, mean - np.percentile(mu0, 16), decimals=2)

lines = ["\\begin{tabular}{|c|" + "c|"*len(header) + "}", "\\hline", " & " + " & ".join(header) + " \\\\ \\hline"]
for dataset in datasets:
    chains = []
    for i in dataset["indices"]:
        chain = utils.load_chain(i, **dataset["load_args"])
        # A few rows of the preliminary no-QSA chains have sigma8 = NaN; drop them
        chain.filter(np.isfinite(chain.samples).all(axis=1))
        chains.append(chain)
    lines.append(f"\\multicolumn{{{len(header)+1}}}{{|c|}}{{{dataset['name']}}} \\\\ \\hline")
    for param, label in params:
        lines.append(f"{label} & " + " & ".join(param_value(chain, param) for chain in chains) + " \\\\ \\hline")
    mu0s = ["$1$" if aktype is None else mu0_value(chain, aktype) for chain, aktype in zip(chains, aktypes)]
    lines.append("$\\mu_0$ & " + " & ".join(mu0s) + " \\\\ \\hline")
lines.append("\\end{tabular}")

table = "\n".join(lines)
print(table)
with open("table2_noqsa.tex", "w") as f:
    f.write(table + "\n")
