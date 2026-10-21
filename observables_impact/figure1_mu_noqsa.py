# Figure 1: mu(z) for the four alpha_K parametrizations, without the quasi-static approximation.
# Without the QSA, mu is scale-dependent (Cataneo+2024 eq. 31a); we show it at a few k and the k -> inf
# limit, which coincides with the QSA mu used in the previous version of the figure (cs2_impact_mg.ipynb).
# Run from this folder: python figure1_mu_noqsa.py
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
h = H0/100
ombh2 = omegab*h*h
omch2 = (omegam - omegab)*h*h

# (w0, wa) chosen per parametrization so that no case has a pole in the scale-dependent mu:
# DESI DR2 BAO + CMB + Pantheon+ mean (phantom crossing) where possible; the cubic Galileon-like case has a pole
# for any DESI-like phantom crossing, and the k-essence-like case for lambda_MG < 3, so they use a non-crossing w0, wa
w0wa_desi = (-0.838, -0.62)
w0wa_no_crossing = (-0.9, -0.1)
# Each figure has one row per (cs2_0, alpha_K_0); with the DESI mean, alpha_K propto Omega_DE has a pole for
# cs2 < 0.73 (with alpha_K_0 = 1) and for alpha_K_0 < 0.73 (with cs2 = 1)
rows_cs2 = [(cs2_0, 1, f"$c_s^2 = {cs2_0}$") for cs2_0 in [0.75, 1, 1.25]]
rows_lambda = [(1, alpha_K_0, f"$\\lambda_\\mathrm{{MG}} = {alpha_K_0}$") for alpha_K_0 in [0.75, 1, 1.25]]

# alpha_K_parametrization: 1 = Omega_DE, 2 = k-essence-like, 3 = cubic Galileon-like, 4 = proportional to alpha_B
parametrizations = [
    (1, "$\\alpha_K \\propto \\Omega_\\mathrm{DE}$", w0wa_desi),
    (4, "$\\alpha_K \\propto \\alpha_B$", w0wa_desi),
    (3, "Cubic Galileon-like", w0wa_no_crossing),
    (2, "K-essence-like", w0wa_no_crossing),
]

ks = [1e-4, 3e-4, 1e-3] # h/Mpc

def get_mu(alpha_K_parametrization, cs2_0, alpha_K_0, w0, wa):
    # Returns z, mu_inf (= QSA mu) and mu(k) for each k in ks; same conventions as solve_mu_k_noqsa in data_analysis/utils.py
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, nnu=3.044, mnu=0.06,
        As=2.1e-9, ns=0.96, tau=0.06, w=w0, wa=wa, dark_energy_model="ppf",
        use_cs2=True, use_qsa=False, alpha_K_parametrization=alpha_K_parametrization, cs2_0=cs2_0, alpha_K_0=alpha_K_0
    )
    results = camb.get_background(cosmo)
    a = np.exp(np.array(results.Params.log_a))
    alpha_B = np.array(results.Params.alpha_B)
    c_sN2 = np.array(results.Params.c_sN2)
    mu_p = np.array(results.Params.mu_p)
    z = 1/a - 1
    aH = a*results.hubble_parameter(z)/299792.458 # 1/Mpc

    safe_c_sN2 = np.where(c_sN2 == 0, 1.0, c_sN2)
    r1 = np.where(c_sN2 == 0, 0.0, mu_p/safe_c_sN2)
    mu_inf = np.where(c_sN2 == 0, np.where(alpha_B == 0, 1.0, 1e20), 1.0 + alpha_B**2/(2.0*safe_c_sN2))
    mu_k = [(r1 + (k*h/aH)**2*mu_inf)/(r1 + (k*h/aH)**2) for k in ks]
    return z, mu_inf, mu_k

# Generated from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=4
styles = [
    {"color": "#fdcc8a", "lw": 3, "ls": "-"},
    {"color": "#fc8d59", "lw": 3, "ls": "--"},
    {"color": "#d7301f", "lw": 3, "ls": "-."},
]
style_inf = {"color": "black", "lw": 3, "ls": ":"}

def plot_mu(rows, output_filename):
    fig, axs = plt.subplots(len(rows), 4, figsize=(16, 3.5*len(rows)), sharex=True, sharey="col", gridspec_kw={"wspace": 0.3, "hspace": 0.05})

    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    for row, (cs2_0, alpha_K_0, row_label) in zip(axs, rows):
        for ax, (alpha_K_parametrization, title, (w0, wa)) in zip(row, parametrizations):
            z, mu_inf, mu_k = get_mu(alpha_K_parametrization, cs2_0, alpha_K_0, w0, wa)
            for k, mu, style in zip(ks, mu_k, styles):
                exponent = int(np.floor(np.log10(k)))
                mantissa = k/10**exponent
                k_label = f"$10^{{{exponent}}}$" if np.isclose(mantissa, 1) else f"${mantissa:.0f}\\times 10^{{{exponent}}}$"
                ax.plot(z, mu, label=f"$k = ${k_label}$\\,h/\\mathrm{{Mpc}}$", **style)
            ax.plot(z, mu_inf, label="$k \\to \\infty$ (QSA)",**style_inf)
            ax.set_xlim([0, 3.2])
            ax.set_xticks([0, 1, 2, 3])
            ax.grid()
            ax.tick_params(axis="both", which="major", labelsize=20)
        row[0].set_ylabel("$\\mu$", fontsize=30)
        row[0].text(-0.45, 0.5, row_label, transform=row[0].transAxes, fontsize=25, ha="right", va="center", rotation=90, bbox=props)

    for ax, (alpha_K_parametrization, title, _) in zip(axs[0], parametrizations):
        ax.set_title(title, fontsize=25)
    for ax in axs[-1]:
        ax.set_xlabel("$z$", fontsize=30)

    handles, labels = axs[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=4, fontsize=20, frameon=False)

    plt.savefig(output_filename, bbox_inches="tight")

plot_mu(rows_cs2, "plots/mg_cases_mu_noqsa.pdf")
plot_mu(rows_lambda, "plots/mg_cases_mu_noqsa_lambda.pdf")
