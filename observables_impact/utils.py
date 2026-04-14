import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
import camb

def get_observables(case, H0, ombh2, omch2, As, ns, tau, w0, wa, dark_energy_type, alpha_k_parametrization, redshifts, ks):
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, As=As, ns=ns, tau=tau,
        num_nu_massless=3.044, omnuh2=0, WantTransfer=True, w=w0, wa=wa,
        dark_energy_model="fluid" if dark_energy_type == "wcdm" else "ppf",
        use_cs2=case["use_cs2"], alpha_K_parametrization=1, cs2_0=case["cs2_0"], alpha_K_0=case["alpha_K_0"]
    )
    cosmo.set_for_lmax(2000, lens_potential_accuracy=1)
    cosmo.set_matter_power(redshifts=redshifts, kmax=20, silent=True, nonlinear=None)
    results = camb.get_results(cosmo)
    P_lin = results.get_matter_power_interpolator(nonlinear=None).P(redshifts, ks)
    cl_tt = results.get_lensed_scalar_cls(CMB_unit='muK')[:,0]
    cl_ee = results.get_lensed_scalar_cls(CMB_unit='muK')[:,1]
    cl_te = results.get_lensed_scalar_cls(CMB_unit='muK')[:,3]
    cl_pp = results.get_lens_potential_cls(lmax=2000)[:,0]
    log_a = results.Params.log_a
    mu = results.Params.mu
    case.update({"P_lin": P_lin, "cl_tt": cl_tt, "cl_ee": cl_ee, "cl_te": cl_te, "cl_pp": cl_pp, "log_a": log_a, "mu": mu})

def plot_pk_differences(cases, base_case, ks, redshifts, plot_name):
    # Precomputing min and max relative errors for unifying the colormaps in each axis
    fig, axs = plt.subplots(len(cases), 1, sharex=True, sharey=True, figsize=(10, 7), gridspec_kw={"hspace": 0.11})
    scale = -np.inf
    for name, case in cases.items():
        relative_errors = 100*np.abs(case["P_lin"]/base_case["P_lin"] - 1)
        scale = max(scale, np.amax(relative_errors))

    norm = Normalize(vmin=-scale, vmax=scale)

    for i, case in enumerate(cases.values()):
        ax = axs[i]
        relative_errors = 100*(case["P_lin"]/base_case["P_lin"] - 1)
        mesh = ax.pcolormesh(ks, redshifts, relative_errors, cmap="seismic", norm=norm)
        ax.set_xscale('log')
        ax.set_ylabel(r"$z$", fontsize=20)
        ax.text(0.83, 0.77, case["label"], transform=ax.transAxes, bbox={'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.5}, fontsize=20)
        ax.set_yticks([0, 1, 2, 3])

    cbar = fig.colorbar(mesh, ax=axs, orientation="vertical", fraction=0.05, pad=0.02)
    cbar.set_label(r"$100 \times \Delta P(k)/P(k)_\mathrm{GR}$", fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    for ax in axs: ax.tick_params(labelsize=17)
    axs[-1].set_xlabel(r"$k \; (h/\mathrm{Mpc})$", fontsize=20)
    plt.savefig(plot_name, bbox_inches="tight")

def plot_cmb_differences(cases, base_case, plot_name):
    fig, axs = plt.subplots(2, 2, figsize=(15, 8), gridspec_kw={'wspace': 0.2, 'hspace': 0.2})
    ells_tt = np.arange(2051)
    ells_pp = np.arange(2001)

    # Generated from https://colorbrewer2.org/#type=sequential&scheme=OrRd&n=4
    styles = [
        {"lw": 2, "color": "#1b9e77"},
        {"lw": 2, "color": "#d95f02"},
        {"lw": 2, "color": "#7570b3"},
    ]

    for style, (name, case) in zip(styles, cases.items()):
        if name == "quint": continue
        axs[0, 0].semilogx(ells_tt, 100*(case["cl_tt"]/base_case["cl_tt"] - 1), label=case["label"], **style)
        axs[0, 1].semilogx(ells_pp, 100*(case["cl_pp"]/base_case["cl_pp"] - 1), label=case["label"], **style)
        axs[1, 0].semilogx(ells_tt, 100*(case["cl_ee"]/base_case["cl_ee"] - 1), label=case["label"], **style)
        axs[1, 1].semilogx(ells_tt, case["cl_te"] - base_case["cl_te"],   label=case["label"], **style)

    axs[0, 0].set_xlabel(r"$\ell$", fontsize=15)
    axs[0, 0].set_ylabel(r"$100 \times \Delta C_\ell^{TT}/C_{\ell, c_s^2=1}^{TT}$", fontsize=15)
    axs[0, 1].set_xlabel(r"$\ell$", fontsize=15)
    axs[0, 1].set_ylabel(r"$100 \times \Delta C_\ell^{\phi\phi}/C_{\ell, c_s^2=1}^{\phi\phi}$", fontsize=15)
    axs[1, 0].set_xlabel(r"$\ell$", fontsize=15)
    axs[1, 0].set_ylabel(r"$100 \times \Delta C_\ell^{EE}/C_{\ell, c_s^2=1}^{EE}$", fontsize=15)
    axs[1, 1].set_xlabel(r"$\ell$", fontsize=15)
    axs[1, 1].set_ylabel(r"$\Delta C_\ell^{TE}$", fontsize=15)
    axs[0, 0].legend(fontsize=15, frameon=True, framealpha=1, edgecolor="black")
    for ax in axs.flatten():
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.grid()
        ax.set_xlim([2, 2000])
    plt.savefig(plot_name, bbox_inches="tight")