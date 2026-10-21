import numpy as np
import getdist

# ----- Matplotlib style and constants -----

import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
mpl.rcParams['mathtext.fontset'] = "stix"
mpl.rcParams['font.family'] = "STIXGeneral"

# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

# ----- GetDist helper functions -----

def load_chain(index, burn_in=0.2, smooth_2D=None, smooth_1D=None):
    settings = {}
    if smooth_2D is not None: settings["smooth_scale_2D"] = smooth_2D
    if smooth_1D is not None: settings["smooth_scale_1D"] = smooth_1D

    chain = getdist.loadMCSamples(f"../chains/MCMC{index}/MCMC{index}", settings=settings)
    chain.removeBurn(burn_in)

    if not chain.paramNames.hasParam("omegam"):
        omegam = (chain["omegach2"]+chain["omegabh2"]+(0.06/93.15))/(chain["H0"]/100)**2
        chain.addDerived(omegam, name="omegam", label="\\Omega_m")
    
    if chain.paramNames.hasParam("alpha_K_0"):
        chain.paramNames.parWithName("alpha_K_0").label = r"\lambda_\mathrm{MG}"
    
    if chain.paramNames.hasParam("cs2_a"):
        chain.paramNames.parWithName("cs2_a").label = r"c_{s,a}^2"

    chain.addDerived(chain["sigma8"]*np.sqrt(chain["omegam"]/0.3), name="S8", label="S_8")
    
    return chain

# ---- Preprocessing `CHAINS.md` for reference -----

from collections import namedtuple

Spec = namedtuple("Spec", ("w_model", "cs2_model", "dataset", "extra_info"))

with open("../CHAINS.md", "r") as f:
    chains_table = f.read().splitlines()

# Chain specs are the settings that are varied in this work (dataset, w model, cs2 model)
specs_by_index = {}
for line in chains_table:
    if line.startswith("#"): continue
    entry = list(map(str.strip, line.split("|")))
    if len(entry) == 4:
        index, w_model, cs2_model, dataset = entry
        extra_info = None
    elif len(entry) == 5:
        index, w_model, cs2_model, dataset, extra_info = entry
    else:
        assert False, f"ERROR: in CHAINS.md, entry {entry} has {len(entry)} fields"
    specs_by_index[int(index)] = Spec(**{"w_model": w_model, "cs2_model": cs2_model, "dataset": dataset, "extra_info": extra_info})

# Reverse table

chains_by_spec = {}

for key, value in specs_by_index.items():
    chains_by_spec[value] = key

dataset_descriptions = {
    "DS1": "Planck 2018 TTTEEE+lowl + DESI DR1 BAO + DESY5 SN",
    "DS2": "Planck 2018 TTTEEE+lowl + DESI DR1 BAO + DESY5 SN + DESY3 Shear"
}

# ---- Running CAMB simulations for extracting \mu and \alpha's from chain samples  -----

import camb

def solve_alpha_B_ft(aktype, alpha_K_0, cs2, cs2_a, ombh2, omch2, w0, wa, H0):
    ombh2 = 0.0122
    h = 0.7
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, nnu=3.044, mnu=0.06,
        As=2.1e-9, ns=0.96, tau=0.06, WantTransfer=True, w=w0, wa=wa,
        dark_energy_model="ppf",
        alpha_K_parametrization=aktype, cs2_0=cs2, cs2_a=cs2_a, use_cs2=True, alpha_K_0=alpha_K_0
    )
    results = camb.get_background(cosmo)
    log_a = results.Params.log_a
    alpha_B = results.Params.alpha_B
    alpha_K = results.Params.alpha_K
    mu = results.Params.mu
    return log_a*np.log10(np.e), alpha_B, alpha_K, mu

def get_mu_alphas_from_chain(chain, aktype, num_samples=1000):
    np.random.shuffle(chain.samples)
    
    thin_factor = len(chain.samples)//num_samples
    thin_samples = chain.samples[::thin_factor]
    
    print(f"Number of samples = {len(chain.samples)}")
    print(f"Number of samples after thinning = {len(thin_samples)}")

    pnames = [line.split("\t")[0] for line in chain.getParamNames().__str__().split("\n")]

    alpha_Bs = []
    alpha_Ks = []
    mus      = []
    
    for sample in thin_samples:
        alpha_K_0 = sample[pnames.index("alpha_K_0")] if "alpha_K_0" in pnames else 1.0
        cs2 = sample[pnames.index("cs2_0")]
        cs2_a = sample[pnames.index("cs2_a")] if "cs2_a" in pnames else 0.0
        w0 = sample[pnames.index("w")] if "w" in pnames else -1.0
        wa = sample[pnames.index("w0pwa")] - w0 if "w0pwa" in pnames else 0.0
        H0 = sample[pnames.index("H0*")]
        ombh2 = sample[pnames.index("omegabh2")]
        omch2 = sample[pnames.index("omegach2")]

        log_a, alpha_B, alpha_K, mu = solve_alpha_B_ft(aktype=aktype, alpha_K_0=alpha_K_0, cs2=cs2, cs2_a=cs2_a, ombh2=ombh2, omch2=omch2, w0=w0, wa=wa, H0=H0)
        
        alpha_Bs.append(alpha_B)
        alpha_Ks.append(alpha_K)
        mus.append(mu)
    
    return log_a, np.array(alpha_Bs), np.array(alpha_Ks), np.array(mus)

# ---- Scale-dependent \mu without the quasi-static approximation (Cataneo+2024 eq. 31a) -----

def solve_mu_k_noqsa(aktype, alpha_K_0, cs2, cs2_a, ombh2, omch2, w0, wa, H0, k_values, use_qsa=False):
    # k_values in h/Mpc. Returns mu_inf (the QSA expression, k -> inf limit) and mu(k) for each k
    # use_qsa=True only skips CAMB's pole check, so it can be used on chains run with the QSA
    cosmo = camb.set_params(
        H0=H0, ombh2=ombh2, omch2=omch2, nnu=3.044, mnu=0.06,
        As=2.1e-9, ns=0.96, tau=0.06, WantTransfer=True, w=w0, wa=wa,
        dark_energy_model="ppf",
        alpha_K_parametrization=aktype, cs2_0=cs2, cs2_a=cs2_a, use_cs2=True, alpha_K_0=alpha_K_0, use_qsa=use_qsa
    )
    results = camb.get_background(cosmo)
    log_a = np.array(results.Params.log_a)
    alpha_B = np.array(results.Params.alpha_B)
    c_sN2 = np.array(results.Params.c_sN2)
    mu_p = np.array(results.Params.mu_p)
    a = np.exp(log_a)
    z = 1/a - 1
    aH = a*results.hubble_parameter(z)/299792.458 # 1/Mpc

    # Same conventions as MG_node_ratios in CAMB
    safe_c_sN2 = np.where(c_sN2 == 0, 1.0, c_sN2)
    r1 = np.where(c_sN2 == 0, 0.0, mu_p/safe_c_sN2)
    mu_inf = np.where(c_sN2 == 0, np.where(alpha_B == 0, 1.0, 1e20), 1.0 + alpha_B**2/(2.0*safe_c_sN2))

    mu_k = []
    for k in k_values:
        x = (k*(H0/100)/aH)**2
        mu_k.append((r1 + x*mu_inf)/(r1 + x))
    return log_a*np.log10(np.e), mu_inf, np.array(mu_k)

def get_mu_k_from_chain(chain, aktype, k_values, num_samples=1000, use_qsa=False, seed=0):
    # Draw samples according to the chain weights (without modifying the chain)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(chain.weights), size=num_samples, p=chain.weights/np.sum(chain.weights))
    thin_samples = chain.samples[idx]

    print(f"Number of samples = {len(chain.samples)}")
    print(f"Number of weighted draws = {len(thin_samples)}")

    pnames = [line.split("\t")[0] for line in chain.getParamNames().__str__().split("\n")]

    mu_infs = []
    mu_ks   = []

    for sample in thin_samples:
        alpha_K_0 = sample[pnames.index("alpha_K_0")] if "alpha_K_0" in pnames else 1.0
        cs2 = sample[pnames.index("cs2_0")]
        cs2_a = sample[pnames.index("cs2_a")] if "cs2_a" in pnames else 0.0
        w0 = sample[pnames.index("w")] if "w" in pnames else -1.0
        wa = sample[pnames.index("w0pwa")] - w0 if "w0pwa" in pnames else 0.0
        H0 = sample[pnames.index("H0*")]
        ombh2 = sample[pnames.index("omegabh2")]
        omch2 = sample[pnames.index("omegach2")]

        log_a, mu_inf, mu_k = solve_mu_k_noqsa(aktype=aktype, alpha_K_0=alpha_K_0, cs2=cs2, cs2_a=cs2_a, ombh2=ombh2, omch2=omch2, w0=w0, wa=wa, H0=H0, k_values=k_values, use_qsa=use_qsa)

        mu_infs.append(mu_inf)
        mu_ks.append(mu_k)

    # mu_ks has shape (num_samples, len(k_values), len(log_a))
    return log_a, np.array(mu_infs), np.array(mu_ks)

# Plot distribution of \mu

def plot_mu_distribution(log_a, mus_all, output_filename, plot_title=None):
    fig, axs = plt.subplots(2, 2, sharey=True, sharex=True, figsize=(9, 9), gridspec_kw={"wspace": 0.05, "hspace": 0.12})
    mean_line = Line2D([], [], ls="-", color="gray", label="Mean")
    error_patch = Patch(color="gray", alpha=0.5, label="95%")
    titles = ["$\\alpha_K \\propto \\Omega_\\mathrm{DE}$", "K-essence-like", "$\\alpha_K \\propto \\alpha_B$", "Cubic Galileon-like"]
    z = 10**-log_a - 1

    for i, color in enumerate(colors[:4]):
        row = i//2
        col = i%2
        axs[row, col].plot(z, np.mean(mus_all[i], axis=0), color=color)
        axs[row, col].fill_between(z, np.percentile(mus_all[i], 2.5, axis=0), np.percentile(mus_all[i], 97.5, axis=0), color=color,  alpha=0.5)

    axs[0, 0].set_ylabel("$\\mu$", fontsize=30)
    axs[1, 0].set_ylabel("$\\mu$", fontsize=30)
    for i, (ax, title) in enumerate(zip(axs.flatten(), titles)):
        row = i//2
        col = i%2
        ax.set_xlim(0, 3)
        ax.grid()
        if row == 1: ax.set_xlabel("$z$", fontsize=30)
        ax.set_title(title, fontsize=25)
        ax.tick_params(axis="both", labelsize=25)
        ax.set_xticks([0, 1, 2, 3])

    fig.legend(handles=[mean_line, error_patch], fontsize=20, bbox_to_anchor=(0.9,0.88))
    if plot_title is not None: fig.suptitle(plot_title, fontsize=35, y=1)

    plt.savefig(output_filename, bbox_inches="tight")

# Analysis of best-fit

def read_minimum(chain_index):
    with open(f"../chains/MCMC{chain_index}/MCMC{chain_index}.minimum.txt", "r") as f:
        contents = f.read().splitlines()
        header, values = contents
        keys = header[1:].split()
        values = list(map(float, values.split()))
        minimum = {key: value for key, value in zip(keys, values)}
    return minimum

def get_chi2(minimum):
    return minimum['chi2__planck_2020_hillipop.TTTEEE']+minimum['chi2__planck_2018_lowl.TT']+minimum['chi2__planck_2020_lollipop.lowlEB']+minimum['chi2__bao.desi_dr2']+minimum['chi2__sn.desy5']

def get_cmb(minimum):
    if "alpha_K_0" not in minimum.keys():
        use_cs2 = False
        alpha_K_0 = 1
    else:
        use_cs2 = True
        alpha_K_0 = minimum['alpha_K_0']

    if "cs2_0" not in minimum.keys():
        cs2_0 = 1.0
    else:
        cs2_0 = minimum['cs2_0']
    
    cosmo = camb.set_params(
        H0=minimum['H0'], ombh2=minimum['omegabh2'], omch2=minimum['omegach2'], nnu=3.044, mnu=0.06,
        As=minimum['As'], ns=minimum['ns'], tau=minimum['tau'], WantTransfer=True, w=minimum['w'], wa=minimum['wa'],
        dark_energy_model="ppf",
        alpha_K_parametrization=minimum['aktype'], cs2_0=cs2_0, cs2_a=0, use_cs2=use_cs2, alpha_K_0=alpha_K_0
    )
    results = camb.get_results(cosmo)
    log_a = results.Params.log_a
    alpha_B = results.Params.alpha_B
    alpha_K = results.Params.alpha_K
    mu = results.Params.mu
    cl_tt = results.get_lensed_scalar_cls(CMB_unit='muK')[:,0]
    cl_ee = results.get_lensed_scalar_cls(CMB_unit='muK')[:,1]
    cl_te = results.get_lensed_scalar_cls(CMB_unit='muK')[:,3]
    cl_pp = results.get_lens_potential_cls(lmax=2000)[:,0]
    return log_a*np.log10(np.e), alpha_B, alpha_K, mu, cl_tt, cl_te, cl_ee, cl_pp