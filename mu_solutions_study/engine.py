"""
    Engine for sound speed modified gravity models
    Author: João Rebouças, December 2025
"""

from enum import IntEnum
import numpy as np

omega_r = 2.5e-5*(1 + 3.044 * 7/8 * (4/11)**(4/3))

class Bg():
    def __init__(self, omega_m, w0, wa):
        self.omega_m = omega_m
        self.w0 = w0
        self.wa = wa
        self.omega_de_0 = 1 - omega_m - omega_r

def rho_de(a, bg):
    return bg.omega_de_0*a**(-3*(1 + bg.w0 + bg.wa))*np.exp(-3*bg.wa*(1-a))

def rho_m(a, bg):
    return bg.omega_m*a**-3

def rho_gamma(a, bg):
    return omega_r*a**-4

def get_bg_funcs(a, bg):
    wde    = bg.w0 + bg.wa*(1-a)
    rhog   = rho_gamma(a, bg)
    rhom   = rho_m(a, bg)
    rhode  = rho_de(a, bg)
    rhotot = rhog + rhode + rhom
    wtot   = (rhog/3 + wde*rhode)/rhotot
    return wde, rhode, wtot, rhotot

class alphaKtype(IntEnum):
    CONST = 1
    OMEGA = 2
    QUINT = 3
    CUGAL = 4
    PROP  = 5
    CUGAL_MOCHI = 6
    DKIN_JOAO = 7

def get_alpha_K(aktype, a, bg, alpha_B, cs2, alpha_K_0):
    match aktype:
        case alphaKtype.CONST: 
            return alpha_K_0
        case alphaKtype.OMEGA:
            rhode  = rho_de(a, bg)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            omega_de = rhode/rhotot
            return alpha_K_0*omega_de/bg.omega_de_0
        case alphaKtype.QUINT:
            rhode  = rho_de(a, bg)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            omega_de = rhode/rhotot
            w_de = bg.w0 + bg.wa*(1 - a)
            return alpha_K_0*omega_de*(1 + w_de)/cs2
        case alphaKtype.CUGAL:
            # NOTE: E^2 = (H/H0)^2 = rhotot/rho_cr_0
            # Our rhotot is already in units of the current critical density
            rhode  = rho_de(a, bg)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            omega_de = rhode/rhotot
            return alpha_K_0*omega_de/rhotot**2 + 6*alpha_B
        case alphaKtype.PROP:
            return alpha_K_0*3*alpha_B
        case alphaKtype.CUGAL_MOCHI:
            # NOTE: E^2 = (H/H0)^2 = rhotot/rho_cr_0
            rhode  = rho_de(a, bg)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            omega_de = rhode/rhotot
            return alpha_K_0*6*omega_de/rhotot**2
        case alphaKtype.DKIN_JOAO:
            # D_kin = lambda*d\ln(H)/d\ln(a) = lambda*(-3*(1 + w_tot)/2)
            rhode  = rho_de(a, bg)
            wde    = bg.w0 + bg.wa*(1 - a)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            Ptot   = rho_gamma(a, bg)/3 + wde*rhode
            w_tot  = Ptot/rhotot
            D_kin  = alpha_K_0*(-1.5)*(1 + w_tot)
            return D_kin - 1.5*alpha_B**2
        case _:
            raise Exception("Unknown alpha_K parametrization")

def deriv(loga, alpha_B, bg, aktype, alpha_K_0, cs2):
    a = np.exp(loga)
    wde, rhode, wtot, rhotot = get_bg_funcs(a, bg)
    d_lnH_d_lna = -1.5*(1 + wtot)
    alpha_K = get_alpha_K(aktype, a, bg, alpha_B, cs2, alpha_K_0)
    return cs2*(alpha_K + 1.5*alpha_B**2) + 0.5*alpha_B**2 - alpha_B*(d_lnH_d_lna + 1) - 3*(1 + wde)*rhode/rhotot

def solve_alpha_B(aktype, omega_m, w0, wa, cs2, alpha_K_0, alpha_B_init=0):
    bg = Bg(omega_m=omega_m, w0=w0, wa=wa)
    
    N       = 200 # Number of steps
    a_ini   = 1e-5
    a_final = 1
    dloga   = (np.log(a_final) - np.log(a_ini))/N
    loga    = np.linspace(np.log(a_ini), np.log(a_final), N+1)
    a       = np.exp(loga)

    alpha_B = np.zeros(N+1)
    alpha_K = np.zeros(N+1)
    mu      = np.zeros(N+1)

    alpha_B[0] = alpha_B_init
    alpha_K[0] = get_alpha_K(aktype, a[0], bg, alpha_B[0], cs2, alpha_K_0)
    D_kin      = alpha_K[0] + 1.5*alpha_B[0]**2
    if alpha_B[0] == 0: mu[0] = 1
    else:
        mu[0]  = 1 + alpha_B[0]**2/(2*cs2*D_kin) if D_kin != 0 else np.inf
    for i in range(N):
        if abs(alpha_B[i]) > 1e6:
            # JVR NOTE: for many cases, \alpha_B just blows up and in Python it becomes \inf.
            # For calculating \mu, this is not a problem, since \mu has a well-defined limit when \alpha_B -> \inf
            # In practice, I enforce this with the threshold defined above in the `if` statement
            # And then I just fill the rest of the arrays with the last values and break out of the integration loop
            for j in range(i, N+1):
                alpha_B[j] = alpha_B[i]
                alpha_K[j] = get_alpha_K(aktype, a[j], bg, alpha_B[j], cs2, alpha_K_0)
                mu[j]      = mu[i]
            break
        alpha_B[i+1] = alpha_B[i] + deriv(loga[i], alpha_B[i], bg, aktype, alpha_K_0, cs2)*dloga
        alpha_K[i+1] = get_alpha_K(aktype, a[i+1], bg, alpha_B[i+1], cs2, alpha_K_0)
        D_kin      = alpha_K[i+1] + 1.5*alpha_B[i+1]**2
        if alpha_B[i+1] == 0: mu[i+1] = 1
        else:
            mu[i+1]  = 1 + alpha_B[i+1]**2/(2*cs2*D_kin) if D_kin != 0 else np.inf
        
    return loga*np.log10(np.e), alpha_B, alpha_K, mu

    
def study_alpha_B(alpha_K, cs2, omega_m, w0, wa):
    # Returns the derivative of alpha_B, evaluated in a grid of alpha_B, log(a)
    # Useful when assessing the stability of the system
    N_loga = 100 # Number of steps
    a_ini = 1e-5
    loga = np.linspace(np.log(a_ini), 0, N_loga)
    N_alpha_B = 100
    alpha_B = np.linspace(-10, 10, N_alpha_B)
    result = np.zeros((N_loga, N_alpha_B))
    for i in range(N_loga):
        for j in range(N_alpha_B):
            result[i, j] = deriv(loga[i], alpha_B[j], alpha_K, cs2, omega_m, w0, wa)
    return result, loga*np.log10(np.e), alpha_B

# ------------------------------------
# Utility functions for model analysis
# ------------------------------------

import matplotlib as mpl
import matplotlib.pyplot as plt

param_name_latex = {
    "omega_m": "$\\Omega_m$",
    "w0"     : "$w_0$",
    "wa"     : "$w_a$",
    "cs2"    : "$c_s^2$",
    "alpha_K_0": "$\\alpha_{K,0}$",
}

def plot_scatter(aktype, priors, fig_path):
    priors_keys  = list(priors.keys())
    priors_array = np.array(list(priors.values()))
    low = priors_array[:, 0]
    high = priors_array[:, 1]

    gr_threshold = 0.05
    extreme_threshold = 0.5
    
    num_samples = 4_000
    samples_params = np.random.rand(num_samples, len(priors)) * (high - low)[None, :] + low[None, :]
    samples = []
    for sample in samples_params:
        params_dict = dict(zip(priors_keys, sample))
        _, _, _, mu = solve_alpha_B(aktype=aktype, **params_dict)
        if   np.any(mu < 0): color = "red" # I don't know if negative mu is a pathology but just in case
        elif np.all(mu < 1 + gr_threshold) and np.all(mu > 1 - gr_threshold): color = "gray" # Mapping points close to GR
        elif np.any(np.abs(mu - 1) > extreme_threshold): color = "orange" # Mapping regions with more extreme variation of mu
        else: color = "blue" # "Typical" points
        samples.append((sample, color))

    fig, axs = plt.subplots(len(priors), len(priors), figsize=(12, 12), gridspec_kw={"hspace": 0, "wspace": 0})
    points = np.array([sample[0] for sample in samples])
    colors = np.array([sample[1] for sample in samples])
    for row in range(len(priors)):
        for col in range(len(priors)):
            ax = axs[row, col]
            if row == len(priors)-1: ax.set_xlabel(list(param_name_latex.values())[col])
            else: ax.set_xticklabels([])
            if col == 0:             ax.set_ylabel(list(param_name_latex.values())[row])
            else: ax.set_yticklabels([])
                
            if col >= row: ax.remove()

            ax.scatter(points[:, col], points[:, row], color=colors, s=1)

    blue_circle = mpl.lines.Line2D([], [], marker="o", color="w", markerfacecolor="blue", label="$0.5 < \\mu < 1.5$ and", )
    gray_circle = mpl.lines.Line2D([], [], marker="o", color="w", markerfacecolor="gray", label=f"$|\\mu - 1| < {gr_threshold}$", )
    orange_circle = mpl.lines.Line2D([], [], marker="o", color="w", markerfacecolor="orange", label=f"$|\\mu - 1|$ above {extreme_threshold}", )
    red_circle = mpl.lines.Line2D([], [], marker="o", color="w", markerfacecolor="red", label=f"$\\mu < 0$", )
    fig.legend(handles=[gray_circle, blue_circle, orange_circle, red_circle], bbox_to_anchor=(0.75, 0.75))
    plt.savefig(fig_path, bbox_inches="tight")
