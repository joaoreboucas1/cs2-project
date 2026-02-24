"""
    Engine for sound speed modified gravity models
    Author: João Rebouças, December 2025
"""

from enum import Enum
import numpy as np
from scipy.integrate import solve_ivp

omega_r = 2.5e-5

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
    w_de     = bg.w0 + bg.wa*(1-a)
    rhog     = rho_gamma(a, bg)
    rhom     = rho_m(a, bg)
    rhode    = rho_de(a, bg)
    rhotot   = rhog + rhode + rhom
    w_no_de  = (rhog/3 + w_de*rhode)/rhotot
    omega_de = rhode/rhotot
    rho_plus_p_no_de_over_rhotot = ((4/3)*rhog + rhom)/rhotot
    return w_de, w_no_de, rhotot, rho_plus_p_no_de_over_rhotot, omega_de

class alphaKtype(Enum):
    CONST = 1
    OMEGA = 2
    QUINT = 3
    CUGAL = 4
    PROP  = 5

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
            return omega_de*(1 + w_de)/cs2
        case alphaKtype.CUGAL:
            # NOTE: E^2 = (H/H0)^2 = rhotot/rho_cr_0
            rhode  = rho_de(a, bg)
            rhotot = rho_m(a, bg) + rho_gamma(a, bg) + rhode
            omega_de = rhode/rhotot
            return alpha_K_0*omega_de/rhotot**2 + 6*alpha_B
        case _:
            raise Exception("Unknown alpha_K parametrization")

def deriv(loga, alpha_B, bg, aktype, alpha_K_0, cs2):
    a = np.exp(loga)
    w_de, w_no_de, rhotot, rho_plus_p_no_de_over_rhotot, omega_de = get_bg_funcs(a, bg)
    alpha_K = get_alpha_K(aktype, a, bg, alpha_B, cs2, alpha_K_0)
    return cs2*(alpha_K + 1.5*alpha_B**2) + (alpha_B - 2)*(1.5*(1 + w_no_de) + 0.5*alpha_B) + 3*rho_plus_p_no_de_over_rhotot

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
    mu[0]      = 1 + alpha_B[0]**2/2/(alpha_K[0] + 1.5*alpha_B[0]**2)/cs2
    for i in range(N):
        if alpha_B[i]**2 > 1e6*abs(alpha_K[i]):
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
        mu[i+1]      = 1 + alpha_B[i+1]**2/2/(alpha_K[i+1] + 1.5*alpha_B[i+1]**2)/cs2
        
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