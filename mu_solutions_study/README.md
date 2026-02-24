# Studying solutions for $\alpha_B$ and $\mu$

See the [Overleaf](https://www.overleaf.com/project/682637cf7d5a109c12b31a0a) for details and references. In our parametrization for modified gravity, the function $\alpha_B$ is given by Equation 2.18, and the $\mu$ and $\Sigma$ functions are given by Equation 2.21. Cosmological perturbations are very sensitive to the $\mu$ and $\Sigma$ functions, so they must behave "nicely" and not deviate extremely from the GR values $\mu = \Sigma = 1$. In this folder, I investigate the solutions for $\mu$ and $\alpha_B$ for different values of the cosmological parameters.

The solutions still depend on the specific form of $\alpha_K$ which is left undetermined. We consider five parametrizations:
- Constant
- Proportional to $\Omega_\mathrm{DE}$
- Quintessence
- Cubic Galileon
- Proportional to $\alpha_B$ (DEPRECATED)

The algorithm that solves Equation 2.18 and calculates $\mu$ is given in `engine.py`.

For each parametrization, I investigate two aspects:
- The form of $\alpha_B$ and $\mu$ as we vary individual cosmological parameters
- Regions of the parameter space where $\mu$ deviates significantly from unity
Both investigations result in Figures which are stored in `plots`.