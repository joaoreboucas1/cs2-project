# Dark Energy Sound Speed

This repository contains analysis files for a research project investigating the dark energy sound speed using fluid and modifed gravity models.

## Quick Links

- Overleaf: https://www.overleaf.com/project/682637cf7d5a109c12b31a0a
- Google Doc: https://docs.google.com/document/d/1-YE5_7iMoF8SzzRpmSuoW7ugUeewO7i0S2rXY5P_Mfk/edit?tab=t.0
- CAMB-cs2: https://github.com/joaoreboucas1/CAMB-cs2

## Contents

- `CHAINS.md` is a list of chains performed within this project and their respective indices.
- The corresponding `.yaml` files are located in `yamls/`.
- `download_chains.py` downloads chains from the CHE cluster.
- `run_mcmc.sh` is a bash script for running chains in the CHE cluster.
- `data_analysis/` contains Jupyter notebooks and plots for the statistical data analysis of the theoretical cs2 models.
- `mu_solutions_study/` contains Jupyter notebooks and plots studying the behavior of the evolution of $\mu$ across the different parametrizations we consider.
- `observables_impact/` contains Jupyter notebooks and plots studying the impact of the fluid and MG models on cosmological power spectra.

## Models

### Background Models
- lcdm
- w0wa

### Stress models
- Fluid
- MG with constant $\alpha_K$
- MG with $\alpha_K$ proportional to $\Omega_\mathrm{DE}$
- MG with $\alpha_K$ like k-essence
- MG with $\alpha_K$ like Cubic Galileon

## Dataset combinations
DS1: P18+DESI+DESY5
DS2: P18+DESI+DESY5+DESY3SHEAR