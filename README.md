# Dark Energy Sound Speed

This repository contains analysis files for a research project investigating the dark energy sound speed.

## Contents

The file `CHAINS.md` enumerates the chains in this project. The chains can be downloaded from CHE using the `download_chains.sh` script. The corresponding `.yaml` files are located in `yamls/`. Once downloaded, the chains are located in `chains/`.

In `analysis/`, we perform the statistical analysis of the chains. The notebook `analysis.ipynb` centralizes the analysis. Plots are stored in `plots/`.

## Models
- w0wa = plain w0wa
- w0waq = w0wa without phantom crossing
- cs2w = sound speed proportional to equation of state
- cs2r = reconstructed (binned) sound speed
- cs2q = constant sound speed equal to 1

## Dataset combinations
DS1: P18+DESI+DESY5
DS2: P18+DESI+DESY5+DESY3SHEAR

## References
- CAMB-cs2: https://github.com/joaoreboucas1/CAMB-cs2