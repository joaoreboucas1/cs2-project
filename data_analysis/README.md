# Data analysis for sound speed models

In this folder, I make statistical analysis of the chains in `../CHAINS.md`. The plots created are stored in `plots/`. The analysis contains several aspects:
- In `triangle_ds1_subluminal.ipynb` I make triangle plots of the chains with subluminal prior and constant sound speed. DS1 means CMB+BAO+SN without CMB lensing.
- In `mu_constraints_ds1_subluminal.ipynb` I take the chains and extract constraints on the $\mu$ function of modified gravity.

The `utils.py` file contains several useful functions used throughout the analysis.