# Data analysis for sound speed models

In this folder, I make statistical analysis of the chains in `../CHAINS.md`. The plots created are stored in `plots/`. The data analysis generates several products:
- I make triangle plots for cosmological parameters of interest using the chains;
- I also extract the correspondent constraints on $\mu$ from the chains.

These two products are generated for several batches. Each batch varies the dataset, the prior on $c_s^2$, whether the sound speed is dynamical and an additional consistency check setting constant dark energy, $w = -1$. Within each batch, we vary the $\alpha_K$ parametrization used.

Tags:
- "DS1" means the combination CMB + BAO + SN;
- "Subluminal" means $c_s^2 < 1$, and "superluminal" means $c_s^2 < 10$;
- "LCDM consistency" means the check we perform setting $w = -1$.

The `utils.py` file contains several useful functions used throughout the analysis.