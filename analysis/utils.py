import numpy as np
import getdist

# ----- Matplotlib style and constants -----

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = "stix"
mpl.rcParams['font.family'] = "STIXGeneral"

# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

# ----- GetDist helper functions -----

DEFAULT_GETDIST_SETTINGS = {
    "ignore_rows": 0.4,
}

def load_chain(index, settings=DEFAULT_GETDIST_SETTINGS):
    chain = getdist.loadMCSamples(f"../chains/MCMC{index}/MCMC{index}", settings=settings)
    if not chain.paramNames.hasParam("omegam"):
        omegam = (chain["omegach2"]+chain["omegabh2"]+0.06/93.15)/(chain["H0"]/100)**2
        chain.addDerived(omegam, name="omegam", label="\Omega_m")
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