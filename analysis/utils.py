import getdist

DEFAULT_GETDIST_SETTINGS = {
    "ignore_rows": 0.4,
}

# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

def load_chain(index, settings=DEFAULT_GETDIST_SETTINGS):
    chain = getdist.loadMCSamples(f"../chains/MCMC{index}/MCMC{index}", settings=settings)
    for i in range(4):
        if chain.paramNames.hasParam(f"log10_cs2_{i}"): chain.paramNames.parWithName(f"log10_cs2_{i}").label = rf"\log_{{10}}(c_{{s,{i}}}^2)"
    return chain

# Parsing CHAINS.md file
with open("../CHAINS.md", "r") as f:
    chains_table = f.read().splitlines()

# Chain specs are the settings that are varied in this work (dataset, w model, cs2 model)
chain_specs_by_index = {}
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
    chain_specs_by_index[int(index)] = {"w_model": w_model, "cs2_model": cs2_model, "dataset": dataset, "extra_info": extra_info}

dataset_descriptions = {
    "DS1": "Planck 2018 TTTEEE+lowl + DESI DR1 BAO + DESY5 SN",
    "DS2": "Planck 2018 TTTEEE+lowl + DESI DR1 BAO + DESY5 SN + DESY3 Shear"
}