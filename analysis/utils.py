import getdist

DEFAULT_GETDIST_SETTINGS = {
    "ignore_rows": 0.4,
}

# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5
colors = ["#1b9e77", "#d95f02", "#7570b3"]

def load_chain(index, settings=DEFAULT_GETDIST_SETTINGS):
    chain = getdist.loadMCSamples(f"../chains/MCMC{index}/MCMC{index}", settings=settings)
    if chain.paramNames.hasParam("log10_cs2_0"):
        chain.paramNames.parWithName("log10_cs2_0").label = r"\log_{10}(c_{s,0}^2)"
    return chain

# Parsing CHAINS.md file
with open("../CHAINS.md", "r") as f:
    chains_table = f.read().splitlines()

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