import yaml
original_indices = range(24, 29)

likelihoods = [
    "bao.sdss_dr16_baoplus_lrg",
    "bao.sdss_dr16_baoplus_elg",
    "bao.sdss_dr16_baoplus_qso",
    "bao.sdss_dr16_baoplus_lyauto",
    "bao.sdss_dr16_baoplus_lyxqso",
]

for i in original_indices:
    with open(f"MCMC{i}.yaml", "r") as f:
        contents = yaml.safe_load(f)
    for like in likelihoods:
        contents["likelihood"][like] = None
    new_i = i + 5
    contents["output"] = f"./projects/cs2-project/chains/MCMC{new_i}/MCMC{new_i}"
    with open(f"MCMC{new_i}.yaml", "w") as f:
        yaml.dump(contents, f)
    print(f"Generated file MCMC{new_i}.yaml")