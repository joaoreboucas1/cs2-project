import re

# Rerun chains with the scale-dependent mu (Cataneo+2024 eq. 31a) instead of the QSA one
original_indices = [56, 57, 58, 59, 79, 80, 81, 82]
new_indices = range(113, 113+8)  # 109-112 reserved for other work

for original_i, new_i in zip(original_indices, new_indices):
    original_filename = f"MCMC{original_i}.yaml"
    new_filename = f"MCMC{new_i}.yaml"
    with open(original_filename, "r") as f: contents = f.read()
    # Add use_qsa: false right after the alpha_K_parametrization line in camb extra_args
    new_contents, n = re.subn(r"(\n( *)alpha_K_parametrization: \d+)", r"\1\n\2use_qsa: false", contents)
    assert n == 1 and "use_qsa" not in contents
    assert f"MCMC{original_i}/MCMC{original_i}" in contents
    new_contents = new_contents.replace(f"MCMC{original_i}/MCMC{original_i}", f"MCMC{new_i}/MCMC{new_i}")
    with open(new_filename, "w") as f: f.write(new_contents)
