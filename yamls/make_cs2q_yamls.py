import os

original_cs2_block = """  cs2_0:
    prior:
      min: 0
      max: 1
    ref:
      min: 0.8
      max: 1
    proposal: 0.1
    latex: c_s^2"""

new_cs2_block = """  cs2_0:
    value: 1.0"""

original_indices = [56, 57, 58, 59, 79, 80, 81, 82]
new_indices = range(85, 85+8)

for original_i, new_i in zip(original_indices, new_indices):
    original_filename = f"MCMC{original_i}.yaml"
    new_filename = f"MCMC{new_i}.yaml"
    with open(original_filename, "r") as f: contents = f.read()
    assert original_cs2_block in contents
    new_contents = contents.replace(original_cs2_block, new_cs2_block)
    new_contents = new_contents.replace(f"MCMC{original_i}/MCMC{original_i}", f"MCMC{new_i}/MCMC{new_i}")
    with open(new_filename, "w") as f: f.write(new_contents)


