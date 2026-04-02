"""
    Simple tool to change the cs2 prior in some yaml files and generate new yaml files
"""

# subluminal_indices = range(56, 60) # This is for constant cs2
# offset = 8 # NOTE: <subluminal chains>|<lcdm chains>|<superluminal chains>

subluminal_indices = range(68, 72)
offset = 4

subluminal_cs2_prior_block = """cs2_0:
    prior:
      min: 0
      max: 1"""
superluminal_cs2_prior_block = """cs2_0:
    prior:
      min: 0
      max: 10"""

for i in subluminal_indices:
    subluminal_file = f"MCMC{i}.yaml"
    with open(subluminal_file, "r") as f: contents = f.read()
    new_i = i + offset
    
    assert(f"MCMC{i}/MCMC{i}" in contents)
    contents = contents.replace(f"MCMC{i}/MCMC{i}", f"MCMC{new_i}/MCMC{new_i}")
    
    assert(subluminal_cs2_prior_block in contents)
    contents = contents.replace(subluminal_cs2_prior_block, superluminal_cs2_prior_block)
    
    superluminal_file = f"MCMC{new_i}.yaml" 
    with open(superluminal_file, "w") as f: f.write(contents)
