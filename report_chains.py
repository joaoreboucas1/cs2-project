"""
    Utility for keeping track of the chains
"""

import os
import time

chain_indices = sorted([folder[4:] for folder in os.listdir("./chains") if folder.startswith("MCMC") and os.path.isdir(f"./chains/{folder}")])

for i in chain_indices:
    checkpoint_file_name = f"./chains/MCMC{i}/MCMC{i}.checkpoint"
    
    with open(checkpoint_file_name, "r") as f: contents = f.read().splitlines()
    
    last_R_minus_one = contents[3].split()[-1]
    print(f"MCMC{i} => |R-1| = {last_R_minus_one} (last update: {time.ctime(os.path.getmtime(checkpoint_file_name))})")
    print("-------")