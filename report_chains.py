"""
    Utility for keeping track of the chains
"""

import os

chain_indices = [folder[4:] for folder in os.listdir("./chains") if folder.startswith("MCMC") and os.path.isdir(f"./chains/{folder}")]
for i in chain_indices:
    # Check contents from MCMC{i}.progress
    with open(f"./chains/MCMC{i}/MCMC{i}.progress", "r") as f: contents = f.read().splitlines()
    last_status = contents[-1].split()
    assert len(last_status) == 5, f"Last status should have 5 entries, found {last_status}"
    num_samples, last_update, acc_rate, rminusone, rminusone_cl = last_status
    print(f"MCMC{i} => num_samples = {round(float(num_samples))}, |R-1| = {rminusone} (updated {last_update})")
    
    # Check contents from cs2_{i}_{jobid}.out
    logs = [log for log in os.listdir("./logs") if log.startswith(f"cs2_{i}") and log.endswith(".out")]
    logs.sort(key=lambda log: int(log.strip(".out").split("_")[-1]))
    last_log = logs[-1]
    with open(f"./logs/{last_log}", "r") as f:
        lines = f.read().splitlines()
    for line in reversed(lines):
        if "R-1" in line: break
    print(f"Last update in {last_log}:", " ".join(line.split()[7:14]))
    print("-------")