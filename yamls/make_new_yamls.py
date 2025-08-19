logcs2_block = r"""  log10_cs2_0:
    prior:
      min: -5
      max: 0
    proposal: 0.2
    latex: \log_10(c_{s,0}^2)
    drop: true"""

cs2_0_block = r"""  cs2_0:
    value: 'lambda log10_cs2_0: 10**log10_cs2_0'"""

cs2_0_new_block = r"""  cs2_0:
    prior:
      min: -0.5
      max: 1.5
    ref:
      dist: norm
      loc: 0.9
      scale: 0.1"""

corresponding = {
    3: 20,
    4: 21,
    7: 22,
    8: 23,
}

targets = [f"MCMC{i}.yaml" for i in corresponding]
for target in targets:
    index = int(target[4:-5])
    new_index = corresponding[index]
    print(f"{index} -> {new_index}")
    with open(target, "r") as f:
        contents = f.read()
    new_contents = contents.replace(f"MCMC{index}", f"MCMC{new_index}")
    assert logcs2_block in new_contents
    assert cs2_0_block in new_contents
    new_contents = new_contents.replace(logcs2_block, "")
    new_contents = new_contents.replace(cs2_0_block, cs2_0_new_block)
    with open(f"MCMC{new_index}.yaml", "w") as f:
        f.write(new_contents)
    