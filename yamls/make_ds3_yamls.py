import yaml

# w0wa free cs2, w0wa cs2=1, lcdm, and the 4 MG cases
original_indices = [47, 48, 49, 56, 57, 58, 59]
new_indices = range(76, 76+len(original_indices))

# Recipe:
# 1. Open template yaml
# 2. Read its likelihood field and params field
# 3. Open original yaml
# 4. Read all fields
# 5. Update likelihood field with template likelihood field
# 6. Update params field to add DES nuisance params
# 7. Change index to new index
# 8. Save file to MCMC{new_index}.yaml

with open("TEMPLATE_SHEAR.yaml", "r") as f:
    template = yaml.safe_load(f)


def add_des_nuisance_params(original_params, template_params):
    """Add DES nuisance params from template to original params without overriding existing keys."""
    added = []
    for k, v in template_params.items():
        if k not in original_params:
            # Heuristic: add parameters that look like DES nuisance (start with DES_)
            if k.startswith("DES_") or k.startswith("DES_A") or k.startswith("DES_BTA"):
                original_params[k] = v
                added.append(k)
    return added


def update_output_path(output, old_index, new_index):
    """Replace occurrences of the old index in the output path with the new index."""
    if not isinstance(output, str):
        return output
    return output.replace(f"MCMC{old_index}", f"MCMC{new_index}")


def main():
    template_likelihood = template.get("likelihood", {})
    template_params = template.get("params", {})

    for orig, new in zip(original_indices, new_indices):
        in_name = f"MCMC{orig}.yaml"
        out_name = f"MCMC{new}.yaml"

        print(f"Processing {in_name} -> {out_name}")

        with open(in_name, "r") as f:
            orig_data = yaml.safe_load(f)

        if orig_data is None:
            print(f"Warning: {in_name} is empty or invalid YAML, skipping.")
            continue

        # 5. Update likelihood field with template likelihood field
        if template_likelihood:
            orig_data["likelihood"].update(template_likelihood)

        # 6. Update params field to add DES nuisance params (without overriding existing ones)
        orig_params = orig_data.get("params", {})
        added = add_des_nuisance_params(orig_params, template_params)
        if added:
            orig_data["params"] = orig_params
            print(f"  Added params: {', '.join(added)}")

        # 7. Change index to new index by updating output path and any explicit mentions
        if "output" in orig_data:
            orig_data["output"] = update_output_path(orig_data["output"], orig, new)

        # Save the updated file
        with open(out_name, "w") as f:
            yaml.safe_dump(orig_data, f, sort_keys=False)

    print("Done.")


if __name__ == "__main__":
    main()