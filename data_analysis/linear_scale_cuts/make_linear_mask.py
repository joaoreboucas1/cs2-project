import numpy as np

data_path = "/home/joao/cosmo/cocoa/Cocoa/projects/des_y3/data"
mask_base = np.loadtxt(f"{data_path}/3x2pt_baseline.mask", unpack=True, usecols=(1,))
halofit_dv = np.loadtxt(f"{data_path}/halofit.datavector", unpack=True, usecols=(1,))
linear_dv = np.loadtxt(f"{data_path}/linear.datavector", unpack=True, usecols=(1,))
cov_data = np.loadtxt(f"{data_path}/des_y3_cov_unblinded_final.txt")

dv_length = 400
halofit_dv = halofit_dv[:dv_length]
linear_dv  = linear_dv[:dv_length]
delta = halofit_dv - linear_dv

# Parsing cov data
cov = np.zeros((dv_length, dv_length))

for line in cov_data:
    i = int(line[0])
    j = int(line[1])

    if i >= dv_length or j >= dv_length: continue

    cov[i,j] = line[2]

    if i != j:
        cov[i,j] *= mask_base[i]*mask_base[j]
        cov[j,i] = cov[i,j]

invcov_base = np.linalg.inv(cov)

def apply_mask(invcov, mask):
    for i in range(dv_length):
        invcov[i,i] *= mask[i]
        for j in range(i, dv_length):
            invcov[i, j] *= mask[i]*mask[j]
            invcov[j, i] = invcov[i, j]

apply_mask(invcov_base, mask_base)
chi2 = delta @ invcov_base @ delta
print(f"Number of unmasked points in fiducial mask = {len(mask_base[mask_base > 0])}")
print(f"chi2 between linear and halofit in fiducial mask = {chi2}")

mask = mask_base.copy()

threshold = 1.0

while chi2 > threshold:
    best_dchi2 = 0.0
    best_idx = None
    for idx in range(len(mask)):
        if mask[idx] == 0.0: continue

        trial_mask = mask.copy()
        trial_mask[idx] = 0.0
        
        masked_invcov = invcov_base.copy()
        apply_mask(masked_invcov, trial_mask)
        
        chi2_trial = delta @ masked_invcov @ delta
        dchi2 = chi2 - chi2_trial
        
        if dchi2 > best_dchi2:
            best_dchi2 = dchi2
            best_idx = idx
    
    if best_idx is None:
        print("Could not find any improvement")
        break
    else:
        print(f"Removing data point {best_idx} with chi2 reduction of {best_dchi2}")
        mask[best_idx] = 0.0
        apply_mask(invcov_base, mask)
        chi2 = delta @ invcov_base @ delta
        print(f"New chi2 = {chi2}")

print(f"Number of unmasked points in linear mask = {len(mask[mask > 0])}")
print(f"chi2 between linear and halofit in linear mask = {chi2}")
print("Saving mask to `linear.mask`")
np.savetxt("linear.mask", mask, fmt="%.1f")