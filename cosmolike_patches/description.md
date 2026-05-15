# Modifications to Cosmolike to be able to pass \mu and \Sigma

- `des_y3/likelihood/_cosmolike_prototype_base.py`:
  1. Add `"theory": {"camb": None}` to the likelihood's requirements so we can access the functions
  2. In `ci.set_cosmology()`, add the `z_MG`, `mu` and `sigma` from CAMB
- `des_y3/interface/interface.cpp`:
  1. L187, `set_cosmology()` function: add mu and sigma as parameters, and then call `set_mg_functions(io_z_2D, mu, sigma)`
- `external_modules/code/cosmolike/generic_interface.cpp`:
  1. Add another function `set_mg_functions(vector io_z, vector io_mu, vector io_sigma)` just like `set_growth(vector io_z, vector io_G)`
- `external_modules/code/cosmolike/structs.h`:
  1. Add two fields to struct `cosmopara`, `size_t MGfuncs_nz` and `double **MGfuncs` such that `MGFuncs[0, :] = z_MG`, `MGFuncs[1, :] = mu` and `MGFuncs[2, :] = sigma`
- `external_modules/code/cosmolike/cosmo2D.c`:
  1. In the function `int_for_C_ss_tomo_limber`, introduce `double mu_a = get_mu(a)` and `double sigma = get_sigma(a)`
  2. In the function `int_for_C_ss_tomo_limber`, within the `switch (nuisance.IA_MODEL)`, in the case `IA_MODEL_NLA` branch, the variable `ans` is a sum of four terms. Multiply each one by the `mu` and `sigma` factors accordingly (i.e.) C_KK goes as sigma^2, C_KI and C_IK goes as mu*sigma, and C_II goes as mu^2.
- `external_modules/code/cosmolike/cosmo3D.c`:
  1. Introduce functions `double get_mu(double a)` and `double get_sigma(double a)` that interpolate your mu and sigma functions at the desired scale factor. The module already imports `gsl_spline.h`.