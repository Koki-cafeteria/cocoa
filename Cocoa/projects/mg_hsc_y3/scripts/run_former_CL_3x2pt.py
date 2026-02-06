import os, sys, numpy as np

GL_ROOT = "/home/tanida/GL"

sys.path.append(os.path.join(GL_ROOT, "self_module"))
import cal_muSigma_test as cal
import set_params

# 出力先（Cocoa側 scripts ディレクトリ）
OUT_DIR = "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/scripts"

# bins
r = np.logspace(-1, 2, 100)
theta = np.loadtxt(os.path.join(GL_ROOT, "covariance", "bin_xi_logcen.dat"))

# lens redshift
z_lz = 0.26982977521734514

# Planck_params (Velocileptors比較と同じ)
cparams = set_params.Planck_params()
cparams["sigma_8"] = 0.8120  # z=0規格化

# nuisance固定
for key in [
    "A_IA","eta_eff","delta_m1","delta_m2","delta_m3","delta_m4",
    "delta_z1","delta_z2","delta_z3","delta_z4","alpha_PSF","beta_PSF"
]:
    cparams[key] = 0.0

b_lz = cparams["b_lz"]

mg = cal.MG({})

for mu0 in [0.0, 1.0]:
    cparams["mu_0"] = mu0
    cparams["Sigma_0"] = 0.0
    mg.set_cosmology(cparams)

    ds = mg.DS(r, z_lz, b_lz)
    wp = mg.W_p(r, z_lz, b_lz)
    xi_p, xi_m = mg.HSC_xi()[0](theta), mg.HSC_xi()[1](theta)

    np.savetxt(
        os.path.join(OUT_DIR, f"former_dswp_mu0_{mu0:.1f}.dat"),
        np.column_stack([r, ds, wp])
    )
    np.savetxt(
        os.path.join(OUT_DIR, f"former_xi_mu0_{mu0:.1f}.dat"),
        np.column_stack([theta, xi_p, xi_m])
    )

print("Step A done")