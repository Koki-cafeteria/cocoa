#!/usr/bin/env python3
"""
Kaiser RSD と Velocileptors で使用している線形P(k)を比較
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# MGCAMBのパスを追加
sys.path.insert(0, '/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/MGCAMB')
import camb
from camb import model

print("=" * 70)
print("線形P(k)の比較: Kaiser RSD vs Velocileptors")
print("=" * 70)

# 宇宙論パラメータ
params = camb.CAMBparams()
params.set_cosmology(
    H0=67.089,
    ombh2=0.0222499167450714,
    omch2=0.11915441939660276,
    mnu=0.06,
    omk=0,
    tau=0.0697186,
    nnu=3.046,
    num_massive_neutrinos=1
)

params.MG_flag = 1
params.pure_MG_flag = 2
params.mu0 = 0.0
params.sigma0 = 0.0
params.InitPower.set_params(As=2.2065e-09, ns=0.9645)
params.set_dark_energy(w=-1.0, dark_energy_model='ppf')

z_eff = 0.5893
params.set_matter_power(redshifts=[z_eff], kmax=5.0, k_per_logint=20)
params.NonLinear = model.NonLinear_none

results = camb.get_results(params)

# k-gridを準備（Kaiser RSDのk-binsを使用）
# ファイルの2列目（k_eff）のみを読み込む
k_kaiser = np.loadtxt('./projects/mg_hsc_y3/data/boss_dr12_cmass_k_bins.txt', usecols=1)

# MGCAMBからP(k)を取得
h = params.H0 / 100.0
kh, z_arr, pk = results.get_matter_power_spectrum(
    minkh=1e-4,
    maxkh=k_kaiser.max(),
    npoints=200
)
pk_z = pk[0, :]  # z=0.5893でのP(k) in (Mpc/h)^3

print(f"\n📊 MGCAMB P(k) (z={z_eff}):")
print(f"   k範囲: {kh.min():.4e} - {kh.max():.4e} h/Mpc")
print(f"   P(k)範囲: {pk_z.min():.4e} - {pk_z.max():.4e} (Mpc/h)³")

# Kaiser RSDのk値でP(k)を補間
pk_kaiser = np.interp(k_kaiser, kh, pk_z)

# Velocileptorsから使用しているP(k)を再現
# （_cosmolike_prototype_base.pyのcompute_velocileptors_pkと同じ処理）
print(f"\n📊 Velocileptors P(k)の取得方法を再現:")

# Velocileptorsのk-bins
k_velo = np.logspace(np.log10(7.424e-3), np.log10(0.395040), 30)
pk_velo = np.interp(k_velo, kh, pk_z)

print(f"   k範囲: {k_velo.min():.4e} - {k_velo.max():.4e} h/Mpc")
print(f"   P(k)範囲: {pk_velo.min():.4e} - {pk_velo.max():.4e} (Mpc/h)³")

# Kaiser RSDのk値でVelociléptorsのP(k)を補間
pk_velo_interp = np.interp(k_kaiser, k_velo, pk_velo)

# 比較
rel_diff = (pk_velo_interp - pk_kaiser) / pk_kaiser * 100

print(f"\n📊 線形P(k)の相対差:")
print(f"   平均: {np.mean(rel_diff):.4f}%")
print(f"   RMS: {np.sqrt(np.mean(rel_diff**2)):.4f}%")
print(f"   最大: {np.max(np.abs(rel_diff)):.4f}%")

# プロット
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# P(k)の比較
ax = axes[0]
ax.loglog(k_kaiser, pk_kaiser, 'b-', linewidth=2, label='Kaiser RSD k-bins')
ax.loglog(k_velo, pk_velo, 'ro', markersize=5, label='Velocileptors k-bins')
ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
ax.set_ylabel(r'$P_{\rm lin}(k)$ [(Mpc/h)$^3$]', fontsize=12)
ax.set_title(f'Linear P(k) at z={z_eff:.4f}', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# 相対差
ax = axes[1]
ax.semilogx(k_kaiser, rel_diff, 'k-', linewidth=2)
ax.axhline(0, color='r', linestyle='--', alpha=0.5)
ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
ax.set_ylabel('Relative Difference [%]', fontsize=12)
ax.set_title(r'(Velo - Kaiser) / Kaiser for $P_{\mathrm{lin}}(k)$', fontsize=14)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('./projects/mg_hsc_y3/data/comparison_linear_pk.png', dpi=150, bbox_inches='tight')
print(f"\n✅ プロット保存: ./projects/mg_hsc_y3/data/comparison_linear_pk.png")
plt.show()

print("\n✅ 完了!")

