#!/usr/bin/env python3
"""
Notebookのコードを完全に再現して、同じパラメータで確認
"""
import numpy as np
import sys

sys.path.append('/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/MGCAMB')
sys.path.append('/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/velocileptors')

import camb
from velocileptors.LPT.moment_expansion_fftw import MomentExpansion

print("=" * 80)
print("Notebookコードの完全再現")
print("=" * 80)

# パラメータ（notebookと同じ）
H0 = 67.089
ombh2 = 0.022242
omch2 = 0.11905
mnu = 0.06
As = 2.2065e-9
ns = 0.9645
w = -1.0

# BOSS DR12 k-bins
k_bins_file = './projects/mg_hsc_y3/data/boss_dr12_cmass_k_bins.txt'
k_bins_data = np.loadtxt(k_bins_file)
k_h_mpc = k_bins_data[:, 1]  # h/Mpc単位

# 赤方偏移とバイアス
z_eff = 0.5893
b1 = 1.2338  # bias（biasのみ、1+biasではない）

print(f"\nパラメータ:")
print(f"   H0 = {H0}")
print(f"   z = {z_eff}")
print(f"   b1 = {b1} (bias only)")
print(f"   (1+b1) = {1+b1}")

# CAMBでP(k)を計算
print(f"\n📊 CAMBでP(k)を計算中...")
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu)
pars.InitPower.set_params(As=As, ns=ns)
pars.set_dark_energy(w=w)

h = H0 / 100.0

# Matter power spectrumを設定
pars.set_matter_power(redshifts=[z_eff], kmax=2.0)
pars.NonLinear = camb.model.NonLinear_none  # 線形のみ
results = camb.get_results(pars)

# sigma8を計算
sigma8_z = results.get_sigma8()[0]
print(f"   sigma8(z={z_eff}) = {sigma8_z:.6f}")

# f*sigma8を取得
fsigma8_z = results.get_fsigma8()[0]
f_z = fsigma8_z / sigma8_z
print(f"   f(z={z_eff}) = {f_z:.6f}")
print(f"   fsigma8(z={z_eff}) = {fsigma8_z:.6f}")

# P(k)を取得（Cobayaと同じ方法）
kh_camb, z_camb, pk_camb = results.get_matter_power_spectrum(
    minkh=k_h_mpc[0]*h, maxkh=k_h_mpc[-1]*h, npoints=200)

# 補間してBOSS DR12のk-binsでの値を取得
from scipy.interpolate import interp1d
pk_interp_mpc3 = interp1d(kh_camb, pk_camb[0, :], kind='cubic', fill_value='extrapolate')

# k値を1/Mpc単位に変換
k_mpc = k_h_mpc * h

# P(k)を取得（Mpc³単位）
plin_mpc3 = pk_interp_mpc3(k_mpc)

print(f"\n最初の3点のP(k) [Mpc³]:")
for i in range(3):
    print(f"   k[{i}] = {k_h_mpc[i]:.6f} h/Mpc → Plin = {plin_mpc3[i]:.2f} Mpc³")

# Velocileptorsの計算（notebookと同じ）
print(f"\n📊 Velocileptorsで計算中...")
klin = k_mpc  # 1/Mpc単位
plin_mpch3 = plin_mpc3 * (h**3)  # (Mpc/h)³単位に変換

print(f"\n最初の3点のP(k) [(Mpc/h)³]:")
for i in range(3):
    print(f"   k[{i}] = {k_h_mpc[i]:.6f} h/Mpc → Plin = {plin_mpch3[i]:.2f} (Mpc/h)³")

# MomentExpansionを初期化
mome = MomentExpansion(klin, plin_mpch3, threads=1,
                      cutoff=10, extrap_min=-4, extrap_max=3, jn=10,
                      nk=len(klin), kmin=k_mpc[0], kmax=k_mpc[-1])

# バイアスパラメータ（線形のみ）
pars_velo = [b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Velocileptorsで計算
kv, p0_velo, p2_velo, p4_velo = mome.compute_redshift_space_power_multipoles(
    pars_velo, f_z, ngauss=4, reduced=True)

print(f"\nVelocileptors結果（最初の3点）:")
for i in range(3):
    print(f"   k[{i}] = {kv[i]:.6f} → P0 = {p0_velo[i]:.2f}, P2 = {p2_velo[i]:.2f}")

# Kaiserで計算（notebookと同じ）
print(f"\n📊 Kaiser RSDで計算中...")
beta = f_z / (1 + b1)
Kaiser_P0 = ((1 + b1)**2) * (1 + 2/3*beta + 1/5*beta**2) * plin_mpch3
Kaiser_P2 = ((1 + b1)**2) * (4/3*beta + 4/7*beta**2) * plin_mpch3

print(f"\nKaiser結果（最初の3点）:")
for i in range(3):
    print(f"   k[{i}] = {k_h_mpc[i]:.6f} → P0 = {Kaiser_P0[i]:.2f}, P2 = {Kaiser_P2[i]:.2f}")

# 比較
print("\n" + "=" * 80)
print("📊 比較（k=0.007424 h/Mpc）:")
print("=" * 80)
print(f"   線形P(k) = {plin_mpch3[0]:.2f} [(Mpc/h)³]")
print(f"")
print(f"   Kaiser P0 = {Kaiser_P0[0]:.2f}")
print(f"   Velo P0 = {p0_velo[0]:.2f}")
print(f"   比率 = {p0_velo[0] / Kaiser_P0[0]:.4f}")
print(f"   差 = {abs(p0_velo[0] - Kaiser_P0[0]) / Kaiser_P0[0] * 100:.1f}%")
print(f"")
print(f"   Kaiser P2 = {Kaiser_P2[0]:.2f}")
print(f"   Velo P2 = {p2_velo[0]:.2f}")
print(f"   比率 = {p2_velo[0] / Kaiser_P2[0]:.4f}")
print(f"   差 = {abs(p2_velo[0] - Kaiser_P2[0]) / Kaiser_P2[0] * 100:.1f}%")

print("\n" + "=" * 80)
print("🎯 結論:")
print("=" * 80)

if abs(p0_velo[0] / Kaiser_P0[0] - 1.0) < 0.05:
    print("✅ NotebookスタイルのKaiserとVelocileptorsはほぼ一致しています！")
    print(f"   P0差: {abs(p0_velo[0] - Kaiser_P0[0]) / Kaiser_P0[0] * 100:.1f}%")
    print(f"   P2差: {abs(p2_velo[0] - Kaiser_P2[0]) / Kaiser_P2[0] * 100:.1f}%")
    
    # 実装との比較
    kaiser_impl = 67515.52
    velo_impl = 119277.25
    
    print(f"\n現在の実装との比較:")
    print(f"   理論Kaiser P0 = {Kaiser_P0[0]:.2f}")
    print(f"   実装Kaiser P0 = {kaiser_impl:.2f}")
    print(f"   比率 = {kaiser_impl / Kaiser_P0[0]:.4f}")
    
    print(f"\n   理論Velo P0 = {p0_velo[0]:.2f}")
    print(f"   実装Velo P0 = {velo_impl:.2f}")
    print(f"   比率 = {velo_impl / p0_velo[0]:.4f}")
else:
    print("❌ NotebookスタイルでもKaiserとVelocileptorsに差があります！")












