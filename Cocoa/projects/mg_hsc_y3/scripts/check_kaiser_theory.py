#!/usr/bin/env python3
"""
Kaiser RSDの理論値と実装値を比較
"""
import numpy as np
import sys

sys.path.append('/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/MGCAMB')
import camb

print("=" * 80)
print("Kaiser RSDの理論値と実装値の比較")
print("=" * 80)

# パラメータ
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
b1 = 1.2338

print(f"\nパラメータ:")
print(f"   z = {z_eff}")
print(f"   b1 = {b1} (bias only)")
print(f"   (1+b1) = {1+b1}")

# CAMBでP(k)を計算
print(f"\n📊 CAMBで線形P(k)を計算中...")
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu)
pars.InitPower.set_params(As=As, ns=ns)
pars.set_dark_energy(w=w)

h = H0 / 100.0

# Matter power spectrumを設定
pars.set_matter_power(redshifts=[z_eff], kmax=2.0)
pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)

# f(z)を計算
sigma8_z = results.get_sigma8()[0]
fsigma8_z = results.get_fsigma8()[0]
f_z = fsigma8_z / sigma8_z

print(f"   sigma8(z={z_eff}) = {sigma8_z:.6f}")
print(f"   f(z={z_eff}) = {f_z:.6f}")
print(f"   fsigma8(z={z_eff}) = {fsigma8_z:.6f}")

# P(k)を取得
kh_camb, z_camb, pk_camb = results.get_matter_power_spectrum(
    minkh=k_h_mpc[0]*h, maxkh=k_h_mpc[-1]*h, npoints=200)

from scipy.interpolate import interp1d
pk_interp_mpc3 = interp1d(kh_camb, pk_camb[0, :], kind='cubic', fill_value='extrapolate')

# k値を1/Mpc単位に変換
k_mpc = k_h_mpc * h

# P(k)を取得（Mpc³単位）
plin_mpc3 = pk_interp_mpc3(k_mpc)

# (Mpc/h)³単位に変換
plin_mpch3 = plin_mpc3 * (h**3)

print(f"\n線形P(k)（最初の3点）:")
for i in range(3):
    print(f"   k = {k_h_mpc[i]:.6f} h/Mpc → Plin = {plin_mpch3[i]:.2f} (Mpc/h)³")

# Kaiser RSDの理論値を計算
beta = f_z / (1 + b1)
Kaiser_P0_theory = ((1 + b1)**2) * (1 + 2/3*beta + 1/5*beta**2) * plin_mpch3
Kaiser_P2_theory = ((1 + b1)**2) * (4/3*beta + 4/7*beta**2) * plin_mpch3

print(f"\n📊 Kaiser RSD理論値（最初の3点）:")
print(f"   β = f/(1+b1) = {beta:.4f}")
print(f"   P0係数 = {((1+b1)**2) * (1 + 2/3*beta + 1/5*beta**2):.4f}")
print(f"   P2係数 = {((1+b1)**2) * (4/3*beta + 4/7*beta**2):.4f}")
print(f"")
for i in range(3):
    print(f"   k = {k_h_mpc[i]:.6f} → P0 = {Kaiser_P0_theory[i]:.2f}, P2 = {Kaiser_P2_theory[i]:.2f}")

# 実装値を読み込む
kaiser_file = './projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector'
kaiser_data = np.loadtxt(kaiser_file)
values = kaiser_data[:, 1]

Nk = 40
Nz = 1
total_size = Nk * Nz

P0_impl = values[:total_size]
P2_impl = values[total_size:2*total_size]

print(f"\n📊 Kaiser RSD実装値（最初の3点）:")
for i in range(3):
    print(f"   k = {k_h_mpc[i]:.6f} → P0 = {P0_impl[i]:.2f}, P2 = {P2_impl[i]:.2f}")

print("\n" + "=" * 80)
print("📊 比較（k=0.007424 h/Mpc）:")
print("=" * 80)
print(f"   線形P(k) = {plin_mpch3[0]:.2f} [(Mpc/h)³]")
print(f"")
print(f"   理論 P0 = {Kaiser_P0_theory[0]:.2f}")
print(f"   実装 P0 = {P0_impl[0]:.2f}")
print(f"   比率 = {P0_impl[0] / Kaiser_P0_theory[0]:.4f}")
print(f"   差 = {abs(P0_impl[0] - Kaiser_P0_theory[0]) / Kaiser_P0_theory[0] * 100:.1f}%")
print(f"")
print(f"   理論 P2 = {Kaiser_P2_theory[0]:.2f}")
print(f"   実装 P2 = {P2_impl[0]:.2f}")
print(f"   比率 = {P2_impl[0] / Kaiser_P2_theory[0]:.4f}")
print(f"   差 = {abs(P2_impl[0] - Kaiser_P2_theory[0]) / Kaiser_P2_theory[0] * 100:.1f}%")

print("\n" + "=" * 80)
print("🎯 結論:")
print("=" * 80)

if abs(P0_impl[0] / Kaiser_P0_theory[0] - 1.0) < 0.05:
    print("✅ Kaiser RSD実装は理論値とほぼ一致しています！")
else:
    print(f"❌ Kaiser RSD実装が理論値と{abs(P0_impl[0] / Kaiser_P0_theory[0] - 1.0) * 100:.1f}%違います！")
    
    # 原因を推測
    ratio_p0 = P0_impl[0] / Kaiser_P0_theory[0]
    ratio_p2 = P2_impl[0] / Kaiser_P2_theory[0]
    
    print(f"\n   P0比率 = {ratio_p0:.4f}")
    print(f"   P2比率 = {ratio_p2:.4f}")
    
    if abs(ratio_p0 - ratio_p2) < 0.01:
        print(f"\n   → P0とP2の比率がほぼ同じなので、線形P(k)の取得に問題があります")
        print(f"   → 実装が使っている線形P(k) = {P0_impl[0] / ((1+b1)**2 * (1 + 2/3*beta + 1/5*beta**2)):.2f}")
        print(f"   → 理論の線形P(k) = {plin_mpch3[0]:.2f}")
        print(f"   → 比率 = {(P0_impl[0] / ((1+b1)**2 * (1 + 2/3*beta + 1/5*beta**2))) / plin_mpch3[0]:.4f}")












