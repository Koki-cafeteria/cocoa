#!/usr/bin/env python3
"""
MGCAMBから直接P(k)を取得して、Kaiser RSDとVelocileptorsの値と比較
"""
import numpy as np
import sys
sys.path.insert(0, '/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/MGCAMB')
import camb
from camb import model

print("=" * 80)
print("MGCAMBから直接P(k)を取得")
print("=" * 80)

# パラメータ（YAMLから）
params = {
    'H0': 67.089,
    'ombh2': 0.022242,  # omegab * h^2
    'omch2': 0.11905,   # omegac * h^2
    'mnu': 0.06,
    'As': 2.2065e-9,
    'ns': 0.9645,
    'w': -1.0,
}

# CAMBパラメータを設定
pars = camb.CAMBparams()
pars.set_cosmology(
    H0=params['H0'],
    ombh2=params['ombh2'],
    omch2=params['omch2'],
    mnu=params['mnu']
)
pars.InitPower.set_params(As=params['As'], ns=params['ns'])
pars.set_dark_energy(w=params['w'])

# k範囲を設定（BOSS DR12のk-bins）
k_bins_file = './projects/mg_hsc_y3/data/boss_dr12_cmass_k_bins.txt'
k_bins_data = np.loadtxt(k_bins_file)
k_bins = k_bins_data[:, 1]  # h/Mpc単位

print(f"\nk-bins（最初の5個）:")
for i in range(5):
    print(f"   k[{i}] = {k_bins[i]:.6f} h/Mpc")

# 赤方偏移
z_eff = 0.5893

# Matter power spectrumを計算
pars.set_matter_power(redshifts=[z_eff], kmax=2.0)
results = camb.get_results(pars)

# 線形P(k)を取得
h = params['H0'] / 100.0

print(f"\n📊 MGCAMB線形P(k) at z={z_eff}:")
print(f"   h = {h:.6f}")

# k値を1/Mpc単位に変換（CAMBは1/Mpc単位を期待）
k_mpc = k_bins * h  # h/Mpc → 1/Mpc

# 線形P(k)を取得（CAMBはMpc³単位で返す）
kh_camb, z_camb, pk_camb = results.get_matter_power_spectrum(
    minkh=k_mpc[0], maxkh=k_mpc[-1], npoints=len(k_bins))

# 線形P(k)をz=z_effで補間
from scipy.interpolate import interp1d
pk_interp = interp1d(kh_camb, pk_camb[0, :], kind='cubic', fill_value='extrapolate')
pk_linear_mpc3 = pk_interp(k_mpc)  # Mpc³単位

# (Mpc/h)³単位に変換
pk_linear_mpch3 = pk_linear_mpc3 * (h**3)

print(f"\n最初の5点:")
print(f"{'k[h/Mpc]':<12} {'k[1/Mpc]':<12} {'Plin[Mpc³]':<15} {'Plin[(Mpc/h)³]':<15}")
print("-" * 70)
for i in range(5):
    print(f"{k_bins[i]:<12.6f} {k_mpc[i]:<12.6f} {pk_linear_mpc3[i]:<15.2f} {pk_linear_mpch3[i]:<15.2f}")

# Kaiser RSDとVelocileptorsの値と比較
kaiser_plin = 10731.28
velo_plin = 18958.56
mgcamb_plin_mpch3 = pk_linear_mpch3[0]
mgcamb_plin_mpc3 = pk_linear_mpc3[0]

print("\n" + "=" * 80)
print("🎯 比較（k=0.007424 h/Mpc）:")
print("=" * 80)
print(f"   Kaiser RSD（逆算）:     Plin = {kaiser_plin:.2f} [(Mpc/h)³]")
print(f"   Velocileptors（逆算）:  Plin = {velo_plin:.2f} [(Mpc/h)³]")
print(f"   MGCAMB直接取得:          Plin = {mgcamb_plin_mpch3:.2f} [(Mpc/h)³]")
print(f"   MGCAMB直接取得:          Plin = {mgcamb_plin_mpc3:.2f} [Mpc³]")

print(f"\n   Kaiser / MGCAMB = {kaiser_plin / mgcamb_plin_mpch3:.4f}")
print(f"   Velo / MGCAMB = {velo_plin / mgcamb_plin_mpch3:.4f}")

if abs(kaiser_plin / mgcamb_plin_mpch3 - 1.0) < 0.01:
    print("\n✅ Kaiser RSDはMGCAMBと一致しています！")
elif abs(velo_plin / mgcamb_plin_mpch3 - 1.0) < 0.01:
    print("\n✅ VelocileptorsはMGCAMBと一致しています！")
else:
    print("\n❌ どちらもMGCAMBと一致していません！")

# 単位変換の確認
print("\n" + "=" * 80)
print("📊 単位変換の確認:")
print("=" * 80)
print(f"   h = {h:.6f}")
print(f"   h³ = {h**3:.6f}")
print(f"   Mpc³ → (Mpc/h)³: ×h³ = ×{h**3:.6f}")
print(f"   (Mpc/h)³ → Mpc³: /h³ = /{h**3:.6f}")

# もしかして、Velocileptorsは単位変換を2回やっている？
velo_plin_corrected = velo_plin / (h**3)
print(f"\n   Velo Plin / h³ = {velo_plin_corrected:.2f} [Mpc³]")
print(f"   MGCAMB Plin = {mgcamb_plin_mpc3:.2f} [Mpc³]")
print(f"   比率 = {velo_plin_corrected / mgcamb_plin_mpc3:.4f}")

if abs(velo_plin_corrected / mgcamb_plin_mpc3 - 1.0) < 0.01:
    print("\n✅ Velocileptorsは単位変換を2回やっている可能性があります！")












