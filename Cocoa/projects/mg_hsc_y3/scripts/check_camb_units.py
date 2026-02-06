#!/usr/bin/env python3
"""
CAMBのPKL.logPが返す単位を確認
"""
import numpy as np
import sys

sys.path.insert(0, '/home/tanida/cocoa_v41/cocoa/Cocoa/external_modules/code/MGCAMB')
import camb

print("="*80)
print("CAMBのPKL.logP単位確認")
print("="*80)

# パラメータ
H0 = 67.089
ombh2 = 0.022242
omch2 = 0.11905
mnu = 0.06
As = 2.2065e-9
ns = 0.9645
w = -1.0

h = H0 / 100.0
z_eff = 0.5893
k_test = 0.007424  # h/Mpc

print(f"\nパラメータ:")
print(f"  H0 = {H0}, h = {h:.6f}")
print(f"  z = {z_eff}")
print(f"  k = {k_test} h/Mpc")

# CAMB設定
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=mnu)
pars.InitPower.set_params(As=As, ns=ns)
pars.set_dark_energy(w=w)
pars.set_matter_power(redshifts=[z_eff], kmax=2.0)
pars.NonLinear = camb.model.NonLinear_none

results = camb.get_results(pars)

# PKLを作成（Cobayaと同じ方法）
PKL = results.get_matter_power_interpolator(
    nonlinear=False, 
    extrap_kmax=3.0,
    hubble_units=False,  # k in 1/Mpc, P in Mpc^3
    k_hunit=False,
    var1='delta_tot',
    var2='delta_tot'
)

# kを1/Mpc単位に変換
k_mpc = k_test * h  # h/Mpc → 1/Mpc

print(f"\n単位変換:")
print(f"  k = {k_test} h/Mpc = {k_mpc:.9f} 1/Mpc")

# P(k)を取得（Cobayaと同じ方法）
pk_mpc3 = PKL.P(z_eff, k_mpc)  # Mpc^3 単位

print(f"\nCAMBから取得した値:")
print(f"  PKL.P(z, k) = {pk_mpc3:.2f} Mpc^3")

# logPを取得（logPメソッドがない場合はnp.logを使う）
logP = np.log(pk_mpc3)
exp_logP = np.exp(logP)

print(f"  logP = log(PKL.P) = {logP:.6f}")
print(f"  exp(logP) = {exp_logP:.2f} Mpc^3")

# (Mpc/h)^3単位に変換
pk_mpch3 = pk_mpc3 * (h**3)
exp_logP_mpch3 = exp_logP * (h**3)

print(f"\n単位変換後:")
print(f"  P × h³ = {pk_mpch3:.2f} (Mpc/h)^3")
print(f"  exp(logP) × h³ = {exp_logP_mpch3:.2f} (Mpc/h)^3")

print(f"\nCobaya側の処理:")
print(f"  lnPL = PKL.logP(...)")
print(f"  lnPL += np.log(h³) = lnPL + {np.log(h**3):.6f}")
print(f"  最終的な lnPL = {logP + np.log(h**3):.9f}")
print(f"  exp(最終的な lnPL) = {np.exp(logP + np.log(h**3)):.2f} (Mpc/h)^3")

print(f"\n比較:")
print(f"  理論値 (CAMB直接) = {pk_mpch3:.2f} (Mpc/h)^3")
print(f"  exp(lnPL + ln(h³)) = {np.exp(logP + np.log(h**3)):.2f} (Mpc/h)^3")
print(f"  比率 = {np.exp(logP + np.log(h**3)) / pk_mpch3:.6f}")

if abs(np.exp(logP + np.log(h**3)) / pk_mpch3 - 1.0) < 0.001:
    print("\n✅ Cobayaの処理は正しいです！")
else:
    print(f"\n❌ 比率が{abs(np.exp(logP + np.log(h**3)) / pk_mpch3 - 1.0) * 100:.2f}%違います！")

print(f"\n実際のC++側で使われている値（デバッグ出力から）:")
print(f"  exp(lnP) = 9553.16 (Mpc/h)^3")
print(f"  理論値 = 2537.0 (Mpc/h)^3")
print(f"  比率 = {9553.16 / 2537.0:.6f}")

# もし、lnPが ln(Mpc^3) なら
lnP_mpc3 = np.log(pk_mpc3)
print(f"\nもし、lnP = ln(Mpc^3) なら:")
print(f"  lnP = {lnP_mpc3:.6f}")
print(f"  exp(lnP) = {pk_mpc3:.2f} Mpc^3")
print(f"  exp(lnP + ln(h³)) = {pk_mpc3 * h**3:.2f} (Mpc/h)^3")
print(f"  理論値との比率 = {(pk_mpc3 * h**3) / pk_mpch3:.6f}")

