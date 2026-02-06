#!/usr/bin/env python3
"""
単位変換のデバッグ
"""
import numpy as np

print("="*80)
print("単位変換の診断")
print("="*80)

H0 = 67.089
h = H0 / 100.0
coverH0 = 2997.92458

print(f"\nパラメータ:")
print(f"  H0 = {H0} km/s/Mpc")
print(f"  h = {h}")
print(f"  coverH0 = {coverH0} Mpc (固定値)")

# c/H0の計算
c = 299792.458  # km/s
c_over_H0 = c / H0  # Mpc

print(f"\n物理量:")
print(f"  c = {c} km/s")
print(f"  c/H0 = {c_over_H0:.4f} Mpc")
print(f"  c/(100 Mpc^-1) = {c/100:.4f} Mpc")

print(f"\n比較:")
print(f"  coverH0 / (c/H0) = {coverH0 / c_over_H0:.6f}")
print(f"  coverH0 × h / (c/H0) = {coverH0 * h / c_over_H0:.6f}")

print(f"\nP(k)の単位変換：")
print(f"  h³ = {h**3:.6f}")
print(f"  1/h³ = {1/h**3:.6f}")
print(f"  coverH0³ = {coverH0**3:.6e}")
print(f"  (c/H0)³ = {c_over_H0**3:.6e}")
print(f"  coverH0³ / (c/H0)³ = {(coverH0/c_over_H0)**3:.6f}")

print(f"\nk値の変換チェック:")
k_test = 0.007424  # h/Mpc (Python側から渡される)
print(f"  k = {k_test} h/Mpc (Python → C++)")
print(f"  k × coverH0 = {k_test * coverH0:.4f} (現在の実装)")
print(f"  k / h × coverH0 = {k_test / h * coverH0:.4f} (正しい？)")

print(f"\n実装/理論の比率: 4.2302")
print(f"  4.2302 ≈ ? を探す")

candidates = [
    ("1/h³", 1/h**3),
    ("(1/h)²", (1/h)**2),
    ("coverH0 / (c/H0)", coverH0 / c_over_H0),
    ("(coverH0 / (c/H0))³", (coverH0 / c_over_H0)**3),
    ("h² × (1/h³)", h**2 / h**3),
    ("(coverH0 / (c/H0))² / h", (coverH0 / c_over_H0)**2 / h),
]

print(f"\n候補の検証:")
for name, value in candidates:
    diff = abs(value - 4.2302)
    print(f"  {name:<30} = {value:.4f}  (差: {diff:.4f})")
    
# 特別な組み合わせ
special = (coverH0 / c_over_H0)**2
print(f"\n⭐ 特別な候補:")
print(f"  (coverH0 / (c/H0))² = ({coverH0} / {c_over_H0:.2f})² = {special:.4f}")
print(f"  差: {abs(special - 4.2302):.6f}")

if abs(special - 4.2302) < 0.01:
    print(f"\n✅ 一致！問題は k 値の変換が2回間違っている可能性があります")
    print(f"   - C++で k × coverH0 するところを k × (coverH0/h) とすべき（1回目のh）")
    print(f"   - p_linで k / coverH0 するところを k / (coverH0/h) とすべき（2回目のh）")
    print(f"   → 合計で h² = 1/{h**2:.4f} = {1/h**2:.4f} の誤差")











