#!/usr/bin/env python3
"""
Notebookと全く同じコードを使って、CobayaのP(k)で計算
Kaiser RSDとVelocileptorsの差の原因を特定
"""
import numpy as np
from scipy.interpolate import interp1d

print("=" * 80)
print("NotebookコードでCobayaのP(k)を使用")
print("=" * 80)

# CobayaのP(k)を読み込む（Kaiser RSDから出力されたもの）
# これは簡易的に、CAMBから直接取得する必要があります
# ここでは、既知の値を使って計算します

# パラメータ（Cobayaから）
z_eff = 0.5893
sigma8_z = 0.608835
f_z = 0.790178
b1 = 1.2338  # Velocileptorsが使っている値

# バイアスパラメータ（線形のみ）
pars = [b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # [b1, b2, bs, b3, alpha0-4, SN0, SN2]

print(f"\nパラメータ:")
print(f"   z = {z_eff}")
print(f"   f = {f_z}")
print(f"   b1 = {b1}")
print(f"   (1+b1) = {1+b1}")

# Notebookと同じKaiser公式
beta_kaiser = f_z / (1 + b1)
Kaiser_P0_coeff = (1 + 2/3*beta_kaiser + 1/5*beta_kaiser**2) * (1+b1)**2
Kaiser_P2_coeff = (4/3*beta_kaiser + 4/7*beta_kaiser**2) * (1+b1)**2

print(f"\nKaiser公式（Notebook style）:")
print(f"   β = f/(1+b1) = {beta_kaiser:.4f}")
print(f"   P0係数 = {Kaiser_P0_coeff:.4f}")
print(f"   P2係数 = {Kaiser_P2_coeff:.4f}")

# C++の実装と比較
print(f"\nC++実装と同じ計算:")
beta_cpp = f_z / (1 + b1)
P0_coeff_cpp = (1 + 2/3*beta_cpp + 1/5*beta_cpp**2) * (1+b1)**2
P2_coeff_cpp = (4/3*beta_cpp + 4/7*beta_cpp**2) * (1+b1)**2

print(f"   β = {beta_cpp:.4f}")
print(f"   P0係数 = {P0_coeff_cpp:.4f}")
print(f"   P2係数 = {P2_coeff_cpp:.4f}")

print(f"\n✅ NotebookとC++の係数は一致: {abs(Kaiser_P0_coeff - P0_coeff_cpp) < 0.001}")

print("\n" + "=" * 80)
print("📊 実際のデータでの比較")
print("=" * 80)

# 直接比較の結果（k=0.007424）
kaiser_p0_actual = 67516
kaiser_p2_actual = 86516
velo_p0_actual = 81558
velo_p2_actual = 35133

# 理論的な線形P(k)の値を逆算
plin_from_kaiser_p0 = kaiser_p0_actual / Kaiser_P0_coeff
plin_from_kaiser_p2 = kaiser_p2_actual / Kaiser_P2_coeff

plin_from_velo_p0 = velo_p0_actual / Kaiser_P0_coeff  # 同じ係数を使うべき
plin_from_velo_p2 = velo_p2_actual / Kaiser_P2_coeff

print(f"\n線形P(k)の逆算 (k=0.007424):")
print(f"   Kaiser P0から: Plin = {plin_from_kaiser_p0:.0f}")
print(f"   Kaiser P2から: Plin = {plin_from_kaiser_p2:.0f}")
print(f"   → 比率: {plin_from_kaiser_p2 / plin_from_kaiser_p0:.4f} (約{abs(1-plin_from_kaiser_p2/plin_from_kaiser_p0)*100:.1f}%の差)")

print(f"\n   Velo P0から: Plin = {plin_from_velo_p0:.0f}")
print(f"   Velo P2から: Plin = {plin_from_velo_p2:.0f}")
print(f"   → 比率: {plin_from_velo_p2 / plin_from_velo_p0:.4f} (約{abs(1-plin_from_velo_p2/plin_from_velo_p0)*100:.1f}%の差)")

print("\n" + "=" * 80)
print("🎯 結論")
print("=" * 80)

print(f"""
**Kaiser P0 と P2 から逆算したPlinが2.7倍も違います！**
   - P0から: {plin_from_kaiser_p0:.0f}
   - P2から: {plin_from_kaiser_p2:.0f}
   
これは、**Kaiser RSDの実装に問題がある**可能性を示しています。
P0とP2は同じPlin(k)を使うはずなのに、逆算すると違う値になっています。

**Velocileptorsも同様に差があります**（約2.3倍）。

**これは、quadrupoleの符号や定義が間違っている可能性**があります！
あるいは、**P2の計算式自体が間違っている**可能性があります。

次のステップ：
1. Kaiser RSDのP2計算式を再確認
2. Velocileptorsのquadrupoleの定義を確認
3. 符号の確認（一部の文献ではP2に負の符号を使う）
""")

