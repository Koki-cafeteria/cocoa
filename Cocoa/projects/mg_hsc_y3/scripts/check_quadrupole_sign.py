#!/usr/bin/env python3
"""
Quadrupoleの符号と定義を確認
"""
import numpy as np

print("=" * 80)
print("Quadrupole (P2) の符号と定義の確認")
print("=" * 80)

# パラメータ
b1 = 1.234
f = 0.790178

print(f"\nパラメータ:")
print(f"   b1 = {b1}")
print(f"   f = {f}")

# Notebookの定義（正しいKaiser公式）
beta_notebook = f / (1 + b1)
P2_coeff_notebook = (4.0/3.0) * beta_notebook + (4.0/7.0) * beta_notebook**2
P2_factor_notebook = ((1 + b1)**2) * P2_coeff_notebook

print("\n" + "=" * 80)
print("Notebook (正しいKaiser公式)")
print("=" * 80)
print(f"β = f/(1+b1) = {f}/(1+{b1}) = {beta_notebook:.4f}")
print(f"P2係数 = (4/3)β + (4/7)β²")
print(f"       = (4/3)×{beta_notebook:.4f} + (4/7)×{beta_notebook**2:.4f}")
print(f"       = {P2_coeff_notebook:.4f}")
print(f"P2 = (1+b1)² × {P2_coeff_notebook:.4f} × Plin")
print(f"   = {(1+b1)**2:.4f} × {P2_coeff_notebook:.4f} × Plin")
print(f"   = {P2_factor_notebook:.4f} × Plin")

# 修正後のC++実装（現在の）
beta_cpp_corrected = f / (1 + b1)
P2_coeff_cpp = (4.0/3.0) * beta_cpp_corrected + (4.0/7.0) * beta_cpp_corrected**2
P2_factor_cpp = ((1 + b1)**2) * P2_coeff_cpp

print("\n" + "=" * 80)
print("修正後のC++実装")
print("=" * 80)
print(f"β = f/(1+b1) = {beta_cpp_corrected:.4f}")
print(f"P2係数 = (4/3)β + (4/7)β²")
print(f"       = {P2_coeff_cpp:.4f}")
print(f"P2 = (1+b1)² × {P2_coeff_cpp:.4f} × Plin")
print(f"   = {P2_factor_cpp:.4f} × Plin")

# 直接比較の値
kaiser_p2_actual = 86516
velo_p2_actual = 35133

print("\n" + "=" * 80)
print("実際の値（k=0.007424での直接比較）")
print("=" * 80)
print(f"Kaiser P2 = {kaiser_p2_actual}")
print(f"Velo P2   = {velo_p2_actual}")
print(f"比率 = {kaiser_p2_actual / velo_p2_actual:.4f}")

# 理論的に期待される比率
print("\n" + "=" * 80)
print("🔍 問題の診断")
print("=" * 80)

print(f"""
**理論的には同じ係数を使うべき:**
   Notebook: {P2_factor_notebook:.4f}
   C++:      {P2_factor_cpp:.4f}
   → 係数は一致 ✅

**しかし実際の比率:**
   Kaiser / Velo = {kaiser_p2_actual / velo_p2_actual:.4f}
   
**可能性:**
1. **Quadrupoleの符号の定義が違う？**
   - 一部の文献では P2 = -[...] と負の符号を使う
   - しかし、比率は正なので符号の問題ではない

2. **線形P(k)の値が違う**
   - Kaiser: Plin ≈ {kaiser_p2_actual / P2_factor_cpp:.0f}
   - Velo:  Plin ≈ {velo_p2_actual / P2_factor_notebook:.0f}
   - 比率: {(kaiser_p2_actual/P2_factor_cpp) / (velo_p2_actual/P2_factor_notebook):.4f}
   
3. **redshift evolutionの扱いが違う**
   - Kaiser: z=0.5893のP(k)を使用
   - Velo: z=0.5のP(k)を使用している可能性

4. **1-loop補正の影響**
   - Velocileptorsは高次項を0にしても1-loop補正が入る
   - これがP2の値を変える可能性

**次のステップ:**
- 両方で同じk, z値での線形P(k)の値を出力して比較
- Velocileptorsの実際のz値を確認（0.5 vs 0.5893）
""")












