#!/usr/bin/env python3
"""
Notebookと同じパラメータでKaiser公式を計算
"""
import numpy as np

print("=" * 80)
print("Notebookと同じパラメータでのKaiser計算")
print("=" * 80)

# Notebookのパラメータ
bias_notebook = 0.7
b1_notebook = bias_notebook + 1.0  # 1.7
f = 0.8076  # notebook uses z=0.8

# 私たちのパラメータ
bias_ours = 1.234
b1_ours = bias_ours + 1.0  # 2.234
f_ours = 0.790178  # z=0.5893

print("\n📊 Notebook (z=0.8):")
print(f"   bias = {bias_notebook}")
print(f"   (1+bias) = {b1_notebook}")
print(f"   f = {f}")

beta_notebook = f / b1_notebook
P0_coeff_notebook = (1 + 2/3*beta_notebook + 1/5*beta_notebook**2) * b1_notebook**2
P2_coeff_notebook = (4/3*beta_notebook + 4/7*beta_notebook**2) * b1_notebook**2

print(f"   β = f/(1+bias) = {beta_notebook:.4f}")
print(f"   P0係数 = {P0_coeff_notebook:.4f}")
print(f"   P2係数 = {P2_coeff_notebook:.4f}")

print("\n📊 私たちの実装 (z=0.5893):")
print(f"   bias = {bias_ours}")
print(f"   (1+bias) = {b1_ours}")
print(f"   f = {f_ours}")

beta_ours = f_ours / b1_ours
P0_coeff_ours = (1 + 2/3*beta_ours + 1/5*beta_ours**2) * b1_ours**2
P2_coeff_ours = (4/3*beta_ours + 4/7*beta_ours**2) * b1_ours**2

print(f"   β = f/(1+bias) = {beta_ours:.4f}")
print(f"   P0係数 = {P0_coeff_ours:.4f}")
print(f"   P2係数 = {P2_coeff_ours:.4f}")

print("\n" + "=" * 80)
print("🎯 結論")
print("=" * 80)

print(f"""
**係数の違い:**
   P0: Notebook={P0_coeff_notebook:.4f}, 私たち={P0_coeff_ours:.4f}, 比率={P0_coeff_ours/P0_coeff_notebook:.4f}
   P2: Notebook={P2_coeff_notebook:.4f}, 私たち={P2_coeff_ours:.4f}, 比率={P2_coeff_ours/P2_coeff_notebook:.4f}

**biasの値が違うため、係数も違います。**

しかし、実際の比較では：
   - 両方とも同じ bias=1.234 を使うべき
   - 両方とも同じ z=0.5893 を使うべき
   - 両方とも同じ線形P(k)を使うべき

**次の確認：**
   Velocileptorsが実際に bias=1.234 を使っているか確認
""")












