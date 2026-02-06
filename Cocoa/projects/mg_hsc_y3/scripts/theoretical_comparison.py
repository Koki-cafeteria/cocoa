#!/usr/bin/env python3
"""
Kaiser RSDとVelociléptorsの理論的な計算を比較
"""
import numpy as np

print("=" * 80)
print("理論的な計算比較")
print("=" * 80)

# 共通パラメータ
b = 1.234
f = 0.790178
beta = f / b

print(f"\n📊 パラメータ:")
print(f"   b1 = {b}")
print(f"   f(z) = {f}")
print(f"   β = f/b = {beta:.4f}")

print("\n" + "=" * 80)
print("Kaiser RSD の理論式")
print("=" * 80)

# Kaiser multipoles (exact)
kaiser_monopole_coeff = 1 + (2.0/3.0) * beta + (1.0/5.0) * beta * beta
kaiser_quadrupole_coeff = (4.0/3.0) * beta + (4.0/7.0) * beta * beta

print(f"\nモノポール (ℓ=0):")
print(f"   P₀(k) = b² · Pₘ(k) · [1 + 2/3·β + 1/5·β²]")
print(f"   係数 = 1 + 2/3·{beta:.4f} + 1/5·{beta**2:.4f}")
print(f"        = 1 + {(2.0/3.0)*beta:.4f} + {(1.0/5.0)*beta**2:.4f}")
print(f"        = {kaiser_monopole_coeff:.4f}")
print(f"   P₀(k) = {b**2:.4f} · {kaiser_monopole_coeff:.4f} · Pₘ(k)")
print(f"   P₀(k) = {b**2 * kaiser_monopole_coeff:.4f} · Pₘ(k)")

print(f"\nクアドルポール (ℓ=2):")
print(f"   P₂(k) = b² · Pₘ(k) · [4/3·β + 4/7·β²]")
print(f"   係数 = 4/3·{beta:.4f} + 4/7·{beta**2:.4f}")
print(f"        = {(4.0/3.0)*beta:.4f} + {(4.0/7.0)*beta**2:.4f}")
print(f"        = {kaiser_quadrupole_coeff:.4f}")
print(f"   P₂(k) = {b**2:.4f} · {kaiser_quadrupole_coeff:.4f} · Pₘ(k)")
print(f"   P₂(k) = {b**2 * kaiser_quadrupole_coeff:.4f} · Pₘ(k)")

print("\n" + "=" * 80)
print("Velocileptors の期待値（線形のみの場合）")
print("=" * 80)

print(f"""
Velocileptorsは EPT (Effective Field Theory) を使用していますが、
高次項を全て0に設定した場合、Kaiser RSDと**ほぼ同じ結果**になるはずです。

ただし、Velocileptorsには追加の項があります：
1. 1-loop補正（小さいが0ではない）
2. IR resummation（赤外発散の繰り込み）
3. 数値的な差異（FFTなど）

これらの効果は通常 **数%以内** です。

しかし、現在の差は **3-20倍** なので、何か根本的な問題があります。
""")

print("=" * 80)
print("可能性のある問題")
print("=" * 80)

print(f"""
1. **Velocileptorsのb1の適用方法が違う可能性**
   - Kaiser: P(k) ∝ b² · Pₘ(k)
   - Velocileptors: P(k) ∝ ??? 
   
2. **線形パワースペクトル Pₘ(k) の値が違う可能性**
   - 同じk, z値で評価しているか？
   - redshift evolutionの扱いが違う？
   
3. **Velocileptorsの"reduced basis"の意味**
   - compute_redshift_space_power_multipoles(pars, f_z, reduced=True)
   - reduced=True が何か特別な正規化をしている？

4. **k-binsの違い**
   - Kaiser: 40点（BOSS DR12のk_bins）
   - Velocileptors: 30点（logspace）
   - 補間の問題？
""")

print("\n" + "=" * 80)
print("📋 次の確認事項")
print("=" * 80)

print("""
1. ✅ f(z)の値 → 同じ (0.790178)
2. ✅ b1の値 → 同じ (1.234)
3. ✅ 単位 → 両方とも (Mpc/h)³
4. ✅ Kaiser計算式 → 正しい
5. ❓ Velocileptorsのb1適用方法 → **要確認**
6. ❓ 線形P(k)の値 → **要確認**
7. ❓ reduced=True の意味 → **要確認**
""")

