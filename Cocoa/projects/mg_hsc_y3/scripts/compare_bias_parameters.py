#!/usr/bin/env python3
"""
Kaiser RSDとVelociléptorsのバイアスパラメータを比較
"""
import numpy as np

print("=" * 80)
print("バイアスパラメータの比較：Kaiser RSD vs Velocileptors")
print("=" * 80)

# 共通パラメータ
sigma8_z = 0.608835  # z=0.5893 での sigma8
f_z = 0.790178       # z=0.5893 での f

print("\n📊 **共通パラメータ (z=0.5893):**")
print(f"   σ8(z) = {sigma8_z}")
print(f"   f(z)  = {f_z}")

print("\n" + "=" * 80)
print("1️⃣  Kaiser RSD")
print("=" * 80)

# Kaiser RSD
b_kaiser = 1.234
print(f"\n🔹 **バイアスパラメータ:**")
print(f"   hsc_B_KAISER_RSD_1 = {b_kaiser}")
print(f"   これは b1 = {b_kaiser} に相当")

print(f"\n🔹 **Kaiser効果の式:**")
print(f"   β = f/b = {f_z}/{b_kaiser} = {f_z/b_kaiser:.4f}")
print(f"""
   モノポール (ℓ=0):
   P₀(k) = b² · Pₘ(k) · [1 + (2/3)β + (1/5)β²]
         = {b_kaiser}² · Pₘ(k) · [1 + (2/3)·{f_z/b_kaiser:.4f} + (1/5)·{(f_z/b_kaiser)**2:.4f}]
         = {b_kaiser**2:.4f} · Pₘ(k) · [{1 + (2/3)*(f_z/b_kaiser) + (1/5)*((f_z/b_kaiser)**2):.4f}]
         = {b_kaiser**2 * (1 + (2/3)*(f_z/b_kaiser) + (1/5)*((f_z/b_kaiser)**2)):.4f} · Pₘ(k)
   
   クアドルポール (ℓ=2):
   P₂(k) = b² · Pₘ(k) · [(4/3)β + (4/7)β²]
         = {b_kaiser}² · Pₘ(k) · [(4/3)·{f_z/b_kaiser:.4f} + (4/7)·{(f_z/b_kaiser)**2:.4f}]
         = {b_kaiser**2:.4f} · Pₘ(k) · [{(4/3)*(f_z/b_kaiser) + (4/7)*((f_z/b_kaiser)**2):.4f}]
         = {b_kaiser**2 * ((4/3)*(f_z/b_kaiser) + (4/7)*((f_z/b_kaiser)**2)):.4f} · Pₘ(k)
""")

print("=" * 80)
print("2️⃣  Velocileptors")
print("=" * 80)

# Velocileptors (DESIスタイルのパラメータ)
b1plus1_sigma8 = 1.36
b2_sigma8_sq = 0.32
bs_sigma8_sq = -0.192
alpha0 = 10.0
alpha2 = 20.0
alpha4 = -60.0
SN0 = 1800.0
SN2 = -1000.0

# 実際のバイアス値に変換
b1 = (b1plus1_sigma8 / sigma8_z) - 1.0
b2 = b2_sigma8_sq / (sigma8_z ** 2)
bs = bs_sigma8_sq / (sigma8_z ** 2)

print(f"\n🔹 **DESIスタイルのパラメータ:**")
print(f"   (1+b1)·σ8 = {b1plus1_sigma8}")
print(f"   b2·σ8²    = {b2_sigma8_sq}")
print(f"   bs·σ8²    = {bs_sigma8_sq}")
print(f"   α₀        = {alpha0}")
print(f"   α₂        = {alpha2}")
print(f"   α₄        = {alpha4}")
print(f"   SN₀       = {SN0}")
print(f"   SN₂       = {SN2}")

print(f"\n🔹 **変換後の実際のバイアス値:**")
print(f"   b1 = ({b1plus1_sigma8} / {sigma8_z}) - 1 = {b1:.4f}")
print(f"   b2 = {b2_sigma8_sq} / {sigma8_z}² = {b2:.4f}")
print(f"   bs = {bs_sigma8_sq} / {sigma8_z}² = {bs:.4f}")

print(f"\n🔹 **Velocileptorsの式（EPT + 1-loop）:**")
print(f"""
   P(k,μ) = 複雑な展開式（EPT + 1-loop補正）
   
   含まれる項：
   - 線形項: b1
   - 2次バイアス: b2, bs (局所的・非局所的)
   - カウンター項: α₀, α₂, α₄ (k² · Pₘ(k) など)
   - ショットノイズ: SN₀, SN₂
   - RSD効果: f(z)を含む複雑な角度依存性
   
   詳細は Velocileptors の LPT_RSD クラスを参照
""")

print("\n" + "=" * 80)
print("3️⃣  比較と分析")
print("=" * 80)

print(f"\n✅ **線形バイアス b1 の比較:**")
print(f"   Kaiser RSD:     b = {b_kaiser:.4f}")
print(f"   Velocileptors: b1 = {b1:.4f}")
print(f"   → **同じ値！** ✅")

print(f"\n⚠️  **高次項の違い:**")
print(f"""
   Kaiser RSD:
   - シンプルなKaiser効果のみ
   - P(k) ∝ b² · (1 + β·μ²)² · Pₘ(k)
   - 高次バイアスなし
   
   Velocileptors:
   - 2次バイアス: b2 = {b2:.4f}, bs = {bs:.4f}
   - カウンター項: α₀ = {alpha0}, α₂ = {alpha2}, α₄ = {alpha4}
   - ショットノイズ: SN₀ = {SN0}, SN₂ = {SN2}
   - EPT + 1-loop補正
""")

print(f"\n🎯 **結論:**")
print(f"""
1. **線形バイアス b1 は同じ** ({b1:.4f}) ← 正しく設定されている ✅

2. **しかし、モデルが根本的に違う:**
   - Kaiser RSD: 線形理論 + シンプルなKaiser効果
   - Velocileptors: EPT + 1-loop + 高次バイアス + カウンター項
   
3. **120%の差の原因:**
   - Velocileptorsの高次項 (b2, bs, α₀, α₂, α₄, SN) が大きく寄与
   - 特にカウンター項とショットノイズが効いている可能性
   
4. **公平な比較のためには:**
   - Velocileptorsの高次項を0に設定: b2=bs=α₀=α₂=α₄=SN₀=SN₂=0
   - または、Kaiser RSDに高次項を追加（現実的ではない）
""")

print("\n" + "=" * 80)
print("📋 次のステップ")
print("=" * 80)
print("""
1. Velocileptorsで高次項を0に設定して再計算
2. 両方を線形Kaiser効果のみで比較
3. どの項が最も効いているか特定
""")












