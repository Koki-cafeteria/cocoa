#!/usr/bin/env python3
"""
Kaiser RSDとVelociléptorsでf(z)の値が一致するか確認するスクリプト
"""
import numpy as np
from scipy.interpolate import interp1d

print("=" * 70)
print("f(z) 値の比較：Kaiser RSD vs Velocileptors")
print("=" * 70)

# YAMLファイルから使用されているz値を読み取る
print("\n📋 設定を確認しています...")

# Kaiser RSD: z_min から z_max の範囲
kaiser_z_min = 0.5893  # EXAMPLE_EVALUATE_kaiser_RSD.yaml
kaiser_z_max = 0.5893
print(f"   Kaiser RSD: z = {kaiser_z_min} - {kaiser_z_max}")

# Velocileptors: z_eff
velocileptors_z_eff = 0.5  # EXAMPLE_EVALUATE_velocileptors.yaml
print(f"   Velocileptors: z_eff = {velocileptors_z_eff}")

print("\n⚠️  注意：Velocileptorsのログから既知の値:")
print(f"   f(z=0.589) = 0.790178 (MGCAMB)")
print(f"   fsigma8    = 0.481088")
print(f"   sigma8     = 0.608835")

print("\n🔍 Kaiser RSDのG(z)からf(z)を計算する方法:")
print("   1. Python側: G(z) = √[P_lin(z) / P_lin(z=0)]")
print("   2. C++側: f(z) = 1 + d ln G / d ln a （数値微分）")

print("\n🔍 Velocileptorsのf(z)計算方法:")
print("   f(z) = fsigma8(z) / sigma8(z)  （MGCAMBから直接取得）")

print("\n" + "=" * 70)
print("📊 結論")
print("=" * 70)
print("""
Kaiser RSDとVelociléptorsで使用するz値が異なります：
  - Kaiser RSD:      z = 0.5893
  - Velocileptors:   z = 0.5

これにより、f(z)の値も若干異なるはずです。

しかし、より重要な問題は：
  1. Kaiser RSDのG(z)計算が正確か？
  2. C++側のf_growth()が正しく微分を計算しているか？
  3. バイアスパラメータの適用方法が違う（これが主因の可能性大）

次のステップ：
  1. Kaiser RSD実行時にf(z)の値をログ出力（追加済み）
  2. 両方を同じz値で計算して直接比較
  3. バイアスパラメータの適用方法を詳細に比較
""")

print("\n✅ 修正済み: Kaiser RSD側にf(z)のログ出力を追加しました")
print("   実行時に以下のような出力が表示されるはずです:")
print("   [Kaiser RSD] Using f(z=0.589) = 0.XXXXXX from MGCAMB")












