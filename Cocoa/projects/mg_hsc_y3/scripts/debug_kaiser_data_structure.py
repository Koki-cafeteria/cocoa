#!/usr/bin/env python3
"""
Kaiser RSDのデータ構造をデバッグ
"""
import numpy as np

print("=" * 80)
print("Kaiser RSDデータ構造のデバッグ")
print("=" * 80)

# Kaiser RSDのmodelベクトルを読み込む
kaiser_file = './projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector'
data = np.loadtxt(kaiser_file)

print(f"\n📊 ファイル: {kaiser_file}")
print(f"   データサイズ: {data.shape}")
print(f"   データタイプ: {data.dtype}")

# 2列形式: [index, value]
indices = data[:, 0]
values = data[:, 1]

print(f"\n最初の20エントリ:")
print(f"{'Index':<10} {'Value':<20} {'推定される多極子':<20}")
print("-" * 80)
for i in range(min(20, len(data))):
    multipole = "P0 (monopole)" if i % 2 == 0 else "P2 (quadrupole)"
    print(f"{int(indices[i]):<10} {values[i]:<20.2f} {multipole:<20}")

# Nk=40, Nz=1 の場合
Nk = 40
Nz = 1
total_size = Nk * Nz  # 40

print(f"\n📊 データ構造の解析:")
print(f"   Nk = {Nk}, Nz = {Nz}")
print(f"   期待されるサイズ: 2 * Nk * Nz = {2 * Nk * Nz}")
print(f"   実際のサイズ: {len(values)}")

# 現在のdirect_comparison.pyの読み方
P0_direct = values[::2]  # 偶数インデックス
P2_direct = values[1::2]  # 奇数インデックス

print(f"\n現在の読み方（direct_comparison.py）:")
print(f"   P0[0] (index=0) = {P0_direct[0]:.2f}")
print(f"   P2[0] (index=1) = {P2_direct[0]:.2f}")
print(f"   P0[1] (index=2) = {P0_direct[1]:.2f}")
print(f"   P2[1] (index=3) = {P2_direct[1]:.2f}")

# C++の意図する構造
# [P_0(k1,z1), ..., P_0(kNk,zNz), P_2(k1,z1), ..., P_2(kNk,zNz)]
P0_intended = values[:total_size]
P2_intended = values[total_size:2*total_size]

print(f"\nC++が意図する構造:")
print(f"   P0[0] (indices 0-39の最初) = {P0_intended[0]:.2f}")
print(f"   P2[0] (indices 40-79の最初) = {P2_intended[0]:.2f}")

print("\n" + "=" * 80)
print("🎯 結論")
print("=" * 80)

if abs(P0_direct[0] - P0_intended[0]) < 1:
    print("✅ P0の読み方は正しい")
else:
    print(f"❌ P0の読み方が間違っている！")
    print(f"   direct: {P0_direct[0]:.2f}")
    print(f"   intended: {P0_intended[0]:.2f}")

if abs(P2_direct[0] - P2_intended[0]) < 1:
    print("✅ P2の読み方は正しい")
else:
    print(f"❌ P2の読み方が間違っている！")
    print(f"   direct: {P2_direct[0]:.2f}")
    print(f"   intended: {P2_intended[0]:.2f}")












