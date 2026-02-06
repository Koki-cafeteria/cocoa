#!/usr/bin/env python3
"""
Kaiser RSDとVelocileptors（線形のみ）の値を直接比較
"""
import numpy as np
import os

print("=" * 80)
print("直接比較：Kaiser RSD vs Velocileptors（線形バイアスのみ）")
print("=" * 80)

# Velocileptors（線形のみ）のデータ
velo_data = np.loadtxt('./projects/mg_hsc_y3/data/velocileptors_pk_multipoles.txt')
k_velo = velo_data[:, 0]
P0_velo = velo_data[:, 1]
P2_velo = velo_data[:, 2]

print("\n📊 Velocileptors（線形のみ）:")
print(f"   k[0] = {k_velo[0]:.6f} h/Mpc")
print(f"   P0[0] = {P0_velo[0]:.2f}")
print(f"   P2[0] = {P2_velo[0]:.2f}")
print(f"\n   k[10] = {k_velo[10]:.6f} h/Mpc")
print(f"   P0[10] = {P0_velo[10]:.2f}")
print(f"   P2[10] = {P2_velo[10]:.2f}")

# Kaiser RSDのデータ（modelve ctor形式）Kaiser RSDのフォーマットを確認
kaiser_file = './projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector'

if os.path.exists(kaiser_file):
    kaiser_data = np.loadtxt(kaiser_file)
    print(f"\n📊 Kaiser RSD:")
    print(f"   データサイズ: {kaiser_data.shape}")
    print(f"   最初の10行:")
    for i in range(min(10, len(kaiser_data))):
        print(f"   {i}: {kaiser_data[i]}")
    
    # modelベクトルは2列形式: [index, value]
    # C++の出力構造: [P0(k1), ..., P0(kNk), P2(k1), ..., P2(kNk)]
    values = kaiser_data[:, 1]  # 2列目が実際の値
    
    # Kaiser RSDの設定（datasetファイルから）
    Nk = 40  # k-binsの数
    Nz = 1   # z-binsの数（BOSS DR12は単一赤方偏移）
    total_size = Nk * Nz
    
    # 正しいデータ構造で読み込む
    P0_kaiser = values[:total_size]  # 最初のNk*Nzエントリ = モノポール
    P2_kaiser = values[total_size:2*total_size]  # 次のNk*Nzエントリ = クアドルポール
    
    print(f"\n   P0[0] = {P0_kaiser[0]:.2f}")
    print(f"   P2[0] = {P2_kaiser[0]:.2f}")
    print(f"\n   P0[10] = {P0_kaiser[10]:.2f}")
    print(f"   P2[10] = {P2_kaiser[10]:.2f}")
    
    print("\n" + "=" * 80)
    print("📊 比較（最初の数点）:")
    print("=" * 80)
    print(f"{'Index':<8} {'P0_Kaiser':<15} {'P0_Velo':<15} {'比率':<10} {'差 [%]':<10}")
    print("-" * 80)
    for i in range(min(15, len(P0_kaiser), len(P0_velo))):
        ratio = P0_velo[i] / P0_kaiser[i] if P0_kaiser[i] != 0 else 0
        diff_pct = abs(P0_velo[i] - P0_kaiser[i]) / P0_kaiser[i] * 100 if P0_kaiser[i] != 0 else 0
        print(f"{i:<8} {P0_kaiser[i]:<15.2f} {P0_velo[i]:<15.2f} {ratio:<10.4f} {diff_pct:<10.1f}")
    
    print("\n" + "=" * 80)
    print("📊 クアドルポール比較:")
    print("=" * 80)
    print(f"{'Index':<8} {'P2_Kaiser':<15} {'P2_Velo':<15} {'比率':<10} {'差 [%]':<10}")
    print("-" * 80)
    for i in range(min(15, len(P2_kaiser), len(P2_velo))):
        ratio = P2_velo[i] / P2_kaiser[i] if P2_kaiser[i] != 0 else 0
        diff_pct = abs(P2_velo[i] - P2_kaiser[i]) / P2_kaiser[i] * 100 if P2_kaiser[i] != 0 else 0
        print(f"{i:<8} {P2_kaiser[i]:<15.2f} {P2_velo[i]:<15.2f} {ratio:<10.4f} {diff_pct:<10.1f}")
    
    print("\n" + "=" * 80)
    print("🎯 統計:")
    print("=" * 80)
    # 共通のk範囲で比較
    n_common = min(len(P0_kaiser), len(P0_velo))
    P0_diff_pct = np.abs(P0_velo[:n_common] - P0_kaiser[:n_common]) / P0_kaiser[:n_common] * 100
    P2_diff_pct = np.abs(P2_velo[:n_common] - P2_kaiser[:n_common]) / P2_kaiser[:n_common] * 100
    
    print(f"P0相対差: 平均 {np.mean(P0_diff_pct):.2f}%, 最大 {np.max(P0_diff_pct):.2f}%")
    print(f"P2相対差: 平均 {np.mean(P2_diff_pct):.2f}%, 最大 {np.max(P2_diff_pct):.2f}%")
    
else:
    print(f"\n⚠️  Kaiser RSDファイルが見つかりません: {kaiser_file}")

