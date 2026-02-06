#!/usr/bin/env python3
"""
Kaiser RSDとVelocileptorsのk-binsを比較
"""
import numpy as np
import os

print("=" * 80)
print("k-binsの比較")
print("=" * 80)

# Kaiser RSDのk-bins
kaiser_dataset = './projects/mg_hsc_y3/data/kaiser_rsd_data.dataset'
with open(kaiser_dataset, 'r') as f:
    for line in f:
        if 'k_bins_file' in line:
            k_bins_file = line.split('=')[1].strip()
            base_dir = os.path.dirname(kaiser_dataset)
            k_bins_file = os.path.join(base_dir, k_bins_file)
            print(f"\nKaiser RSDのk-binsファイル: {k_bins_file}")
            
            k_bins_kaiser = np.loadtxt(k_bins_file)
            print(f"   データ形状: {k_bins_kaiser.shape}")
            
            if k_bins_kaiser.ndim == 2:
                # 2列形式の場合、最後の列を使用
                k_bins_kaiser = k_bins_kaiser[:, -1]
            
            print(f"   k-bins数: {len(k_bins_kaiser)}")
            print(f"   最初の10個:")
            for i in range(min(10, len(k_bins_kaiser))):
                print(f"      k[{i}] = {k_bins_kaiser[i]:.6f} h/Mpc")
            
            break

# Velocileptorsのk-bins
velo_data = np.loadtxt('./projects/mg_hsc_y3/data/velocileptors_pk_multipoles.txt')
k_bins_velo = velo_data[:, 0]

print(f"\nVelocileptorsのk-bins:")
print(f"   k-bins数: {len(k_bins_velo)}")
print(f"   最初の10個:")
for i in range(min(10, len(k_bins_velo))):
    print(f"      k[{i}] = {k_bins_velo[i]:.6f} h/Mpc")

print("\n" + "=" * 80)
print("📊 比較:")
print("=" * 80)

n_common = min(len(k_bins_kaiser), len(k_bins_velo))
print(f"\n共通点数: {n_common}")
print(f"\n{'Index':<8} {'Kaiser k':<15} {'Velo k':<15} {'差':<15} {'相対差 [%]':<15}")
print("-" * 80)

for i in range(min(15, n_common)):
    diff = k_bins_velo[i] - k_bins_kaiser[i]
    rel_diff = abs(diff) / k_bins_kaiser[i] * 100 if k_bins_kaiser[i] != 0 else 0
    print(f"{i:<8} {k_bins_kaiser[i]:<15.6f} {k_bins_velo[i]:<15.6f} {diff:<15.6e} {rel_diff:<15.2f}")

print("\n" + "=" * 80)
print("🎯 結論:")
print("=" * 80)

if n_common > 0:
    max_rel_diff = max([abs(k_bins_velo[i] - k_bins_kaiser[i]) / k_bins_kaiser[i] * 100 
                        for i in range(n_common) if k_bins_kaiser[i] != 0])
    print(f"最大相対差: {max_rel_diff:.4f}%")
    
    if max_rel_diff < 0.01:
        print("✅ k-binsはほぼ完全に一致しています")
    elif max_rel_diff < 1.0:
        print(f"⚠️  k-binsに{max_rel_diff:.2f}%の差があります（許容範囲内かもしれません）")
    else:
        print(f"❌ k-binsに{max_rel_diff:.2f}%の大きな差があります！")
else:
    print("❌ 共通点がありません")












