#!/usr/bin/env python3
"""
Kaiser RSD vs Velocileptors 比較（単位系修正後）
"""
import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("Kaiser RSD vs Velocileptors 比較（単位系修正後）")
print("="*80)

# ============================================================================
# 1. データを読み込む
# ============================================================================

# Kaiser RSD
kaiser_file = '/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector'
kaiser_data = np.loadtxt(kaiser_file)
values = kaiser_data[:, 1]

# BOSS DR12のデータ構造: Nk=40, Nz=1
Nk = 40
Nz = 1
total_size = Nk * Nz

# P0とP2を正しく分割
P0_kaiser = values[:total_size]
P2_kaiser = values[total_size:2*total_size]

# k-bins
k_bins_file = '/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/boss_dr12_cmass_k_bins.txt'
k_bins_data = np.loadtxt(k_bins_file)
k_h_mpc = k_bins_data[:, 1]  # h/Mpc 単位

print(f"\n✅ Kaiser RSD loaded:")
print(f"   Nk = {Nk}, Nz = {Nz}")
print(f"   k range: {k_h_mpc[0]:.6f} - {k_h_mpc[-1]:.6f} h/Mpc")
print(f"   P0[0] = {P0_kaiser[0]:.2f}")
print(f"   P2[0] = {P2_kaiser[0]:.2f}")

# Velocileptors
velo_file = '/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/velocileptors_pk_multipoles.txt'
velo_data = np.loadtxt(velo_file)
k_velo = velo_data[:, 0]
P0_velo = velo_data[:, 1]
P2_velo = velo_data[:, 2]

print(f"\n✅ Velocileptors loaded:")
print(f"   Nk = {len(k_velo)}")
print(f"   k range: {k_velo[0]:.6f} - {k_velo[-1]:.6f} h/Mpc")
print(f"   P0[0] = {P0_velo[0]:.2f}")
print(f"   P2[0] = {P2_velo[0]:.2f}")

# ============================================================================
# 2. 比較統計
# ============================================================================

# k範囲を合わせる（最小の共通範囲）
k_min_common = max(k_h_mpc[0], k_velo[0])
k_max_common = min(k_h_mpc[-1], k_velo[-1])

mask_kaiser = (k_h_mpc >= k_min_common) & (k_h_mpc <= k_max_common)
mask_velo = (k_velo >= k_min_common) & (k_velo <= k_max_common)

# 補間して同じk点で比較
from scipy.interpolate import interp1d
P0_kaiser_interp = interp1d(k_h_mpc, P0_kaiser, kind='cubic', fill_value='extrapolate')
P2_kaiser_interp = interp1d(k_h_mpc, P2_kaiser, kind='cubic', fill_value='extrapolate')

P0_kaiser_at_velo_k = P0_kaiser_interp(k_velo)
P2_kaiser_at_velo_k = P2_kaiser_interp(k_velo)

# 相対差
rel_diff_P0 = (P0_kaiser_at_velo_k - P0_velo) / P0_velo * 100
rel_diff_P2 = (P2_kaiser_at_velo_k - P2_velo) / P2_velo * 100

print(f"\n{'='*80}")
print(f"統計情報（単位系修正後）")
print(f"{'='*80}")
print(f"\nk·P0 相対差:")
print(f"  平均: {np.mean(rel_diff_P0):.2f}%")
print(f"  RMS:  {np.sqrt(np.mean(rel_diff_P0**2)):.2f}%")
print(f"  最小: {np.min(rel_diff_P0):.2f}%")
print(f"  最大: {np.max(rel_diff_P0):.2f}%")

print(f"\nk·P2 相対差:")
print(f"  平均: {np.mean(rel_diff_P2):.2f}%")
print(f"  RMS:  {np.sqrt(np.mean(rel_diff_P2**2)):.2f}%")
print(f"  最小: {np.min(rel_diff_P2):.2f}%")
print(f"  最大: {np.max(rel_diff_P2):.2f}%")

# ============================================================================
# 3. プロット
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# P0の比較
ax = axes[0, 0]
ax.plot(k_h_mpc, k_h_mpc * P0_kaiser, 'b-', label='Kaiser RSD', linewidth=2)
ax.plot(k_velo, k_velo * P0_velo, 'r--', label='Velocileptors', linewidth=2)
ax.set_xlabel('k [h/Mpc]', fontsize=12)
ax.set_ylabel('k·P₀(k) [h⁻²Mpc²]', fontsize=12)
ax.set_title('Monopole (P₀)', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# P2の比較
ax = axes[0, 1]
ax.plot(k_h_mpc, k_h_mpc * P2_kaiser, 'b-', label='Kaiser RSD', linewidth=2)
ax.plot(k_velo, k_velo * P2_velo, 'r--', label='Velocileptors', linewidth=2)
ax.set_xlabel('k [h/Mpc]', fontsize=12)
ax.set_ylabel('k·P₂(k) [h⁻²Mpc²]', fontsize=12)
ax.set_title('Quadrupole (P₂)', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# P0の相対差
ax = axes[1, 0]
ax.plot(k_velo, rel_diff_P0, 'g-', linewidth=2)
ax.axhline(0, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('k [h/Mpc]', fontsize=12)
ax.set_ylabel('相対差 [%]', fontsize=12)
ax.set_title(f'P₀ 相対差 (平均: {np.mean(rel_diff_P0):.2f}%)', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

# P2の相対差
ax = axes[1, 1]
ax.plot(k_velo, rel_diff_P2, 'g-', linewidth=2)
ax.axhline(0, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('k [h/Mpc]', fontsize=12)
ax.set_ylabel('相対差 [%]', fontsize=12)
ax.set_title(f'P₂ 相対差 (平均: {np.mean(rel_diff_P2):.2f}%)', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = '/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/plots/kaiser_velocileptors_comparison_unit_fix.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\n✅ プロットを保存しました: {output_file}")

plt.show()

print(f"\n{'='*80}")
print("完了！")
print(f"{'='*80}")







