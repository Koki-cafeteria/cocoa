"""
COCOAで計算したmonopoleとquadrupoleをBOSS DR12の観測データと比較してプロット
Notebookにコピペして使用
"""

import numpy as np
import matplotlib.pyplot as plt

# ===== 設定 =====
# ファイルパス（必要に応じて変更）
data_dir = "./projects/mg_hsc_y3/data/"
theory_file = data_dir + "mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector"
obs_file = data_dir + "boss_dr12_cmass_datavector.txt"
boss_original_file = data_dir + "Power_Spectrum_cmass_ngc_v5.txt"  # kの値を取得するため

# データ読み込み
# COCOA理論値（chunked形式: 最初の40個がmonopole、次の40個がquadrupole）
theory_data = np.loadtxt(theory_file)
theory_indices = theory_data[:, 0].astype(int)
theory_values = theory_data[:, 1]

# BOSS DR12観測データ（chunked形式: 最初の40個がmonopole、次の40個がquadrupole）
obs_data = np.loadtxt(obs_file)
obs_indices = obs_data[:, 0].astype(int)
obs_values = obs_data[:, 1]

# kの値を取得（元のBOSS DR12ファイルから）
k_eff = []
with open(boss_original_file, 'r') as f:
    for line in f:
        if line.startswith('#') or len(line.strip()) == 0:
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                k_eff.append(float(parts[1]))  # k-eff列
            except ValueError:
                continue

k_eff = np.array(k_eff[:40])  # 最初の40個のkビン

# データをmonopoleとquadrupoleに分離
Nk = 40
theory_p0 = theory_values[:Nk]
theory_p2 = theory_values[Nk:]
obs_p0 = obs_values[:Nk]
obs_p2 = obs_values[Nk:]

# プロット
fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Monopole (P0)
ax = axes[0]
ax.plot(k_eff, theory_p0, 'b-', label='COCOA Theory', linewidth=2, marker='o', markersize=4)
ax.plot(k_eff, obs_p0, 'r--', label='BOSS DR12 Obs', linewidth=2, marker='s', markersize=4, alpha=0.7)
ax.set_ylabel(r'$P_0(k)$ [$(h^{-1}\mathrm{Mpc})^3$]', fontsize=12)
ax.set_title('Monopole (ℓ=0) Comparison', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# Quadrupole (P2)
ax = axes[1]
ax.plot(k_eff, theory_p2, 'b-', label='COCOA Theory', linewidth=2, marker='o', markersize=4)
ax.plot(k_eff, obs_p2, 'r--', label='BOSS DR12 Obs', linewidth=2, marker='s', markersize=4, alpha=0.7)
ax.set_xlabel(r'$k$ [$h\,\mathrm{Mpc}^{-1}$]', fontsize=12)
ax.set_ylabel(r'$P_2(k)$ [$(h^{-1}\mathrm{Mpc})^3$]', fontsize=12)
ax.set_title('Quadrupole (ℓ=2) Comparison', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

plt.tight_layout()
plt.show()

# 残差プロット
fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Monopole残差
ax = axes[0]
residual_p0 = (theory_p0 - obs_p0) / obs_p0 * 100  # パーセント残差
ax.plot(k_eff, residual_p0, 'b-', linewidth=2, marker='o', markersize=4)
ax.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax.set_ylabel(r'Residual [%]', fontsize=12)
ax.set_title('Monopole (ℓ=0) Residual: (Theory - Obs) / Obs × 100', fontsize=14)
ax.grid(True, alpha=0.3)

# Quadrupole残差
ax = axes[1]
residual_p2 = (theory_p2 - obs_p2) / obs_p2 * 100  # パーセント残差
ax.plot(k_eff, residual_p2, 'b-', linewidth=2, marker='o', markersize=4)
ax.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel(r'$k$ [$h\,\mathrm{Mpc}^{-1}$]', fontsize=12)
ax.set_ylabel(r'Residual [%]', fontsize=12)
ax.set_title('Quadrupole (ℓ=2) Residual: (Theory - Obs) / Obs × 100', fontsize=14)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# 統計情報を表示
print("=== 統計情報 ===")
print(f"\nMonopole (P0):")
print(f"  Theory range: [{theory_p0.min():.2e}, {theory_p0.max():.2e}]")
print(f"  Obs range: [{obs_p0.min():.2e}, {obs_p0.max():.2e}]")
print(f"  Mean residual: {residual_p0.mean():.2f}%")
print(f"  RMS residual: {np.sqrt(np.mean(residual_p0**2)):.2f}%")

print(f"\nQuadrupole (P2):")
print(f"  Theory range: [{theory_p2.min():.2e}, {theory_p2.max():.2e}]")
print(f"  Obs range: [{obs_p2.min():.2e}, {obs_p2.max():.2e}]")
print(f"  Mean residual: {residual_p2.mean():.2f}%")
print(f"  RMS residual: {np.sqrt(np.mean(residual_p2**2)):.2f}%")

