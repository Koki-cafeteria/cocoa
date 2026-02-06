#!/usr/bin/env python3
"""
Kaiser RSD vs Velocileptors の比較プロット

使用方法:
1. Kaiser RSDを実行: cobaya-run EXAMPLE_EVALUATE_kaiser_RSD.yaml
2. Velocileptorsを実行: cobaya-run EXAMPLE_EVALUATE_velocileptors.yaml
3. このスクリプトを実行: python scripts/compare_kaiser_velocileptors.py
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from getdist import IniFile

def read_dataset_config(dataset_file):
    """データセットファイルから設定を読み込む"""
    config = {}
    with open(dataset_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                try:
                    if '.' in value:
                        config[key] = float(value)
                    else:
                        config[key] = int(value)
                except ValueError:
                    config[key] = value
    return config

def load_kaiser_data(kaiser_file, dataset_file):
    """
    Kaiser RSDのデータを読み込む
    
    データ構造:
    [P0(k1,z1), P0(k1,z2), ..., P0(kNk,zNz), 
     P2(k1,z1), P2(k1,z2), ..., P2(kNk,zNz)]
    """
    # データセット設定を読み込む
    config = read_dataset_config(dataset_file)
    
    # k-binsを読み込む
    if 'k_bins_file' in config:
        k_bins_file = config['k_bins_file']
        if not os.path.isabs(k_bins_file):
            base_dir = os.path.dirname(dataset_file)
            k_bins_file = os.path.join(base_dir, k_bins_file)
        k_bins = np.loadtxt(k_bins_file)
        if k_bins.ndim == 2:
            k_bins = k_bins[:, -1]  # 最後の列をk配列として採用
        Nk = len(k_bins)
    else:
        Nk = config.get('Nk', 100)
        k_min = config.get('k_min', 0.01)
        k_max = config.get('k_max', 1.0)
        k_bins = np.logspace(np.log10(k_min), np.log10(k_max), Nk)
    
    Nz = config.get('Nz', 10)
    z_min = config.get('z_min', 0.1)
    z_max = config.get('z_max', 1.0)
    z_bins = np.linspace(z_min, z_max, Nz)
    
    # Kaiser RSDデータを読み込む
    data = np.loadtxt(kaiser_file)
    if data.ndim == 2:
        data_vector = data[:, 1]  # 2列目がデータ
    else:
        data_vector = data
    
    # データを分割: P0とP2
    total_size = Nk * Nz
    P0_flat = data_vector[:total_size]
    P2_flat = data_vector[total_size:2*total_size]
    
    # 2次元配列に変換 (Nk x Nz)
    P0 = P0_flat.reshape(Nk, Nz)
    P2 = P2_flat.reshape(Nk, Nz)
    
    return k_bins, z_bins, P0, P2

def load_velocileptors_data(velo_file):
    """
    Velocileptorsのデータを読み込む
    
    データ構造: [k, P0, P2, P4]
    """
    data = np.loadtxt(velo_file)
    k = data[:, 0]
    P0 = data[:, 1]
    P2 = data[:, 2]
    P4 = data[:, 3]
    
    return k, P0, P2, P4

def plot_comparison_single_z(k_kaiser, P0_kaiser, P2_kaiser, 
                             k_velo, P0_velo, P2_velo, P4_velo, 
                             z_value, output_file=None):
    """
    Kaiser RSD vs Velocileptors の比較プロット（単一のz値）
    k*P(k) をプロット
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # P0 (Monopole) の比較: k*P0(k) - LINEAR SCALE
    ax = axes[0, 0]
    ax.plot(k_kaiser, k_kaiser * P0_kaiser, 'b-', linewidth=2, label='Kaiser RSD')
    ax.plot(k_velo, k_velo * P0_velo, 'r--', linewidth=2, label='Velocileptors')
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_ylabel(r'$k \cdot P_0(k)$ [$h^{-2}$ Mpc$^2$]', fontsize=12)
    ax.set_title(f'Monopole (k·P0) at z={z_value:.2f}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.25)  # 大スケール（小さなk）重視
    
    # P2 (Quadrupole) の比較: k*P2(k) - LINEAR SCALE
    ax = axes[0, 1]
    ax.plot(k_kaiser, k_kaiser * P2_kaiser, 'b-', linewidth=2, label='Kaiser RSD')
    ax.plot(k_velo, k_velo * P2_velo, 'r--', linewidth=2, label='Velocileptors')
    ax.axhline(0, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_ylabel(r'$k \cdot P_2(k)$ [$h^{-2}$ Mpc$^2$]', fontsize=12)
    ax.set_title(f'Quadrupole (k·P2) at z={z_value:.2f}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.25)  # 大スケール（小さなk）重視
    
    # P4 (Hexadecapole) - Velocileptorsのみ: k*P4(k) - LINEAR SCALE
    ax = axes[1, 0]
    ax.plot(k_velo, k_velo * P4_velo, 'r-', linewidth=2, label='Velocileptors')
    ax.axhline(0, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_ylabel(r'$k \cdot P_4(k)$ [$h^{-2}$ Mpc$^2$]', fontsize=12)
    ax.set_title(f'Hexadecapole (k·P4) at z={z_value:.2f} (Velocileptors only)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.25)  # 大スケール（小さなk）重視
    
    # 相対差: (Velo - Kaiser) / Kaiser (k*P(k)ベース)
    ax = axes[1, 1]
    # k値を補間して比較（k範囲が異なる場合）
    kP0_velo_interp = np.interp(k_kaiser, k_velo, k_velo * P0_velo)
    kP2_velo_interp = np.interp(k_kaiser, k_velo, k_velo * P2_velo)
    
    kP0_kaiser = k_kaiser * P0_kaiser
    kP2_kaiser = k_kaiser * P2_kaiser
    
    rel_diff_P0 = (kP0_velo_interp - kP0_kaiser) / kP0_kaiser * 100
    rel_diff_P2 = (kP2_velo_interp - kP2_kaiser) / np.abs(kP2_kaiser) * 100
    
    ax.plot(k_kaiser, rel_diff_P0, 'b-', linewidth=2, label='k·P0 relative diff')
    ax.plot(k_kaiser, rel_diff_P2, 'r-', linewidth=2, label='k·P2 relative diff')
    ax.axhline(0, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_ylabel('Relative Difference [%]', fontsize=12)
    ax.set_title(f'(Velocileptors - Kaiser) / Kaiser at z={z_value:.2f}', fontsize=14)
    ax.set_xlim(0, 0.25)  # 大スケール重視
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Kaiser RSD vs Velocileptors Comparison (k·P(k))', fontsize=16, y=0.995)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✅ プロットを保存しました: {output_file}")
    
    plt.show()

def plot_multiple_z(k_kaiser, z_bins, P0_kaiser, P2_kaiser,
                    k_velo, P0_velo, P2_velo, output_file=None):
    """
    Kaiser RSD の複数のz値でのプロット（Velocileptorsは単一z）
    k*P(k) をプロット
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(z_bins)))
    
    # P0 (Monopole): k*P0(k) - LINEAR SCALE
    ax = axes[0]
    for i, z in enumerate(z_bins):
        ax.plot(k_kaiser, k_kaiser * P0_kaiser[:, i], color=colors[i], 
                 linewidth=1.5, alpha=0.7, label=f'Kaiser z={z:.2f}')
    ax.plot(k_velo, k_velo * P0_velo, 'r-', linewidth=3, label='Velocileptors', zorder=100)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_ylabel(r'$k \cdot P_0(k)$ [$h^{-2}$ Mpc$^2$]', fontsize=12)
    ax.set_title('Monopole (k·P0)', fontsize=14)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 0.25)  # 大スケール重視
    
    # P2 (Quadrupole): k*P2(k) - LINEAR SCALE
    ax = axes[1]
    for i, z in enumerate(z_bins):
        ax.plot(k_kaiser, k_kaiser * P2_kaiser[:, i], color=colors[i], 
                   linewidth=1.5, alpha=0.7, label=f'Kaiser z={z:.2f}')
    ax.plot(k_velo, k_velo * P2_velo, 'r-', linewidth=3, label='Velocileptors', zorder=100)
    ax.axhline(0, color='k', linestyle=':', alpha=0.5)
    ax.set_xlabel(r'$k$ [$h$ Mpc$^{-1}$]', fontsize=12)
    ax.set_xlim(0, 0.25)  # 大スケール重視
    ax.set_ylabel(r'$k \cdot P_2(k)$ [$h^{-2}$ Mpc$^2$]', fontsize=12)
    ax.set_title('Quadrupole (k·P2)', fontsize=14)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Kaiser RSD (multiple z) vs Velocileptors (single z) - k·P(k)', fontsize=16)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"✅ プロットを保存しました: {output_file}")
    
    plt.show()

def main():
    """メイン処理"""
    print("=" * 70)
    print("Kaiser RSD vs Velocileptors 比較プロット")
    print("=" * 70)
    
    # ファイルパスを設定
    base_dir = "./projects/mg_hsc_y3/data"
    kaiser_file = os.path.join(base_dir, "mg_hsc_y3_theory_kaiser_rsd_boss_dr12.modelvector")
    kaiser_dataset = os.path.join(base_dir, "boss_dr12_cmass_kaiser_rsd.dataset")
    velo_file = os.path.join(base_dir, "velocileptors_pk_multipoles.txt")  # GR版を使用
    
    # ファイルの存在確認
    if not os.path.exists(kaiser_file):
        print(f"❌ Kaiser RSDファイルが見つかりません: {kaiser_file}")
        print("   → cobaya-run EXAMPLE_EVALUATE_kaiser_RSD.yaml を実行してください")
        return
    
    if not os.path.exists(velo_file):
        print(f"❌ Velocileptorsファイルが見つかりません: {velo_file}")
        print("   → cobaya-run EXAMPLE_EVALUATE_velocileptors.yaml を実行してください")
        return
    
    # データを読み込む
    print("\n📊 データを読み込んでいます...")
    k_kaiser, z_bins, P0_kaiser, P2_kaiser = load_kaiser_data(kaiser_file, kaiser_dataset)
    k_velo, P0_velo, P2_velo, P4_velo = load_velocileptors_data(velo_file)
    
    print(f"   Kaiser RSD: Nk={len(k_kaiser)}, Nz={len(z_bins)}")
    print(f"   Kaiser k range: {k_kaiser[0]:.4f} - {k_kaiser[-1]:.4f} h/Mpc")
    print(f"   Kaiser z range: {z_bins[0]:.4f} - {z_bins[-1]:.4f}")
    print(f"   Velocileptors: Nk={len(k_velo)}")
    print(f"   Velocileptors k range: {k_velo[0]:.4f} - {k_velo[-1]:.4f} h/Mpc")
    
    # Velocileptorsのz_effに最も近いz binを見つける
    # EXAMPLE_EVALUATE_velocileptors.yamlのz_effを確認
    velo_z_eff = 0.5893  # デフォルト値（YAMLファイルから読み取る必要がある場合は変更）
    z_idx = np.argmin(np.abs(z_bins - velo_z_eff))
    z_match = z_bins[z_idx]
    
    print(f"\n🎯 比較に使用するz値: Kaiser z={z_match:.4f} (index={z_idx}) ≈ Velocileptors z_eff={velo_z_eff}")
    
    # 単一z値での比較プロット
    print("\n📈 単一z値での比較プロットを作成しています...")
    output_file1 = os.path.join(base_dir, "comparison_kaiser_velocileptors_single_z.png")
    plot_comparison_single_z(
        k_kaiser, P0_kaiser[:, z_idx], P2_kaiser[:, z_idx],
        k_velo, P0_velo, P2_velo, P4_velo,
        z_match, output_file=output_file1
    )
    
    # 複数z値でのプロット
    print("\n📈 複数z値でのプロットを作成しています...")
    output_file2 = os.path.join(base_dir, "comparison_kaiser_velocileptors_multiple_z.png")
    plot_multiple_z(
        k_kaiser, z_bins, P0_kaiser, P2_kaiser,
        k_velo, P0_velo, P2_velo,
        output_file=output_file2
    )
    
    # 統計情報 (k*P(k) ベース)
    print("\n📊 統計情報 (k·P(k) ベース):")
    kP0_velo_interp = np.interp(k_kaiser, k_velo, k_velo * P0_velo)
    kP2_velo_interp = np.interp(k_kaiser, k_velo, k_velo * P2_velo)
    
    kP0_kaiser = k_kaiser * P0_kaiser[:, z_idx]
    kP2_kaiser = k_kaiser * P2_kaiser[:, z_idx]
    
    rel_diff_P0 = (kP0_velo_interp - kP0_kaiser) / kP0_kaiser * 100
    rel_diff_P2 = (kP2_velo_interp - kP2_kaiser) / np.abs(kP2_kaiser) * 100
    
    print(f"   k·P0相対差: 平均 {np.mean(rel_diff_P0):.2f}%, RMS {np.sqrt(np.mean(rel_diff_P0**2)):.2f}%")
    print(f"   k·P2相対差: 平均 {np.mean(rel_diff_P2):.2f}%, RMS {np.sqrt(np.mean(rel_diff_P2**2)):.2f}%")
    
    print("\n✅ 比較完了!")

if __name__ == "__main__":
    main()

